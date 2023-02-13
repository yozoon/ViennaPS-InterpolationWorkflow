#ifndef FEATURE_EXTRACTION_HPP
#define FEATURE_EXTRACTION_HPP

#include <algorithm>
#include <array>
#include <vector>

#include <lsToDiskMesh.hpp>

#include <psDomain.hpp>
#include <psKDTree.hpp>

#define SYMMETRIC

template <typename NumericType, int D> class FeatureExtraction {
public:
  using ConstPtr = psSmartPointer<const std::vector<NumericType>>;

  FeatureExtraction() {}

  FeatureExtraction(psSmartPointer<psDomain<NumericType, D>> passedDomain)
      : domain(passedDomain) {}

  void setDomain(psSmartPointer<psDomain<NumericType, D>> passedDomain) {
    domain = passedDomain;
  }

  void setNumberOfSamples(int passedNumberOfSamples, bool passedClosed = true) {
    assert(numberOfSamples > 1);
    numberOfSamples = passedNumberOfSamples;
    closed = passedClosed;
    sampleLocations.clear();
  }

  void setEdgeAffinity(NumericType passedEdgeAffinity) {
    edgeAffinity = passedEdgeAffinity;
    sampleLocations.clear();
  }

  ConstPtr getFeatures() const { return ConstPtr::New(features); }

  ConstPtr getSampleLocations() const { return ConstPtr::New(sampleLocations); }

  void setOrigin(const std::array<NumericType, 3> &passedOrigin) {
    origin = passedOrigin;
  }

  void apply() {
    initializeSampleLocations();

    if (!domain)
      return;

    // Re-initialize the feature vector with a value of zero.
    features.clear();
#ifdef SYMMETRIC
    features.resize(sampleLocations.size(), 0.);
#else
    features.resize(2 * sampleLocations.size() - 1, 0.);
#endif

    // Convert the geometry to a disk mesh and extract the nodes
    auto mesh = psSmartPointer<lsMesh<>>::New();
    lsToDiskMesh<NumericType, D>(domain->getLevelSets()->back(), mesh).apply();
    auto nodes = mesh->getNodes();

    // Extract the depth of the trench
    auto [min, max] = getMinMax(nodes, verticalDir);
    NumericType depth = max - min;
    features[0] = depth;

    // Only use the vertical trench axis for the kDtree
    std::vector<std::vector<NumericType>> vertical;
    vertical.reserve(nodes.size());
    std::transform(nodes.begin(), nodes.end(), std::back_inserter(vertical),
                   [=](auto &node) {
                     return std::vector<NumericType>{node[verticalDir]};
                   });

    // Initialize a 1D KDTree (~= multiset data structure with convencience
    // functions)
    psKDTree<NumericType> tree;
    tree.setPoints(vertical);
    tree.build();

    // The extract the diameters along its depth at the relative coordinates
    // given by depths
    const NumericType gridDelta = domain->getGrid().getGridDelta();
#ifdef SYMMETRIC
    for (unsigned i = 1; i < sampleLocations.size(); ++i) {
      std::vector<NumericType> loc = {max - depth * sampleLocations[i]};

      // Extraction that assumes that the trench is symmetric
      auto neighborsOpt = tree.findNearestWithinRadius(loc, gridDelta / 2);
      if (!neighborsOpt)
        continue;

      auto neighbors = neighborsOpt.value();

      // Here we assume that the trench is centered at the origin and symmetric.
      // with its vertical axis being the axis of symmetry int idxL = -1;
      int idxL = -1;
      int idxR = -1;
      for (auto &nb : neighbors) {
        // if the point is on the left trench sidewall
        if (idxL < 0 &&
            nodes[nb.first][horizontalDir] - origin[horizontalDir] < 0) {
          idxL = nb.first;
        }

        // if the point is on the right trench sidewall
        if (idxR < 0 &&
            nodes[nb.first][horizontalDir] - origin[horizontalDir] >= 0) {
          idxR = nb.first;
        }

        // if both indices were found
        if (idxL >= 0 && idxR >= 0)
          break;
      }

      if (idxL >= 0 && idxR >= 0) {
        const auto d =
            std::abs(nodes[idxL][horizontalDir] - origin[horizontalDir]) +
            std::abs(nodes[idxR][horizontalDir] - origin[horizontalDir]);
        features[i] = d;
      }
    }
#else
    for (unsigned i = 1; i < sampleLocations.size(); ++i) {
      std::vector<NumericType> loc = {max - depth * sampleLocations[i]};

      // Extraction that assumes that the trench is symmetric
      auto neighborsOpt = tree.findNearestWithinRadius(loc, gridDelta / 2);
      if (!neighborsOpt)
        continue;

      auto neighbors = neighborsOpt.value();

      // Here we assume that the trench is centered at the origin and symmetric.
      // with its vertical axis being the axis of symmetry int idxL = -1;
      int idxL = -1;
      int idxR = -1;
      for (auto &nb : neighbors) {
        // if the point is on the left trench sidewall
        if (idxL < 0 &&
            nodes[nb.first][horizontalDir] - origin[horizontalDir] < 0) {
          idxL = nb.first;
        }

        // if the point is on the right trench sidewall
        if (idxR < 0 &&
            nodes[nb.first][horizontalDir] - origin[horizontalDir] >= 0) {
          idxR = nb.first;
        }

        // if both indices were found
        if (idxL >= 0 && idxR >= 0)
          break;
      }

      if (idxL >= 0 && idxR >= 0) {
        features[2 * i - 1] =
            std::abs(nodes[idxL][horizontalDir] - origin[horizontalDir]);
        features[2 * i] =
            std::abs(nodes[idxR][horizontalDir] - origin[horizontalDir]);
      }
    }
#endif
  }

  // Distribute n points in the range from 0 to 1 and place them closer to the
  // edges or closer to the center based on the edgeAffinity parameter.
  void
  distributeSampleLocations(typename std::vector<NumericType>::iterator start,
                            typename std::vector<NumericType>::iterator end,
                            NumericType edgeAffinity = 0.0,
                            bool closed = true) const {
    auto n = std::distance(start, end);
    if (n < 1)
      return;

    // Generate n evenly spaced points in the interval [-1, 1] (or (-1, 1) if
    // closed==false)
    if (closed) {
      std::generate(start, end, [=, i = 0]() mutable {
        return -1.0 + 2.0 * i++ / (n - 1);
      });
    } else {
      std::generate(start, end, [=, i = 0]() mutable {
        return -1.0 + 2.0 * (i++ + 0.5) / n;
      });
    }

    NumericType maxVal = 1.0;
    if (edgeAffinity != 0.0) {
      // Spread the points. Increase density around zero (if edgeAffinity < 0)
      // or increase density at the boundaries (if edgeAffinity > 0)
      std::transform(start, end, start, [edgeAffinity](NumericType xi) {
        return xi < 0 ? 1.0 - std::exp(edgeAffinity * xi)
                      : std::exp(-edgeAffinity * xi) - 1.0;
      });
      maxVal = std::abs(std::exp(-edgeAffinity) - 1.);
    }
    // Now transform the points back into the interval [0,1] (or (0,1) if
    // closed==False)
    std::transform(start, end, start,
                   [=](NumericType xi) { return (xi / maxVal + 1.0) / 2.0; });

    if (edgeAffinity < 0)
      std::reverse(start, end);
  }

  std::tuple<NumericType, NumericType>
  getMinMax(const std::vector<std::array<NumericType, 3>> &nodes,
            int axis) const {
    if (nodes.empty())
      return {};

    const auto [minposIter, maxposIter] = std::minmax_element(
        nodes.cbegin(), nodes.cend(),
        [&](const auto &a, const auto &b) { return a.at(axis) < b.at(axis); });

    return {minposIter->at(axis), maxposIter->at(axis)};
  }

  void initializeSampleLocations() {
    if (!sampleLocations.empty())
      return;

    // The first sample location is negative - indicating that this is the
    // height measurement
    sampleLocations.clear();
    sampleLocations.resize(numberOfSamples);
    sampleLocations[0] = -1.;

    // The remaining sample locations are distributed in the range 0 to 1.
    distributeSampleLocations(std::next(sampleLocations.begin(), 1),
                              sampleLocations.end(), edgeAffinity, closed);
  }

private:
  psSmartPointer<psDomain<NumericType, D>> domain = nullptr;

  std::array<NumericType, 3> origin{0.};
  int numberOfSamples = 32;
  NumericType edgeAffinity = 0.;
  bool closed = false;

  int verticalDir = D - 1;
  int horizontalDir = 0;

  std::vector<NumericType> sampleLocations;
  std::vector<NumericType> features;
};
#endif