#ifndef FEATURE_EXTRACTION_HPP
#define FEATURE_EXTRACTION_HPP

#include <algorithm>
#include <array>
#include <vector>

#include <lsCalculateNormalVectors.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>

#include <psDomain.hpp>
#include <psKDTree.hpp>

template <typename NumericType, int D> class FeatureExtraction {
  enum FeatureLabelEnum : unsigned {
    NONE = 0,
    LEFT_SIDEWALL = 1,
    RIGHT_SIDEWALL = 2,
    TRENCH_BOTTOM = 3,
  };

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

  void setOrigin(const std::array<NumericType, 3> &passedOrigin) {
    origin = passedOrigin;
  }

  ConstPtr getFeatures() const { return ConstPtr::New(features); }

  ConstPtr getSampleLocations() const { return ConstPtr::New(sampleLocations); }

  void apply() {
    initializeSampleLocations();

    if (!domain)
      return;

    const NumericType gridDelta = domain->getGrid().getGridDelta();

    // Re-initialize the feature vector with a value of zero.
    features.clear();
    features.resize(sampleLocations.size(), 0.);

    auto &levelset = domain->getLevelSets()->back();

    // lsCalculateNormalVectors<NumericType, D>(levelset).apply();

    // Check if the trench is pinched off
    lsMarkVoidPoints<NumericType, D>(levelset).apply();

    lsExpand<NumericType, D>(levelset, 4).apply();

    // Convert the geometry to a surface mesh and extract the nodes as well as
    // the void points
    auto mesh = psSmartPointer<lsMesh<>>::New();
    lsToDiskMesh<NumericType, D>(levelset, mesh).apply();
    auto nodes = mesh->getNodes();

    auto normalsPtr = mesh->getCellData().getVectorData(
        lsCalculateNormalVectors<NumericType, D>::normalVectorsLabel);
    if (normalsPtr == nullptr) {
      return;
    }

    auto normals = *normalsPtr;

    std::vector<NumericType> featureLabels(nodes.size(), 0.);
    for (unsigned i = 0; i < normals.size(); ++i) {
      auto &normal = normals[i];
      auto nx = normal[horizontalDir];
      NumericType threshold = 0.2;
      NumericType label = FeatureLabelEnum::NONE;
      if (nx >= threshold) {
        label = FeatureLabelEnum::LEFT_SIDEWALL;
      } else if (nx <= -threshold) {
        label = FeatureLabelEnum::RIGHT_SIDEWALL;
      } else if (nodes[i][verticalDir] < origin[verticalDir] - gridDelta) {
        label = FeatureLabelEnum::TRENCH_BOTTOM;
      }
      featureLabels[i] = label;
    }

#ifndef NDEBUG
    mesh->getPointData().insertNextScalarData(featureLabels, "feature");
    lsVTKWriter<NumericType>(mesh, "FeatureExtractionLabels.vtp").apply();
#endif

    // Extract the depth of the trench
    auto [min, max] = getMinMax(nodes, verticalDir);

    // // Determine the height of the pinchoff (if one occured)
    // auto voidPoints = mesh->getPointData().getScalarData(
    //     lsMarkVoidPoints<NumericType, D>::voidPointLabel);
    // bool isClosed = std::any_of(voidPoints->begin(), voidPoints->end(),
    //                             [](int p) { return p == 1; });
    // NumericType depthOffset = 0.;
    // if (isClosed) {
    //   NumericType pinchoffHeight = min;
    //   for (unsigned i = 0; i < nodes.size(); ++i) {
    //     if (voidPoints->at(i) != 0 && nodes[i][verticalDir] > pinchoffHeight)
    //       pinchoffHeight = nodes[i][verticalDir];
    //   }
    //   depthOffset = max - pinchoffHeight;
    //   // max = pinchoffHeight;
    //   std::cout << "Pinchoff at: " << pinchoffHeight << '\n';
    // }

    NumericType depth = max - min;
    features[0] = depth;
    // features[1] = depthOffset;

    // Project all points onto the hyperplane spanned by all axis except the
    // one of the direction of the trench diameter. (i.e. project on the
    // vertical axis in 2D and in 3D project onto the symmetry plane spanned by
    // the vertical axis and the extrusion axis).
    std::vector<std::vector<NumericType>> vertical;
    vertical.reserve(nodes.size());
    std::transform(
        nodes.begin(), nodes.end(), std::back_inserter(vertical),
        [=](auto &node) {
          if constexpr (D == 3) {
            return std::vector<NumericType>{node[verticalDir], node[D - 2]};
          } else {
            return std::vector<NumericType>{node[verticalDir]};
          }
        });

    // Initialize a 1D KDTree (~= multiset data structure with convencience
    // functions)
    psKDTree<NumericType> tree;
    tree.setPoints(vertical);
    tree.build();

    // The extract the diameters along its depth at the relative coordinates
    // given by depths
    for (unsigned i = 1; i < sampleLocations.size(); ++i) {
      std::vector<NumericType> loc = {min + depth * sampleLocations[i]};
      if constexpr (D == 3) {
        loc.push_back(origin[D - 2]);
      }

      auto neighborsOpt = tree.findNearestWithinRadius(loc, gridDelta / 2);
      if (!neighborsOpt)
        continue;

      auto neighbors = neighborsOpt.value();

      if (neighbors.empty())
        continue;

      if (i < numSamplesRight + 1) { // Right sidewall
        int idx = -1;
        for (auto &nb : neighbors) {
          NumericType label = featureLabels[nb.first];
          // NumericType horizontalDistance =
          //     std::abs(nodes[nb.first][horizontalDir] -
          //     origin[horizontalDir]);
          if (label == FeatureLabelEnum::RIGHT_SIDEWALL) {
            idx = nb.first;
            break;
          }
        }
        if (idx >= 0)
          features[i] = nodes[idx][horizontalDir] - origin[horizontalDir];
      } else { // Left sidewall
        int idx = -1;
        for (auto &nb : neighbors) {
          NumericType label = featureLabels[nb.first];
          // NumericType horizontalDistance =
          //     std::abs(origin[horizontalDir] -
          //     nodes[nb.first][horizontalDir]);
          if (label == FeatureLabelEnum::LEFT_SIDEWALL) {
            idx = nb.first;
            break;
          }
        }
        if (idx >= 0)
          features[i] = nodes[idx][horizontalDir] - origin[horizontalDir];
      }
    }
  }

  // Distribute n points in the range from 0 to 1 and place them closer to the
  // edges or closer to the center based on the edgeAffinity parameter.
  void
  distributeSampleLocations(typename std::vector<NumericType>::iterator start,
                            typename std::vector<NumericType>::iterator end,
                            NumericType edgeAffinity = 0.0,
                            bool ascending = true, bool closed = true) const {
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
      // or increase density at the edges (if edgeAffinity > 0)
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

    if (ascending) {
      if (edgeAffinity > 0)
        std::transform(start, end, start, [](auto &v) { return 1.0 - v; });
    } else {
      if (edgeAffinity < 0)
        std::transform(start, end, start, [](auto &v) { return 1.0 - v; });
    }
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
    // Sample locations are in the range 0 (bottom) ... 1 (top)
    if (!sampleLocations.empty())
      return;

    // The first sample location is negative - indicating that this is the
    // height measurement
    sampleLocations.clear();
    sampleLocations.resize(numberOfSamples, 0.0);
    sampleLocations[0] = -1.;

    // Ensure that the number of samples of the right sidewall is greater than
    // or equal to the number of samples of the left sidewall.
    numSamplesRight =
        static_cast<unsigned>(std::ceil(1.0 * (numberOfSamples - 1) / 2));

    // The remaining sample locations are distributed in the range 0 to 1.
    // Left sidewall sample locations (top to bottom -> descending)
    distributeSampleLocations(
        std::next(sampleLocations.begin(), 1),
        std::next(sampleLocations.begin(), numSamplesRight + 1), edgeAffinity,
        /*descending*/ false, closed);

    // Right sidewall sample locations (top to bottom -> descending)
    distributeSampleLocations(
        std::next(sampleLocations.begin(), numSamplesRight + 1),
        sampleLocations.end(), edgeAffinity, /*descending*/ false, closed);
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
  unsigned numSamplesRight;
};
#endif