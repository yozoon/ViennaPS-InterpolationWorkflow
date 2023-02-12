#ifndef FEATURE_EXTRACTION_HPP
#define FEATURE_EXTRACTION_HPP

#include <algorithm>
#include <vector>

#include <lsToDiskMesh.hpp>

#include <psDomain.hpp>
#include <psKDTree.hpp>

template <typename NumericType, int D> class FeatureExtraction {
public:
  using ConstPtr = psSmartPointer<const std::vector<NumericType>>;

  FeatureExtraction()
      : sampleLocations(
            distributeSampleLocations(numberOfSamples, -edgeAffinity)) {}

  FeatureExtraction(psSmartPointer<psDomain<NumericType, D>> passedDomain)
      : domain(passedDomain), sampleLocations(distributeSampleLocations(
                                  numberOfSamples, -edgeAffinity)) {}

  void setDomain(psSmartPointer<psDomain<NumericType, D>> passedDomain) {
    domain = passedDomain;
  }

  void setNumberOfSamples(int passedNumberOfSamples) {
    assert(numberOfSamples > 1);
    numberOfSamples = passedNumberOfSamples;
    // The first sample location is negative - indicating that this is the
    // height measurement
    sampleLocations = {-1.};
    // The remaining sample locations are distributed in the range 0 to 1.
    auto sampleDistrinution =
        distributeSampleLocations(numberOfSamples - 1, -edgeAffinity);
    std::copy(sampleDistrinution.begin(), sampleDistrinution.end(),
              std::back_inserter(sampleLocations));
  }

  void setEdgeAffinity(NumericType passedEdgeAffinity) {
    edgeAffinity = passedEdgeAffinity;
  }

  ConstPtr getFeatures() const { return ConstPtr::New(features); }

  ConstPtr getSampleLocations() const { return ConstPtr::New(sampleLocations); }

  void apply() {
    if (!domain)
      return;

    // Re-initialize the feature vector with a value of zero.
    features.clear();
    features.resize(sampleLocations.size(), 0.);

    // Convert the geometry to a disk mesh and extract the nodes
    auto mesh = psSmartPointer<lsMesh<>>::New();
    lsToDiskMesh<NumericType, D>(domain->getLevelSets()->back(), mesh).apply();
    auto nodes = mesh->getNodes();

    // Extract the depth of the trench
    auto [min, max] = getMinMax(nodes, 1 /* axis */);
    NumericType depth = max - min;
    features[0] = depth;

    // Only use the vertical trench axis for the kDtree
    std::vector<std::vector<NumericType>> ys;
    ys.reserve(nodes.size());
    std::transform(
        nodes.begin(), nodes.end(), std::back_inserter(ys),
        [=](auto &node) { return std::vector<NumericType>{node[D - 1]}; });

    // Initialize a 1D KDTree (~ multiset data structure with convencience
    // functions)
    psKDTree<NumericType> tree;
    tree.setPoints(ys);
    tree.build();

    // The extract the diameters along its depth at the relative coordinates
    // given by depths
    int i = 1;
    const NumericType gridDelta = domain->getGrid().getGridDelta();
    for (auto sl : sampleLocations) {
      if (sl < 0)
        continue;

      std::vector<NumericType> loc = {max - depth * sl};

      auto neighborsOpt = tree.findNearestWithinRadius(loc, gridDelta / 2);
      if (!neighborsOpt) {
        ++i;
        continue;
      }
      auto neighbors = neighborsOpt.value();

      // auto [minIt, maxIt] =
      //     std::minmax_element(neighbors->begin(), neighbors->end(),
      //                         [&](const auto &a, const auto &b) {
      //                           return nodes[a.first][0] <
      //                           nodes[b.first][0];
      //                         });
      // if (minIt == neighbors->end() || maxIt == neighbors->end()) {
      //   ++i;
      //   continue;
      // }

      // int idxL = minIt->first;
      // int idxR = maxIt->first;
      // const auto d = std::abs(nodes[idxL][0]) + std::abs(nodes[idxR][0]);
      // dimensions->at(i) = d;

      // Here we assume that the trench is centered at zero and symmetric with
      // its vertical axis being the axis of symmetry int idxL = -1;
      int idxL = -1;
      int idxR = -1;
      for (auto &nb : neighbors) {
        // if the point is on the left trench sidewall
        if (idxL < 0 && nodes[nb.first][0] < 0) {
          idxL = nb.first;
        }

        // if the point is on the right trench sidewall
        if (idxR < 0 && nodes[nb.first][0] >= 0) {
          idxR = nb.first;
        }

        // if both indices were found
        if (idxL >= 0 && idxR >= 0)
          break;
      }

      if (idxL >= 0 && idxR >= 0) {
        const auto d = std::abs(nodes[idxL][0]) + std::abs(nodes[idxR][0]);
        features[i] = d;
      }

      ++i;
    }
  }

private:
  static std::vector<NumericType>
  distributeSampleLocations(int n, NumericType alpha = 1.0) {
    std::vector<NumericType> x;
    x.reserve(n);
    for (int i = 0; i < n; ++i)
      x.emplace_back(-1.0 + 2.0 * i / (n - 1));

    if (alpha != 0) {
      std::transform(x.begin(), x.end(), x.begin(), [alpha](NumericType xi) {
        return xi < 0 ? 1 - std::exp(-alpha * xi) : std::exp(alpha * xi) - 1;
      });
      std::transform(
          x.begin(), x.end(), x.begin(),
          [maxVal = x.back()](NumericType xi) { return xi / maxVal; });
    }
    std::transform(x.begin(), x.end(), x.begin(),
                   [](NumericType xi) { return (xi + 1) / 2; });
    if (alpha < 0) {
      std::reverse(x.begin(), x.end());

      return x;
    } else {
      return x;
    }
  }

  std::tuple<NumericType, NumericType>
  getMinMax(const std::vector<std::array<NumericType, 3>> &nodes, int axis) {
    if (nodes.empty())
      return {};

    const auto [minposIter, maxposIter] = std::minmax_element(
        nodes.cbegin(), nodes.cend(),
        [&](const auto &a, const auto &b) { return a.at(axis) < b.at(axis); });

    return {minposIter->at(axis), maxposIter->at(axis)};
  }

  psSmartPointer<psDomain<NumericType, D>> domain = nullptr;

  int numberOfSamples = 10;
  NumericType edgeAffinity = 3.;

  std::vector<NumericType> features;
  std::vector<NumericType> sampleLocations;
};
#endif