#ifndef FEATURE_EXTRACTION_HPP
#define FEATURE_EXTRACTION_HPP

#include <algorithm>
#include <array>
#include <vector>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>

#ifndef NDEBUG
#include <lsVTKWriter.hpp>
#endif

#include "KDTree.hpp"
#include "span.hpp"

template <typename NumericType, int D> class FeatureExtraction {
  enum FeatureLabelEnum : unsigned {
    NONE = 0,
    LEFT_SIDEWALL = 1,
    RIGHT_SIDEWALL = 2,
    TRENCH_BOTTOM = 3,
  };

public:
  using ConstPtr = lsSmartPointer<const std::vector<NumericType>>;
  FeatureExtraction() {}

  FeatureExtraction(lsSmartPointer<lsDomain<NumericType, D>> passedDomain)
      : levelset(passedDomain) {}

  void setDomain(lsSmartPointer<lsDomain<NumericType, D>> passedDomain) {
    levelset = passedDomain;
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

    if (!levelset)
      return;

    const NumericType gridDelta = levelset->getGrid().getGridDelta();

    // Re-initialize the feature vector with a value of zero.
    features.clear();
    features.resize(sampleLocations.size(), 0.);

    // Convert the geometry to a surface mesh and extract the nodes as well as
    // the void points
    auto mesh = lsSmartPointer<lsMesh<>>::New();
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
      const auto &normal = normals[i];
      auto nx = normal[horizontalDir];
      NumericType threshold = 0.1;
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
      std::vector<NumericType> loc = {min + gridDelta / 2 +
                                      (depth - gridDelta) * sampleLocations[i]};
      if constexpr (D == 3) {
        loc.push_back(origin[D - 2]);
      }

      auto neighborsOpt = tree.findNearestWithinRadius(loc, 2.0 * gridDelta);
      if (!neighborsOpt)
        continue;

      auto neighbors = neighborsOpt.value();

      if (neighbors.empty())
        continue;

      FeatureLabelEnum matchLabel;
      if (i < numSamplesRight + 1) {
        matchLabel = FeatureLabelEnum::RIGHT_SIDEWALL;
      } else {
        matchLabel = FeatureLabelEnum::LEFT_SIDEWALL;
      }

      int upperIdx = -1;
      int lowerIdx = -1;
      NumericType upperDistance = std::numeric_limits<NumericType>::max();
      NumericType lowerDistance = std::numeric_limits<NumericType>::max();
      NumericType upperWidth = 0.;
      NumericType lowerWidth = 0.;
      for (const auto &nb : neighbors) {
        NumericType label = featureLabels[nb.first];
        NumericType nodeZ = nodes[nb.first][verticalDir];
        NumericType nodeX = nodes[nb.first][horizontalDir];
        if (label == matchLabel) {
          if (lowerIdx == -1 && nodeZ < loc[0]) {
            lowerIdx = nb.first;
            lowerDistance = loc[0] - nodeZ;
            lowerWidth = nodeX - origin[horizontalDir];
          } else if (upperIdx == -1 && nodeZ >= loc[0]) {
            upperIdx = nb.first;
            upperDistance = nodeZ - loc[0];
            upperWidth = nodeX - origin[horizontalDir];
          }
        }
        if (upperIdx != -1 && lowerIdx != -1)
          break;
      }

      if (upperIdx != 0 && upperDistance < 1e-4) {
        // If the vertical position of the upper point coincides with the sample
        // position up to a certain epsilon, use its width.
        features[i] = upperWidth;
      } else if (lowerIdx != 0 && lowerDistance < 1e-4) {
        // If the vertical position of the lower point coincides with the sample
        // position up to a certain epsilon, use its width.
        features[i] = lowerWidth;
      } else if (upperIdx != -1 && lowerIdx != -1) {
        // Otherwise linearly interpolate between the two widths based on the
        // offset to to sample position.
        NumericType totalDistance = lowerDistance + upperDistance;
        features[i] =
            (lowerDistance * upperWidth + upperDistance * lowerWidth) /
            totalDistance;
      }
    }
  }

  // Populate the given range with a sequence of values in the range from 0
  // to 1 and place them closer to the edges or closer to the center based on
  // the edgeAffinity parameter.
  void distributeSampleLocations(nonstd::span<NumericType> range,
                                 NumericType edgeAffinity = 0.0,
                                 bool ascending = true,
                                 bool closed = true) const {
    auto n = range.size();
    if (n < 1)
      return;

    // Generate n evenly spaced points in the interval [-1, 1] (or (-1, 1) if
    // closed==false)
    if (closed) {
      std::generate(range.begin(), range.end(), [=, i = 0]() mutable {
        return -1.0 + 2.0 * i++ / (n - 1);
      });
    } else {
      std::generate(range.begin(), range.end(), [=, i = 0]() mutable {
        return -1.0 + 2.0 * (i++ + 0.5) / n;
      });
    }

    NumericType maxVal = 1.0;
    if (edgeAffinity != 0.0) {
      // Spread the points. Increase density around zero (if edgeAffinity < 0)
      // or increase density at the edges (if edgeAffinity > 0)
      std::transform(range.begin(), range.end(), range.begin(),
                     [edgeAffinity](NumericType xi) {
                       return xi < 0 ? 1.0 - std::exp(edgeAffinity * xi)
                                     : std::expm1(-edgeAffinity * xi);
                     });
      maxVal = std::abs(std::expm1(-edgeAffinity));
    }
    // Now transform the points back into the interval [0,1] (or (0,1) if
    // closed==False)
    std::transform(range.begin(), range.end(), range.begin(),
                   [=](NumericType xi) { return (xi / maxVal + 1.0) / 2.0; });

    if (ascending) {
      if (edgeAffinity > 0)
        std::transform(range.begin(), range.end(), range.begin(),
                       [](const auto &v) { return 1.0 - v; });
    } else {
      if (edgeAffinity <= 0)
        std::transform(range.begin(), range.end(), range.begin(),
                       [](const auto &v) { return 1.0 - v; });
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
        nonstd::span(std::next(sampleLocations.begin(), 1),
                     std::next(sampleLocations.begin(), numSamplesRight + 1)),
        edgeAffinity,
        /*descending*/ false, closed);

    // Right sidewall sample locations (top to bottom -> descending)
    distributeSampleLocations(
        nonstd::span(std::next(sampleLocations.begin(), numSamplesRight + 1),
                     sampleLocations.end()),
        edgeAffinity, /*descending*/ false, closed);
  }

private:
  lsSmartPointer<lsDomain<NumericType, D>> levelset = nullptr;

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