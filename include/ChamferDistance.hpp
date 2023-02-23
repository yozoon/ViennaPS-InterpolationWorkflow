#ifndef CHAMFER_DISTANCE_HPP
#define CHAMFER_DISTANCE_HPP

#include <algorithm>
#include <array>

#include <psKDTree.hpp>

/**
 * Based on the chamfer distance score introduced in "A Point Set Generation
 * Network for 3D Object Reconstruction from a Single Image" by H. Fan et. al.
 * (https://arxiv.org/pdf/1612.00603.pdf)
 */
template <typename NumericType, class VectorType = std::array<NumericType, 3>,
          class LocatorType = psKDTree<NumericType>>
class ChamferDistanceScore {
  using PCType = std::vector<VectorType>;
  using PCVec = std::vector<std::vector<NumericType>>;

  const PCType &firstPointCloud;
  const PCType &secondPointCloud;

public:
  ChamferDistanceScore(const PCType &passedFirstPointCloud,
                       const PCType &passedSecondPointCloud)
      : firstPointCloud(passedFirstPointCloud),
        secondPointCloud(passedSecondPointCloud) {}

  NumericType calculate() {
    if (firstPointCloud.empty() || secondPointCloud.empty()) {
      lsMessage::getInstance()
          .addWarning("At least one of the meshes provided to "
                      "ChamferDistance is empty.")
          .print();
      return std::numeric_limits<NumericType>::infinity();
    }

    PCVec firstPointVec;
    firstPointVec.reserve(firstPointCloud.size());
    std::transform(firstPointCloud.cbegin(), firstPointCloud.cend(),
                   std::back_inserter(firstPointVec), [](const auto &pt) {
                     return std::vector<NumericType>(pt.begin(), pt.end());
                   });

    LocatorType firstLocator(firstPointVec);
    firstLocator.build();

    PCVec secondPointVec;
    secondPointVec.reserve(secondPointCloud.size());
    std::transform(secondPointCloud.cbegin(), secondPointCloud.cend(),
                   std::back_inserter(secondPointVec), [](const auto &pt) {
                     return std::vector<NumericType>(pt.begin(), pt.end());
                   });

    LocatorType secondLocator(secondPointVec);
    secondLocator.build();

    NumericType sum1 = 0.;

#pragma omp parallel for default(shared) reduction(+ : sum1)
    for (const auto &node : secondPointVec) {
      auto nearestOpt = firstLocator.findNearest(node);
      sum1 += nearestOpt ? nearestOpt.value().second : 0.;
    }
    sum1 /= secondPointVec.size();

    NumericType sum2 = 0.;
#pragma omp parallel for default(shared) reduction(+ : sum2)
    for (const auto &node : firstPointVec) {
      auto nearestOpt = secondLocator.findNearest(node);
      sum2 += nearestOpt ? nearestOpt.value().second : 0.;
    }
    sum2 /= firstPointVec.size();

    return sum1 + sum2;
  }
};
#endif