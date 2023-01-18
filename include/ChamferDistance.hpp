#ifndef CHAMFER_DISTANCE_HPP
#define CHAMFER_DISTANCE_HPP

#include <array>

#include <psKDTree.hpp>
#include <psPointLocator.hpp>
#include <psSmartPointer.hpp>

/**
 * Based on the chamfer distance score introduced in "A Point Set Generation
 * Network for 3D Object Reconstruction from a Single Image" by H. Fan et. al.
 * (https://arxiv.org/pdf/1612.00603.pdf)
 */
template <typename NumericType, class VectorType = std::array<NumericType, 3>,
          class LocatorType =
              psKDTree<NumericType, std::tuple_size_v<VectorType>>>
class ChamferDistanceScore {
  static_assert(
      std::is_base_of_v<
          psPointLocator<NumericType, std::tuple_size_v<VectorType>>,
          LocatorType>,
      "The passed point locator is not a subclass of psPointLocator.");

  psSmartPointer<std::vector<VectorType>> firstPointCloud = nullptr;
  psSmartPointer<LocatorType> firstLocator = nullptr;
  bool buildLocator = true;

public:
  ChamferDistanceScore() {}
  ChamferDistanceScore(
      psSmartPointer<std::vector<VectorType>> passedFirstPointCloud)
      : firstPointCloud(passedFirstPointCloud) {}

  void setFirstPointCloud(
      psSmartPointer<std::vector<VectorType>> passedFirstPointCloud) {
    firstPointCloud = passedFirstPointCloud;
    buildLocator = true;
  }

  NumericType
  calculate(psSmartPointer<std::vector<VectorType>> secondPointCloud) {
    if (firstPointCloud == nullptr || secondPointCloud == nullptr) {
      lsMessage::getInstance()
          .addWarning("At least one of the meshes provided to "
                      "cmChamferDistance is a nullptr.")
          .print();
      return std::numeric_limits<NumericType>::infinity();
    }

    // Only build the first locator if it hasn't been build yet
    if (buildLocator) {
      buildLocator = false;
      firstLocator = psSmartPointer<LocatorType>::New(*firstPointCloud);
      firstLocator->build();
    }

    auto secondLocator = psSmartPointer<LocatorType>::New(*secondPointCloud);
    secondLocator->build();

    NumericType sum1 = 0.;

#pragma omp parallel for default(shared) reduction(+ : sum1)
    for (const auto &node : *secondPointCloud) {
      auto nearestOpt = firstLocator->findNearest(node);
      sum1 += nearestOpt.has_value() ? nearestOpt.value().second : 0.;
    }
    sum1 /= secondPointCloud->size();

    NumericType sum2 = 0.;
#pragma omp parallel for default(shared) reduction(+ : sum2)
    for (const auto &node : *firstPointCloud) {
      auto nearestOpt = secondLocator->findNearest(node);
      sum2 += nearestOpt.has_value() ? nearestOpt.value().second : 0.;
    }
    sum2 /= firstPointCloud->size();

    return sum1 + sum2;
  }
};
#endif