#ifndef FEATURE_RECONSTRUCTION_HPP
#define FEATURE_RECONSTRUCTION_HPP

#include <array>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsSmartPointer.hpp>

#include "MakeTrenchStamp.hpp"

template <typename NumericType, int D> class FeatureReconstruction {
  lsSmartPointer<lsDomain<NumericType, D>> levelset = nullptr;

  const std::array<NumericType, 3> &origin;
  const std::vector<NumericType> &sampleLocations;
  const std::vector<NumericType> &features;

  int verticalDir = D - 1;
  int horizontalDir = 0;

public:
  FeatureReconstruction(lsSmartPointer<lsDomain<NumericType, D>> passedLevelset,
                        const std::array<NumericType, 3> &passedOrigin,
                        const std::vector<NumericType> &passedSampleLocations,
                        const std::vector<NumericType> &passedFeatures)
      : levelset(passedLevelset), origin(passedOrigin),
        sampleLocations(passedSampleLocations), features(passedFeatures) {}

  void apply() {
    if (!levelset) {
      std::cout << "No levelset provided!\n";
      return;
    }

    auto stamp =
        MakeTrenchStamp(levelset->getGrid(), origin, sampleLocations, features);
    if (!stamp) {
      return;
    }

    // First: generate an initial plane from which we will remove the trench
    // geometry later on
    {
      NumericType normal[D] = {0.};
      normal[verticalDir] = 1.;

      auto plane = lsSmartPointer<lsPlane<NumericType, D>>::New(origin.data(), normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    // Second: remove the negative geometry
    lsBooleanOperation<NumericType, D>(
        levelset, stamp, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
};
#endif