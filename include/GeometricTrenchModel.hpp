#ifndef GEOMETRIC_TRENCH_DEPOSITION_MODEL_HPP
#define GEOMETRIC_TRENCH_DEPOSITION_MODEL_HPP

#include <array>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include "MakeTrenchStamp.hpp"
#include "SplineGridInterpolation.hpp"

template <typename NumericType, int D> class GeometricTrenchDepositionModel {
  using LSPtrType = lsSmartPointer<lsDomain<NumericType, D>>;

  LSPtrType levelset = nullptr;

  SplineGridInterpolation<NumericType> gridInterpolation;

  std::vector<NumericType> sampleLocations;
  std::array<NumericType, 3> origin;

  NumericType extractionInterval;

  bool advanceTopSurface = true;

  int horizontalDir = 0;
  int verticalDir = D - 1;
  bool initialized = false;

public:
  GeometricTrenchDepositionModel(
      LSPtrType passedLevelset,
      const std::vector<NumericType> &passedSampleLocations,
      const std::array<NumericType, 3> &passedOrigin,
      NumericType passedExtractionInterval = 1.0,
      bool passedAdvanceTopSurface = true)
      : levelset(passedLevelset), sampleLocations(passedSampleLocations),
        origin(passedOrigin), extractionInterval(passedExtractionInterval),
        advanceTopSurface(passedAdvanceTopSurface) {}

  void setLevelset(LSPtrType passedLevelset) { levelset = passedLevelset; }

  void
  setSampleLocations(const std::vector<NumericType> &passedSampleLocations) {
    sampleLocations = passedSampleLocations;
  }

  void setOrigin(const std::array<NumericType, 3> &passedOrigin) {
    origin = passedOrigin;
  }

  void setAdvanceTopSurface(bool passedAdvanceTopSurface) {
    advanceTopSurface = passedAdvanceTopSurface;
  }

  void setData(lsSmartPointer<const std::vector<std::vector<NumericType>>> data,
               int inputDim) {
    if (!data) {
      lsMessage::getInstance().addError("No data provided!").print();
      return;
    }
    if (data->empty()) {
      lsMessage::getInstance().addError("Provided data is empty!").print();
      return;
    }

    auto numFeatures = data->at(0).size() - inputDim;

    if (numFeatures < 1) {
      lsMessage::getInstance()
          .addError("Invalid input data dimension!")
          .print();
      return;
    }

    gridInterpolation.setDataDimensions(inputDim, numFeatures);
    gridInterpolation.setData(data);
    gridInterpolation.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    gridInterpolation.initialize();

    initialized = true;
  }

  void apply(const NumericType taperAngle,
             const NumericType stickingProbability,
             const NumericType processTime) {
    if (!levelset) {
      lsMessage::getInstance().addError("No levelset provided.").print();
      return;
    }

    if (!initialized) {
      lsMessage::getInstance()
          .addError("Interpolation does not have up-to-date data. Call "
                    "`setData` first.")
          .print();
    }

    std::vector<NumericType> evaluationPoint = {
        taperAngle, stickingProbability, processTime / extractionInterval};

    auto estimationOpt = gridInterpolation.estimate(evaluationPoint);
    if (!estimationOpt) {
      lsMessage::getInstance().addError("Value estimation failed.").print();
      return;
    }

    auto [estimatedFeatures, isInside] = estimationOpt.value();

    lsMessage::getInstance()
        .addDebug(
            std::string("Evaluation point within data grid boundaries: ") +
            (isInside ? std::string("true") : std::string("false")))
        .print();

    if (sampleLocations.size() != estimatedFeatures.size()) {
      lsMessage::getInstance()
          .addError("Mismatch of feature dimensions!")
          .print();
      return;
    }

    // First: generate an initial plane from which we will remove the trench
    // geometry later on
    if (advanceTopSurface) {
      NumericType normal[D] = {0.};
      normal[verticalDir] = 1.;

      auto plane =
          lsSmartPointer<lsPlane<NumericType, D>>::New(origin.data(), normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    // Second: Generate the trench stamp based on the interpolated features
    auto stamp = MakeTrenchStamp(levelset->getGrid(), origin, sampleLocations,
                                 estimatedFeatures);
    if (!stamp) {
      lsMessage::getInstance()
          .addError("Error while generating the trench stamp.")
          .print();
      return;
    }

    // Third: stamp out the trench
    lsBooleanOperation<NumericType, D>(
        levelset, stamp, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
};
#endif