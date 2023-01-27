#ifndef ADVECTION_CALLBACK_HPP
#define ADVECTION_CALLBACK_HPP

#include <iomanip>
#include <iostream>

#include <psAdvectionCallback.hpp>
#include <psCSVWriter.hpp>
#include <psSmartPointer.hpp>

template <typename NumericType, int D, typename FeatureExtractionType>
class AdvectionCallback : public psAdvectionCalback<NumericType, D> {
protected:
  using psAdvectionCalback<NumericType, D>::domain;

public:
  AdvectionCallback() : deltaT(0.5) {}

  AdvectionCallback(NumericType passedDeltaT) : deltaT(passedDeltaT) {}

  void setFeatureExtraction(
      psSmartPointer<FeatureExtractionType> passedFeatureExtraction) {
    featureExtraction = passedFeatureExtraction;
  }

  void setDataPtr(psSmartPointer<std::vector<NumericType>> passedDataPtr) {
    dataPtr = passedDataPtr;
  }

  void apply() {
    if (!featureExtraction)
      return;

    if (dataPtr) {
      featureExtraction->setDomain(domain);
      featureExtraction->apply();

      auto features = featureExtraction->getFeatures();
      if (features) {
        dataPtr->push_back(processTime / deltaT);
        std::copy(features->begin(), features->end(),
                  std::back_inserter(*dataPtr));
      }
    }
  }

  void applyPreAdvect(const NumericType passedProcessTime) override {
    if (passedProcessTime == 0.) {
      apply();
      ++counter;
      lastUpdateTime = 0.;
    }
    processTime = passedProcessTime;
  }

  void applyPostAdvect(const NumericType advectionTime) override {
    processTime += advectionTime;
    if (processTime - lastUpdateTime >= deltaT) {
      apply();
      lastUpdateTime = counter * deltaT;
      ++counter;
    }
  }

private:
  NumericType deltaT = 0.5;

  NumericType processTime = 0.0;
  NumericType lastUpdateTime = 0.0;
  size_t counter = 0;

  psSmartPointer<FeatureExtractionType> featureExtraction = nullptr;
  psSmartPointer<std::vector<NumericType>> dataPtr = nullptr;
};
#endif