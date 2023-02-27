#ifndef ADVECTION_CALLBACK_HPP
#define ADVECTION_CALLBACK_HPP

#include <algorithm>
#include <vector>

#include <psAdvectionCallback.hpp>
#include <psCSVWriter.hpp>
#include <psSmartPointer.hpp>

template <typename NumericType, int D, typename FeatureExtractionType,
          typename WriterType>
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

  void setWriter(psSmartPointer<WriterType> passedWriter) {
    writer = passedWriter;
  }

  void setPrefixData(const std::vector<NumericType> &passedPrefixData) {
    prefixData = passedPrefixData;
  }

  void setDataPtr(psSmartPointer<std::vector<NumericType>> passedDataPtr) {
    dataPtr = passedDataPtr;
  }

  void apply() {
    if (!featureExtraction)
      return;

    featureExtraction->setDomain(domain->getLevelSets()->back());
    featureExtraction->apply();

    auto features = featureExtraction->getFeatures();
    if (features) {
      std::vector<NumericType> row(prefixData.begin(), prefixData.end());
      row.push_back(counter);
      std::copy(features->begin(), features->end(), std::back_inserter(row));
      if (writer)
        writer->writeRow(row);
      if (dataPtr)
        std::copy(row.begin(), row.end(), std::back_inserter(*dataPtr));
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
  std::vector<NumericType> prefixData;

  psSmartPointer<FeatureExtractionType> featureExtraction = nullptr;
  psSmartPointer<WriterType> writer = nullptr;
  psSmartPointer<std::vector<NumericType>> dataPtr = nullptr;
};
#endif