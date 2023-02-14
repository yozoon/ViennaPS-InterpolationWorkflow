#include <fmt/color.h>
#include <fmt/core.h>

#include "FeatureExtraction.hpp"
#include "Parameters.hpp"
#include "TrenchDeposition.hpp"

int main(int argc, const char *const *const) {
  using NumericType = double;
  static constexpr int D = 2;
  static constexpr int numberOfSamples = 60;

  // Whether to concatenate data of different timesteps into one row, or record
  // it in separate rows
  bool concatenateData = false;
  if (argc > 1)
    concatenateData = true;

  // How long to run the process and at which intervals to do the extraction
  static constexpr NumericType processDuration = 5.0;
  static constexpr NumericType extractionInterval = 1.0;

  // The parameters we are interested in
  std::vector<NumericType> stickingProbabilities = {1., 0.7, 0.4, 0.1};
  std::vector<NumericType> taperAngles = {-15., -10, -5., 0., 5., 10., 15.};

  int concatenatedDataDimension =
      (numberOfSamples + 1) *
      (static_cast<int>(std::ceil(processDuration)) + 1);

  // The data we are going to store consists of stickingProbability,
  // taperAngle, time and the sampled geometry descriptors as provided by the
  // feature extraction.
  int InputDimension = concatenateData ? 2 : 3;

  // Instantiate the featureExtraction
  auto featureExtraction =
      psSmartPointer<FeatureExtraction<NumericType, D>>::New();
  featureExtraction->setNumberOfSamples(numberOfSamples,
                                        false /* => open interval*/);
  // Uniformly spaced sample points
  featureExtraction->setEdgeAffinity(2.0);
  featureExtraction->initializeSampleLocations();

  // The locations at which the features are extracted (normalized to
  // the trench depth at each timestep)
  auto sampleLocations = featureExtraction->getSampleLocations();

  std::string filename = concatenateData ? "data_concat.csv" : "data.csv";
  auto writer = psSmartPointer<psCSVWriter<NumericType>>::New(filename);

  // Creation of a descriptive/ detailed header
  std::string header = "taperAngle,stickingProbability,";
  if (concatenateData) {
    header += "{extractionStep,depth,diameters...}";
  } else {
    header += "extractionStep,depth";
    for (unsigned i = 0; i < numberOfSamples - 1; ++i)
      header += ",diameter_" + std::to_string(i);
  }

  header += "\nData generated by simple trench deposition example.\nRelative "
            "locations of diameter measurements:";
  header += "\n!" + join(sampleLocations->begin(), sampleLocations->end());
  header += "\n!InputDimension=" + std::to_string(InputDimension) +
            ",ExtractionInterval=" + std::to_string(extractionInterval);

  writer->setHeader(header);
  writer->initialize();

  for (auto taperAngle : taperAngles) {
    Parameters<NumericType> params;
    params.taperAngle = taperAngle;
    params.processTime = processDuration;
    const NumericType offset =
        std::tan(params.taperAngle * rayInternal::PI / 180.) *
        params.trenchHeight;

    if (params.trenchWidth / 2. + offset <= params.gridDelta) {
      fmt::print(
          fmt::fg(fmt::color::yellow),
          "Warning: a trench with the provided height {:.3f} and taper angle "
          "{:+.3f} would have an initial opening size of less than one grid "
          "delta! Skipping.\n",
          params.trenchHeight, params.taperAngle);
      continue;
    }

    // Make sure that the trench sidewalls stay inside the simulation domain,
    // even if they are tapered.
    params.xExtent =
        2. * std::max(params.trenchWidth / 2 - offset + params.gridDelta,
                      params.xExtent / 2);

    for (auto stickingProbability : stickingProbabilities) {
      params.stickingProbability = stickingProbability;
      fmt::print("AR={:.3f}, taperAngle={:+.3f}, sticking probability={:.3f}\n",
                 params.trenchHeight / params.trenchWidth, params.taperAngle,
                 params.stickingProbability);

      // Using the advection callback, we can run the extraction at
      // certain pre-defined advection timesteps.
      NumericType deltaT = extractionInterval / params.stickingProbability;

      auto advectionCallback = psSmartPointer<AdvectionCallback<
          NumericType, D, decltype(featureExtraction)::element_type,
          decltype(writer)::element_type>>::New(deltaT);

      advectionCallback->setFeatureExtraction(featureExtraction);

      psSmartPointer<std::vector<NumericType>> dataPtr = nullptr;
      if (concatenateData) {
        dataPtr = decltype(dataPtr)::New();
        dataPtr->reserve(InputDimension + concatenatedDataDimension);
        dataPtr->push_back(params.taperAngle);
        dataPtr->push_back(params.stickingProbability);
        advectionCallback->setDataPtr(dataPtr);
      } else {
        advectionCallback->setPrefixData(std::vector<NumericType>{
            params.taperAngle, params.stickingProbability});
        advectionCallback->setWriter(writer);
      }

      auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
      psMakeTrench<NumericType, D>(geometry, params.gridDelta, params.xExtent,
                                   params.yExtent, params.trenchWidth,
                                   params.trenchHeight, params.taperAngle)
          .apply();

      executeProcess(geometry, params, advectionCallback);

      if (concatenateData) {
        if (dataPtr)
          writer->writeRow(*dataPtr);
      }

      writer->flush();
    }
  }

  return EXIT_SUCCESS;
}