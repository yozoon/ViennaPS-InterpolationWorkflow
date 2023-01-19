#include <chrono>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>

#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsToDiskMesh.hpp>

#include <psCSVDataSource.hpp>
#include <psDataScaler.hpp>
#include <psMakeTrench.hpp>
#include <psNearestNeighborsInterpolation.hpp>

#include "ChamferDistance.hpp"
#include "GeometryReconstruction.hpp"
#include "NaturalCubicSplineInterpolation.hpp"
#include "Parameters.hpp"
#include "TrenchDeposition.hpp"

namespace fs = std::filesystem;

template <typename NumericType, int D>
psSmartPointer<lsDomain<NumericType, D>>
createEmptyLevelset(const Parameters<NumericType> &params) {

  double bounds[2 * D];
  bounds[0] = -params.xExtent / 2.;
  bounds[1] = params.xExtent / 2.;

  if constexpr (D == 3) {
    bounds[2] = -params.yExtent / 2.;
    bounds[3] = params.yExtent / 2.;
    bounds[4] = -params.gridDelta;
    bounds[5] = params.gridDelta;
  } else {
    bounds[2] = -params.gridDelta;
    bounds[3] = params.gridDelta;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];

  for (int i = 0; i < D - 1; i++)
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  return psSmartPointer<lsDomain<NumericType, D>>::New(bounds, boundaryCons,
                                                       params.gridDelta);
}

int main(int argc, char *argv[]) {
  using NumericType = double;
  static constexpr int D = 2;

  using Clock = std::chrono::high_resolution_clock;
  using Duration = std::chrono::duration<double>;

  fs::path dataFile = "./data.csv";
  if (argc > 1)
    dataFile = fs::path{argv[1]};

  Parameters<NumericType> params;
  if (argc > 2) {
    auto config = psUtils::readConfigFile(argv[2]);
    if (config.empty()) {
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }

  // Interpolation based on previous simulation results
  auto start = Clock::now();
  auto interpolatedGeometry = psSmartPointer<psDomain<NumericType, D>>::New();
  {
    psCSVDataSource<NumericType> dataSource;
    dataSource.setFilename(dataFile.string());

    // Get a copy of the data from the data source
    auto data = dataSource.get();

    // Also read the data that was stored alongside the actual data describing
    // where it was sampled from (at which relative height)
    auto sampleLocations = dataSource.getPositionalParameters();

    // The number of samples is the number of sample locations plus one, since
    // we also have the total depth of the resulting trench in the dataset.
    int numberOfSamples = sampleLocations.size() + 1;

    // The input dimension is provided by the csv file named parameter
    // 'InputDim' In our case it is 2 since we use taper angle and sticking
    // probability as input parameters.
    int InputDim;
    auto namedParams = dataSource.getNamedParameters();
    if (auto id = namedParams.find("InputDimension"); id != namedParams.end()) {
      InputDim = std::round(id->second);
    } else {
      std::cout << "'InputDimension' not found in provided data CSV file.\n";
      return EXIT_FAILURE;
    }

    // The dimension of the data that is to be interpolated. In our case this
    // are the extracted dimensions at different timesteps (the timesteps are
    // also included in the data itself)
    int TargetDim = data->at(0).size() - InputDim;

    // The number of timesteps at which the dimensions were extracted.
    // Since the data is stored in the pattern (time+samples), we have to divide
    // by number of samples plus one.
    int numberOfTimesteps = TargetDim / (numberOfSamples + 1);

    // ================ Step 1: input interpolation ================ //
    // This is the point in the Input parameter space where we want to get an
    // estimate of the geometric dimensions.
    std::vector<NumericType> x = {params.taperAngle,
                                  params.stickingProbability};

    // Interpolate the samples based on the provided input. For this we use
    // nearest neighbors interpolation.
    int numberOfNeighbors = 3;
    NumericType distanceExponent = 2.;
    psNearestNeighborsInterpolation<NumericType,
                                    psMedianDistanceScaler<NumericType>>
        estimator;

    estimator.setNumberOfNeighbors(numberOfNeighbors);
    estimator.setDistanceExponent(distanceExponent);
    estimator.setDataDimensions(InputDim, TargetDim);
    estimator.setData(data);

    if (!estimator.initialize())
      return EXIT_FAILURE;

    auto estimateOpt = estimator.estimate(x);
    if (!estimateOpt) {
      std::cout << "Error during estimation.\n";
      return EXIT_FAILURE;
    }

    auto [result, distance] = estimateOpt.value();

    std::cout << std::setw(40) << "Distance to nearest data point: ";
    std::cout << distance << std::endl;

    // ================ Step 2: Time interpolation ================ //
    int stepSize = (numberOfSamples + 1);

    std::vector<NumericType> timesteps;
    timesteps.reserve(numberOfTimesteps);
    std::copy_if(result.begin(), result.end(), std::back_inserter(timesteps),
                 [=, i = 0](auto) mutable { return (i++ % stepSize) == 0; });

#ifndef LINEAR
    std::vector<std::vector<NumericType>> dimensionsOverTime;
    dimensionsOverTime.reserve(numberOfTimesteps);

    for (int i = 0; i < numberOfTimesteps; ++i) {
      std::vector<NumericType> vec(numberOfSamples);
      for (int j = 0; j < numberOfSamples; ++j)
        vec[j] = result[i * stepSize + 1 + j];

      dimensionsOverTime.push_back(vec);
    }

    NaturalCubicSplineInterpolation<NumericType> spline(timesteps,
                                                        dimensionsOverTime);

    auto interpolated = spline(params.processTime);

    std::vector<NumericType> dimensions(interpolated.begin(),
                                        interpolated.end());

#else
    // Since we use regular intervals we can simply take the difference between
    // the first two timesteps
    int extractionInterval = timesteps.at(1) - timesteps.at(0);

    // Now determine which two timesteps we should consider for interpolating
    // along the time axis
    NumericType extractionStep = params.processTime / extractionInterval;

    int lowerIdx = std::clamp(static_cast<int>(std::floor(extractionStep)) *
                                  (numberOfSamples + 1),
                              0, TargetDim - stepSize - 1);

    int upperIdx = std::clamp(static_cast<int>(std::ceil(extractionStep)) *
                                  (numberOfSamples + 1),
                              0, TargetDim - stepSize - 1);
    NumericType distanceToLower = 1.0 * lowerIdx - extractionStep;
    NumericType distanceToUpper = 1.0 * extractionStep - upperIdx;
    NumericType totalDistance = distanceToUpper + distanceToLower;

    std::vector<NumericType> dimensions;

    // Copy the data corresponding to the dimensions of the lower timestep
    auto lowerDimensions =
        psSmartPointer<std::vector<NumericType>>::New(numberOfSamples);
    for (int i = 0; i < numberOfSamples; ++i)
      lowerDimensions->at(i) = result.at(lowerIdx + 1 + i);

    if (totalDistance > 0) {
      // Copy the data corresponding to the dimensions of the upper timestep
      auto upperDimensions =
          psSmartPointer<std::vector<NumericType>>::New(numberOfSamples);
      for (int i = 0; i < numberOfSamples; ++i)
        upperDimensions->at(i) = result.at(upperIdx + 1 + i);

      // Now for each individual dimension do linear interpolation between upper
      // and lower value based on the relative distance
      for (unsigned i = 0; i < lowerDimensions->size(); ++i)
        dimensions.emplace_back((distanceToUpper * lowerDimensions->at(i) +
                                 distanceToLower * upperDimensions->at(i)) /
                                totalDistance);
    } else {
      std::copy(lowerDimensions->begin(), lowerDimensions->end(),
                std::back_inserter(dimensions));
    }
#endif
    NumericType origin[D] = {0.};
    origin[D - 1] = params.processTime + params.trenchHeight;

    auto substrate = createEmptyLevelset<NumericType, D>(params);
    interpolatedGeometry->insertNextLevelSet(substrate);
    auto geometry = psSmartPointer<lsDomain<NumericType, D>>::New(substrate);

    GeometryReconstruction<NumericType, D>(geometry, origin, sampleLocations,
                                           dimensions)
        .apply();

    interpolatedGeometry->insertNextLevelSet(geometry);
#ifndef LINEAR
    interpolatedGeometry->printSurface("interpolated_spline.vtp");
#else
    interpolatedGeometry->printSurface("interpolated_linear.vtp");
#endif
  }
  auto stop = Clock::now();
  auto interpDuration = Duration(stop - start).count();
  std::cout << std::setw(40) << "Interpolation and reconstruction took: ";
  std::cout << std::scientific << interpDuration << "s\n";

  // Actual simulation for reference
  start = Clock::now();
  auto referenceGeometry = psSmartPointer<psDomain<NumericType, D>>::New();
  {
    psMakeTrench<NumericType, D>(
        referenceGeometry, params.gridDelta /* grid delta */,
        params.xExtent /*x extent*/, params.yExtent /*y extent*/,
        params.trenchWidth /*trench width*/,
        params.trenchHeight /*trench height*/,
        params.taperAngle /* tapering angle */)
        .apply();

    referenceGeometry->printSurface("initial.vtp");

    params.processTime = params.processTime / params.stickingProbability;

    executeProcess<NumericType, D>(referenceGeometry, params);
    referenceGeometry->printSurface("reference.vtp");
  }
  stop = Clock::now();
  auto simDuration = Duration(stop - start).count();
  std::cout << std::setw(40) << "Physical simulation took: ";
  std::cout << std::scientific << simDuration << "s\n";

  std::cout << std::setw(40) << "Speedup: ";
  std::cout << std::fixed << simDuration / interpDuration << '\n';

  // Now compare the two geometries by using the chamfer distance as a measure
  // of similarity
  auto mesh1 = psSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(referenceGeometry->getLevelSets()->back(), mesh1)
      .apply();
  auto nodes1 = mesh1->getNodes();

  auto mesh2 = psSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(interpolatedGeometry->getLevelSets()->back(),
                               mesh2)
      .apply();
  auto nodes2 = mesh2->getNodes();

  auto chamferDistance =
      ChamferDistanceScore<NumericType>(nodes1, nodes2).calculate();

  std::cout << std::setw(40) << "Chamfer distance (0==perfect match): ";
  std::cout << std::fixed << chamferDistance << std::endl;
  return EXIT_SUCCESS;
}