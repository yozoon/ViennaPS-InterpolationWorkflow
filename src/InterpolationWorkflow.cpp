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

#include "ChamferDistance.hpp"
#include "FeatureReconstruction.hpp"
#include "Parameters.hpp"
#include "SplineGridInterpolation.hpp"
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
    auto data = dataSource.getData();

    // Also read the data that was stored alongside the actual data describing
    // where it was sampled from (at which relative height)
    auto sampleLocations = dataSource.getPositionalParameters();

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
    NumericType extractionInterval;
    if (auto id = namedParams.find("ExtractionInterval");
        id != namedParams.end()) {
      extractionInterval = std::round(id->second);
    } else {
      std::cout
          << "'ExtractionInterval' not found in provided data CSV file.\n";
      return EXIT_FAILURE;
    }

    // The dimension of the data that is to be interpolated. In our case this
    // are the extracted dimensions at different timesteps (the timesteps are
    // also included in the data itself)
    int numFeatures = data->at(0).size() - InputDim;

    SplineGridInterpolation<NumericType> gridInterpolation;
    gridInterpolation.setDataDimensions(InputDim, numFeatures);
    gridInterpolation.setData(data);
    gridInterpolation.setBCType(SplineBoundaryConditionType::NATURAL);
    gridInterpolation.initialize();

    std::vector<NumericType> evaluationPoint = {
        params.taperAngle, params.stickingProbability,
        params.processTime / extractionInterval};

    for (auto e : evaluationPoint)
      std::cout << e << ", ";
    std::cout << "\n";

    auto estimationOpt = gridInterpolation.estimate(evaluationPoint);
    if (!estimationOpt)
      return EXIT_FAILURE;

    auto [estimatedFeatures, isInside] = estimationOpt.value();

    std::cout << "Inside: " << isInside << '\n';
    std::cout << estimatedFeatures.at(0) << "\n";

    NumericType origin[D] = {0.};
    origin[D - 1] = params.processTime + params.trenchHeight;

    auto substrate = createEmptyLevelset<NumericType, D>(params);
    interpolatedGeometry->insertNextLevelSet(substrate);
    auto geometry = psSmartPointer<lsDomain<NumericType, D>>::New(substrate);

    FeatureReconstruction<NumericType, D>(geometry, origin, sampleLocations,
                                          estimatedFeatures)
        .apply();

    interpolatedGeometry->insertNextLevelSet(geometry);

    interpolatedGeometry->printSurface("interpolated.vtp");
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