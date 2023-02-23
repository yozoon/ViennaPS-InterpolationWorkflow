#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>

#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <psCSVDataSource.hpp>

#include "GeometricTrenchModel.hpp"
#include "MakeTrench.hpp"
#include "Parameters.hpp"

namespace fs = std::filesystem;

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
    params.taperAngle = std::stof(argv[2]);
  }
  if (argc > 3) {
    NumericType tmp = std::stof(argv[3]);
    if (tmp > 0 && tmp <= 1.0)
      params.stickingProbability = tmp;
  }
  if (argc > 4) {
    NumericType tmp = std::stof(argv[4]);
    if (tmp >= 0.0)
      params.processTime = tmp;
  }

  // Interpolation based on previous simulation results
  auto start = Clock::now();

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
    std::cout << "'ExtractionInterval' not found in provided data CSV file.\n";
    return EXIT_FAILURE;
  }

  // The dimension of the data that is to be interpolated. In our case this
  // are the extracted dimensions at different timesteps (the timesteps are
  // also included in the data itself)
  int numFeatures = data->at(0).size() - InputDim;

  if (numFeatures != static_cast<int>(sampleLocations.size())) {
    lsMessage::getInstance().addError("Invalid feature dimension.").print();
    return EXIT_FAILURE;
  }

  std::array<NumericType, 3> origin{0.};
  origin[D - 1] = params.processTime + params.trenchHeight;

  auto substrate = MakeTrench<NumericType, D>(
      params.gridDelta, params.xExtent, params.yExtent, params.trenchWidth,
      params.trenchHeight, params.taperAngle, params.taperAngle);

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "IW_initial.vtp").apply();
  }

  GeometricTrenchDepositionModel<NumericType, D> geometricModel(
      substrate, sampleLocations, origin, extractionInterval, true);
  geometricModel.setData(data, InputDim);

  geometricModel.apply(params.taperAngle, params.stickingProbability,
                       params.processTime);

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "IW_interpolated_spline.vtp").apply();
  }

  auto stop = Clock::now();
  std::cout << std::setw(40) << "Interpolation and reconstruction took: ";
  std::cout << std::scientific << Duration(stop - start).count() << "s\n";
  return EXIT_SUCCESS;
}