#include <iostream>

#include <psSmartPointer.hpp>

#include "AdvectionCallback.hpp"
#include "ChamferDistance.hpp"
#include "FeatureExtraction.hpp"
#include "FeatureReconstruction.hpp"
#include "Parameters.hpp"
#include "TrenchDeposition.hpp"

int main() {
  using NumericType = double;
  static constexpr int D = 2;

  static constexpr int numberOfSamples = 30;

  Parameters<NumericType> params;

  // Generate the initial trench geometry
  auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  psMakeTrench<NumericType, D>(
      geometry, params.gridDelta /* grid delta */, params.xExtent /*x extent*/,
      params.yExtent /*y extent*/, params.trenchWidth /*trench width*/,
      params.trenchHeight /*trench height*/,
      params.taperAngle /* tapering angle */)
      .apply();

  // Run a physical deposition simulation
  executeProcess<NumericType, D>(geometry, params);
  geometry->printSurface("GR_simulation.vtp");

  // Extract features from the geometry
  FeatureExtraction<NumericType, D> extraction;
  extraction.setDomain(geometry);
  extraction.setNumberOfSamples(numberOfSamples, false /* closed */);
  extraction.setEdgeAffinity(2.);
  extraction.apply();

  auto sampleLocations = extraction.getSampleLocations();
  auto features = extraction.getFeatures();

  assert(sampleLocations->size() == features->size());

  // Now reconstruct the geometry based on the extracted features
  NumericType origin[D] = {0.};
  origin[D - 1] = params.processTime + params.trenchHeight;

  auto ls = psSmartPointer<lsDomain<NumericType, D>>::New(
      geometry->getLevelSets()->back()->getGrid());

  FeatureReconstruction<NumericType, D>(ls, origin, *sampleLocations, *features)
      .apply();

  auto reconstruction = psSmartPointer<psDomain<NumericType, D>>::New();
  reconstruction->insertNextLevelSet(ls);
  reconstruction->printSurface("GR_reconstruction.vtp");

  // Now compare the two geometries by using the chamfer distance as a measure
  // of similarity
  auto mesh1 = psSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(ls, mesh1).apply();
  auto nodes1 = mesh1->getNodes();

  auto mesh2 = psSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(geometry->getLevelSets()->back(), mesh2).apply();
  auto nodes2 = mesh2->getNodes();

  auto chamferDistance =
      ChamferDistanceScore<NumericType>(nodes1, nodes2).calculate();

  std::cout << "Chamfer distance (0==perfect match): " << chamferDistance
            << std::endl;
}