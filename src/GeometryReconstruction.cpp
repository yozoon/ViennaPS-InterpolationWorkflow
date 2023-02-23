#include <iostream>

#include <lsSmartPointer.hpp>

#include "AdvectionCallback.hpp"
#include "ChamferDistance.hpp"
#include "FeatureExtraction.hpp"
#include "FeatureReconstruction.hpp"
#include "MakeTrench.hpp"
#include "Parameters.hpp"
#include "TrenchDeposition.hpp"

int main() {
  using NumericType = double;
  static constexpr int D = 2;

  static constexpr int numberOfSamples = 101;

  Parameters<NumericType> params;
  params.taperAngle = -5.;
  params.processTime = 4.;
  params.stickingProbability = .4;

  // Generate the initial trench geometry
  auto geometry = lsSmartPointer<psDomain<NumericType, D>>::New();
  MakeTrench<NumericType, D>(geometry, params.gridDelta, params.xExtent,
                             params.yExtent, params.trenchWidth,
                             params.trenchHeight, params.taperAngle, 0.)
      .apply();
  geometry->printSurface("GR_initial.vtp");

  // Run a physical deposition simulation
  executeProcess<NumericType, D>(geometry, params);
  geometry->printSurface("GR_simulation.vtp");

  // Extract features from the geometry
  FeatureExtraction<NumericType, D> extraction;
  extraction.setDomain(geometry->getLevelSets()->back());
  extraction.setNumberOfSamples(numberOfSamples, false /* open */);
  extraction.setEdgeAffinity(0.0);
  extraction.setOrigin(std::array<NumericType, 3>{0., params.trenchHeight, 0.});
  extraction.apply();

  auto sampleLocations = extraction.getSampleLocations();
  auto features = extraction.getFeatures();

#ifndef NDEBUG
  std::cout << "Number of features=" << features->size() << std::endl;
  for (unsigned i = 0; i < sampleLocations->size(); ++i) {
    std::cout << i << ": " << std::setprecision(4) << sampleLocations->at(i)
              << ", " << features->at(i) << '\n';
  }
#endif

  assert(sampleLocations->size() == features->size());

  // Now reconstruct the geometry based on the extracted features
  NumericType origin[D] = {0.};
  origin[D - 1] = params.processTime + params.trenchHeight;

  auto ls = lsSmartPointer<lsDomain<NumericType, D>>::New(
      geometry->getLevelSets()->back()->getGrid());

  FeatureReconstruction<NumericType, D>(ls, origin, *sampleLocations, *features)
      .apply();

  auto reconstruction = lsSmartPointer<psDomain<NumericType, D>>::New();
  reconstruction->insertNextLevelSet(ls);
  reconstruction->printSurface("GR_reconstruction.vtp");

  // Now compare the two geometries by using the chamfer distance as a measure
  // of similarity
  auto mesh1 = lsSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(ls, mesh1).apply();
  auto nodes1 = mesh1->getNodes();

  auto mesh2 = lsSmartPointer<lsMesh<>>::New();
  lsToDiskMesh<NumericType, D>(geometry->getLevelSets()->back(), mesh2).apply();
  auto nodes2 = mesh2->getNodes();

  auto chamferDistance =
      ChamferDistanceScore<NumericType>(nodes1, nodes2).calculate();

  std::cout << "Chamfer distance (0==perfect match): " << chamferDistance
            << std::endl;
}