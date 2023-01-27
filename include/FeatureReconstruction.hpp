#ifndef FEATURE_RECONSTRUCTION_HPP
#define FEATURE_RECONSTRUCTION_HPP

#include <lsBooleanOperation.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsVTKWriter.hpp>

template <typename NumericType, int D> class FeatureReconstruction {

  psSmartPointer<lsDomain<NumericType, D>> &levelset;
  const NumericType (&origin)[D];
  const std::vector<NumericType> &sampleLocations;
  const std::vector<NumericType> &features;
  NumericType eps;

public:
  FeatureReconstruction(
      psSmartPointer<lsDomain<NumericType, D>> &passedLevelset,
      const NumericType (&passedOrigin)[D],
      const std::vector<NumericType> &passedSampleLocations,
      const std::vector<NumericType> &passedFeatures)
      : levelset(passedLevelset), origin(passedOrigin),
        sampleLocations(passedSampleLocations), features(passedFeatures),
        eps(1e-4) {}

  void apply() {
    // First generate an initial plane from which we will remove the trench
    // geometry later on
    {
      NumericType normal[D] = {0.};
      normal[D - 1] = 1.;

      auto plane = psSmartPointer<lsPlane<NumericType, D>>::New(origin, normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    // Now create a mesh that reconstructs the trench profile by using the
    // extracted features. Use this mesh to generate a new levelset, which
    // will be subtracted from the plane.
    {
      NumericType depth = features.at(0);

      // Manually create a surface mesh based on the extracted features
      auto mesh = psSmartPointer<lsMesh<>>::New();

      NumericType gridDelta = levelset->getGrid().getGridDelta();

      for (unsigned i = 1; i < features.size(); ++i) {
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[0] -= std::max(features.at(i) / 2, eps);
        point[D - 1] -= depth * sampleLocations.at(i - 1);

        mesh->insertNextNode(point);
      }

      {
        // Add one point on top of the geometry, so that we avoid potential
        // sharp corners
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[0] -= std::max(features.back() / 2, eps);
        point[D - 1] += 2 * gridDelta;

        mesh->insertNextNode(point);

        // Now also do the same thing for the right side
        point[0] = origin[0];
        point[0] += std::max(features.back() / 2, eps);

        mesh->insertNextNode(point);
      }

      for (unsigned i = features.size() - 1; i >= 1; --i) {
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[0] += std::max(features.at(i) / 2, eps);
        point[D - 1] -= depth * sampleLocations.at(i - 1);

        mesh->insertNextNode(point);
      }

      for (unsigned i = 0; i < mesh->nodes.size() - 1; ++i)
        mesh->lines.emplace_back(std::array<unsigned, 2>{i, i + 1});

      mesh->lines.emplace_back(std::array<unsigned, 2>{
          static_cast<unsigned>(mesh->lines.size()), 0U});

      lsVTKWriter<NumericType>(mesh, "hullMesh.vtp").apply();

      // Create the new levelset based on the mesh and substract it from the
      // plane
      auto hull =
          psSmartPointer<lsDomain<NumericType, D>>::New(levelset->getGrid());

      lsFromSurfaceMesh<NumericType, D>(hull, mesh).apply();

      lsBooleanOperation<NumericType, D>(
          levelset, hull, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();
    }
  }
};
#endif