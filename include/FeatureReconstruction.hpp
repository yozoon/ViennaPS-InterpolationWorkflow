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

  int verticalDir = D - 1;
  int horizontalDir = 0;

  psSmartPointer<lsMesh<>>
  pointsToMesh(const std::vector<NumericType> &x,
               const std::vector<NumericType> &y) const {}

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
      normal[verticalDir] = 1.;

      auto plane = psSmartPointer<lsPlane<NumericType, D>>::New(origin, normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    if (features.size() != sampleLocations.size()) {
      std::cout << "Not enough features provided for reconstruction!\n";
      return;
    }

    // INFO: Mesh point numbering CW: solid is enclosed inside points

    // Now create a mesh that reconstructs the trench profile by using the
    // extracted features. Use this mesh to generate a new levelset, which
    // will be subtracted from the plane.
    {
      NumericType depth = features[0];

      // Manually create a surface mesh based on the extracted features
      auto mesh = psSmartPointer<lsMesh<>>::New();

      NumericType gridDelta = levelset->getGrid().getGridDelta();

      unsigned samplesPerSide = (sampleLocations.size() - 1) / 2;
      std::cout << samplesPerSide << std::endl;

      if (std::abs(features[1]) > eps) {
        // Add one point on top of the geometry, so that we avoid potential
        // sharp corners
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[horizontalDir] += features[1]; // std::max(features[2], eps);
        point[verticalDir] += 2 * gridDelta;

        mesh->insertNextNode(point);
      }

      for (unsigned i = 1; i < samplesPerSide; ++i) {
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[horizontalDir] += features[i]; // std::max(features[i], eps);
        point[verticalDir] -= depth * sampleLocations[i];

        mesh->insertNextNode(point);

        if (features[i] == 0.0 && features[i + samplesPerSide] == 0.0) {
          std::cout << "Pinchoff point detected\n";
        }
      }

      for (unsigned i = samplesPerSide; i < sampleLocations.size(); ++i) {
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[horizontalDir] -= features[i]; // std::max(features[i], eps);
        point[verticalDir] -= depth * sampleLocations[i];

        mesh->insertNextNode(point);
      }

      if (std::abs(features.back()) > eps) {
        // Add one point on top of the geometry, so that we avoid potential
        // sharp corners
        std::array<NumericType, 3> point{0.};
        point[0] = origin[0];
        point[1] = origin[1];
        if constexpr (D == 3)
          point[2] = origin[2];

        point[horizontalDir] -=
            features.back(); // std::max(features.back(), eps);
        point[verticalDir] += 2 * gridDelta;

        mesh->insertNextNode(point);
      }

      for (unsigned i = 0; i < mesh->nodes.size() - 1; ++i)
        mesh->lines.emplace_back(std::array<unsigned, 2>{i, i + 1});

      mesh->lines.emplace_back(std::array<unsigned, 2>{
          static_cast<unsigned>(mesh->lines.size()), 0U});

#ifndef NDEBUG
      // Print the created hull mesh
      lsVTKWriter<NumericType>(mesh, "hullMesh.vtp").apply();
#endif

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