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
    if (features.size() != sampleLocations.size()) {
      std::cout << "Not enough features provided for reconstruction!\n";
      return;
    }

    // First generate an initial plane from which we will remove the trench
    // geometry later on
    {
      NumericType normal[D] = {0.};
      normal[verticalDir] = 1.;

      auto plane = psSmartPointer<lsPlane<NumericType, D>>::New(origin, normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    // INFO: Mesh point numbering CW: solid is enclosed inside points

    // Now create a mesh that reconstructs the trench profile by using the
    // extracted features. Use this mesh to generate a new levelset, which
    // will be subtracted from the plane.
    {
      NumericType depth = features[0];

      NumericType gridDelta = levelset->getGrid().getGridDelta();

      auto numSamplesRight = static_cast<unsigned>(
          std::ceil(1.0 * (sampleLocations.size() - 1) / 2));

      // std::cout << numSamplesRight << std::endl;

#ifndef NDEBUG
      int meshCount = 0;
#endif

      unsigned j = 0;
      while (j < numSamplesRight) {
        // Manually create a surface mesh based on the extracted features
        auto mesh = psSmartPointer<lsMesh<>>::New();

        if (j == 0) {
          // Add one point on top of the geometry, so that we avoid potential
          // sharp corners
          std::array<NumericType, 3> point{0.};
          point[0] = origin[0];
          point[1] = origin[1];
          if constexpr (D == 3)
            point[2] = origin[2];

          point[horizontalDir] += features[1];
          point[verticalDir] += 2 * gridDelta;

          mesh->insertNextNode(point);
        }

        unsigned nextJ = j;
        for (unsigned i = 1 + j; i < numSamplesRight + 1; ++i) {
          nextJ = i - 1;

          if (std::abs(features[i] - features[i + numSamplesRight]) <
              gridDelta / 2) {
#ifndef NDEBUG
            std::cout << "Pinchoff point detected!\n";
#endif
            break;
          }

          std::array<NumericType, 3> point{0.};
          point[0] = origin[0];
          point[1] = origin[1];
          if constexpr (D == 3)
            point[2] = origin[2];

          point[horizontalDir] += features[i];
          point[verticalDir] += depth * (sampleLocations[i] - 1.0);

          mesh->insertNextNode(point);
        }

        for (unsigned i = numSamplesRight + nextJ + 1; i > numSamplesRight + j;
             --i) {
          std::array<NumericType, 3> point{0.};
          point[0] = origin[0];
          point[1] = origin[1];
          if constexpr (D == 3)
            point[2] = origin[2];

          point[horizontalDir] += features[i];
          point[verticalDir] += depth * (sampleLocations[i] - 1.0);

          mesh->insertNextNode(point);
        }

        if (j == 0) {
          // Add one point on top of the geometry, so that we avoid potential
          // sharp corners
          std::array<NumericType, 3> point{0.};
          point[0] = origin[0];
          point[1] = origin[1];
          if constexpr (D == 3)
            point[2] = origin[2];

          point[horizontalDir] += features[numSamplesRight + 1];
          point[verticalDir] += 2 * gridDelta;

          mesh->insertNextNode(point);
        }

        j = nextJ + 1;

        // If there are not enough points for a triangle skip the boolean
        // operation.
        if (mesh->nodes.size() < 3) {
#ifndef NDEBUG
          std::cout << "Mesh is just a line. Skipping.\n";
#endif
          continue;
        }

        // Connect all nodes of the mesh with lines
        for (unsigned i = 0; i < mesh->nodes.size() - 1; ++i)
          mesh->lines.emplace_back(std::array<unsigned, 2>{i, i + 1});

        // Close the hull mesh
        mesh->lines.emplace_back(std::array<unsigned, 2>{
            static_cast<unsigned>(mesh->lines.size()), 0U});

#ifndef NDEBUG
        // Print the created hull mesh
        lsVTKWriter<NumericType>(mesh, "hullMesh" +
                                           std::to_string(meshCount++) + ".vtp")
            .apply();
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
  }
};
#endif