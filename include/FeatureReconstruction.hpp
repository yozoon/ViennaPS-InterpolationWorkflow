#ifndef FEATURE_RECONSTRUCTION_HPP
#define FEATURE_RECONSTRUCTION_HPP

#include <lsBooleanOperation.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsSmartPointer.hpp>

#include "span.hpp"

template <typename NumericType, int D> class FeatureReconstruction {

  lsSmartPointer<lsDomain<NumericType, D>> levelset;
  const NumericType (&origin)[D];
  const std::vector<NumericType> &sampleLocations;
  const std::vector<NumericType> &features;
  NumericType eps;

  int verticalDir = D - 1;
  int horizontalDir = 0;

public:
  FeatureReconstruction(lsSmartPointer<lsDomain<NumericType, D>> passedLevelset,
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

      auto plane = lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal);
      lsMakeGeometry<NumericType, D>(levelset, plane).apply();
    }

    // INFO: Mesh point numbering CW: solid is enclosed inside points

    // Now create a mesh that reconstructs the trench profile by using the
    // extracted features. Use this mesh to generate a new levelset, which
    // will be subtracted from the plane.
    {
      NumericType gridDelta = levelset->getGrid().getGridDelta();

      unsigned numSamplesRight = static_cast<unsigned>(
          std::ceil(1.0 * (sampleLocations.size() - 1) / 2));
      unsigned numSamplesLeft = sampleLocations.size() - 1 - numSamplesRight;

      NumericType depth = features.front();

      // Create spans for specific ranges in the feature and location vectors
      // (to simplify indexing and make the algorithm more readable)
      const auto rightFeatures =
          nonstd::span(std::next(features.begin(), 1),
                       std::next(features.begin(), numSamplesRight + 1));
      const auto rightLocations =
          nonstd::span(std::next(sampleLocations.begin(), 1),
                       std::next(sampleLocations.begin(), numSamplesRight + 1));

      const auto leftFeatures = nonstd::span(
          std::next(features.begin(), numSamplesLeft + 1), features.end());
      const auto leftLocations =
          nonstd::span(std::next(sampleLocations.begin(), numSamplesLeft + 1),
                       sampleLocations.end());

#ifndef NDEBUG
      int meshCount = 0;
#endif

      unsigned j = 0;
      while (j < std::min(numSamplesRight, numSamplesLeft)) {
        // Manually create a surface mesh based on the extracted features
        auto mesh = lsSmartPointer<lsMesh<>>::New();

        if (j == 0) {
          // Add one point on top of the geometry, in order to avoid potential
          // sharp corners
          std::array<NumericType, 3> point{0.};
          std::copy(std::begin(origin), std::end(origin), point.begin());

          point[horizontalDir] += rightFeatures.front();
          point[verticalDir] += 2 * gridDelta;

          mesh->insertNextNode(point);
        }

        unsigned nextJ = j;

        // Add points of the right sidewall to the mesh until we encounter a
        // pinch-off (the distance between the left and right sidwall features
        // at that location falling below a certain threshold)
        for (unsigned i = nextJ; i < rightFeatures.size(); ++i) {
          nextJ = i;

          if (rightFeatures[i] - leftFeatures[std::min(numSamplesLeft - 1, i)] <
              gridDelta / 5) {
#ifndef NDEBUG
            std::cout << "Pinchoff point detected!\n";
#endif
            break;
          }

          std::array<NumericType, 3> point{0.};
          std::copy(std::begin(origin), std::end(origin), point.begin());

          point[horizontalDir] += rightFeatures[i];
          point[verticalDir] += depth * (rightLocations[i] - 1.0);

          mesh->insertNextNode(point);
        }

        // Add the points of the left sidewall in reverse order to the mesh,
        // starting from the last, or the one at the previously determined
        // pinch-off location.
        for (int i = std::min(numSamplesLeft - 1, nextJ);
             i > static_cast<int>(j); --i) {
          std::array<NumericType, 3> point{0.};
          std::copy(std::begin(origin), std::end(origin), point.begin());

          point[horizontalDir] += leftFeatures[i];
          point[verticalDir] += depth * (leftLocations[i] - 1.0);

          mesh->insertNextNode(point);
        }

        if (j == 0) {
          // Add one point on top of the geometry, in order to avoid potential
          // sharp corners
          std::array<NumericType, 3> point{0.};
          std::copy(std::begin(origin), std::end(origin), point.begin());

          point[horizontalDir] += leftFeatures.back();
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
            lsSmartPointer<lsDomain<NumericType, D>>::New(levelset->getGrid());

        lsFromSurfaceMesh<NumericType, D>(hull, mesh).apply();

        lsBooleanOperation<NumericType, D>(
            levelset, hull, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
            .apply();
      }
    }
    //
  }
};
#endif