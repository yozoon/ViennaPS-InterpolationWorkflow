#pragma once

#include <SimpleDeposition.hpp>
#include <psCSVWriter.hpp>
#include <psMakeTrench.hpp>
#include <psProcess.hpp>
#include <psToSurfaceMesh.hpp>
#include <psUtils.hpp>

#include "AdvectionCallback.hpp"
#include "Parameters.hpp"

template <typename NumericType, int D, typename AdvectionCallbackType>
void executeProcess(lsSmartPointer<psDomain<NumericType, D>> geometry,
                    const Parameters<NumericType> &params,
                    lsSmartPointer<AdvectionCallbackType> advectionCallback) {
  // copy top layer to capture deposition
  auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(
      geometry->getLevelSets()->back());
  geometry->insertNextLevelSet(depoLayer);

  auto processModel =
      SimpleDeposition<NumericType, D>(
          params.stickingProbability /*particle sticking probability*/,
          params.sourcePower /*particle source power*/)
          .getProcessModel();

  processModel->setAdvectionCallback(advectionCallback);

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(processModel);
  process.setNumberOfRaysPerPoint(2000);
  process.setProcessDuration(params.processTime / params.stickingProbability);

  process.apply();
}

template <typename NumericType, int D>
void executeProcess(lsSmartPointer<psDomain<NumericType, D>> geometry,
                    const Parameters<NumericType> &params) {
  // copy top layer to capture deposition
  auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(
      geometry->getLevelSets()->back());
  geometry->insertNextLevelSet(depoLayer);

  auto processModel =
      SimpleDeposition<NumericType, D>(
          params.stickingProbability /*particle sticking probability*/,
          params.sourcePower /*particel source power*/)
          .getProcessModel();

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(processModel);
  process.setNumberOfRaysPerPoint(2000);
  process.setProcessDuration(params.processTime / params.stickingProbability);

  process.apply();
}