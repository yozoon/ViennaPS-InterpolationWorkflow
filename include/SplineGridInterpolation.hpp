#ifndef NATURAL_CUBIC_SPLINE_GRID_INTERPOLATION_HPP
#define NATURAL_CUBIC_SPLINE_GRID_INTERPOLATION_HPP

#include <algorithm>
#include <array>
#include <iostream>
#include <optional>
#include <set>
#include <vector>

#include "CubicSplineInterpolation.hpp"

// Class providing cubic spline interpolation on rectilinear data grids
template <typename NumericType> class SplineGridInterpolation {
  using SizeType = size_t;
  using ItemType = std::vector<NumericType>;
  using VectorType = std::vector<ItemType>;
  using VectorPtr = psSmartPointer<std::vector<ItemType>>;
  using ConstPtr = psSmartPointer<const std::vector<ItemType>>;

public:
  // Future improvement: parallelize the recursive sorting using OpenMP taks
  static bool rearrange(typename VectorType::iterator start,
                        typename VectorType::iterator end, SizeType axis,
                        std::vector<std::set<NumericType>> &uniqueValues,
                        SizeType inputDim, bool capture) {
    bool equalSize = true;

    // We only reorder based on the input dimension, not the output dimension
    if (axis >= inputDim)
      return equalSize;

    SizeType size = end - start;
    if (size > 1) {
      // Now sort the data along the given axis
      std::sort(start, end, [&](const auto &a, const auto &b) {
        return a[axis] < b[axis];
      });

      // Record the indices that separate ranges of the same value
      SizeType rangeSize = 0;
      std::vector<SizeType> rangeBreaks;
      rangeBreaks.push_back(0);
      if (capture)
        uniqueValues[axis].insert(start->at(axis));

      for (unsigned i = 1; i < size; ++i)
        if ((start + i - 1)->at(axis) != (start + i)->at(axis)) {
          if (rangeSize == 0)
            rangeSize = i;

          SizeType tmp = rangeBreaks.back();

          rangeBreaks.push_back(i);
          uniqueValues[axis].insert((start + i)->at(axis));

          if (rangeSize != i - tmp) {
            std::cout << "Data is not arranged in a rectilinear grid!\n";
            equalSize = false;
            return false;
          }
        }

      rangeBreaks.push_back(size);

      // Recursively launch sorting tasks for each of the separated ranges. Only
      // the first leaf in each level of the tree is instructed to capture the
      // unique values.
      for (unsigned i = 1; i < rangeBreaks.size(); ++i)
        equalSize =
            equalSize &&
            rearrange(start + rangeBreaks[i - 1], start + rangeBreaks[i],
                      axis + 1, uniqueValues, inputDim, capture && (i == 1));
    }

    return equalSize;
  }

private:
  std::vector<SplineBoundaryConditionType> bcTypes;

  SizeType inputDim{0};
  SizeType outputDim{0};

  ConstPtr data = nullptr;

  bool dataChanged = true;

  // The unique values along each dimension
  std::vector<std::set<NumericType>> uniqueValues;

  // A local copy of the provided data
  VectorType localData;

public:
  SplineGridInterpolation() {}

  void
  setBCTypes(const std::vector<SplineBoundaryConditionType> &passedBCTypes) {
    assert(passedBCTypes.size() == inputDim &&
           "Not enough boundary conditions passed.");
    bcTypes = passedBCTypes;
  }

  std::vector<std::set<NumericType>> getUniqueValues() const {
    return uniqueValues;
  }

  ConstPtr getSortedData() const { return ConstPtr::New(localData); }

  void setBCType(SplineBoundaryConditionType passedBCType) {
    bcTypes = std::vector<SplineBoundaryConditionType>(inputDim, passedBCType);
  }

  void setDataDimensions(SizeType passedInputDim, SizeType passedOutputDim) {
    inputDim = passedInputDim;
    outputDim = passedOutputDim;
    if (bcTypes.empty()) {
      bcTypes.resize(inputDim, SplineBoundaryConditionType::NOT_A_KNOT);
    }
  }

  void setData(ConstPtr passedData) {
    data = passedData;
    dataChanged = true;
  }

  bool initialize() {
    if (!data || (data && data->empty())) {
      std::cout << "SplineGridInterpolation: the provided data is "
                   "empty.\n";
      return false;
    }

    if (data->at(0).size() != inputDim + outputDim) {
      std::cout << "SplineGridInterpolation: the sum of the provided "
                   "InputDimension and OutputDimension does not match the "
                   "dimension of the provided data.\n";
      return false;
    }

    if (bcTypes.size() != inputDim) {
      std::cout
          << "SplineGridInterpolation: the provided boundary condition vector "
             "has a different size than the input dimensions require.\n";
      return false;
    }

    localData.clear();
    localData.reserve(data->size());
    std::copy(data->begin(), data->end(), std::back_inserter(localData));

    uniqueValues.resize(inputDim);

    auto equalSize = rearrange(localData.begin(), localData.end(), 0,
                               uniqueValues, inputDim, true);

    if (!equalSize) {
      std::cout << "SplineGridInterpolation: Data is not arranged "
                   "in a rectilinear grid!\n";
      return false;
    }

    for (SizeType i = 0; i < inputDim; ++i)
      if (uniqueValues[i].empty()) {
        std::cout << "SplineGridInterpolation: The grid has no "
                     "values along dimension "
                  << i << std::endl;
        return false;
      }

    dataChanged = false;
    return true;
  }

  std::optional<std::tuple<std::vector<NumericType>, bool>>
  estimate(const std::vector<NumericType> &input) {
    if (dataChanged)
      if (!initialize())
        return {};

    bool isInside = true;
    for (SizeType i = 0; i < inputDim; ++i) {
      if (!uniqueValues[i].empty()) {
        // Check if the input lies within the bounds of our data grid
        if (input[i] < *(uniqueValues[i].begin()) ||
            input[i] > *(uniqueValues[i].rbegin())) {
          isInside = false;
          break;
        }
      } else {
        return {};
      }
    }

    // Interpolate in each dimension
    SizeType numPoints = localData.size();

    // Copy only the output dimensions of the data into a new temporary data
    // vector
    std::vector<std::vector<NumericType>> tmpData;
    tmpData.reserve(numPoints);
    for (auto &ld : localData) {
      tmpData.emplace_back(std::next(ld.begin(), inputDim), ld.end());
    }

    SizeType numberOfSplines = numPoints;
    SizeType numUniqueAlongPreviousAxis = 1;
    for (int i = inputDim - 1; i >= 0; --i) {
      // The knots used by the spline interpolation are just the unique values
      // along the current axis
      std::vector<NumericType> x(uniqueValues[i].begin(),
                                 uniqueValues[i].end());
      SizeType numUniqueAlongAxis = x.size();
      numberOfSplines /= numUniqueAlongAxis;
      for (SizeType j = 0; j < numberOfSplines; ++j) {
        // Copy the appropriate data points from the temporary data vector
        std::vector<std::vector<NumericType>> y;
        y.reserve(numUniqueAlongAxis);
        for (SizeType k = 0; k < numUniqueAlongAxis; ++k)
          y.push_back(tmpData.at(j * numUniqueAlongAxis +
                                 k * numUniqueAlongPreviousAxis));

        // Instantiate the spline interpolation
        CubicSplineInterpolation<NumericType> interpolation(x, y, bcTypes[i]);
        // And evaluate the interpolation function at the location of the input.
        // The output overwrites part of the temporary data vector, so that we
        // can use it as input in the next interpolation iteration.
        tmpData[j] = interpolation(input[i]);
      }
      numUniqueAlongAxis = numUniqueAlongAxis;
    }

    return {{tmpData[0], isInside}};
  }
};

#endif