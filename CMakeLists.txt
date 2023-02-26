cmake_minimum_required(VERSION 3.20)

project("InterpolationWorkflow")

set(VIENNAPS_VERBOSE OFF)

# Check dependencies
find_package(ViennaLS REQUIRED)
find_package(LAPACK REQUIRED)

# Get the pre-c++20 span implementation library
include(FetchContent)
FetchContent_Declare(
  spanlite
  GIT_REPOSITORY https://github.com/martinmoene/span-lite.git
)
FetchContent_MakeAvailable(spanlite)

# Set default build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

# Enable more warning levels
if (MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Add include directories
list(APPEND PROJECT_INCLUDE_DIRS ${VIENNALS_INCLUDE_DIRS} include/ ${spanlite_SOURCE_DIR}/include/nonstd)
# Add libraries
list(APPEND PROJECT_LIBRARIES ${VIENNALS_LIBRARIES} ${LAPACK_LIBRARIES})

# The main executable
set(WORKFLOW InterpolationWorkflow)
add_executable(${WORKFLOW} src/${WORKFLOW}.cpp)
target_include_directories(${WORKFLOW} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${WORKFLOW} PRIVATE ${PROJECT_LIBRARIES})

# Optional executables for generating the dataset as well as demonstrating the extraction and reconstruction flow
find_package(ViennaPS)
if(ViennaPS_FOUND)
  list(APPEND PROJECT_INCLUDE_DIRS ${VIENNAPS_INCLUDE_DIRS})
  list(APPEND PROJECT_LIBRARIES ${VIENNAPS_LIBRARIES})

  target_compile_definitions(${WORKFLOW} PRIVATE USE_VIENNAPS)

  set(RECONSTRUCT GeometryReconstruction)
  add_executable(${RECONSTRUCT} src/${RECONSTRUCT}.cpp)
  target_include_directories(${RECONSTRUCT} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${RECONSTRUCT} PRIVATE ${PROJECT_LIBRARIES})

  set(CREATE_DATA CreateData)
  add_executable(${CREATE_DATA} src/${CREATE_DATA}.cpp)
  target_include_directories(${CREATE_DATA} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${CREATE_DATA} PRIVATE ${PROJECT_LIBRARIES})
endif()

# Copy data and config files to the build directory
configure_file(data/data.csv ${CMAKE_CURRENT_BINARY_DIR}/data.csv COPYONLY)
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
