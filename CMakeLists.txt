cmake_minimum_required(VERSION 3.20)

project("InterpolationWorkflow")

set(VIENNAPS_VERBOSE OFF)

# Check dependencies
find_package(ViennaPS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(fmt)

# Get the pre-c++20 span implementation library
include(FetchContent)
FetchContent_Declare(
  spanlite
  GIT_REPOSITORY https://github.com/martinmoene/span-lite.git
)
FetchContent_MakeAvailable(spanlite)

if (MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

list(APPEND PROJECT_INCLUDE_DIRS ${VIENNAPS_INCLUDE_DIRS} include/ ${spanlite_SOURCE_DIR}/include/nonstd)
list(APPEND PROJECT_LIBRARIES ${VIENNAPS_LIBRARIES} ${LAPACK_LIBRARIES})

set(WORKFLOW InterpolationWorkflow)
add_executable(${WORKFLOW} src/${WORKFLOW}.cpp)
target_include_directories(${WORKFLOW} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${WORKFLOW} PRIVATE ${PROJECT_LIBRARIES})

set(RECONSTRUCT GeometryReconstruction)
add_executable(${RECONSTRUCT} src/${RECONSTRUCT}.cpp)
target_include_directories(${RECONSTRUCT} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${RECONSTRUCT} PRIVATE ${PROJECT_LIBRARIES})

set(CREATE_DATA CreateData)
add_executable(${CREATE_DATA} src/${CREATE_DATA}.cpp)
target_include_directories(${CREATE_DATA} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries(${CREATE_DATA} PRIVATE ${PROJECT_LIBRARIES})

if(fmt_FOUND)
  list(APPEND PROJECT_LIBRARIES fmt::fmt)

  set(SPLINE TestSpline)
  add_executable(${SPLINE} src/${SPLINE}.cpp)
  target_include_directories(${SPLINE} PUBLIC ${PROJECT_INCLUDE_DIRS})
  target_link_libraries(${SPLINE} PRIVATE ${PROJECT_LIBRARIES})
endif()

configure_file(data/data.csv ${CMAKE_CURRENT_BINARY_DIR}/data.csv COPYONLY)
configure_file(data/data_pivot.csv ${CMAKE_CURRENT_BINARY_DIR}/data_pivot.csv COPYONLY)
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
