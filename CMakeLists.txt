cmake_minimum_required(VERSION 3.20)

project("InterpolationWorkflow")

set(VIENNAPS_VERBOSE OFF)

find_package(ViennaPS REQUIRED)
find_package(fmt REQUIRED) 
find_package(LAPACK REQUIRED)

# Get the pre-c++20 span implementation library
include(FetchContent)
FetchContent_Declare(
  spanlite
  GIT_REPOSITORY https://github.com/martinmoene/span-lite.git
)
FetchContent_MakeAvailable(spanlite)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

include_directories(include)
include_directories(${spanlite_SOURCE_DIR}/include/nonstd)

list(APPEND VIENNAPS_LIBRARIES fmt::fmt ${LAPACK_LIBRARIES})

set(WORKFLOW InterpolationWorkflow)
add_executable(${WORKFLOW} src/${WORKFLOW}.cpp)
target_include_directories(${WORKFLOW} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${WORKFLOW} PRIVATE ${VIENNAPS_LIBRARIES})

set(CREATE_DATA CreateData)
add_executable(${CREATE_DATA} src/${CREATE_DATA}.cpp)
target_include_directories(${CREATE_DATA} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${CREATE_DATA} PRIVATE ${VIENNAPS_LIBRARIES})

set(RECONSTRUCT GeometryReconstruction)
add_executable(${RECONSTRUCT} src/${RECONSTRUCT}.cpp)
target_include_directories(${RECONSTRUCT} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${RECONSTRUCT} PRIVATE ${VIENNAPS_LIBRARIES})


set(SPLINE TestSpline)
add_executable(${SPLINE} src/${SPLINE}.cpp)
target_include_directories(${SPLINE} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${SPLINE} PRIVATE ${VIENNAPS_LIBRARIES})

configure_file(data/data.csv ${CMAKE_CURRENT_BINARY_DIR}/data.csv COPYONLY)
configure_file(data/data_concat.csv ${CMAKE_CURRENT_BINARY_DIR}/data_concat.csv COPYONLY)
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
