cmake_minimum_required(VERSION 3.20)

project("InterpolationWorkflow")

set(VIENNAPS_VERBOSE OFF)

find_package(ViennaPS REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED) 

if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # lots of warnings
    add_compile_options(-Wall -Wextra -pedantic)
endif()

include_directories(include)

list(APPEND VIENNAPS_LIBRARIES Eigen3::Eigen fmt::fmt)

set(WORKFLOW InterpolationWorkflow)
add_executable(${WORKFLOW} src/${WORKFLOW}.cpp)
target_include_directories(${WORKFLOW} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
target_link_libraries(${WORKFLOW} PRIVATE ${VIENNAPS_LIBRARIES})

# set(WORKFLOW_LIN InterpolationWorkflowLinear)
# add_executable(${WORKFLOW_LIN} src/${WORKFLOW}.cpp)
# target_include_directories(${WORKFLOW_LIN} PUBLIC ${VIENNAPS_INCLUDE_DIRS})
# target_link_libraries(${WORKFLOW_LIN} PRIVATE ${VIENNAPS_LIBRARIES})
# target_compile_definitions(${WORKFLOW_LIN} PUBLIC LINEAR)

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
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
