cmake_minimum_required(VERSION 3.20)

project("InterpolationWorkflow")

set(VIENNAPS_VERBOSE OFF)

find_package(ViennaPS REQUIRED)
find_package(fmt REQUIRED) 
find_package(LAPACK REQUIRED)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

include_directories(include)

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
configure_file(data/data_consol.csv ${CMAKE_CURRENT_BINARY_DIR}/data_consol.csv COPYONLY)
configure_file(data/config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
