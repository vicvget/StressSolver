# CMake project file
cmake_minimum_required(VERSION 2.6) 
set (PROJECT test) 
 
project (${PROJECT}) 

# утстановка пути до директории с исходными файлами
set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (EXECUTABLE_OUTPUT_PATH ${SRC_DIR}/../_Debug/bin/)

set (SOURCES
${SRC_DIR}/BlenderExporter.cpp
${SRC_DIR}/ChartsExporter.cpp
${SRC_DIR}/CommonSolvers.cpp
${SRC_DIR}/CommonSolversTest.cpp
${SRC_DIR}/ComparativeTest.cpp
${SRC_DIR}/DummyExporter.cpp
${SRC_DIR}/FrameChartsExporter.cpp
${SRC_DIR}/main.cpp
${SRC_DIR}/MprExporter.cpp
${SRC_DIR}/StressSolver.cpp
${SRC_DIR}/StressSolverTest.cpp
${SRC_DIR}/Test10x1x1.cpp
${SRC_DIR}/Test10x1x3.cpp
${SRC_DIR}/Test10x3x1.cpp
${SRC_DIR}/Test10x3x3.cpp
${SRC_DIR}/Test10x5x5.cpp
${SRC_DIR}/Test10x5x5_2.cpp
${SRC_DIR}/Test10x7x7.cpp
${SRC_DIR}/Test1x100x100a.cpp
${SRC_DIR}/Test1x10x10.cpp
${SRC_DIR}/Test1x11x11.cpp
${SRC_DIR}/Test1x1x3.cpp
${SRC_DIR}/Test1x22x22a.cpp
${SRC_DIR}/Test1x44x44a.cpp
${SRC_DIR}/Test1x51x51.cpp
${SRC_DIR}/Test1x66x66a.cpp
${SRC_DIR}/Test50x5x5.cpp
${SRC_DIR}/Test50x5x5_2.cpp
${SRC_DIR}/TestFactory.cpp
${SRC_DIR}/TestPlate1xXxXa.cpp
)

link_directories (${SRC_DIR}/../_Debug/bin/)
add_executable(${PROJECT} ${SOURCES})
target_link_libraries (${PROJECT} MathHelpers  GridProvider ResProvider stress)
