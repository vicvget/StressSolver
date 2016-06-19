# CMake project file for avmodel lib
cmake_minimum_required(VERSION 2.6) 
set (PROJECT test) 
 
project (${PROJECT}) 

# утстановка пути до директории с исходными файлами
set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (EXECUTABLE_OUTPUT_PATH ${SRC_DIR}/../_Debug/bin/)

set (SOURCES
${SRC_DIR}/BlenderExporter.cpp
${SRC_DIR}/CommonSolvers.cpp
${SRC_DIR}/CommonSolversTest.cpp
${SRC_DIR}/ComparativeTest.cpp
${SRC_DIR}/main.cpp
${SRC_DIR}/StressSolver.cpp
${SRC_DIR}/StressSolverTest.cpp
)

link_directories (${SRC_DIR}/../_Debug/bin/)
add_executable(${PROJECT} ${SOURCES})
target_link_libraries (${PROJECT} MathHelpers  GridProvider ResProvider stress)
