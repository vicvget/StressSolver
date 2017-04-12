# CMake project file
cmake_minimum_required(VERSION 2.6) 
set (PROJECT stress) 
 
project (${PROJECT}) 

# утстановка пути до директории с исходными файлами
set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (LIBRARY_OUTPUT_PATH ${SRC_DIR}/../../_Debug/bin/)


set_property(TARGET MathHelpers PROPERTY
             IMPORTED_LOCATION ${SRC_DIR}/../../_Debug/bin/libMathHelpers.a)
set (SOURCES
${SRC_DIR}/BoundaryParams.cpp
${SRC_DIR}/RotationSolver.cpp
${SRC_DIR}/StressStrainCppIterativeSolver.cpp
${SRC_DIR}/StressStrainCppIterativeSolverAVX.cpp
${SRC_DIR}/StressStrainCppIterativeSolverFMA.cpp
${SRC_DIR}/StressStrainCppIterativeSolverKNC.cpp
${SRC_DIR}/StressStrainCppIterativeSolverKNC2.cpp
${SRC_DIR}/StressStrainCppSolver.cpp
${SRC_DIR}/StressStrainSolver.cpp
${SRC_DIR}/StressStrainSolverExports.cpp
)

add_library (${PROJECT} STATIC ${SOURCES})


