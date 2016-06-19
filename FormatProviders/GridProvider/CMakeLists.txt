# CMake project file for avmodel lib
cmake_minimum_required(VERSION 2.6) 
set (PROJECT GridProvider) 
 
project (${PROJECT}) 

set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (LIBRARY_OUTPUT_PATH ${SRC_DIR}/../../_Debug/bin/)

set (SOURCES

${SRC_DIR}/ArrayVertexToGridFormatter.cpp
${SRC_DIR}/ArrayVertexToGridFormatter.h
${SRC_DIR}/BoundaryCondition.cpp
${SRC_DIR}/BoundaryCondition.h
${SRC_DIR}/BoundaryNormal.cpp
${SRC_DIR}/BoundaryNormal.h
${SRC_DIR}/ControlReader.cpp
${SRC_DIR}/ControlReader.h
${SRC_DIR}/MeshDataCleaner.cpp
${SRC_DIR}/MeshDataCleaner.h
${SRC_DIR}/MeshDataGraphFormatter.cpp
${SRC_DIR}/MeshDataGraphFormatter.h
${SRC_DIR}/MeshDataProvider.cpp
${SRC_DIR}/MeshDataProvider.h
${SRC_DIR}/MeshDataSolverFormatter.cpp
${SRC_DIR}/MeshDataSolverFormatter.h
${SRC_DIR}/MeshLayer.cpp
${SRC_DIR}/MeshLayer.h
${SRC_DIR}/MeshRow.cpp
${SRC_DIR}/MeshRow.h
${SRC_DIR}/NodeLinks.cpp
${SRC_DIR}/NodeLinks.h
${SRC_DIR}/OccRectilinearGrid.cpp
${SRC_DIR}/OccRectilinearGrid.h
${SRC_DIR}/Point3D.cpp
${SRC_DIR}/Point3D.h
${SRC_DIR}/RLCControlReader.cpp
${SRC_DIR}/RLCControlReader.h
${SRC_DIR}/RLCControlWriter.cpp
${SRC_DIR}/RLCControlWriter.h
${SRC_DIR}/RLCHeader.cpp
${SRC_DIR}/RLCHeader.h
${SRC_DIR}/SurfaceParam.cpp
${SRC_DIR}/SurfaceParam.h
${SRC_DIR}/UsingPolicy.cpp
${SRC_DIR}/UsingPolicy.h
${SRC_DIR}/Vertex.h
)

add_library (${PROJECT} STATIC ${SOURCES})
