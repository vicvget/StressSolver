# CMake project file for avmodel lib
cmake_minimum_required(VERSION 2.6) 
set (PROJECT ResProvider) 
 
project (${PROJECT}) 

set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (LIBRARY_OUTPUT_PATH ${SRC_DIR}/../../_Debug/bin/)

set (SOURCES
#${SRC_DIR}/AbstractProviderRezr.cpp
#${SRC_DIR}/AbstractProviderRezr.h
${SRC_DIR}/Buffer.cpp
${SRC_DIR}/Buffer.h
#${SRC_DIR}/FileStream.cpp
#${SRC_DIR}/FileStream.h
#${SRC_DIR}/MechanicalObjects.h
#${SRC_DIR}/MechanicalObjects.cpp
#${SRC_DIR}/OldAddresses.cpp
#${SRC_DIR}/OldAddresses.h
${SRC_DIR}/ProviderMpr.cpp
${SRC_DIR}/ProviderMpr.h
#${SRC_DIR}/ProviderRezr.cpp
#${SRC_DIR}/ProviderRezr.h
${SRC_DIR}/ResultsFactory.cpp
${SRC_DIR}/ResultsFactory.h
${SRC_DIR}/ResultsRawReader.cpp
${SRC_DIR}/ResultsRawReader.h
${SRC_DIR}/ResultsReader.cpp
${SRC_DIR}/ResultsReader.h
${SRC_DIR}/ResultsWriter.cpp
${SRC_DIR}/ResultsWriter.h
${SRC_DIR}/Types.h
)

add_library (${PROJECT} STATIC ${SOURCES})
