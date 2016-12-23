# CMake project file
cmake_minimum_required(VERSION 2.6) 
set (PROJECT MathHelpers) 
 
project (${PROJECT}) 

# утстановка пути до директории с исходными файлами
set (SRC_DIR ${PROJECT_SOURCE_DIR} )

set (LIBRARY_OUTPUT_PATH ${SRC_DIR}/../../_Debug/bin/)

set (SOURCES
	${SRC_DIR}/MathFunctions.cpp
	${SRC_DIR}/Matrix3x3.cpp
	${SRC_DIR}/Matrix3x4.cpp
	${SRC_DIR}/Vector3.cpp
	${SRC_DIR}/Vector3Data.cpp
	${SRC_DIR}/Vector3Ref.cpp
)

add_library (${PROJECT} STATIC ${SOURCES})
