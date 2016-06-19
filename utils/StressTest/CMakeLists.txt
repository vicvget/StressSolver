# CMake project file for avmodel lib
cmake_minimum_required(VERSION 2.6)

include(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")

set (PROJECT Stress) 

project (${PROJECT}) 

#SET(CMAKE_BUILD_TYPE release)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE release CACHE STRING
		"Choose the type of build, options are: None Debug Release."
		FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# FFLAGS depend on the compiler

if(MSVC) # MSVC compiler (Win32 only)
	# Display more warnings
	set(CMAKE_CXX_FLAGS "/W3")
elseif(UNIX OR CMAKE_COMPILER_IS_GNUCXX) # Clang OR Gcc (Linux, Mac OS or Win32 with MingW)
	# Enable C++11 and displays all warnings
	#set(CMAKE_CXX_FLAGS "-Wall -std=c++11 -Wextra -o3")
	set(CMAKE_CXX_FLAGS "-std=c++11 -openmp -fPIC -mmic -Wall -Wextra -pedantic -Wno-variadic-macros -Wno-unknown-pragmas -Wno-sign-compare -Wno-switch -O3")
endif()

add_definitions(-DNOMPI -DGNUCPP -DGFORTRAN -DGNU_CPP -DUSE_KNC)

#include_directories($(PROJECT_SOURCE_DIR)/fcore/wrappers)

# установка пути до директории с исходными файлами
set (SRC_DIR ${PROJECT_SOURCE_DIR} )

#set (LIBRARY_OUTPUT_PATH ${SRC_DIR}_Debug/bin/)
#set(CMAKE_DEBUG_POSTFIX "d")

add_subdirectory (AdditionalModules/fmath)
add_subdirectory (FormatProviders/GridProvider)
add_subdirectory (FormatProviders/ResProvider)

add_subdirectory (Solvers/Stress)
add_subdirectory (Test)

#target_link_libraries (${PROJECT} ИМЯ_БИБЛИОТЕКИ)
