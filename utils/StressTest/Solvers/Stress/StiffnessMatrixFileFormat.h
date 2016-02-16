#ifndef STIFFNESS_MATRIX_FILE_FORMAT_H

#define STIFFNESS_MATRIX_FILE_FORMAT_H


#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"

#include <cstdint>


// целочисленный тип, предназначенный для хранения списка форматов файлов матрицы жесткости
using StiffnessMatrixFileFormats = uint32_t;


// формат файла матрицы жесткости в качестве элемента вспомогательного перечисления
// форматов файлов матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT(fileFormatName, _) \
	fileFormatName,

// определелить макрос, отвечающий за развертывание формата файла матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT

// вспомогательное перечисление форматов файлов матрицы жесткости
enum class StiffnessMatrixFileFormat_ : size_t
{

#	include "StiffnessMatrixFileFormatList.h"

	FormatsCount

};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT


// константа, равная количеству форматов файлов матрицы жесткости
const size_t StiffnessMatrixFileFormatsCount = TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat_::FormatsCount);


// формат файла матрицы жесткости в качестве элемента перечисления форматов файлов матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT(fileFormatName, _) \
	fileFormatName = 1 << TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat_::fileFormatName),

// определелить макрос, отвечающий за развертывания формата файла матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT

// перечисление форматов файлов матрицы жесткости
enum class StiffnessMatrixFileFormat : StiffnessMatrixFileFormats
{
#	include "StiffnessMatrixFileFormatList.h"

	// конец списка форматов файлов матриц жесткости
	EndFormat = 1 << StiffnessMatrixFileFormatsCount

};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT




#endif