#ifndef STIFFNESS_MATRIX_FILE_FORMAT_H

#define STIFFNESS_MATRIX_FILE_FORMAT_H


#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"

#include <cstdint>


// ������������� ���, ��������������� ��� �������� ������ �������� ������ ������� ���������
using StiffnessMatrixFileFormats = uint32_t;


// ������ ����� ������� ��������� � �������� �������� ���������������� ������������
// �������� ������ ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT(fileFormatName, _) \
	fileFormatName,

// ������������ ������, ���������� �� ������������� ������� ����� ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT

// ��������������� ������������ �������� ������ ������� ���������
enum class StiffnessMatrixFileFormat_ : size_t
{

#	include "StiffnessMatrixFileFormatList.h"

	FormatsCount

};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT


// ���������, ������ ���������� �������� ������ ������� ���������
const size_t StiffnessMatrixFileFormatsCount = TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat_::FormatsCount);


// ������ ����� ������� ��������� � �������� �������� ������������ �������� ������ ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT(fileFormatName, _) \
	fileFormatName = 1 << TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat_::fileFormatName),

// ������������ ������, ���������� �� ������������� ������� ����� ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT

// ������������ �������� ������ ������� ���������
enum class StiffnessMatrixFileFormat : StiffnessMatrixFileFormats
{
#	include "StiffnessMatrixFileFormatList.h"

	// ����� ������ �������� ������ ������ ���������
	EndFormat = 1 << StiffnessMatrixFileFormatsCount

};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_ENUM_ELEMENT




#endif