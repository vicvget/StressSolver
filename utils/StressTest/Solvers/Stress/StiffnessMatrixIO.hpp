#ifndef STIFFNESS_MATRIX_IO_HPP

#define STIFFNESS_MATRIX_IO_HPP


#include "StiffnessMatrixIO.h"


/**
* �������� ������� ��������� � ����
* @param outputFileName ������������ �����, � ������� ������������ ������ ������� ���������
* @param fileFormats ����� �������� ������, � ������� ������ ������������� ������
*/
template
	<
		class Matrix // ��� ������� ���������
	>
void StiffnessMatrixIO<Matrix>::WriteToFile
	(
		const std::string& outputFileName,
		StiffnessMatrixFileFormats fileFormats
	)	const
{
	const Matrix& thisMatrix = static_cast<const Matrix&>(*this);
	const std::string fullOutputFileName = outputFileName + Matrix::_fileNamePostfix;
	StiffnessMatrixFileFormats currentFormat = 1;

	for
		(
			size_t formatIndex = 0;
			formatIndex < StiffnessMatrixFileFormatsCount;
			++formatIndex,
			currentFormat <<= 1
		)
	{
		const std::string outputFileNameWithExtension =
			fullOutputFileName + '.' + _fileFormatsExtensions[formatIndex];

		if (currentFormat & fileFormats)
		{
			(thisMatrix.*_writeMembers[formatIndex])(outputFileNameWithExtension);
		}
	}
}


// ������ ����� ������� ��������� � �������� �������� ������ �������������
// ������� �������-������ ������ ������� ���������, ��������������� ��� �� ������
#define STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT(fileFormatName, _) \
	&Matrix::WriteTo##fileFormatName,

// ������������ ������, ���������� �� ������������� ������� ����� ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT

// ������ �������-������ ������ ������� ���������, ��������������� ��� �� ������
template
	<
		class Matrix // ��� ������� ���������
	>
// static
const typename StiffnessMatrixIO<Matrix>::WriteMember StiffnessMatrixIO<Matrix>::_writeMembers[] =
	{
#		include "StiffnessMatrixFileFormatList.h"
	};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT


// ������ ����� ������� ��������� � �������� �������� ������ �������������
// ������� �����, ����������� ���������� ������ ������ ���������
#define STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT(_, extension) \
	extension,

// ������������ ������, ���������� �� ������������� ������� ����� ������� ���������
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT

// ������ �����, ���������� ���������� ������ ������ ���������
template
	<
		class Matrix // ��� ������� ���������
	>
// static
const std::string StiffnessMatrixIO<Matrix>::_fileFormatsExtensions[] =
	{
#		include "StiffnessMatrixFileFormatList.h"
	};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT


#endif // STIFFNESS_MATRIX_IO_HPP