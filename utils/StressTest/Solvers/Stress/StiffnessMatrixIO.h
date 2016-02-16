#ifndef STIFFNESS_MATRIX_IO_H

#define STIFFNESS_MATRIX_IO_H


#include "StiffnessMatrixFileFormat.h"

#include <string>


/**
* ����� ��� ����������� �����/������ ������� ��������� �/�� �����/-�� ��������� ��������
*/
template
	<
		class Matrix // ��� ������� ���������
	>
class StiffnessMatrixIO
{
public:

	// ������ � ���� / ������ �� �����

	/**
	* �������� ������� ��������� � ����
	* @param outputFileName ������������ �����, � ������� ������������ ������ ������� ���������
	* @param fileFormats ����� �������� ������, � ������� ������ ������������� ������
	*/
	void WriteToFile
		(
			const std::string& outputFileName,
			StiffnessMatrixFileFormats fileFormats
		)	const;


private:

	// ���������� ����

	// ������� ��� ������ ��� ����������� �����/������ ������� ���������
	using ThisType = StiffnessMatrixIO<Matrix>;

	// ��� �������-����� ������ ������� ���������, ���������������� ��� �� ������
	typedef
	void (Matrix::*WriteMember)
		(
			const std::string& outputFileName // ������������ �����, � ������� ��������� ������� ���������
		)	const;


	// ����������� ����������-�����

	// ������ �������-������ ������ ������� ���������, ��������������� ��� �� ������
	static
	const WriteMember _writeMembers[];

	// ������ �����, ���������� ���������� ������ ������ ���������
	static
	const std::string _fileFormatsExtensions[];

};


// �������� ����� ������� ���������
#define DECLARE_STIFFNESS_MATRIX(Matrix) \
 \
class Matrix final \
	: \
		public StiffnessMatrixIO<Matrix>


#include "StiffnessMatrixIO.hpp"


#endif // STIFFNESS_MATRIX_IO_H