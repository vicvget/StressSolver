#ifndef STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H

#define STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H


#include "StiffnessMatrixIO.h"

#include <vector>
#include <string>


/**
* ������� ������� ���������, ����������� � ������������ �������
*/
struct MatrixElement final
{

	// ������ ������ ������� ���������, � ������� ���������� ������ �������
	int rowIndex;

	// ������ ������� ������� ���������, � ������� ���������� ������ �������
	int columnIndex;

	// �������� �������� ������� ��������� � ������� ���������
	double elementValue;

};


/**
* ����� ��� �������� ������� ���������, ����������� � ������������ �������
* (��� ������� �������� ������ ������� �������� ��� ������� � ��������)
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrixPackedInCoordinateFormat)
{
public:

	// ����������� ���� ������

	// �������, ����������� � ����� �����, � ������� ������������ ��� �� �������� �������� ������� ���������
	static
	const std::string _fileNamePostfix;


	// ������������ � ����������

	/**
	* �����������
	* @param nodesCount ���������� ����� ��������� ������������� ����
	*/
	StiffnessMatrixPackedInCoordinateFormat
		(
			size_t nodesCount
		);


	// ������ � ���������� �������

	/**
	* ��������� ������� ������� ���������
	* @param rowIndex ������ ������ ������� ���������
	* @param columnIndex ������ ������� ������� ���������
	* @param elementValue �������� ������������ �������� ������� ���������
	*/
	void FillElement
		(
			size_t rowIndex,
			size_t columnIndex,
			double elementValue
		);

	/**
	* ��������� ������� ������� ���������
	* @param rowIndex ������ ���� ������� ���������, ���������� �� ��������� 6 x 6
	* @param rowDofIndex ������ ������ ���������� 6 x 6
	* @param columnIndex ������ ������� ������� ���������, ��������� �� ��������� 6 x 6
	* @param columnDofIndex ������ ������� ���������� 6 x 6
	* @param elementValue �������� ������������ �������� ������� ���������
	*/
	void FillElement
		(
			size_t rowIndex,
			size_t rowDofIndex,
			size_t columnIndex,
			size_t columnDofIndex,
			double elementValue
		);


	// ������ � ���� / ������ �� �����

	/**
	* �������� ������� ��������� � ��������� ����
	* @param textFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
	*/
	void WriteToTextFile
		(
			const std::string& textFileName
		)	const;

	/**
	* �������� ������� ��������� � ��������� ���� ������� GNU Octave
	* @param octaveTextFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
	*/
	void WriteToOctaveTextFile
		(
			const std::string& octaveTextFileName
		)	const;

	/**
	* �������� ������� ��������� � �������� ����
	* @param binaryFileName ������������ ��������� �����, � ������� ������������ ������ ������� ���������
	*/
	void WriteToBinaryFile
		(
			const std::string& binaryFileName
		)	const;

	/**
	* ������� ������� ��������� �� ��������� �����
	* @param binaryFileName ������������ ��������� �����, �� �������� ������������ ���������� ������� ���������
	*/
	void ReadFromBinaryFile
		(
			const std::string& binaryFileName
		);


private:

	// ������ ������ (�������������) ������� ���������
	size_t _matrixSize{};

	// �������� ����������� ������� ���������
	std::vector<MatrixElement> _matrixElements;

};


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrixPackedInCoordinateFormat
extern
template class StiffnessMatrixIO<StiffnessMatrixPackedInCoordinateFormat>;


#endif // STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H