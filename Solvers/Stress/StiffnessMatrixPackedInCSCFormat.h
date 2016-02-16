#ifndef STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H

#define STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H


#include "StiffnessMatrixIO.h"

#include <vector>


/**
* ����� ��� �������� ������� ���������, ����������� � CSC (compressed sparse column) �������
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrixPackedInCSCFormat)
{
public:

	// ����������� ���� ������

	// �������, ����������� � ����� �����, � ������� ������������ ��� �� �������� �������� ������� ���������
	static
	const std::string _fileNamePostfix;


	// ������������ � ����������

	/**
	* ����������� �� ���������
	*/
	StiffnessMatrixPackedInCSCFormat() = default;

	/**
	* �����������
	* @param nodesCount ���������� ����� ��������� ������������� ����
	*/
	StiffnessMatrixPackedInCSCFormat
		(
			size_t nodesCount
		);


	// ���������

	/**
	* �������� ������ ������ (�������������) ������� ���������
	* @return ������ ������ (�������������) ������� ���������
	*/
	size_t GetMatrixSize() const;

	/**
	* �������� ���������� ��������� ��������� � ������� ���������
	* @return ���������� ��������� ��������� � ������� ���������
	*/
	size_t GetMatrixNonZeroElementsCount() const;

	/**
	* �������� ��������� �� ������ ��������� ������� ���������
	* @return ��������� �� ������ ��������� ������� ���������
	*/
	const double* GetElements() const;

	/**
	* ����������� �������� ������� ��������� � ������ ������
	* @param elements - ������, � ������� ������������ ����������� ��������� ������� ���������
	*/
	void CopyElements
		(
			double* elements
		)	const;

	/**
	* ����������� ������� ����� ��������� � ������� ��������� � ������ ������
	* @param rowIds - ������, � ������� ������������ ����������� �������� ����� ��������� � ������� ���������
	*/
	void CopyRowIds
		(
			int* rowIds
		)	const;

	/**
	* ����������� ������� ������ ��������� ���������� ������� � ���������� ������� � ������ ������
	* @param columnsFirstIds - ������, � ������� ������������ ����������� �������� ������ ���������
	* ���������� ������� � ���������� ������� � ������ ������
	*/
	void CopyColumnFirstElements
		(
			int* columnsFirstIds
		)	const;


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

	// �������� ��������� ����������� ������� ���������
	std::vector<double> _values;

	// ������� ����� ��������� � ������� ���������
	std::vector<size_t> _rows;

	// ������� ������ ��������� ���������� ������� � ���������� �������
	std::vector<size_t> _columnBeginnings = std::vector<size_t>{0};

};


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrixPackedInCSCFormat
extern
template class StiffnessMatrixIO<StiffnessMatrixPackedInCSCFormat>;


#endif // STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H