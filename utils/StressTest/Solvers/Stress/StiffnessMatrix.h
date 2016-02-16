#ifndef STIFFNESS_MATRIX_H

#define STIFFNESS_MATRIX_H


#include "AuxiliaryStressStuff.h"
#include "StiffnessMatrixIO.h"


#include <memory>
#include <string>


// ������� 6 x 6, ������� ������������ ��� �������� ��������� ������� ���������
using Matrix6x6 = double[FreedomsCount][FreedomsCount];

// ��������� �� ������� 6 x 6
using PMatrix6x6 = std::unique_ptr<Matrix6x6[]>;

// ������� ��������� �� ������� 6 x 6
using PPMatrix6x6 = std::unique_ptr<PMatrix6x6[]>;


/**
* ������� ���������, ������������ �� ��������� ������� 6 x 6
* (���������� ����� ��������� ����� _nodesCount x _nodesCount)
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrix)
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
	StiffnessMatrix
		(
			size_t nodesCount
		);

	/**
	* ����������� �����������
	* @param otherStiffnessMatrix - ������� ���������, �� ������� ������������ �����������
	*/
	StiffnessMatrix
		(
			StiffnessMatrix&& otherStiffnessMatrix
		);


	// ������ � ���������� �������

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


	// ������ � ��������� ������� ���������
	
	/**
	* �������� ��������� �� ��� ������� ��������� � �������� rowIndex (��� ������� �� ��������� 6 x 6)
	* @param rowIndex ������ ���� ������� ���������
	* @return ��������� �� ��� ������� ��������� � �������� rowIndex
	*/
	PMatrix6x6& operator[]
		(
			size_t rowIndex
		);


private:

	// ���������� ����� � �������� ������������� ����
	size_t _nodesCount{};

	// ������� ��������� �� ������� 6 x 6, �� �������� �������� ������� ���������
	PPMatrix6x6 _matrix{};


	// ��������������� �������

	/**
	* �������� ������ ��� �������� ������� ���������
	* (�� ���������� ����� ����� � �������� ������������� _nodesCount)
	*/
	void Allocate();

};


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrix
extern
template class StiffnessMatrixIO<StiffnessMatrix>;


#endif // STIFFNESS_MATRIX_H