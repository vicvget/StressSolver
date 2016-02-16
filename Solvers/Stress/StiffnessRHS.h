#ifndef STIFFNESS_RHS_H

#define STIFFNESS_RHS_H


#include "AuxiliaryStressStuff.h"

#include <vector>
#include <string>


// ������ ������ ������ ���� �� ������ ������� ��������� (RHS - right hand side)
using StiffnessRHSVector = std::vector<double>;


/**
* ��������� ������ ������� ������ ������������ �������� ������� ��� ��������� ��������� ������������� ����
* � ������� ������ ������ ���� �� ������ ������� ��������� (��������������� �������� ������� ����� ��������)
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
* ��������� ��������� �������������
*/
void ApplySealedElementsToStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector,
		const IsSealedFlagsList& isSealedFlagsList
	);

/**
* ������������� ������ ������ ������ ���� �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void InvertStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector
	);

/**
* �������� ������ ������ ������ ���� �� ������ ������� ��������� � ��������� ����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void WriteStiffnessRHSVectorToTextFile
	(
		const std::string& stiffnessRHSTextFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	);

/**
* �������� ������ ������ ������ ���� �� ������ ������� ��������� � �������� ����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void WriteStiffnessRHSVectorToBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	);

/**
* ������� ������ ������ ������ ���� �� ������ ������� ��������� �� ��������� �����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @return ��������� ������ ������ ������ ���� �� ������ ������� ���������
*/
StiffnessRHSVector ReadStiffnessRHSVectorFromBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName
	);


#endif // STIFFNESS_RHS_H