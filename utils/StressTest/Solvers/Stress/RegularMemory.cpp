#if 0
#include "RegularMemory.h"
#include "RegularMemoryAllocator.h"
#include <cmath>

RegularMemory::RegularMemory(const int nNodes)
		_nElements(nNodes)
		_dataSize(nNodes * 4),
{
	_data = new float[_dataSize];	// ��������� ������:  3x-������, 1x-��������� ��������
	_dataInternal = RegularMemoryAllocator::AllocVector(nNodes * 3); // 3 ������� ��� ���������, ���������, �����������, �� 2 ������� (�������������� � ������������)
	_stress = RegularMemoryAllocator::AllocVector(nNodes);

	_vecStride = RegularMemoryAllocator::GetVectorStride();
	_matStride = vecStride * 3;
	_varStride = vecStride * 2;
	_dataStride = RegularMemoryAllocator::GetVectorStride()*2*nNodes;
}

RegularMemory::~RegularMemory()
{
	RegularMemoryAllocator::Free(_data);
	RegularMemoryAllocator::Free(_dataInternal);
	RegularMemoryAllocator::Free(_stress);
}

float* RegularMemory::GetMemoryPointer() const
{
	return _data;
}

unsigned int RegularMemory::GetMemorySize() const
{
	return _dataSize * sizeof(float);
}

double* RegularMemory::GetDataInternal(DataType dataType) const
{
	return _dataInternal + _dataStride * (int)dataType;
}

double RegularMemory::GetData(int nNode, int nDof, DataType dataType)
{
	return _dataInternal + _dataStride * (int)dataType + nNode * _varStride + nDof;
}

/** �������� ��������� ��������
* @param data - ������ ��� ������ ���������� ���������
*/
// virtual
void RegularMemory::GetScalarParameter
	(
		float* data
	)
{
	GetStressesByVonMises(data);
}
#endif