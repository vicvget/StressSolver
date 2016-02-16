#pragma once
#include "IMemory.h"

class RegularMemory: public IMemory
{
protected:
	size_t _vecStride;	// �������� ��������
	size_t _matStride;	// �������� ������ ��������
	size_t _varStride;	// �������� �������� ����������� ��� ����
	size_t _dataStride	// �������� �� ���������

public:
#pragma region overriden

float* GetMemoryPointer() const;
int GetMemorySize() const;
double GetData(int nNode, int nDof, DataType dataType);
double* GetDataInternal(DataType dataType) const;

#pragma endregion overriden
};