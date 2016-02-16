#pragma once
#include "IMemory.h"

class RegularMemory: public IMemory
{
protected:
	size_t _vecStride;	// смещение векторов
	size_t _matStride;	// смещения матриц поворота
	size_t _varStride;	// смещение векторов неизвестных для узла
	size_t _dataStride	// смещение до скоростей

public:
#pragma region overriden

float* GetMemoryPointer() const;
int GetMemorySize() const;
double GetData(int nNode, int nDof, DataType dataType);
double* GetDataInternal(DataType dataType) const;

#pragma endregion overriden
};