#pragma once

class IMemory
{
	virtual float* GetMemoryPointer() const = 0;
	virtual unsigned int GetMemorySize() const = 0;
	virtual double GetData(int nNode, int nDof, DataType dataType) = 0;
	virtual double* GetDataInternal(DataType dataType) const = 0;
};