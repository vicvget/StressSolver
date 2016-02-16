#pragma once
//#include

template<class T>
class IMemoryAllocator
{
public:
	static virtual T* AllocVector(std::size_t size) = 0;
	static virtual T* AllocMatrix(std::size_t size) = 0;
	static virtual void Free(T* data) = 0;
	static unsigned int GetVectorStride() = 0;
	static unsigned int GetMatrixStride() = 0;
};