#pragma once

template<class T>
class RegularMemoryAllocator: public IMemoryAllocator
{
public:
	static T* AllocVector(std::size_t size) = 0;
	static T* AllocMatrix(std::size_t size) = 0;
	static void Free(T* data) = 0;
	static unsigned int GetVectorStride() = 0;
	static unsigned int GetMatrixStride() = 0;
};

template<class T>
int GetVectorStride()
{
	return 3;
}

template<class T>
int GetMatrixStride()
{
	return 9;
}

template<class T>
void Free(T* data)
{
	if(data != NULL)
		delete [] data;
}


template<class T>
T* RegularMemoryAllocator::AllocVector(std::size_t size)
{
	T* res = new T[size*GetVectorStride()*2];
	memset(res, 0, size * sizeof(T));
	return res;
}

template<class T>
T* RegularMemoryAllocator::AllocMatrix(std::size_t size)
{
	T* res = new T[size*GetMatrixStride()];
	memset(res, 0, size * sizeof(T));
	return res;
}
