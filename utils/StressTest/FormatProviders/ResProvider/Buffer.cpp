#include "Buffer.h"

#include <cstring>


//класс буфера
//class Buffer

//конструктор и деструктор

//конструктор
Buffer::Buffer
	(
		Integer size
	)
	:
		_buffer(nullptr),
		_size(size),
		_contentsLength(0),
		_shift(0)
{
	if (_size > 0)
	{
		_buffer = new char[_size];
	}
	else
	{
		_size = 0;
	}
}

//деструктор
Buffer::~Buffer()
{
	delete[] _buffer;
}

//селекторы

//получение указателя на буфер
char* Buffer::GetBuffer() const
{
	return _buffer;
}

//получение размера буфера
Integer Buffer::GetSize() const
{
	return _size;
}

//получение размера содержимого буфера
Integer Buffer::GetContentsSize() const
{
	return _contentsLength;
}

//возвращает истину, если буфер заполнен
bool Buffer::IsFilled() const
{
	return (_contentsLength >= _size);
}

//возвращает истину, если буфер пуст
bool Buffer::IsEmpty() const
{
	return (_shift >= _contentsLength);
}

//модификаторы

//устанавливает новое значение размера содержимого
void Buffer::SetContentsLength
	(
		Integer contentsLength
	)
{
	_contentsLength = contentsLength;
}

//функции для заполнения/опустошения буфера

//извлекает из буфера еще один блок
bool Buffer::Pop
	(
		void* block,
		Integer blockSize
	)
{
	if (_shift + blockSize > _contentsLength)
	{
		return false;
	}
	std::memcpy(block, _buffer + _shift, blockSize);
	_shift += blockSize;

	return true;
}

//заносит в буфер еще один блок
bool Buffer::Push
	(
		const void* block,
		Integer blockSize
	)
{
	if (_contentsLength + blockSize > _size)
	{
		return false;
	}
	std::memcpy(_buffer + _contentsLength, block, blockSize);
	_contentsLength += blockSize;

	return true;
}

//опустошает буфер
void Buffer::Flush()
{
	_contentsLength = 0;
	_shift = 0;
}