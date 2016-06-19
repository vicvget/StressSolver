#include "ResultsWriter.h"

//#include "../../Fcore/Exceptions/fcExceptions.h"
#include <algorithm>

ResultsWriter::ResultsWriter
	(
		const string& filename
	)
{
	Open(filename);
}

ResultsWriter::~ResultsWriter()
{
	Close();
}

void ResultsWriter::AllocateFloatBuffer
	(
		size_t bufferSize
	)
{
	_buf = std::unique_ptr<float[]>(new float[bufferSize]);
}

void ResultsWriter::Open
	(
		const string& filename
	)
{
	_ofs.open(filename, ofstream::binary | ofstream::out | ofstream::trunc);
	if(!_ofs.is_open())
	{
		throw "FileNotFound";
		//exceptions::ThrowFileNotFound(filename);
	}
}

void ResultsWriter::Close()
{
	if(_ofs.is_open())
	{
		_ofs.close();
	}
}

void ResultsWriter::WriteBuffer
	(
		const void* buffer,
		size_t bufferSize
	)
{
	if (_ofs.is_open())
	{
		_ofs.flush();
		_ofs.write(static_cast<const char*>(buffer), bufferSize);
	}
}

void ResultsWriter::WriteFloatBuffer(const double* dbuf, float* fbuf, size_t bufferSize)
{
	std::transform(dbuf, dbuf + bufferSize, fbuf,
		[](const double element)
	{
		return static_cast<float>(element);
	});
	WriteBuffer(fbuf, bufferSize * sizeof(float));
}

/** Записать буфер double с преобразованием во float
* @param buffer - буфер
* @param bufferSize - размер буфера в элементах
*/
void ResultsWriter::WriteBufferToFloat
	(
		const double* buffer,
		size_t bufferSize
	)
{
	auto tmp = std::unique_ptr<float[]>(new float[bufferSize]);
	float* pbuf = _buf == nullptr ? tmp.get() : _buf.get();
	WriteFloatBuffer(buffer, pbuf, bufferSize);
}