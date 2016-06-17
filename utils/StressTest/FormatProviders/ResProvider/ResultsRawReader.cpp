#include "ResultsRawReader.h"

//#include "../../Fcore/Exceptions/fcExceptions.h"
#include <algorithm>

ResultsRawReader::ResultsRawReader
	(
		const string& filename
	)
{
	Open(filename);
}

ResultsRawReader::~ResultsRawReader()
{
	Close();
}

void ResultsRawReader::Open
	(
		const string& filename
	)
{
	_ifs.open(filename, ifstream::binary | ifstream::in);
	if(!_ifs.is_open())
	{
		throw "FileNotFound";
		//exceptions::ThrowFileNotFound(filename);
	}
}

void ResultsRawReader::Close()
{
	if(_ifs.is_open())
	{
		_ifs.close();
	}
}

bool ResultsRawReader::ReadBuffer
	(
		void* buffer,
		size_t bufferSize
	)
{
	if (_ifs.is_open())
	{
		_ifs.read(static_cast<char*>(buffer), bufferSize);
		return !(_ifs.eof());
	}
	return false;
}

