#include "ResultsReader.h"


ResultsReader::ResultsReader
	(
		const string& rezrName,
		const string& fadresName
	)
	:
		AbstractProviderRezr(rezrName, fadresName)
{
}

// virtual
ResultsReader::~ResultsReader() // override
{
	if (_resultsFile.is_open())
	{
		_resultsFile.close();
	}
}

// virtual
bool ResultsReader::OpenStream() // override
{
	_resultsFile.open(_rezrName.c_str(), ifstream::binary);
	return _resultsFile.is_open();
}

// virtual
void ResultsReader::CloseStream() // override
{
	_resultsFile.close();
}

// virtual
void ResultsReader::ReadHeader() // override
{
	_resultsFile.read((char*)&_header, sizeof(RezrHeader));
	CalcFramesCount();
}

// virtual
bool ResultsReader::TryReadFrame() // override
{
	if (_resultsFile.eof())
	{
		return false;
	}
	if (_currentFrameNumber >= _framesCount)
	{
		return false;
	}
	_resultsFile.read((char*)&_currentTime, sizeof(float));
	if (_resultsFile.eof())
	{
		return false;
	}
	_resultsFile.read((char*)_frameBuffer.data(), _header.kgr * sizeof(float));
	if (_resultsFile.eof())
	{
		return false;
	}
	++_currentFrameNumber;

	return true;
}

/**
* Перейти к кадру с данным номером
* @param frameNumber - номер кадра
* @return признак успешного (true) или неуспешного (false) перехода
*/
// virtual
bool ResultsReader::GoToFrame
	(
		int frameNumber
	)	// override
{
	const auto currentPosition = _resultsFile.tellg();

	_resultsFile.seekg(0, std::ios::end);

	const auto fileSize = _resultsFile.tellg();
	const size_t currentFramePosition = sizeof(RezrHeader) + frameNumber * (_header.kgr + 1) * sizeof(float);
	const bool result = currentFramePosition < fileSize;

	if (result)
	{
		_currentFrameNumber = frameNumber;
		_resultsFile.seekg(currentFramePosition);
	}
	else
	{
		_resultsFile.seekg(currentPosition);
	}

	return result;
}