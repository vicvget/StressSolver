#include "ProviderRezr.h"

//#include "../../Fcore/Exceptions/fcExceptions.h"


/**
* Конструктор
* @param rezrName - имя файла с результатами
* @param fadresName - имя файла со структурой результатов
*/
ProviderRezr::ProviderRezr
	(
		const string& rezrName,
		const string& fadresName
	)
	:
		AbstractProviderRezr(rezrName, fadresName),
		_powerStationFileStream(nullptr)
{
	
}

/** Деструктор
*/
//virtual
ProviderRezr::~ProviderRezr() // override
{
	if(_powerStationFileStream != nullptr)
	{
		delete _powerStationFileStream;
	}
}

/**
* Открытие файла результатов
* @return false в случае ошибки
*/
// virtual
bool ProviderRezr::OpenStream() // override
{
	if (_powerStationFileStream != nullptr)
	{
		delete _powerStationFileStream;
	}
	_powerStationFileStream = new PowerStationFileStream(_rezrName);

	return _powerStationFileStream->Open(FileStream::OpenMode::ReadOnly);
}

/**
* Закрыть поток результатов
*/
// virtual
void ProviderRezr::CloseStream() // override
{
	_powerStationFileStream->Close();
}

/**
* Считывание заголовка
*/
// virtual
void ProviderRezr::ReadHeader() // override
{
	*_powerStationFileStream
		>> _header.simv
		>> _header.v
		>> _header.hz
		>> _header.tk
		>> _header.nvar
		>> _header.kgr
		>> _header.so;
	CalcFramesCount();
}

/**
* Чтение кадра
* @return false, если кадр не прочитан целиком или конец потока
*/
// virtual
bool ProviderRezr::TryReadFrame() // override
{
	if(_currentFrameNumber >= _framesCount)
	{
		return false;
	}
	*_powerStationFileStream >> _currentTime;
	// TODO: здесь оптимизировать на чтение куска памяти
	for (int i = 0; i < _header.kgr; ++i)
	{
		*_powerStationFileStream >>_frameBuffer[i];
	}
	if (_powerStationFileStream->IsEof())
	{
		return false;
	}
	++_currentFrameNumber;

	return true;
}