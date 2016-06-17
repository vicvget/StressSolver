#include "AbstractProviderRezr.h"

//#include "../../Fcore/Exceptions/fcExceptions.h"

#include <cmath>
#include <cstring>


// struct RezrBodyPosition

/**
* Оператор сравнения двух положений тел
* @param other - другое положение тела
* @return признак, идентичны ли два положения тел
*/
bool RezrBodyPosition::operator ==
	(
		const RezrBodyPosition& other
	)	const 
{
	for (int i = 0; i < 6; ++i)
	{
		if (fabs(Coordinates[i] - other.Coordinates[i]) >= 1e-6)
		{
			return false;
		}
	}
	for (int i = 0; i < 9; ++i)
	{
		if (fabs(Rotation[i] - other.Rotation[i]) >= 1e-6)
		{
			return false;
		}
	}

	return true;
}

/**
* Конструктор
* @param rezrName - имя файла с результатами
* @param fadresName - имя файла со структурой результатов
*/
AbstractProviderRezr::AbstractProviderRezr
	(
		const string& rezrName,
		const string& fadresName
	)
	:
		_rezrName(rezrName),
		_oldAddresses(std::make_unique<OldAddresses>())
{
	_oldAddresses->ReadFadres(fadresName);
}

/**
* Перейти к кадру с данным номером
* @param frameNumber - номер кадра
* @return признак успешного (true) или неуспешного (false) перехода
*/
// virtual
bool AbstractProviderRezr::GoToFrame
	(
		int frameNumber
	)
{
	while (GetCurrentFrameNumber() < frameNumber)
	{
		if (!TryReadFrame())
		{
			return false;
		}
	}

	return true;
}

/**
* Получить данные заголовка
*/
RezrHeader AbstractProviderRezr::GetHeader() const
{
	return _header;
}

/**
* Получить по номеру тела данные о его положении
* @param bodyNumber - номер тела
* @return данные о положении тела в текущем кадре
*/
RezrBodyPosition AbstractProviderRezr::GetRezrBodyPosition
	(
		size_t bodyNumber
	)	const
{
	RezrBodyPosition rezrBodyPosition;
	BodiesDofOffsetsCache::const_iterator it;

	it = _oldAddresses->GetBodiesDofOffsetsCache().find(bodyNumber);
	if (it != _oldAddresses->GetBodiesDofOffsetsCache().end())
	{
		for (int i = 0; i < 6; ++i)
		{
			if (it->second.PositionOffsets[i] != -1)
			{
				rezrBodyPosition.Coordinates[i] = _frameBuffer[it->second.PositionOffsets[i]];
			}
			else 
			{
				rezrBodyPosition.Coordinates[i] = 0.0f;
			}
		}
		if (it->second.RotationOffset != -1)
		{
			std::memcpy
				(
					&(rezrBodyPosition.Rotation[0]),
					&_frameBuffer[it->second.RotationOffset],
					9 * sizeof(float)
				);
		}
		else
		{
			memset(rezrBodyPosition.Rotation, 0, 9 * sizeof(float));
			rezrBodyPosition.Rotation[0] = rezrBodyPosition.Rotation[4] = rezrBodyPosition.Rotation[8] = 1.0f;
		}
	}
	else
	{
		//exceptions::ThrowBodyNotFoundException(bodyNumber);
	}

	return rezrBodyPosition;
}

/**
* Изменить данные о положении тела по его номеру
* @param bodyNumber - номер тела
* @param bodyPosition - данные о положении тела в текущем кадре
* @return признак того, что тело с данным номером было найдено
*/
bool AbstractProviderRezr::GetRezrBodyPosition
	(
		size_t bodyNumber,
		RezrBodyPosition& rezrBodyPosition
	)	const
{
	BodiesDofOffsetsCache::const_iterator it;
	
	it = _oldAddresses->GetBodiesDofOffsetsCache().find(bodyNumber);
	if (it != _oldAddresses->GetBodiesDofOffsetsCache().end())
	{
		for (int i = 0; i < 6; ++i)
		{
			if (it->second.PositionOffsets[i] != -1)
			{
				rezrBodyPosition.Coordinates[i] = _frameBuffer[it->second.PositionOffsets[i]];
			}
		}
		if (it->second.RotationOffset != -1)
		{
			std::memcpy
				(
					&(rezrBodyPosition.Rotation[0]),
					&_frameBuffer[it->second.RotationOffset],
					9 * sizeof(float)
				);
		}
	}
	else
	{
		return false;
	}

	return true;
}

/**
* Получить список номеров тел
* @return список номеров тел
*/
NumbersOfBodies AbstractProviderRezr::GetNumbersOfBodies() const
{
	NumbersOfBodies numbersOfBodies;

	for
		(
			Bodies::const_iterator bt = _oldAddresses->GetBodies().begin();
			bt != _oldAddresses->GetBodies().end();
			++bt
		)
	{
		numbersOfBodies.push_back(bt->first);
	}

	return numbersOfBodies;
}

/**
* Получить текущее время (время, соответствующее текущему кадру)
* @return текущее время
*/
float AbstractProviderRezr::GetCurrentTime() const
{
	return _currentTime;
}

/**
* Получить номер текущего кадра (нумерация начианется с нуля)
* @return номер текущего кадра
*/
int AbstractProviderRezr::GetCurrentFrameNumber() const
{
	return _currentFrameNumber;
}

/**
* Получить списки смещений различных объектов
* @return списки смещений различных объектов
*/
const OldAddresses& AbstractProviderRezr::GetOldAddresses() const
{
	return *_oldAddresses;
}

/**
* Вычислить количество кадров и выделить память
*/
void AbstractProviderRezr::CalcFramesCount()
{
	_frameBuffer.resize(_header.kgr);

	// TODO: Check It!
	_framesCount = (int)(_header.tk / _header.hz - _header.so / _header.v);
}
