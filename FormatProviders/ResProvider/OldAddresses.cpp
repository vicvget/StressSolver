#include "OldAddresses.h"

#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"

#include <cstdio>


//класс списков смещений различных объектов
// class OldAddresses

/**
* Конструктор по умолчанию
*/
OldAddresses::OldAddresses()
	:
		_ksts(0),
		_kgr(0)
{

}

/**
* Заполнение списков из файла fadres.dat
* @param fadresName - имя файла fadres (возможно, с путем к нему)
*/
void OldAddresses::ReadFadres
	(
		const string& fadresName
	)
{
	FILE* f; //считываемый файл
	int objectType; //тип элемента
	int bodyNumber; //номер тела
	int connectionElementNumber; //номер соединительного элемента
	int rotationMatrixNumber; //номер тела для соеответствующей матрицы поворота
	int freedom; //номер степени свободы
	int address; //адрес данной степени свободы в файле "rezr"
	int emptyField; //пустое поле
	
	Bodies::iterator bt; //итератор, по которому находится очередное тело
	ConnectionElements::iterator ct; //итератор, по которому находится очередной с.э.
	RotationMatrixes::iterator mt; //итератор, по которому находится очередная матрица поворота
	
	f = fopen(fadresName.c_str(), "rt"); //открываем файл (для чтения в текстовом режиме)
	if (f != nullptr) //если файл успешно открылся
	{
		_bodies.clear();
		_connectionElements.clear();
		_rotationMatrixes.clear();
		fscanf(f, "%d", &_ksts); //считываем число степеней свободы
		_kgr = 3 * _ksts;
		while (!feof(f)) //пока не достигнут конец файла
		{
			//считываем тип элемента
			fscanf(f, "%d", &objectType);
			switch (static_cast<ObjectType>(objectType))
			{
			case ObjectType::Body:
				fscanf(f, "%d%d%d", &bodyNumber, &freedom, &address);
				// находим этот тело
				// (если массив еще не содержит данного тела, оно будет добавлено)
				// и вставляем в него новую пару "степень свободы - смещение"
				_bodies[bodyNumber].emplace
					(
						static_cast<FreedomType>(freedom),
						address
					);
				break;

			case ObjectType::ConnectionElement:
				fscanf(f, "%d%d%d", &connectionElementNumber, &freedom, &address);
				// находим этот соединительный элемент
				// (если массив еще не содержит данного элемента, он будет добавлен)
				// и вставляем в него новую пару "степень свободы - смещение"
				_connectionElements[connectionElementNumber].emplace
					(
						static_cast<FreedomType>(freedom),
						address
					);
				if (address + 3 > _kgr)
				{
					_kgr = address + 3;
				}
				break;

			case ObjectType::RotationMatrix:
				fscanf(f, "%d%d%d", &rotationMatrixNumber, &address, &emptyField);
				// находим эту матрицу поворота
				// (если массив еще не содержит данной матрицы поворота, она будет добавлена)
				// и вставляем в нее новый адрес
				_rotationMatrixes[rotationMatrixNumber] = address;
				if (address + 8 > _kgr)
				{
					_kgr = address + 8;
				}
				break;
			}
		}
		fclose(f); //закрываем файл
		BuildBodiesDofOffsetsCache();
		++_kgr;
	}
}

/**
* Получить значение числа степеней свободы для тел
* @return ksts - число степеней свободы для тел
*/
int OldAddresses::GetKsts() const
{
	return _ksts;
}

/**
* Получить значение общего числа степеней свободы
* @return kgr - общее число степеней свободы
*/
int OldAddresses::GetKgr() const
{
	return _kgr;
}

/**
* Получение списка смещений для соединительных элементов
* @return список смещений для соединительных элементов
*/
const ConnectionElements& OldAddresses::GetConnectionElements() const
{
	return _connectionElements;
}

/**
* Получение списка смещений для матриц поворота
* @return список смещений для матриц поворота
*/
const RotationMatrixes& OldAddresses::GetRotationMatrixes() const
{
	return _rotationMatrixes;
}

/**
* Проверить, есть ли данная степень свободы у тела с данным номером
* @param bodyNumber - номер тела
* @param freedom - тип степени свободы
* @return признак, есть ли данная степень свободы у тела с данным номером (true) или ее нет (false)
*/
bool OldAddresses::CheckBodyFreedom
	(
		Integer bodyNumber,
		FreedomType freedom
	)	const
{
	Bodies::const_iterator bodyIterator = _bodies.find(bodyNumber);

	if (bodyIterator == _bodies.end())
	{
		return false;
	}

	const Freedoms& freedoms = bodyIterator->second;
	Freedoms::const_iterator freedomIterator = freedoms.find(freedom);

	return freedomIterator != freedoms.end();
}

/**
* Проверить, есть ли данная степень свободы у тела с данным номером
* @param bodyNumber - номер тела
* @param motion - тип движения
* @param direction - тип направления движения
* @return признак, есть ли данная степень свободы у тела с данным номером (true) или ее нет (false)
*/
bool OldAddresses::CheckBodyFreedom
	(
		Integer bodyNumber,
		MotionType motion,
		DirectionType direction
	)	const
{
	FreedomType freedom = ConstructFreedom(motion, direction);

	return CheckBodyFreedom(bodyNumber, freedom);
}

/**
* Проверить, есть ли данная степень свободы у соединительного элемента с данным номером
* @param connectionElementNumber - номер соединительного элемента
* @param freedom - тип степени свободы
* @return смещение
*/
bool OldAddresses::CheckConnectionElementFreedom
	(
		Integer connectionElementNumber,
		FreedomType freedom
	)	const
{
	ConnectionElements::const_iterator connectionElementIterator =
		_connectionElements.find(connectionElementNumber);

	if (connectionElementIterator == _connectionElements.end())
	{
		return false;
	}

	const Freedoms& freedoms = connectionElementIterator->second;
	Freedoms::const_iterator freedomIterator = freedoms.find(freedom);

	return freedomIterator != freedoms.end();
}

/**
* Проверить, есть ли данная степень свободы у соединительного элемента с данным номером
* @param connectionElementNumber - номер соединительного элемента
* @param motion - тип движения
* @param direction - тип направления движения
* @return смещение
*/
bool OldAddresses::CheckConnectionElementFreedom
	(
		Integer connectionElementNumber,
		MotionType motion,
		DirectionType direction
	)	const
{
	FreedomType freedom = ConstructFreedom(motion, direction);

	return CheckConnectionElementFreedom(connectionElementNumber, freedom);
}

/**
* Получение точного смещения для тела
* @param bodyNumber - номер тела
* @param freedom - тип степени свободы
* @param value - величина
* @return смещение
*/
Integer OldAddresses::GetBodyAddress
	(
		Integer bodyNumber,
		FreedomType freedom,
		Value value
	)	const
{
	Bodies::const_iterator it;

	it = _bodies.find(bodyNumber);
	if (it != _bodies.end())
	{
		const Freedoms& freedoms = it->second;

		Freedoms::const_iterator jt;

		jt = freedoms.find(freedom);
		if (jt != freedoms.end())
		{
			return jt->second + ToIntegralType(value) * _ksts - 1;
		}
	}
	return -1;
}

/**
* Получение точного смещения для тела
* @param bodyNumber - номер тела
* @param motion - тип движения
* @param direction - тип направления движения
* @param value - величина
* @return смещение
*/
Integer OldAddresses::GetBodyAddress
	(
		Integer bodyNumber,
		MotionType motion,
		DirectionType direction,
		Value value
	)	const
{
	FreedomType freedom = ConstructFreedom(motion, direction);

	return GetBodyAddress(bodyNumber, freedom, value);
}

/**
* Получение точного смещения для соединительного элемента
* @param connectionElementNumber - номер соединительного элемента
* @param freedom - тип степени свободы
* @param component - компонента характеристики силы
* @return смещение
*/
Integer OldAddresses::GetConnectionElementAddress
	(
		Integer connectionElementNumber,
		FreedomType freedom,
		CharacteristicComponent component
	)	const
{
	ConnectionElements::const_iterator it;

	it = _connectionElements.find(connectionElementNumber);
	if (it != _connectionElements.end())
	{
		const Freedoms& freedoms = it->second;

		Freedoms::const_iterator jt;

		jt = freedoms.find(freedom);
		if (jt != freedoms.end())
		{
			return jt->second + ToIntegralType(component) - 1;
		}
	}

	return -1;
}

/**
* Получение точного смещения для соединительного элемента
* @param connectionElementNumber - номер соединительного элемента
* @param motion - тип движения
* @param direction - тип направления движения
* @param component - компонента характеристики силы
* @return смещение
*/
Integer OldAddresses::GetConnectionElementAddress
	(
		Integer connectionElementNumber,
		MotionType motion,
		DirectionType direction,
		CharacteristicComponent component
	)	const
{
	FreedomType freedom = ConstructFreedom(motion, direction);

	return GetConnectionElementAddress(connectionElementNumber, freedom, component);
}

/**
* Получение точного значения смещения для матрицы поворота
* @param rotationMatrixNumber - номер матрицы поворота
* @return смещение
*/
Integer OldAddresses::GetRotationMatrixAddress
	(
		Integer rotationMatrixNumber
	)	const
{
	RotationMatrixes::const_iterator it;

	it = _rotationMatrixes.find(rotationMatrixNumber);
	if (it != _rotationMatrixes.end())
	{
		return it->second - 1;
	}
	return -1;
}

/**
* Получение списка смещений для степеней свободы тел
* @return список смещений для степеней свободы тел
*/
const Bodies& OldAddresses::GetBodies() const
{
	return _bodies;
}

/**
* Получение списка кэшированных степеней свободы для тел
* @return кэш смещений для тел
*/
const BodiesDofOffsetsCache& OldAddresses::GetBodiesDofOffsetsCache() const
{
	return _bodiesDofOffsetsCache;
}

/**
* Заполнение кэша степеней свободы для тел (-1, если отсутствует)
*/
void OldAddresses::BuildBodiesDofOffsetsCache()
{
	_bodiesDofOffsetsCache.clear();
	for (const auto& body : _bodies)
	{
		int bodyNumber = body.first;
		BodyDofOffsetsCache offsetsCache;

		for (int positionOffset = 0; positionOffset < 6; ++positionOffset)
		{
			offsetsCache.PositionOffsets[positionOffset] = GetBodyAddress
				(
					bodyNumber,
					static_cast<FreedomType>(positionOffset + 1),
					Value::Shift
				);
		}
		offsetsCache.RotationOffset = GetRotationMatrixAddress(bodyNumber);
		_bodiesDofOffsetsCache.emplace(bodyNumber, offsetsCache);
	}
}
