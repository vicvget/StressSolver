#include "NodeLinks.h"


// union LinkDirection

/**
* Создать объект по умолчанию
* @return объект по умолчанию
*/
// static
LinkDirection LinkDirection::Create()
{
	LinkDirection linkDirection;

	linkDirection._plusLink = NO_LINK;
	linkDirection._minusLink = NO_LINK;

	return linkDirection;
}

/**
* Создать объект с заданными значениями полей
* @param plusLink - ссылка на соседний узел по направлению
* @param minusLink - ссылка на соседний узел против направления
* @return созданный объект с заданными значениями полей
*/
// static
LinkDirection LinkDirection::Create
	(
		int plusLink,
		int minusLink
	)
{
	LinkDirection linkDirection;

	linkDirection._plusLink = plusLink;
	linkDirection._minusLink = minusLink;

	return linkDirection;
}

/**
* Создать объект с заданными значениями полей
* @param links - массив ссылок на соседние узлы вдоль данного направления
* @return созданный объект с заданными значениями полей
*/
// static
LinkDirection LinkDirection::Create
	(
		const int links[]
	)
{
	LinkDirection linkDirection;

	for (int variant = 0; variant < VariantsCount; variant++)
	{
		linkDirection._links[variant] = links[variant];
	}

	return linkDirection;
}

/**
* Получить ссылку на соседний узел по направлению
* @return ссылка на соседний узел по направлению
*/
int LinkDirection::Plus() const
{
	return _plusLink;
}

/**
* Получить доступ к ссылке на соседний узел по направлению
* @return ссылка на соседний узел по направлению
*/
int& LinkDirection::Plus()
{
	return _plusLink;
}

/**
* Получить ссылку на соседний узел против направления
* @return ссылка на соседний узел против направления
*/
int LinkDirection::Minus() const
{
	return _minusLink;
}

/**
* Получить доступ к ссылке на соседний узел против направления
* @return ссылка на соседний узел против направления
*/
int& LinkDirection::Minus()
{
	return _minusLink;
}

/**
* Получить ссылку на соседний узел по данному варианту связи вдоль направления
* @param direction - вариант связи вдоль направления
* @return ссылка на соседний узел по данному варианту связи вдоль направления
*/
int LinkDirection::Link
	(
		DirectionVariantEnumeration variant
	)	const
{
	return _links[variant];
}

/**
* Получить доступ к ссылке на соседний узел по данному варианту связи вдоль направления
* @param direction - вариант связи вдоль направления
* @return ссылка на соседний узел по данному варианту связи вдоль направления
*/
int& LinkDirection::Link
	(
		DirectionVariantEnumeration variant
	)
{
	return _links[variant];
}

/**
* Оператор получения ссылки на соседний узел по данному варианту связи вдоль направления
* @param direction - вариант связи вдоль направления
* @return ссылка на соседний узел по данному варианту связи вдоль направления
*/
int LinkDirection::operator[]
	(
		DirectionVariantEnumeration variant
	)	const
{
	return _links[variant];
}

/**
* Оператор получения доступа к ссылке на соседний узел по данному варианту связи вдоль направления
* @param direction - вариант связи вдоль направления
* @return ссылка на соседний узел по данному варианту связи вдоль направления
*/
int& LinkDirection::operator[]
	(
		DirectionVariantEnumeration variant
	)
{
	return _links[variant];
}


// union NodeLinks

/**
* Конструктор по умолчанию
*/
NodeLinks::NodeLinks()
	:
		_xDirection(LinkDirection::Create()),
		_yDirection(LinkDirection::Create()),
		_zDirection(LinkDirection::Create())
{
}

/**
* Конструктор
* @param xDirection - ссылки на соседние узлы вдоль оси Ox
* @param yDirection - ссылки на соседние узлы вдоль оси Oy
* @param zDirection - ссылки на соседние узлы вдоль оси Oz
*/
NodeLinks::NodeLinks
	(
		const LinkDirection& xDirection,
		const LinkDirection& yDirection,
		const LinkDirection& zDirection
	)
	:
		_xDirection(xDirection),
		_yDirection(yDirection),
		_zDirection(zDirection)
{
	
}

/**
* Конструктор
* @param directions - массив ссылок на соседние узлы по направлениям
*/
NodeLinks::NodeLinks
	(
		const LinkDirection directions[]
	)
{
	for (int direction = 0; direction < DirectionsCount; direction++)
	{
		_directions[direction] = directions[direction];
	}
}

/**
* Конструктор
* @param links - список ссылок на связанные узлы в виде массива
*/
NodeLinks::NodeLinks
	(
		const int links[]
	)
{
	for (int linkIndex = 0; linkIndex < _linksCount; linkIndex++)
	{
		_links[linkIndex] = links[linkIndex];
	}
}

/**
* Получить ссылки на соседние узлы по направлению вдоль оси Ox
* @return набор ссылок на соседние узлы по направлению вдоль оси Ox
*/
const LinkDirection& NodeLinks::X() const
{
	return _xDirection;
}

/**
* Получить доступ к ссылкам на соседние узлы по направлению вдоль оси Ox
* @return набор ссылок на соседние узлы по направлению вдоль оси Ox
*/
LinkDirection& NodeLinks::X()
{
	return _xDirection;
}

/**
* Получить ссылки на соседние узлы по направлению вдоль оси Oy
* @return набор ссылок на соседние узлы по направлению вдоль оси Oy
*/
const LinkDirection& NodeLinks::Y() const
{
	return _yDirection;
}

/**
* Получить доступ к ссылкам на соседние узлы по направлению вдоль оси Oy
* @return набор ссылок на соседние узлы по направлению вдоль оси Oy
*/
LinkDirection& NodeLinks::Y()
{
	return _yDirection;
}

/**
* Получить ссылки на соседние узлы по направлению вдоль оси Oz
* @return набор ссылок на соседние узлы по направлению вдоль оси Oz
*/
const LinkDirection& NodeLinks::Z() const
{
	return _zDirection;
}

/**
* Получить доступ к ссылкам на соседние узлы по направлению вдоль оси Oz
* @return набор ссылок на соседние узлы по направлению вдоль оси Oz
*/
LinkDirection& NodeLinks::Z()
{
	return _zDirection;
}

/**
* Получить ссылки на соседние узлы по данному направлению
* @param direction - направление для связи с соседними узлами
* @return набор ссылок на соседние узлы по данному направлению
*/
const LinkDirection& NodeLinks::Direction
	(
		DirectionEnumeration direction
	)	const
{
	return _directions[direction];
}

/**
* Получить доступ к ссылкам на соседние узлы по данному направлению
* @param direction - направление для связи с соседними узлами
* @return набор ссылок на соседние узлы по данному направлению
*/
LinkDirection& NodeLinks::Direction
	(
		DirectionEnumeration direction
	)
{
	return _directions[direction];
}

/**
* Оператор получения ссылки на соседние узлы по данному направлению
* @param direction - направление для связи с соседними узлами
* @return набор ссылок на соседние узлы по данному направлению
*/
const LinkDirection& NodeLinks::operator[]
	(
		DirectionEnumeration direction
	)	const
{
	return _directions[direction];
}

/**
* Оператор доступа к ссылкам на соседние узлы по данному направлению
* @param direction - направление для связи с соседними узлами
* @return набор ссылок на соседние узлы по данному направлению
*/
LinkDirection& NodeLinks::operator[]
	(
		DirectionEnumeration direction
	)
{
	return _directions[direction];
}

/**
* Получить ссылку на соседний узел по индексу связи
* @param linkIndex - индекс связи
* @return ссылка на соседний узел по данному индексу связи
*/
int NodeLinks::Link
	(
		int linkIndex
	)	const
{
	return _links[linkIndex];
}

/**
* Получить доступ к ссылке на соседний узел по индексу связи
* @param linkIndex - индекс связи
* @return ссылка на соседний узел по данному индексу связи
*/
int& NodeLinks::Link
	(
		int linkIndex
	)
{
	return _links[linkIndex];
}

/**
* Оператор получения ссылки на соседний узел по индексу связи
* @param linkIndex - индекс связи
* @return ссылка на соседний узел по данному индексу связи
*/
int NodeLinks::operator[]
	(
		int linkIndex
	)	const
{
	return _links[linkIndex];
}

/**
* Оператор получения доступа к ссылке на соседний узел по индексу связи
* @param linkIndex - индекс связи
* @return ссылка на соседний узел по данному индексу связи
*/
int& NodeLinks::operator[]
	(
		int linkIndex
	)
{
	return _links[linkIndex];
}

/**
* Получить общее максимально возможное количество связанных с данным узлом узлов сетки
* @return общее максимально возможное количество связанных с данным узлом узлов сетки
*/
// static
int NodeLinks::LinksCount()
{
	return _linksCount;
}