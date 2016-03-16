#include "NodeLinks.h"


// union LinkDirection

/**
* ������� ������ �� ���������
* @return ������ �� ���������
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
* ������� ������ � ��������� ���������� �����
* @param plusLink - ������ �� �������� ���� �� �����������
* @param minusLink - ������ �� �������� ���� ������ �����������
* @return ��������� ������ � ��������� ���������� �����
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
* ������� ������ � ��������� ���������� �����
* @param links - ������ ������ �� �������� ���� ����� ������� �����������
* @return ��������� ������ � ��������� ���������� �����
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
* �������� ������ �� �������� ���� �� �����������
* @return ������ �� �������� ���� �� �����������
*/
int LinkDirection::Plus() const
{
	return _plusLink;
}

/**
* �������� ������ � ������ �� �������� ���� �� �����������
* @return ������ �� �������� ���� �� �����������
*/
int& LinkDirection::Plus()
{
	return _plusLink;
}

/**
* �������� ������ �� �������� ���� ������ �����������
* @return ������ �� �������� ���� ������ �����������
*/
int LinkDirection::Minus() const
{
	return _minusLink;
}

/**
* �������� ������ � ������ �� �������� ���� ������ �����������
* @return ������ �� �������� ���� ������ �����������
*/
int& LinkDirection::Minus()
{
	return _minusLink;
}

/**
* �������� ������ �� �������� ���� �� ������� �������� ����� ����� �����������
* @param direction - ������� ����� ����� �����������
* @return ������ �� �������� ���� �� ������� �������� ����� ����� �����������
*/
int LinkDirection::Link
	(
		DirectionVariantEnumeration variant
	)	const
{
	return _links[variant];
}

/**
* �������� ������ � ������ �� �������� ���� �� ������� �������� ����� ����� �����������
* @param direction - ������� ����� ����� �����������
* @return ������ �� �������� ���� �� ������� �������� ����� ����� �����������
*/
int& LinkDirection::Link
	(
		DirectionVariantEnumeration variant
	)
{
	return _links[variant];
}

/**
* �������� ��������� ������ �� �������� ���� �� ������� �������� ����� ����� �����������
* @param direction - ������� ����� ����� �����������
* @return ������ �� �������� ���� �� ������� �������� ����� ����� �����������
*/
int LinkDirection::operator[]
	(
		DirectionVariantEnumeration variant
	)	const
{
	return _links[variant];
}

/**
* �������� ��������� ������� � ������ �� �������� ���� �� ������� �������� ����� ����� �����������
* @param direction - ������� ����� ����� �����������
* @return ������ �� �������� ���� �� ������� �������� ����� ����� �����������
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
* ����������� �� ���������
*/
NodeLinks::NodeLinks()
	:
		_xDirection(LinkDirection::Create()),
		_yDirection(LinkDirection::Create()),
		_zDirection(LinkDirection::Create())
{
}

/**
* �����������
* @param xDirection - ������ �� �������� ���� ����� ��� Ox
* @param yDirection - ������ �� �������� ���� ����� ��� Oy
* @param zDirection - ������ �� �������� ���� ����� ��� Oz
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
* �����������
* @param directions - ������ ������ �� �������� ���� �� ������������
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
* �����������
* @param links - ������ ������ �� ��������� ���� � ���� �������
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
* �������� ������ �� �������� ���� �� ����������� ����� ��� Ox
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Ox
*/
const LinkDirection& NodeLinks::X() const
{
	return _xDirection;
}

/**
* �������� ������ � ������� �� �������� ���� �� ����������� ����� ��� Ox
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Ox
*/
LinkDirection& NodeLinks::X()
{
	return _xDirection;
}

/**
* �������� ������ �� �������� ���� �� ����������� ����� ��� Oy
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Oy
*/
const LinkDirection& NodeLinks::Y() const
{
	return _yDirection;
}

/**
* �������� ������ � ������� �� �������� ���� �� ����������� ����� ��� Oy
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Oy
*/
LinkDirection& NodeLinks::Y()
{
	return _yDirection;
}

/**
* �������� ������ �� �������� ���� �� ����������� ����� ��� Oz
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Oz
*/
const LinkDirection& NodeLinks::Z() const
{
	return _zDirection;
}

/**
* �������� ������ � ������� �� �������� ���� �� ����������� ����� ��� Oz
* @return ����� ������ �� �������� ���� �� ����������� ����� ��� Oz
*/
LinkDirection& NodeLinks::Z()
{
	return _zDirection;
}

/**
* �������� ������ �� �������� ���� �� ������� �����������
* @param direction - ����������� ��� ����� � ��������� ������
* @return ����� ������ �� �������� ���� �� ������� �����������
*/
const LinkDirection& NodeLinks::Direction
	(
		DirectionEnumeration direction
	)	const
{
	return _directions[direction];
}

/**
* �������� ������ � ������� �� �������� ���� �� ������� �����������
* @param direction - ����������� ��� ����� � ��������� ������
* @return ����� ������ �� �������� ���� �� ������� �����������
*/
LinkDirection& NodeLinks::Direction
	(
		DirectionEnumeration direction
	)
{
	return _directions[direction];
}

/**
* �������� ��������� ������ �� �������� ���� �� ������� �����������
* @param direction - ����������� ��� ����� � ��������� ������
* @return ����� ������ �� �������� ���� �� ������� �����������
*/
const LinkDirection& NodeLinks::operator[]
	(
		DirectionEnumeration direction
	)	const
{
	return _directions[direction];
}

/**
* �������� ������� � ������� �� �������� ���� �� ������� �����������
* @param direction - ����������� ��� ����� � ��������� ������
* @return ����� ������ �� �������� ���� �� ������� �����������
*/
LinkDirection& NodeLinks::operator[]
	(
		DirectionEnumeration direction
	)
{
	return _directions[direction];
}

/**
* �������� ������ �� �������� ���� �� ������� �����
* @param linkIndex - ������ �����
* @return ������ �� �������� ���� �� ������� ������� �����
*/
int NodeLinks::Link
	(
		int linkIndex
	)	const
{
	return _links[linkIndex];
}

/**
* �������� ������ � ������ �� �������� ���� �� ������� �����
* @param linkIndex - ������ �����
* @return ������ �� �������� ���� �� ������� ������� �����
*/
int& NodeLinks::Link
	(
		int linkIndex
	)
{
	return _links[linkIndex];
}

/**
* �������� ��������� ������ �� �������� ���� �� ������� �����
* @param linkIndex - ������ �����
* @return ������ �� �������� ���� �� ������� ������� �����
*/
int NodeLinks::operator[]
	(
		int linkIndex
	)	const
{
	return _links[linkIndex];
}

/**
* �������� ��������� ������� � ������ �� �������� ���� �� ������� �����
* @param linkIndex - ������ �����
* @return ������ �� �������� ���� �� ������� ������� �����
*/
int& NodeLinks::operator[]
	(
		int linkIndex
	)
{
	return _links[linkIndex];
}

/**
* �������� ����� ����������� ��������� ���������� ��������� � ������ ����� ����� �����
* @return ����� ����������� ��������� ���������� ��������� � ������ ����� ����� �����
*/
// static
int NodeLinks::LinksCount()
{
	return _linksCount;
}