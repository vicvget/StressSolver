#include "MechanicalObjects.h"

#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"


/**
* ������������ ������� ������� �� ���� �������� � �����������
* @param motion - ��� ��������
* @param direction - ��� ����������� ��������
* @return �������������� ������� �������
*/
FreedomType ConstructFreedom
	(
		MotionType motion,
		DirectionType direction
	)
{
	return static_cast<FreedomType>((motion == MotionType::Linear ? 1 : 4) + ToIntegralType(direction));
}