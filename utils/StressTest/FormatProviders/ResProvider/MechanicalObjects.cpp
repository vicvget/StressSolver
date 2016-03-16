#include "MechanicalObjects.h"

#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"


/**
* Сформировать степень свободы по типу движения и направлению
* @param motion - тип движения
* @param direction - тип направления движения
* @return сформированная степень свободы
*/
FreedomType ConstructFreedom
	(
		MotionType motion,
		DirectionType direction
	)
{
	return static_cast<FreedomType>((motion == MotionType::Linear ? 1 : 4) + ToIntegralType(direction));
}