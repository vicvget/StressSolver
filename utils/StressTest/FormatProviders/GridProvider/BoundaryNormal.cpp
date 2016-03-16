#include "BoundaryNormal.h"


/**
* Конструктор
*/
BoundaryNormal::BoundaryNormal
	(
		unsigned long pointNumber,
		const VertexPoint& unitaryNormal
	)
	:
		_pointNumber(pointNumber),
		_unitaryNormal(unitaryNormal)
{
}