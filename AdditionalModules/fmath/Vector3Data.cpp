#include "Vector3Data.h"

#include <algorithm>


namespace MathHelpers
{

	// class Vector3Data

	Vector3Data::Vector3Data()
		:
			x(),
			y(),
			z()
	{
	}

	Vector3Data::Vector3Data
		(
			const double* data
		)
	{
		std::copy(data, data + STRIDE3, _data);
	}

}