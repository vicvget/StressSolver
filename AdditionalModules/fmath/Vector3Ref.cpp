#include "Vector3Ref.h"


namespace MathHelpers
{

	// class Vector3Reference

	template
		<
			typename Real // тип вещественного числа
		>
	Vector3Reference<Real>::Vector3Reference()
		:
			_data(nullptr)
	{

	}

	template
		<
			typename Real // тип вещественного числа
		>
	Vector3Reference<Real>::Vector3Reference
		(
			Real* data
		)
		:
			_data(data)
	{
	}


	template
	class Vector3Reference<double>;

	template
	class Vector3Reference<const double>;

}
