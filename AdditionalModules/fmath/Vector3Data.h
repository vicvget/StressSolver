#ifndef VECTOR_3_DATA_H

#define VECTOR_3_DATA_H


#include "Strides.h"


namespace MathHelpers
{

	/**
	* Класс для хранения компонентов вектора 3 x 1
	*/
	class Vector3Data
	{
	public:

		Vector3Data();
		
		Vector3Data
			(
				const double* data
			);


		union
		{

			struct
			{

				double x;
				double y;
				double z;

			};

			double _data[STRIDE3];

		};

	};

}


#endif // VECTOR_3_DATA_H