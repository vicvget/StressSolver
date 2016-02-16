#ifndef VECTOR_3_REF_H

#define VECTOR_3_REF_H


namespace MathHelpers
{

	/**
	* Класс для хранения ссылки на массив double* или const double*,
	* интерпретируемый как последовательные компоненты вектора 3 x 1
	*/
	template
		<
			typename Real // тип вещественного числа
		>
	class Vector3Reference
	{
	public:

		Vector3Reference();

		Vector3Reference
			(
				Real* data
			);


		Real* _data;

	};


	// ссылка на неконстантный массив
	typedef Vector3Reference<double> Vector3Ref;

	// ссылка на константный массив
	typedef Vector3Reference<const double> Vector3CRef;

}


#endif // VECTOR_3_REF_H
