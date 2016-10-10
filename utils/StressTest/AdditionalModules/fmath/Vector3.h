#ifndef VECTOR_3_H
#define VECTOR_3_H


#include "Vector3Data.h"
#include "Vector3Ref.h"


namespace MathHelpers
{

	/**
	* Класс для вектора 3 x 1
	*/
	template
		<
			typename DataType
		>
	class Vector3;


	// Вектор 3 x 1, который сам хранит свои компоненты
	typedef Vector3<Vector3Data> Vec3Data;

	// Вектор 3 x 1, который хранит ссылку на массив double*,
	// интерпретируемый как последовательные компоненты векторы
	typedef Vector3<Vector3Ref> Vec3Ref;

	// Вектор 3 x 1, который хранит ссылку на массив const double*,
	// интерпретируемый как последовательные компоненты векторы
	typedef Vector3<Vector3CRef> Vec3CRef;

	// Вектор 3 x 1 по умолчанию
	typedef Vec3Data Vec3;


	/**
	* Класс для вектора 3 x 1
	*/
	template<typename DataType>
	class Vector3
	{
	public:

		Vector3() = default;

		Vector3
			(
				const Vector3<DataType>& other
			)
			= default;

		template
			<
				typename DataType2
			>
		Vector3
			(
				const Vector3<DataType2>& other
			);


		template
			<
				typename DataType2
			>
		Vector3<DataType>& operator =
			(
				const Vector3<DataType2>& copy
			);

		Vector3<DataType>& operator /= (double rhs);
		
		template<typename DataType2>
		//void operator += (const Vector3<DataType2>& rhs);
		Vector3<DataType>& operator += (const Vector3<DataType2>& rhs);
		
		template<typename DataType2>
		//void operator -= (const Vector3<DataType2>& rhs);
		Vector3<DataType>& operator -= (const Vector3<DataType2>& rhs);

		void Init
			(
				double x,
				double y,
				double z
			);

		void Export
			(
				double* data
			)	const;

		void Export
			(
				float* data
			)	const;


		double X() const;

		double Y() const;

		double Z() const;


		double operator []
			(
				unsigned int id
			)	const;

		double& operator []
			(
				unsigned int id
			);

		/* Инверсия компонентов
		* @return вектор с компонентами противоположного знака
		*/
		Vec3Data operator - () const;
	
		/* Модуль компонентов
		* @return вектор с компонентами, взятыми по модулю
		*/
		Vec3Data Abs() const;

		/* Нормализация вектора
		* @return вектор единичной длины
		*/
		Vec3Data Normalize() const;

		/* Нормальный вектор
		* @return нормальный вектор
		*/
		Vec3Data Normal() const;

		/* Длина вектора
		* @return длина вектора
		*/
		double Magnitude() const;

		template
			<
				typename DataType2
			>
		Vec3Data Cross
			(
				const Vector3<DataType2>& vec
			)	const;
		
		template
			<
				typename DataType2
			>
		double Angle
			(
				const Vector3<DataType2>& vec
			)	const;

		template
			<
				typename DataType2
			>
		friend
		class Vector3;
		
		friend
		Vector3<Vector3Data> MakeVec3();

		friend
		Vector3<Vector3Data> MakeVec3
			(
				double x,
				double y,
				double z
			);

		friend
		Vector3<Vector3Ref> MakeVec3
			(
				double* data
			);

		friend
		Vector3<Vector3CRef> MakeVec3
			(
				const double* data
			);

		friend
		Vector3<Vector3Ref> AssignRef
			(
				Vector3<Vector3Data>& vec
			);

	private:

		DataType _rep;

		Vector3
			(
				double* data
			);

		template
			<
				typename = void
			>
		Vector3
			(
				const double* data
			);

		Vector3
			(
				double x,
				double y,
				double z
			);

	};


	template <>
	template <>
	Vector3<Vec3Ref>::Vector3
		(
			const Vector3<Vec3Data>& other
		)
		= delete;

	template <>
	template <>
	Vector3<Vec3Ref>::Vector3
		(
			const Vector3<Vec3CRef>& other
		)
		= delete;


	template
		<
			class DataType1,
			class DataType2
		>
	Vec3Data operator +
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);
	
	template
		<
			class DataType1,
			class DataType2
		>
	Vec3Data operator -
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);
	
	template
		<
			class DataType
		>
	Vec3Data operator *
		(
			const Vector3<DataType>& vec,
			double scalar
		);
	
	template
		<
			class DataType
		>
	Vec3Data operator *
		(
			double scalar,
			const Vector3<DataType>& vec
		);

	template
		<
			class DataType
		>
	Vec3Data operator /
		(
			const Vector3<DataType>& vec,
			double scalar
		);

	template
		<
			class DataType
		>
	Vec3Data operator /
		(
			double scalar,
			const Vector3<DataType>& vec
		);

	template
		<
			class DataType1,
			class DataType2
		>
	Vec3Data Cross
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);
	
	template
		<
			class DataType1,
			class DataType2
		>
	double operator *
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);

	Vec3Data MakeVec3
		(
			double x,
			double y,
			double z
		);

	Vec3Ref MakeVec3
		(
			double* x
		);

	Vec3CRef MakeVec3
		(
			const double* x
		);

	template
		<
			class DataType1,
			class DataType2
		>
	Vec3Data MakeVec3
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);

	Vec3Ref AssignRef
		(
			Vec3Data& vec
		);

	template
		<
			class DataType1,
			class DataType2
		>
	bool operator ==
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		);

}


#endif // VECTOR_3_H
