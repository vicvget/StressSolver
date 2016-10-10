#include "Vector3.h"

//#include "../RealsComparing/RealComparing.h"
#include "../Macros/AdditionalMacros.h"

#include <cmath>
#include <algorithm>


namespace MathHelpers
{

	// template class Vector3

	template
		<
			typename DataType
		>
	template
		<
			typename DataType2
		>
	Vector3<DataType>::Vector3
		(
			const Vector3<DataType2>& other
		)
		:
			_rep(other._rep._data)
	{
	}


	template
	Vector3<Vector3Data>::Vector3
		(
			const Vector3<Vector3Ref>& other
		);

	template
	Vector3<Vector3Data>::Vector3
		(
			const Vector3<Vector3CRef>& other
		);

	template
	Vector3<Vector3CRef>::Vector3
		(
			const Vector3<Vector3Data>& other
		);

	template
	Vector3<Vector3CRef>::Vector3
		(
			const Vector3<Vector3Ref>& other
		);


	template
		<
			typename DataType
		>
	Vector3<DataType>::Vector3
		(
			double* data
		)
		:
			_rep(data)
	{
	}

	template
		<
			typename DataType
		>
	template
		<
			typename
		>
	Vector3<DataType>::Vector3
		(
			const double* data
		)
		:
			_rep(data)
	{
	}

	template
	Vector3<Vector3Data>::Vector3
		(
			const double* data
		);

	template
	Vector3<Vector3CRef>::Vector3
		(
			const double* data
		);


	template <>
	Vector3<Vector3Data>::Vector3
		(
			double x,
			double y,
			double z
		)
	{
		_rep.x = x;
		_rep.y = y;
		_rep.z = z;
	}


	template
		<
			typename DataType
		>
	template
		<
			typename DataType2
		>
	Vector3<DataType>& Vector3<DataType>::operator = 
		(
			const Vector3<DataType2>& copy
		)
	{
		std::copy(copy._rep._data, copy._rep._data + STRIDE3, _rep._data);

		return *this;
	}

	template<typename DataType>
	Vector3<DataType>& Vector3<DataType>::operator /= (double rhs)
	{
		for (size_t coord = 0; coord < 3; coord++)
			this->_rep._data[coord] /= rhs;
		return *this;
	}

	template
		<
			typename DataType
		>
	template
		<
			typename DataType2
		>
	Vector3<DataType>& Vector3<DataType>::operator += (const Vector3<DataType2>& rhs)
	{
		for (size_t coord = 0; coord < 3; coord++)
			this->_rep._data[coord] += rhs[coord];
		return *this;
	}

	template
		<
			typename DataType
		>
	template
		<
			typename DataType2
		>
	Vector3<DataType>& Vector3<DataType>::operator -= (const Vector3<DataType2>& rhs)
	{
		for (size_t coord = 0; coord < 3; coord++)
			this->_rep._data[coord] -= rhs[coord];
		return *this;
	}	
	
	template
		<
			typename DataType
		>
	void Vector3<DataType>::Init
		(
			double x,
			double y,
			double z
		)
	{
		_rep._data[0] = x;
		_rep._data[1] = y;
		_rep._data[2] = z;
	}

	template
		<
			typename DataType
		>
	void Vector3<DataType>::Export
		(
			double* data
		)	const
	{
		std::copy(_rep._data, _rep._data + STRIDE3, data);
	}

	template
		<
			typename DataType
		>
	void Vector3<DataType>::Export
		(
			float* data
		)	const
	{
		std::transform
			(
				_rep._data,
				_rep._data + STRIDE3,
				data,
				[](double val)
				{
					return static_cast<float>(val);
				}
			);
	}


	template
		<
			typename DataType
		>
	double Vector3<DataType>::X() const
	{
		return _rep._data[0];
	}

	template
		<
			typename DataType
		>
	double Vector3<DataType>::Y() const
	{
		return _rep._data[1];
	}

	template
		<
			typename DataType
		>
	double Vector3<DataType>::Z() const
	{
		return _rep._data[2];
	}

	template
		<
			typename DataType
		>
	double Vector3<DataType>::operator []
		(
			unsigned int id
		)	const
	{
		return _rep._data[id];
	}

	template
		<
			typename DataType
		>
	double& Vector3<DataType>::operator []
		(
			unsigned int id
		)
	{
		return _rep._data[id];
	}

	template
		<
			typename DataType
		>
	Vec3Data Vector3<DataType>::operator - () const
	{
		Vec3Data vecNeg;

		for (int i = 0; i < STRIDE3; ++i)
		{
			vecNeg[i] = -_rep._data[i];
		}

		return vecNeg;
	}

	template
		<
			typename DataType
		>
	template
		<
			typename DataType2
		>
	Vec3Data Vector3<DataType>::Cross
		(
			const Vector3<DataType2>& vec
		)	const
	{
		Vec3Data vecCross;

		vecCross[0] =  _rep._data[1] * vec[2] - _rep._data[2] * vec[1];
		vecCross[1] = -_rep._data[0] * vec[2] + _rep._data[2] * vec[0];
		vecCross[2] =  _rep._data[0] * vec[1] - _rep._data[1] * vec[0];

		return vecCross;
	}


#	define METHOD_INSTANTIATION(ReturnType, Method, DataType1, DataType2) \
		template \
		ReturnType Vector3<DataType1>::Method \
			( \
				const Vector3<DataType2>& vec \
			)	const;

#	define SAME_METHOD_INSTANTIATION(Method, DataType1, DataType2) \
		template \
		Vector3<DataType1>& Vector3<DataType1>::Method \
			( \
				const Vector3<DataType2>& vec \
			);

#	define BINARY_VV_INSTANTIATION(ReturnType, Method, DataType1, DataType2) \
		template \
		ReturnType Method \
			( \
				const Vector3<DataType1>& vec1, \
				const Vector3<DataType2>& vec2 \
			);

#	define BINARY_SV_INSTANTIATION(ReturnType, Method, ScalarDataType, DataType) \
		template \
		ReturnType Method \
			( \
				const Vector3<DataType>& vec, \
				ScalarDataType scalar \
			); \
		 \
		template \
		ReturnType Method \
			( \
				ScalarDataType scalar, \
				const Vector3<DataType>& vec \
			);


#	ifndef MS_VC_COMPILER

#		define UNARY_FUNCTON_ISTANTIATION(Macro, ...) \
			Macro(__VA_ARGS__, Vector3Data) \
			Macro(__VA_ARGS__, Vector3Ref) \
			Macro(__VA_ARGS__, Vector3CRef)

#		define BINARY_FUNCTON_ISTANTIATION(Macro, ...) \
			Macro(__VA_ARGS__, Vector3Data, Vector3Data) \
			Macro(__VA_ARGS__, Vector3Data, Vector3Ref) \
			Macro(__VA_ARGS__, Vector3Data, Vector3CRef) \
			Macro(__VA_ARGS__, Vector3Ref, Vector3Data) \
			Macro(__VA_ARGS__, Vector3Ref, Vector3Ref) \
			Macro(__VA_ARGS__, Vector3Ref, Vector3CRef) \
			Macro(__VA_ARGS__, Vector3CRef, Vector3Data) \
			Macro(__VA_ARGS__, Vector3CRef, Vector3Ref) \
			Macro(__VA_ARGS__, Vector3CRef, Vector3CRef)

#	else

		// обход бага MS Visual C++ с variadic macros

#		define UNARY_FUNCTON_ISTANTIATION(Macro, ...) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3CRef))

#		define BINARY_FUNCTON_ISTANTIATION(Macro, ...) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3CRef)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3CRef)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3CRef, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3CRef, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3CRef, Vector3CRef))

#		define BINARY_SAME_FUNCTON_ISTANTIATION(Macro, ...) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Data, Vector3CRef)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3Data)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3Ref)) \
			COMBINE(Macro, (__VA_ARGS__, Vector3Ref, Vector3CRef)) 

#	endif


#	define FULL_METHOD_INSTANTIATION(ReturnType, Method) \
		BINARY_FUNCTON_ISTANTIATION(METHOD_INSTANTIATION, ReturnType, Method)

#	define FULL_SAME_METHOD_INSTANTIATION(Method) \
		BINARY_SAME_FUNCTON_ISTANTIATION(SAME_METHOD_INSTANTIATION, Method)


#	define FULL_BINARY_VV_INSTANTIATION(ReturnType, Method) \
		BINARY_FUNCTON_ISTANTIATION(BINARY_VV_INSTANTIATION, ReturnType, Method)

#	define FULL_BINARY_SV_INSTANTIATION(ReturnType, Method) \
		UNARY_FUNCTON_ISTANTIATION(BINARY_SV_INSTANTIATION, ReturnType, Method, double)


	FULL_METHOD_INSTANTIATION(Vec3Data, Cross)
	FULL_SAME_METHOD_INSTANTIATION(operator +=)		
	FULL_SAME_METHOD_INSTANTIATION(operator -=)

	FULL_METHOD_INSTANTIATION(double, Angle)

	FULL_BINARY_VV_INSTANTIATION(bool, operator ==)
	FULL_BINARY_VV_INSTANTIATION(double, operator *)
	FULL_BINARY_VV_INSTANTIATION(Vec3Data, operator +)
	FULL_BINARY_VV_INSTANTIATION(Vec3Data, operator -)
	FULL_BINARY_VV_INSTANTIATION(Vec3Data, Cross)
	FULL_BINARY_VV_INSTANTIATION(Vec3Data, MakeVec3)

	FULL_BINARY_SV_INSTANTIATION(Vec3Data, operator *)
	FULL_BINARY_SV_INSTANTIATION(Vec3Data, operator /)


	template
		<
			typename DataType1,
			typename DataType2
		>
	bool operator ==
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		)
	{
		double magnitude = (vec1 - vec2).Magnitude();

		return DoubleComparing::IsZero(magnitude);
	}

	template
		<
			typename DataType
		>
	Vec3Data Vector3<DataType>::Abs() const
	{
		Vec3Data vecAbs = MakeVec3(_rep._data[0],_rep._data[1],_rep._data[2]);

		std::transform
			(
				vecAbs._rep._data,
				vecAbs._rep._data + STRIDE3,
				vecAbs._rep._data,
				[](double val)
				{
					return std::abs(val);
				}
			);

		return vecAbs;
	}

	template
		<
			typename DataType
		>
	Vec3Data Vector3<DataType>::Normalize() const
	{
		double magnitude = Magnitude();
		if(DoubleComparing::IsZero(magnitude))
			return MakeVec3(0., 0., 0.);
		return *this / magnitude;
	}


	template
		<
			typename DataType
		>
	Vec3Data Vector3<DataType>::Normal() const
	{
		Vec3Data vec1 =  MakeVec3(1., 0., 0.);
		if (DoubleComparing::IsZero(Cross(vec1).Magnitude()))
		{
			vec1[0] = 0.;
			vec1[1] = 1.;
		}
		return Cross(vec1);
	}

	template
		<
			typename DataType
		>
	double Vector3<DataType>::Magnitude() const
	{
		return sqrt(X()*X() + Y()*Y() + Z()*Z());
	}


	template
	class Vector3<Vector3Ref>;

	template
	class Vector3<Vector3Data>;


	Vec3Data MakeVec3(double x, double y, double z)
	{
		return {x, y, z};
	}

	Vec3Ref MakeVec3
		(
			double* x
		)
	{
		return x;
	}

	Vec3CRef MakeVec3
		(
			const double* x
		)
	{
		return x;
	}

	template
		<
			class DataType1,
			class DataType2
		>
	Vec3Data MakeVec3
		(
			const Vector3<DataType1>& vec1,
			const Vector3<DataType2>& vec2
		)
	{
		return vec2 - vec1;
	}

	Vec3Ref AssignRef
		(
			Vec3Data& vec
		)
	{
		return MakeVec3(vec._rep._data);
	}


	template <typename DataType1>
	template <typename DataType2>
	double Vector3<DataType1>::Angle(const Vector3<DataType2>& vec) const
	{
		double denominator = Magnitude() * vec.Magnitude();
		if (DoubleComparing::IsZero(denominator))
		{
			return 0.0;
		}

		double cosAngle = *this * vec / denominator;

		if (cosAngle < -1.0) cosAngle = -1.0;
		else if (cosAngle > 1.0) cosAngle = 1.0;
		return acos(cosAngle);
	}

	template <class DataType1, class DataType2>
	double operator * (const Vector3<DataType1>& vec1, const Vector3<DataType2>& vec2)
	{
		double res = 0.0;

		for (int i = 0; i < STRIDE3; ++i)
		{
			res += vec1[i] * vec2[i];
		}

		return res;
	}

	template <class DataType1, class DataType2>
	Vec3Data Cross(const Vector3<DataType1>& vec1, const Vector3<DataType2>& vec2)
	{
		return vec1.Cross(vec2);
	}

	template <class DataType1, class DataType2>
	Vec3Data operator + (const Vector3<DataType1>& vec1, const Vector3<DataType2>& vec2)
	{
		return MakeVec3(vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]);
	}

	template <class DataType1, class DataType2>
	Vec3Data operator - (const Vector3<DataType1>& vec1, const Vector3<DataType2>& vec2)
	{
		return MakeVec3(vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]);
	}

	template <class DataType>
	Vec3Data operator * (const Vector3<DataType>& vec, double scalar)
	{
		Vec3Data res;
		for (int i = 0; i < 3; i++)
		{
			res[i] = vec[i] * scalar;
		}
		return res;
	}

	template <class DataType>
	Vec3Data operator * (double scalar, const Vector3<DataType>& vec)
	{
		return vec * scalar;
	}

	template <class DataType>
	Vec3Data operator / (const Vector3<DataType>& vec, double scalar)
	{
		Vec3Data res;
		for (int i = 0; i < 3; i++)
		{
			res[i] = vec[i] / scalar;
		}
		return res;
	}

	template <class DataType>
	Vec3Data operator / (double scalar, const Vector3<DataType>& vec)
	{
		return vec / scalar;
	}


}
