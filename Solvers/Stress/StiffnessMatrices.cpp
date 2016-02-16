#include "StiffnessMatrices.h"

#include "../../AdditionalModules/AuxiliaryModules/ForEachTupleElement.h"


// template class StiffnessMatricesTuple

/**
* �����������
* @param nodesCount ���������� ����� ��������� ������������� ����
*/
template
	<
		class... Matrices // ������ ������ ���������
	>
// explicit
StiffnessMatricesTuple<Matrices...>::StiffnessMatricesTuple
	(
		size_t nodesCount
	)
	:
		_tuple(Accumulate<Matrices>(nodesCount)...)
{
}

/**
* ����������� �����������
* @param otherStiffnessMatricesTuple - ������ ������ ������ ���������, �� �������� ������������ �����������
*/
template
	<
		class... Matrices // ������ ������ ���������
	>
StiffnessMatricesTuple<Matrices...>::StiffnessMatricesTuple
	(
		StiffnessMatricesTuple<Matrices...>&& otherStiffnessMatricesTuple
	)
	:
		_tuple(std::move(otherStiffnessMatricesTuple._tuple))
{
}

/**
* �������� ������ ������ ���������
* @return ������ ������ ���������
*/
template
	<
		class... Matrices // ������ ������ ���������
	>
const typename StiffnessMatricesTuple<Matrices...>::Tuple& StiffnessMatricesTuple<Matrices...>::GetTuple() const
{
	return _tuple;
}

/**
* �������� ������ ������ ���������
* @return ������ ������ ���������
*/
template
	<
		class... Matrices // ������ ������ ���������
	>
typename StiffnessMatricesTuple<Matrices...>::Tuple& StiffnessMatricesTuple<Matrices...>::GetTuple()
{
	return _tuple;
}


// ���������� ������� ������� ������ ���������

// ������������ ������������ ��������, ��������������� ��� ������ ������� ������� - �������� �������
#define FUNCTOR_NAME(Function) \
	Function##Functor


// ���������� ������ ������� ������ ���������
#define STIFFNESS_MATRICES_FUNCTION_IMPLEMENTATION(Function, ReturnType, ...) \
 \
struct FUNCTOR_NAME(Function) \
{ \
 \
	CREATE_FIELD_LIST(__VA_ARGS__) \
 \
	template \
		< \
			class StiffnessMatrix \
		> \
	void operator() \
		( \
			StiffnessMatrix& stiffnessMatrix \
		) \
	{ \
		stiffnessMatrix.Function \
			( \
				CREATE_ARGUMENTS_LIST(__VA_ARGS__) \
			); \
	} \
 \
}; \
 \
 \
template \
	< \
		class... Matrices \
	> \
ReturnType StiffnessMatricesTuple<Matrices...>::Function \
	( \
		CREATE_PARAMETERS_LIST(__VA_ARGS__) \
	) \
{ \
	ForEachTupleElement \
		( \
			_tuple, \
			FUNCTOR_NAME(Function) \
				{ \
					CREATE_ARGUMENTS_LIST(__VA_ARGS__) \
				} \
		); \
}

#define STIFFNESS_MATRICES_FUNCTION STIFFNESS_MATRICES_FUNCTION_IMPLEMENTATION

#include "StiffnessMatricesFunctions.h"

#undef STIFFNESS_MATRICES_FUNCTION

#undef STIFFNESS_MATRICES_FUNCTION_IMPLEMENTATION

#undef FUNCTOR_NAME


// ����� ������������� StiffnessMatrices

// ��� ����� StiffnessMatrix � StiffnessMatrixPackedInCoordinateFormat
template class StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;