#include "StiffnessMatrices.h"

#include "../../AdditionalModules/AuxiliaryModules/ForEachTupleElement.h"


// template class StiffnessMatricesTuple

/**
*  онструктор
* @param nodesCount количество узлов сеточного представлени€ тела
*/
template
	<
		class... Matrices // классы матриц жесткости
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
*  онструктор перемещени€
* @param otherStiffnessMatricesTuple - другой кортеж матриц жесткости, из которого производитс€ перемещение
*/
template
	<
		class... Matrices // классы матриц жесткости
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
* ѕолучить кортеж матриц жесткости
* @return кортеж матриц жесткости
*/
template
	<
		class... Matrices // классы матриц жесткости
	>
const typename StiffnessMatricesTuple<Matrices...>::Tuple& StiffnessMatricesTuple<Matrices...>::GetTuple() const
{
	return _tuple;
}

/**
* ѕолучить кортеж матриц жесткости
* @return кортеж матриц жесткости
*/
template
	<
		class... Matrices // классы матриц жесткости
	>
typename StiffnessMatricesTuple<Matrices...>::Tuple& StiffnessMatricesTuple<Matrices...>::GetTuple()
{
	return _tuple;
}


// –еализаци€ методов кортежа матриц жесткости

// —формировать наименование функтора, использующегос€ дл€ вызова функции матрицы - элемента кортежа
#define FUNCTOR_NAME(Function) \
	Function##Functor


// –еализаци€ метода кортежа матриц жесткости
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


// явна€ конкретизаци€ StiffnessMatrices

// дл€ типов StiffnessMatrix и StiffnessMatrixPackedInCoordinateFormat
template class StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;