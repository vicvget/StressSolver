#ifndef FOR_EACH_TUPLE_ELEMENT_H

#define FOR_EACH_TUPLE_ELEMENT_H


#include "IntegerSequence.h"
#include "ForEachParameterInPack.h"

#include <tuple>


/**
* ���������� ������� ���������� �������� ��� ������� �������� �������
* @param tuple ������, ��� ���������� �������� ����������� �������
* @param functor �������, ������� ����������� ��� ���������� �������
* @param indexSequence ��������������� ������,
* ����������� ��� ��������� ������������������ �������� ��������� �������
*/
template
	<
		class Tuple, // ��� �������
		typename Functor, // ��� ��������
		size_t... Indices // ������������������ �������� ��������� �������
	>
void ForEachTupleElement_
	(
		Tuple& tuple,
		Functor&& functor,
		IndexSequence<Indices...>&& indexSequence
	)
{
	Do{(std::forward<Functor>(functor)(std::get<Indices>(tuple)), void(), 0)...};
}


/**
* ��������� ������� ��� ������� �������� �������
* @param tuple ������, ��� ���������� �������� ����������� �������
* @param functor �������, ������� ����������� ��� ���������� �������
*/
template
	<
		class Tuple, // ��� �������
		typename Functor // ��� ��������
	>
void ForEachTupleElement
	(
		Tuple& tuple,
		Functor&& functor
	)
{
	ForEachTupleElement_
		(
			tuple,
			functor,
			IncreasingIndexSequence<std::tuple_size<Tuple>::value>{}
		);
}


#endif // FOR_EACH_TUPLE_ELEMENT_H