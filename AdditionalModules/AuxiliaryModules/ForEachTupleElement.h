#ifndef FOR_EACH_TUPLE_ELEMENT_H

#define FOR_EACH_TUPLE_ELEMENT_H


#include "IntegerSequence.h"
#include "ForEachParameterInPack.h"

#include <tuple>


/**
* Реализация функции выполнения функтора для каждого элемента кортежа
* @param tuple кортеж, над элементами которого выполняется функтор
* @param functor функтор, который выполняется над элементами кортежа
* @param indexSequence вспомогательный объект,
* необходимый для получения последовательности индексов элементов кортежа
*/
template
	<
		class Tuple, // тип кортежа
		typename Functor, // тип функтора
		size_t... Indices // последовательность индексов элементов кортежа
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
* Выполнить функтор для каждого элемента кортежа
* @param tuple кортеж, над элементами которого выполняется функтор
* @param functor функтор, который выполняется над элементами кортежа
*/
template
	<
		class Tuple, // тип кортежа
		typename Functor // тип функтора
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