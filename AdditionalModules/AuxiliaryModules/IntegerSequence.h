#ifndef INTEGER_SEQUENCE_H

#define INTEGER_SEQUENCE_H


/**
* Шаблонный класс, параметрами которого являются целочисленный тип
* и последовательность целочисленных значений
*/
template
	<
		typename Integer, // целочисленный тип
		Integer... Numbers // последовательность значений этого типа
	>
struct IntegerSequence
{
};


/**
* Последовательность индексов (индексы имеют тип size_t)
*/
template
	<
		size_t... Numbers // последовательность индексов
	>
using IndexSequence = IntegerSequence<size_t, Numbers...>;


/**
* Реализация генератора возрастающей целочисленной последовательности
*/
template
	<
		typename Integer, // целочисленный тип
		Integer CurrentNumber, // текущее целое число, достигнутое последовательностью
		Integer EndNumber, // конечное целое число, которое может быть достигнуто последовательностью
		Integer... Numbers // хвостовая последовательность значений этого типа
	>
struct IntegerSequenceGenerator_
	:
		IntegerSequenceGenerator_<Integer, CurrentNumber + 1, EndNumber, Numbers..., CurrentNumber + 1>
{
};


/**
* Реализация генератора возрастающей целочисленной последовательности
* (специализация окончания генерации)
*/
template
	<
		typename Integer, // целочисленный тип
		Integer EndNumber, // конечное целое число, которое может быть достигнуто последовательностью
		Integer... Numbers // хвостовая последовательность значений этого типа
	>
struct IntegerSequenceGenerator_<Integer, EndNumber, EndNumber, Numbers...>
{

	// тип целочисленной последовательности, сгенерированной генератором
	using type = IntegerSequence<Integer, Numbers...>;

};


/**
* Генератор возрастающей целочисленной последовательности
*/
template
	<
		typename Integer, // целочисленный тип
		Integer BeginNumber, // целое число, с которого должна начаться генерация последовательности
		Integer EndNumber // конечное целое число, которое может быть достигнуто последовательностью
	>
struct IntegerSequenceGenerator
	:
		IntegerSequenceGenerator_<Integer, BeginNumber, EndNumber, BeginNumber>
{
};


/**
* Возрастающая целочисленная последовательность
*/
template
	<
		typename Integer, // целочисленный тип
		Integer BeginNumber, // целое число, с которого должна начаться генерация последовательности
		Integer EndNumber // конечное целое число, которое может быть достигнуто последовательностью
	>
using IncreasingIntegerSequence = typename IntegerSequenceGenerator<Integer, BeginNumber, EndNumber>::type;


/**
* Возрастающая последовательность индексов (индексы имеют тип size_t, их нумерация начинается с 0)
*/
template
	<
		size_t Number // количество чисел в последовательности
	>
using IncreasingIndexSequence = IncreasingIntegerSequence<size_t, 0, Number - 1>;


#endif // INTEGER_SEQUENCE_H