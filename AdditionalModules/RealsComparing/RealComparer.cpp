#include "RealComparer.h"

#include <cmath>
#include <algorithm>

#include "RealErrors.h"
#include "../AuxiliaryModules/Enumerations.h"


// template class RealComparer

/**
* Конструктор по умолчанию
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
RealComparer<Real, UseReferences>::RealComparer()
{
	FillErrors();
}

/**
* Конструктор
* @param absoluteError - минимальная абсолютная ошибка сравнений
* @param relativeError - минимальная относительная ошибка сравнений
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
RealComparer<Real, UseReferences>::RealComparer
	(
		RealArgument absoluteError,
		RealArgument relativeError
	)
{
	// TODO: обработать результат некорректной инициализации
	Init(absoluteError, relativeError);
}

/**
* Инициализация компаратора значениями ошибок по умолчанию
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
void RealComparer<Real, UseReferences>::Init()
{
	FillErrors();
}

/**
* Инициализация компаратора
* @param absoluteError - минимальная абсолютная ошибка сравнений
* @param relativeError - минимальная относительная ошибка сравнений
* @return признак успешной инициализации
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::Init
	(
		RealArgument absoluteError,
		RealArgument relativeError
	)
{
	if (!SetAbsoluteError(absoluteError))
	{
		return false;
	}
	if (!SetRelativeError(relativeError))
	{
		return false;
	}
	FillBoundaryValue();

	return true;
}

/**
* Получить минимальную абсолютную ошибку сравнений
* @return - минимальная абсолютная ошибка сравнений
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
typename RealComparer<Real, UseReferences>::RealArgument
RealComparer<Real, UseReferences>::GetAbsoluteError() const
{
	return _absoluteError;
}

/**
* Получить минимальную относительную ошибку сравнений
* @return - минимальная относительная ошибка сравнений
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
typename RealComparer<Real, UseReferences>::RealArgument
RealComparer<Real, UseReferences>::GetRelativeError() const
{
	return _relativeError;
}

/**
* Установить минимальную абсолютную ошибку сравнений
* @param absoluteError - минимальная абсолютная ошибка сравнений
* @return признак успешной установки
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::SetAbsoluteError
	(
		RealArgument absoluteError
	)
{
	return SetError(_absoluteError, absoluteError);
}

/**
* Установить минимальную относительную ошибку сравнений
* @param relativeError - минимальная относительная ошибка сравнений
* @return признак успешной установки
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::SetRelativeError
	(
		RealArgument relativeError
	)
{
	if (relativeError >= static_cast<Real>(1.0L))
	{
		return false;
	}

	return SetError(_relativeError, relativeError);
}

/**
* Сравнение на равенство с нулем
* @param number - сравниваемое число
* @return результат сравнения (true - число равно 0, false - число не равно 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::IsZero
	(
		RealArgument number
	)	const
{
	if (number == static_cast<Real>(0.0L))
	{
		return true;
	}

	return std::fabs(number) < _absoluteError;
}

/**
* Сравнение на неравенство с нулем
* @param number - сравниваемое число
* @return результат сравнения (true - число не равно 0, false - число равно 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::IsNotZero
	(
		RealArgument number
	)	const
{
	return !IsZero(number);
}

/**
* Сравнение, определяющее, является ли число-аргумент меньшим нуля
* @param number - сравниваемое число
* @return результат сравнения (true - число меньше 0, false - число не меньше 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::LessThanZero
	(
		RealArgument number
	)	const
{
	return number < -_absoluteError;
}

/**
* Сравнение, определяющее, является ли число-аргумент большим нуля
* @param number - сравниваемое число
* @return результат сравнения (true - число больше 0, false - число не больше 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::GreaterThanZero
	(
		RealArgument number
	)	const
{
	return number > _absoluteError;
}

/**
* Сравнение, определяющее, является ли число-аргумент меньшим или равным нулю
* @param number - сравниваемое число
* @return результат сравнения (true - число меньше или равно 0, false - число больше 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::LessThanOrEqualToZero
	(
		RealArgument number
	)	const
{
	return !GreaterThanZero(number);
}

/**
* Сравнение, определяющее, является ли число-аргумент большим или равным нулю
* @param number - сравниваемое число
* @return результат сравнения (true - число больше или равно 0, false - число меньше 0)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::GreaterThanOrEqualToZero
	(
		RealArgument number
	)	const
{
	return !LessThanZero(number);
}

/**
* Сравнение на равенство
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - числа равны, false - числа неравны)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::Equal
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	if (number1 == number2)
	{
		return true;
	}

	Real maximalNumber;

	maximalNumber = std::max
		(
			std::fabs(number1),
			std::fabs(number2)
		);
	if (maximalNumber < _boundaryValue)
	{
		return std::fabs(number1 - number2) < _absoluteError;
	}
	else
	{
		return std::fabs(number1 - number2) < _relativeError * maximalNumber;
	}
}

/**
* Сравнение на неравенство
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - числа неравны, false - числа равны)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::Unequal
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	return !Equal(number1, number2);
}

/**
* Сравнение на "меньше"
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - первое число меньше второго, false - второе)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::Less
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	return Greater(number2, number1);
}

/**
* Сравнение на "больше"
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - первое число больше второго, false - второе)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::Greater
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	if (number1 == number2)
	{
		return false;
	}
	
	Real maximalNumber;

	maximalNumber = std::max
		(
			std::fabs(number1),
			std::fabs(number2)
		);
	if (maximalNumber < _boundaryValue)
	{
		return number1 > number2 + _absoluteError;
	}
	else
	{
		return number1 > number2 + _relativeError * maximalNumber;
	}
}

/**
* Сравнение на "меньше или равно"
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - первое число меньше или равно второму, false - второе)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::LessOrEqual
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	return !Greater(number1, number2);
}

/**
* Сравнение на "больше или равно"
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @return результат сравнения (true - первое число больше или равно второму, false - второе)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::GreaterOrEqual
	(
		RealArgument number1,
		RealArgument number2
	)	const
{
	return !Less(number1, number2);
}

/**
* Обертка над функциями сравнения, использующая теги для определения типа сравнения
* @param number1 - первое сравниваемое число
* @param number2 - второе сравниваемое число
* @param comparisonTag - тег, используемый для определения типа сравнения
* @return результат сравнения (true - первое число больше или равно второму, false - второе)
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::operator()
	(
		RealArgument number1,
		RealArgument number2,
		ComparisonTag comparisonTag
	)	const
{
	const std::size_t comparisonIndex = ToIntegralType(comparisonTag);
	const ComparisonMember comparisonMember = _comparisonMembers[comparisonIndex];
	
	return (this->*comparisonMember)(number1, number2);
}

// тег сравнения в качестве элемента списка инициализации
#define COMPARISON_TAG_INITIALIZER_ELEMENT(tagName) \
	&RealComparer<Real, UseReferences>::tagName,

// определелить макрос, отвечающий за развертывание тега сравнения
#define COMPARISON_TAG COMPARISON_TAG_INITIALIZER_ELEMENT

// массив функций-членов сравнения вещественных чисел
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
// static
const typename RealComparer<Real, UseReferences>::ComparisonMember
RealComparer<Real, UseReferences>::_comparisonMembers[] =
	{
#		include "ComparisonTagList.h"
	};

#undef COMPARISON_TAG

#undef COMPARISON_TAG_ENUM_ELEMENT

/**
* Заполнить значения минимальных абсолютной и относительной ошибок сравнения
* данными класса RealErrors
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
void RealComparer<Real, UseReferences>::FillErrors()
{
	_absoluteError = RealErrors<Real>::absoluteError;
	_relativeError = RealErrors<Real>::relativeError;
	FillBoundaryValue();
}

/**
* Заполнить граничное значение при переходе от одной ошибки к другой
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
void RealComparer<Real, UseReferences>::FillBoundaryValue()
{
	_boundaryValue = _absoluteError / _relativeError;
}

/**
* Установить минимальную ошибку сравнений
* @param error - минимальная ошибка сравнений
* @param errorValue - новое значение ошибки сравнений
* @return признак успешной установки
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
bool RealComparer<Real, UseReferences>::SetError
	(
		Real& error,
		RealArgument errorValue
	)
{
	if (errorValue <= 0)
	{
		return false;
	}
	error = errorValue;
	FillBoundaryValue();

	return true;
}


// Явное инстанцирование класса RealComparer для вещественнозначных типов

// типы, использующие передачу параметров по ссылке в функциях сравнения

// для типа float
template
class RealComparer<float, true>;

// для типа double
template
class RealComparer<double, true>;

// для типа long double
template
class RealComparer<long double, true>;


// типы, использующие передачу параметров по значению в функциях сравнения

// для типа float
template
class RealComparer<float, false>;

// для типа double
template
class RealComparer<double, false>;

// для типа long double
template
class RealComparer<long double, false>;