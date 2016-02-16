#ifndef REAL_COMPARER_H

#define REAL_COMPARER_H


#include <type_traits>

#include "ComparisonTag.h"


/**
* Класс, предназначенный для сравнения вещественных чисел разного типа (компаратор).
* Осуществляет сравнение чисел с учетом абсолютной и/или относительной погрешностей
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences = true // признак, передавать ли параметры в функции сравнения по ссылке
	>
class RealComparer
{
public:

	// Типы

	// тип параметра, передаваемого в функции сравнения
	typedef
		typename std::conditional
			<
				UseReferences,
				const Real&,
				Real
			>
		::type
		RealArgument;


	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	RealComparer();

	/**
	* Конструктор
	* @param absoluteError - минимальная абсолютная ошибка сравнений
	* @param relativeError - минимальная относительная ошибка сравнений
	*/
	RealComparer
		(
			RealArgument absoluteError,
			RealArgument relativeError
		);


	// Инициализация компаратора

	/**
	* Инициализация компаратора значениями ошибок по умолчанию
	*/
	void Init();

	/**
	* Инициализация компаратора
	* @param absoluteError - минимальная абсолютная ошибка сравнений
	* @param relativeError - минимальная относительная ошибка сравнений
	* @return признак успешной инициализации
	*/
	bool Init
		(
			RealArgument absoluteError,
			RealArgument relativeError
		);


	// Селекторы

	/**
	* Получить минимальную абсолютную ошибку сравнений
	* @return - минимальная абсолютная ошибка сравнений
	*/
	RealArgument GetAbsoluteError() const;

	/**
	* Получить минимальную относительную ошибку сравнений
	* @return - минимальная относительная ошибка сравнений
	*/
	RealArgument GetRelativeError() const;


	// Модификаторы

	/**
	* Установить минимальную абсолютную ошибку сравнений
	* @param absoluteError - минимальная абсолютная ошибка сравнений
	* @return признак успешной установки
	*/
	bool SetAbsoluteError
		(
			RealArgument absoluteError
		);

	/**
	* Установить минимальную относительную ошибку сравнений
	* @param relativeError - минимальная относительная ошибка сравнений
	* @return признак успешной установки
	*/
	bool SetRelativeError
		(
			RealArgument relativeError
		);


	// Функции сравнения

	/**
	* Сравнение на равенство с нулем
	* @param number - сравниваемое число
	* @return результат сравнения (true - число равно 0, false - число не равно 0)
	*/
	bool IsZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение на неравенство с нулем
	* @param number - сравниваемое число
	* @return результат сравнения (true - число не равно 0, false - число равно 0)
	*/
	bool IsNotZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение, определяющее, является ли число-аргумент меньшим нуля
	* @param number - сравниваемое число
	* @return результат сравнения (true - число меньше 0, false - число не меньше 0)
	*/
	bool LessThanZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение, определяющее, является ли число-аргумент большим нуля
	* @param number - сравниваемое число
	* @return результат сравнения (true - число больше 0, false - число не больше 0)
	*/
	bool GreaterThanZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение, определяющее, является ли число-аргумент меньшим или равным нулю
	* @param number - сравниваемое число
	* @return результат сравнения (true - число меньше или равно 0, false - число больше 0)
	*/
	bool LessThanOrEqualToZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение, определяющее, является ли число-аргумент большим или равным нулю
	* @param number - сравниваемое число
	* @return результат сравнения (true - число больше или равно 0, false - число меньше 0)
	*/
	bool GreaterThanOrEqualToZero
		(
			RealArgument number
		)	const;

	/**
	* Сравнение на равенство
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - числа равны, false - числа неравны)
	*/
	bool Equal
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Сравнение на неравенство
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - числа неравны, false - числа равны)
	*/
	bool Unequal
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Сравнение на "меньше"
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - первое число меньше второго, false - второе)
	*/
	bool Less
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Сравнение на "больше"
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - первое число больше второго, false - второе)
	*/
	bool Greater
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Сравнение на "меньше или равно"
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - первое число меньше или равно второму, false - второе)
	*/
	bool LessOrEqual
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Сравнение на "больше или равно"
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @return результат сравнения (true - первое число больше или равно второму, false - второе)
	*/
	bool GreaterOrEqual
		(
			RealArgument number1,
			RealArgument number2
		)	const;

	/**
	* Обертка над функциями сравнения, использующая теги для определения типа сравнения
	* @param number1 - первое сравниваемое число
	* @param number2 - второе сравниваемое число
	* @param comparisonTag - тег, используемый для определения типа сравнения
	* @return результат сравнения (true - первое число больше или равно второму, false - второе)
	*/
	bool operator()
		(
			RealArgument number1,
			RealArgument number2,
			ComparisonTag comparisonTag = ComparisonTag::Less
		)	const;


protected:

	// внутренние типы

	// текущий тип сравнивателя вещественных чисел
	typedef RealComparer<Real, UseReferences> ThisType;

	// тип функции-члена, предназначенной для сравнения вещественных чисел
	typedef
	bool (ThisType::*ComparisonMember)
		(
			RealArgument number1,
			RealArgument number2
		)	const;


	// статические переменные члены

	// массив функций-членов сравнения вещественных чисел
	static
	const ComparisonMember _comparisonMembers[];


	// переменные-члены

	// минимальная абсолютная ошибка сравнений
	Real _absoluteError;

	// минимальная относительная ошибка сравнений
	Real _relativeError;

	// граничное значение при переходе
	// от сравнения с использованием абсолютной ошибки
	// к сравнению с использованием относительной ошибки
	Real _boundaryValue;


	// Вспомогательные функции

	/**
	* Заполнить значения минимальных абсолютной и относительной ошибок сравнения
	* данными класса RealErrors
	*/
	void FillErrors();

	/**
	* Заполнить граничное значение при переходе от одной ошибки к другой
	*/
	void FillBoundaryValue();

	/**
	* Установить минимальную ошибку сравнений
	* @param error - минимальная ошибка сравнений
	* @param errorValue - новое значение ошибки сравнений
	* @return признак успешной установки
	*/
	bool SetError
		(
			Real& error,
			RealArgument errorValue
		);


};


// класс, предназначенный для сравнения вещественных чисел типа float
typedef RealComparer<float, false> FloatComparer;

// класс, предназначенный для сравнения вещественных чисел типа double
typedef RealComparer<double, false> DoubleComparer;

// класс, предназначенный для сравнения вещественных чисел типа long double
typedef RealComparer<double, true> LongDoubleComparer;


#endif