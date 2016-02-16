// Функции компаратора


// Инициализация компаратора значениями ошибок по умолчанию
REAL_COMPARING_FUNCTION
	(
		Init,
		void
	)

// Инициализация компаратора
REAL_COMPARING_FUNCTION
	(
		Init,
		bool,
		RealArgument,
		RealArgument
	)

// Получить минимальную абсолютную ошибку сравнений
REAL_COMPARING_FUNCTION
	(
		GetAbsoluteError,
		EXPAND_MULTIPLE(typename RealComparing<Real, UseReferences>::RealArgument)
	)

// Получить минимальную относительную ошибку сравнений
REAL_COMPARING_FUNCTION
	(
		GetRelativeError,
		EXPAND_MULTIPLE(typename RealComparing<Real, UseReferences>::RealArgument)
	)

// Установить минимальную абсолютную ошибку сравнений
REAL_COMPARING_FUNCTION
	(
		SetAbsoluteError,
		bool,
		RealArgument
	)

// Установить минимальную относительную ошибку сравнений
REAL_COMPARING_FUNCTION
	(
		SetRelativeError,
		bool,
		RealArgument
	)

// Сравнение на равенство с нулем
REAL_COMPARING_FUNCTION
	(
		IsZero,
		bool,
		RealArgument
	)

// Сравнение на неравенство с нулем
REAL_COMPARING_FUNCTION
	(
		IsNotZero,
		bool,
		RealArgument
	)

// Сравнение, определяющее, является ли число-аргумент меньшим нуля
REAL_COMPARING_FUNCTION
	(
		LessThanZero,
		bool,
		RealArgument
	)

// Сравнение, определяющее, является ли число-аргумент большим нуля
REAL_COMPARING_FUNCTION
	(
		GreaterThanZero,
		bool,
		RealArgument
	)

// Сравнение, определяющее, является ли число-аргумент меньшим или равным нулю
REAL_COMPARING_FUNCTION
	(
		LessThanOrEqualToZero,
		bool,
		RealArgument
	)

// Сравнение, определяющее, является ли число-аргумент большим или равным нулю
REAL_COMPARING_FUNCTION
	(
		GreaterThanOrEqualToZero,
		bool,
		RealArgument
	)

// Сравнение на равенство
REAL_COMPARING_FUNCTION
	(
		Equal,
		bool,
		RealArgument,
		RealArgument
	)

// Сравнение на неравенство
REAL_COMPARING_FUNCTION
	(
		Unequal,
		bool,
		RealArgument,
		RealArgument
	)

// Сравнение на "больше"
REAL_COMPARING_FUNCTION
	(
		Greater,
		bool,
		RealArgument,
		RealArgument
	)

// Сравнение на "меньше"
REAL_COMPARING_FUNCTION
	(
		Less,
		bool,
		RealArgument,
		RealArgument
	)

// Сравнение на "больше или равно"
REAL_COMPARING_FUNCTION
	(
		GreaterOrEqual,
		bool,
		RealArgument,
		RealArgument
	)

// Сравнение на "меньше или равно"
REAL_COMPARING_FUNCTION
	(
		LessOrEqual,
		bool,
		RealArgument,
		RealArgument
	)