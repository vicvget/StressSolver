// Функции кортежа матриц жесткости


// Заполнить элемент матрицы жесткости
STIFFNESS_MATRICES_FUNCTION
	(
		FillElement,
		void,
		size_t,
		size_t,
		size_t,
		size_t,
		double
	)

// Записать матрицу жесткости в текстовый файл
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToTextFile,
		void,
		const std::string&
	)

// Записать матрицу жесткости в двоичный файл
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToBinaryFile,
		void,
		const std::string&
	)

//  Записать матрицу жесткости в файл
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToFile,
		void,
		const std::string&,
		StiffnessMatrixFileFormats
	)