// ������� ������� ������ ���������


// ��������� ������� ������� ���������
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

// �������� ������� ��������� � ��������� ����
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToTextFile,
		void,
		const std::string&
	)

// �������� ������� ��������� � �������� ����
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToBinaryFile,
		void,
		const std::string&
	)

//  �������� ������� ��������� � ����
STIFFNESS_MATRICES_FUNCTION
	(
		WriteToFile,
		void,
		const std::string&,
		StiffnessMatrixFileFormats
	)