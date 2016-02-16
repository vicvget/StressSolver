#ifndef CsrSymmetricMatrixH
#define CsrSymmetricMatrixH

struct CsrSymmetricMatrix
{
	int* ia;
	int* ja;
	double* a;
	double* rhs;
	double* x;
	int _size;
	int _memsize;

	void Sort();
	void RemoveNearestToZero();
	void Print() const;
	void PrintPortrait() const;
	void Allocate();
	void Allocate1d(int nNodes);
	void Allocate3d(int nNodes);
	void ApplyRhs(int varId, double value);
	void Fix3dNode(int nodeNumber);
	void Fix3dNode(int nodeNumber, int direction);
	void Fix3dNodeRotation(int nodeNumber);
	void ApplyRhsTo3dNode(int nodeNumber, int direction, double value);
	void FixVar(int varId);
	CsrSymmetricMatrix();
	CsrSymmetricMatrix(int size);
	~CsrSymmetricMatrix();

private:

	// запрещаем использование конструктор копирования
	CsrSymmetricMatrix(const CsrSymmetricMatrix&);

};
#endif