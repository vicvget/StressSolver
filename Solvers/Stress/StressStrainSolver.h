#pragma once
#include <iostream>

using std::cout;
using std::string;

enum DataType
{
	DT_Shifts=0,
	DT_Speeds=1,
	DT_Accelerations=2
};

#define ALIGNED_MEM
namespace stress
{
#define MAX(x, y) ((x) > (y) ? x : y)
	double DMOD(const double d1, const double d2);

//#ifdef ALIGNED_MEM
//	const size_t vecStride = 4;				// смещение векторов
//#else
//	const size_t vecStride = 3;				// смещение векторов
//#endif
//
//	const size_t matStride = vecStride * 3; // смещения матриц поворота (3x3, 3x4)
//	const size_t vecStride2 = vecStride * 2; // смещение векторов неизвестных для узла((3+0)X, (3+0)R)
//	const size_t alignment = 32;
}

#define MAX(x, y) ((x) > (y) ? x : y)
double DMOD(const double d1, const double d2);

//#ifdef ALIGNED_MEM
//const size_t vecStride = 4;				// смещение векторов
//#else
//const size_t vecStride = 3;				// смещение векторов
//#endif

//const size_t matStride = vecStride * 3; // смещения матриц поворота (3x3, 3x4)
//const size_t vecStride2 = vecStride * 2; // смещение векторов неизвестных для узла((3+0)X, (3+0)R)
const size_t alignment = 32;


//typedef T double;
class StressStrainSolver
{
	//size_t _stride;

protected:
	string _fileName;
	bool _readIco;
	bool _writeIco;
	int _nWriteIteration;
public:

	size_t vecStride;
	size_t vecStride2;
	size_t matStride;

	StressStrainSolver
		(
			int const nNodes,
			int stride
		);


	void SetZeroVelocities();

	virtual
	~StressStrainSolver();
	
	void InitIco ( const string& fileName,
		bool readIco,
		bool writeIco,
		int nWriteIteration);

	bool ReadIco(const char* fileName);
	void WriteIco(const char* fileName) const;

	float* GetMemoryPointer() const;

	int GetMemorySize() const;

	/** Возвращает указатель на массив внутренних данных
	* @param dataType - тип данных (перемещения, скорости, ускорения)
	* @return - указатель на массив данных
	*/
	double* GetDataInternal(DataType dataType) const;
	
	/** Возвращает указатель на массив матриц поворота
	* @return - указатель на массив матриц поворота
	*/
	double* GetDataRotaionMtx() const;

	/** Добавление связей
	* @param links - массив 6*_nElements с признаками связей (1 - есть связь, 0 - нет связи)
	*/
	virtual
	void AddLinks
		(
			const int* links
		)
		= 0;

	/** Добавить граничное условие подобласти
	* @param boundaryNodesIndices - индексы граничных точек
	* @param numberOfNodesInPartialBoundary - количество граничных точек в частиной границе
	* @param numberOfNodesInFullBoundary - количество граничных точек в полной границе
	* @param bcKind - тип граничных словий (=1 - сила, равномерно разделенная по точкам, 
	* задается вектором из 6 компонент - 3 поступательных компоненты вектора силы 
	* по координатным осям и три момента относительно осей)
	* @param bcParams - переменные граничных условий
	* =1: 6 параметров - компоненты силы Fx,Fy,Fz,Mx,My,Mz
	*/
	virtual
	void AddPartialBoundary
		(
			int* boundaryNodesIndices, 
			int numberOfNodesInPartialBoundary ,
			int numberOfNodesInFullBoundary ,
			int bcKind,
			double* bcParams
		)
		= 0;

	/** Добавить граничное условие
	* @param boundaryNodesIndices - индексы граничных точек
	* @param numberOfBoundaryNodes - количество граничных точек
	* @param bcKind - тип граничных словий (=1 - сила, равномерно разделенная по точкам, 
	* задается вектором из 6 компонент - 3 поступательных компоненты вектора силы 
	* по координатным осям и три момента относительно осей)
	* @param bcParams - переменные граничных условий
	* =1: 6 параметров - компоненты силы Fx,Fy,Fz,Mx,My,Mz
	*/
	virtual
	void AddBoundary
		(
			int* boundaryNodesIndices, 
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
		)
		= 0;

	/** Изменить парметр граничного условия
	* @param bcNumber - номер граничного условия (в порядке добавления),
	* @param bcParamNumber - номер параметра в параметрах граничных условий
	* @param bcParamValue - новое значение параметра
	*/
	virtual
	void ChangeBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber,
			const double bcParamValue
		)
		= 0;
	
	/** Получить параметр по номеру граничного условия и по индексу
	* @param bcNumber - номер граничного условия
	* @param bcParamNumber - номер параметра
	*/
	virtual
	double GetBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber
		)
		= 0;

	virtual
	void UpdateBuffer
		(
			double scale = 1
		)
		= 0;

	/** Расчет нескольких итераций
	* @param nIteratons - количество итераций
	*/
	virtual
	void Solve
		(
			const int nIteratons
		)
		= 0;

	/**
	* Расчет первой стадии метода Рунге-Кутты
	*/
	virtual
	void Solve1()
		= 0;

	/**
	* Расчет второй стадии метода Рунге-Кутты
	*/
	virtual
	void Solve2()
		= 0;

	/**
	* Расчет третьей стадии метода Рунге-Кутты
	*/
	virtual
	void Solve3()
		= 0;

	/**
	* Расчет четвертой стадии метода Рунге-Кутты
	*/
	virtual
	void Solve4()
		= 0;

	/**
	* Расчет пятой стадии метода Рунге-Кутты
	*/
	virtual
	void Solve5()
		= 0;

	double GetData
		(
			int nNode,
			int dir,
			int type
		);

	/** Получить скалярный параметр
	* @param data - массив для записи скалярного параметра
	*/
	virtual
	void GetScalarParameter
		(
			float* data
		);

	/** Получить смещения
	* @param data - массив для записи смещений как скалярного параметра
	*/
	virtual
	void GetDisplacement
		(
			float* data
		)
		= 0;

	/** Получить напряжения по первой теории прочности
	* @param data - массив для записи напряжений как скалярного параметра
	*/
	virtual
	void GetStressesByFirstTheoryOfStrength
		(
			float* data
		)
		= 0;

	/** Получить напряжения по von Mises
	* @param data - массив для записи напряжений как скалярного параметра
	*/
	virtual
	void GetStressesByVonMises
		(
			float* data
		)
		= 0;


public:

	float*	_data;				// неизвестные для сброса в файл (x,y,z,p)
	double* _dataInternal;		// неизвестные для расчета (X,DX,DDX)
	double* _dataRotationMtx;	// матрицы поворота для расчета

	double* _elements;			// исходные координаты элементов

	int _dataSize;
	double* _stress;			// напряжения
	int _nElements;				// число элементов


	virtual double* GetElementGridCoordinates(int elementId) const;;
	virtual double* GetElementStress(int elementId) const;
	virtual double* GetElementStressAngular(int elementId) const;
	virtual double* GetElementShift(int elementId) const;
	virtual double* GetElementVelocity(int elementId) const;
	virtual double* GetElementAcceleration(int elementId) const;
	virtual double* GetElementShiftAngular(int elementId) const;
	virtual double* GetElementVelocityAngular(int elementId) const;
	virtual double* GetElementAccelerationAngular(int elementId) const;
	virtual double* GetRotationMatrix(int elementId) const;
};

