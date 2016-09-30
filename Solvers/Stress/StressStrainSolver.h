#pragma once
#include <iostream>

using std::cout;
using std::string;

#define SQR(x) ((x) * (x))
#define MAX(x, y) ((x) > (y) ? x : y)

#define MeasuredRun(TIMER, COMMAND) \
	_testTimer.Start(TIMER); \
	COMMAND; \
	_testTimer.Stop(TIMER);

#define ALIGNMENT 64 // KNC

enum DataType
{
	DT_Shifts=0,
	DT_Speeds=1,
	DT_Accelerations=2
};

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
	float UpdateBuffer
		(
			double scale = 1
		)
		= 0;

	/** Расчет начального времени
	*/
	virtual
	void InitialSolve
		(
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

	/** Получить напряжения по X
	* @param data - массив для записи напряжений как скалярного параметра
	*/
	virtual
		void GetStressesX
		(
		float* data
		)
		= 0;

public:

	float*	_data;				// неизвестные для сброса в файл (x,y,z,p)
	double* _dataInternal;		// неизвестные для расчета (X,DX,DDX)
	double* _dataRotationMtx;	// матрицы поворота для расчета
	double* _buffer;			// массив для KNC

	double* _elements;			// исходные координаты элементов

	int _dataSize;
	double* _stress;			// напряжения
	int _nElements;				// число элементов


	// Методы для получения характеристик элемента
	virtual double  GetElementDisplacement(size_t elementId) const;
	virtual double* GetElementGridCoordinates(size_t elementId) const;
	virtual double* GetElementStress				(size_t elementId) const;
	virtual double* GetElementStressAngular			(size_t elementId) const;
	virtual double* GetElementShift					(size_t elementId) const;
	virtual double* GetElementVelocity				(size_t elementId) const;
	virtual double* GetElementAcceleration			(size_t elementId) const;
	virtual double* GetElementShiftAngular			(size_t elementId) const;
	virtual double* GetElementVelocityAngular		(size_t elementId) const;
	virtual double* GetElementAccelerationAngular	(size_t elementId) const;
	virtual double* GetRotationMatrix				(size_t elementId) const;

	string _uid;
	virtual void SetUid(const string& uid);
};

