#pragma once
#include "IMemory.h"

/** Интерфейсный класс для решателя НДС
*/
class IStressStrainSolver: public IMemory
{
public:
	virtual ~IStressStrainSolver() = 0;

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

	/** Обновление выходных данных с заданным масштабированием смещений
	* @param scale - масштаб смещений
	*/
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


protected:

	float* _data;
	double* _dataInternal;
	int _dataSize;
	double* _stress;
	int _nNodes;

};

#endif