#pragma once

#ifndef STANDALONE
#define STANDALONE
#endif

#ifdef STANDALONE
#define EXTERN
#else
#define EXTERN extern "C"
#endif


#ifndef STANDALONE
#ifdef _USRDLL
#define DLL_FUNCTION EXTERN __declspec(dllexport)
#else
#define DLL_FUNCTION EXTERN __declspec(dllimport)
#endif
#else
#define DLL_FUNCTION
#endif
#include <string>


namespace Stress
{

	/** Инициализация
	* @param params - параметры
	* @param links - связи
	* @param nLinks - количество связей
	* @param nodes - координаты узлов (тройки чисел)
	* @param nNodes - количество узлов
	* @param solverType - тип решателя 0-iterative, 1-statics
	* @return указатель на решатель
	*/
	DLL_FUNCTION
	void* Init
		(
			const std::string& solverUid,
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
			double gridStep,
			double timeStep,
			int procType = 0,
			int numThreads = 0,
			bool useFloat = false,
			int solverType = 0
		);

	DLL_FUNCTION
	void InitIco
		(
			void* hsolver,
			const char* fileName,
			bool readIco,
			bool writeIco
		);

	DLL_FUNCTION
	bool ReadIco
		(
			void* hsolver,
			const char* fileName
		);

	DLL_FUNCTION
	void WriteIco
		(
			void* hsolver,
			const char* fileName
		);

	/** Инициализация без связей
	* @param params - параметры
	* @param nNodes - количество узлов
	* @param nodes - координаты узлов (тройки чисел)
	*
	* @return указатель на решатель
	*/
	/*
	DLL_FUNCTION
	void* Init
	(
	double* params,
	float *nodes,
	int nNodes,
	double gridStep,
	double timeStep,
	int procType = 0,
	bool useFloat = false
	);
	*/

	/** Удалить решатель
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void ReleaseMemory
		(
			void* &hsolver
		);

	/** Получить адрес указателя на память
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	float* GetMemoryPointer
		(
			const void* hsolver
		);

	/** Получить размер памятм
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	int GetMemorySize
		(
			const void* hsolver
		);

	/** Получить адрес указателя на память
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void* GetAccelerations(const void* hsolver);

	/** Добавление связей
	* @param hsolver - дескриптор решателя
	* @param links - массив 6*_nElements с признаками связей (1 - есть связь, 0 - нет связи)
	*/
	DLL_FUNCTION
	void AddLinks
		(
			const void* hsolver,
			int* links
		);

	/** Добавить граничное условие
	* @param hsolver - дескриптор решателя
	* @param boundaryNodesIndices - индексы граничных точек
	* @param numberOfBoundaryNodes - количество граничных точек
	* @param bcKind - тип граничных словий (=1 - сила, равномерно разделенная по точкам, 
	* задается вектором из 6 компонент - 3 поступательных компоненты вектора силы 
	* по координатным осям и три момента относительно осей)
	* @param bcParams - переменные граничных условий
	* =1: 6 параметров - компоненты силы Fx,Fy,Fz,Mx,My,Mz
	*/
	DLL_FUNCTION
	void AddBoundary
		(
			const void *hsolver,
			int* boundaryNodesIndices, 
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
		);

	/** Изменить парметр граничного условия
	* @param hsolver - дескриптор решателя
	* @param bcNumber - номер граничного условия (в порядке добавления),
	* @param bcParamNumber - номер параметра в параметрах граничных условий
	* @param bcParamValue - новое значение параметра
	*/
	DLL_FUNCTION
	void ChangeBoundaryParam
		(
			const void *hsolver,
			int bcNumber,
			int bcParamNumber,
			double bcParamValue
		);

	/** Получить параметр по номеру граничного условия и по индексу
	* @param hsolver - дескриптор решателя
	* @param bcNumber - номер граничного условия
	* @param bcParamNumber - номер параметра
	*/
	DLL_FUNCTION
	double GetBoundaryParam
		(
			const void *hsolver,
			int bcNumber,
			int bcParamNumber
		);

	/**
	* Обновить буфер
	* @param hsolver - дескриптор решателя
	* @param scaleFactor - масштабный коэффициент
	*/
	DLL_FUNCTION
	float UpdateBuffer
		(
			const void* hsolver,
			double scaleFactor = 1.0
		);

	DLL_FUNCTION
		void InitialSolve
		(
		void* hsolver
		);

	/** Расчет
	* @param hsolver - дескриптор решателя
	* @param nIterations - количество итераций
	*/
	DLL_FUNCTION
	void Solve
		(
			void* hsolver,
			int nIterations
		);

	/**
	* Расчет первой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void Solve1
		(
			void* hsolver
		);

	/**
	* Расчет второй стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void Solve2
		(
			void* hsolver
		);

	/**
	* Расчет третьей стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void Solve3
		(
			void* hsolver
		);

	/**
	* Расчет четвертой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void Solve4
		(
			void* hsolver
		);

	/**
	* Расчет пятой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
	void Solve5
		(
			void* hsolver
		);

	DLL_FUNCTION
	double GetData
		(
			const void* hsolver,
			int nNode,
			int dir,
			int type
		);

	//DLL_FUNCTION
	//void SolveSystemOfLinearEquationsForStiffness();

}

