#pragma once

//#define _USE_MATH_DEFINES
//#include <cmath>
#include <vector>

//#include "StressStrainSolver.h"


using std::string;
using std::vector;

namespace Stress
{

	class RotationSolver
	{
		double* _varR;	// текущие переменные
		double* _varDR;	// текущие производные по времени
		double* _initR;	// значения переменных предыдущего шага
		double* _hDR1;	// вспомогательные переменные (RK4)
		double* _hDR2; // вспомогательные переменные (RK4)
		double* _hDR3; // вспомогательные переменные (RK4)

// Кэш синусов и косинусов для вычисления через регистры
#ifdef USE_SVML
		double* _sincache;
		double* _coscache;
#endif

		// матрица для изменения системы кординат отсчета углов при приближении к сингулярности
		// для выбранной системы углов (x,y,z) при abs(y % pi) близком к pi/2 матрица 
		// присваивается текущей матрице поворота, а углы обнуляются
		double* _rframeMtx;	

		double* _wPointer;
		double* _mtxPointer;
		double _timeStep;

		size_t _vecStride; 
		size_t _vecStride2;
		size_t _matStride;

		bool _isValid;
		int _maxRegSize;
	public:

		bool IsSingularityAngle(size_t elementId) const;

		bool IsValid() const;

		void Update(size_t elementId, int stageRK = 0);
		void MakeZeroVectors(size_t elementId) const;

		void InitIteration() const;
		void InitialSolve();

		void ReadIco(std::ifstream& ifs);
		void WriteIco(std::ofstream& ofs);

		void Solve1();
		void Solve2();
		void Solve3();
		void Solve4();

		RotationSolver
			(
				int nElements, 
				int stride,
				double timeStep,
				double* wPointer,
				double* mtxPointer
			);
		~RotationSolver();

	protected:
		size_t _nElements;
		size_t _nRVariables;
		void CalculateRHS();		
		
		// обновление правых частей
		bool UpdateRHS(size_t elementId) const;
		
		// обновление матриц
		void UpdateMtx(size_t elementId) const;
		void UpdateMtxs() const;

		//void UpdateR2(const int offet, const int stageRK) const;

		// Методы для получения характеристик элемента
		double* GetRframeMtx(size_t elementId) const;
		double* GetRotationMtx(size_t elementId) const;
		double* GetAngles(size_t elementId) const;
		double* GetDerivatives(size_t elementId) const;
		double* GetAngularVelocity(size_t elementId) const;

#ifdef USE_SVML
		void FillSinCosCaches();
		double* GetCos(size_t elementId) const;
		double* GetSin(size_t elementId) const;
#endif
	};

};