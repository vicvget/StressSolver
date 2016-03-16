#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "StressStrainSolver.h"


using std::string;
using std::vector;

namespace Stress
{

	class RotationSolver
	{
		double* _varX;	// текущие переменные
		double* _varDX;	// текущие производные по времени
		double* _initX;	// значения переменных предыдущего шага
		double* _hDDX1;	// вспомогательные переменные (RK4)
		double* _hDDX2; // вспомогательные переменные (RK4)
		double* _hDDX3; // вспомогательные переменные (RK4)

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

	public:

		double* GetRframeMtx(size_t elementId) const;
		double* GetRotationMtx(size_t elementId) const;
		double* GetAngles(size_t elementId) const;
		double* GetAngularVelocity(size_t elementId) const;

		bool IsSingularityAngle(size_t elementId) const;

		bool IsValid() const;

		void Update(size_t elementId, int stageRK = 0);
		void MakeZeroVectors(size_t elementId) const;

		void InitIteration() const;
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
		size_t _nVariables;
		void CalculateRHS();		
		bool UpdateRHS(int elementId) const;
		void UpdateMtx(int elementId) const;
		void UpdateMtxs() const;

		void UpdateR2(const int offet, const int stageRK) const;

	};


	//struct Copym
	//{
	//	bool* isFirstCopy;
	//	//double T;
	//	//COMMON/UGR/ KGR1,KGR2,T
	//	double* TM;
	//	double* AZ;
	//	double T0;

	//	int vecStride, vecStride2, matStride;
	//	Copym(const int nNodes, int stride) :
	//		vecStride(stride),
	//		vecStride2(stride * 2),
	//		matStride(stride * 3)
	//	{
	//		isFirstCopy = new bool[nNodes];
	//		for (int i = 0; i < nNodes; i++)
	//			isFirstCopy[i] = true;
	//		TM = new double[nNodes];
	//		AZ = new double[nNodes*matStride];
	//		T0 = 0;
	//	}

	//	~Copym()
	//	{
	//		delete[] isFirstCopy;
	//		delete[] TM;
	//		delete[] AZ;
	//	}

	//	void Copy(double* A, double* A1, const int N, const double T)
	//	{
	//		if (isFirstCopy[N] || T == T0)
	//		{
	//			memcpy(A1, A, matStride*sizeof(double));
	//			isFirstCopy[N] = false;
	//			TM[N] = T;
	//		}
	//		else
	//		{
	//			const int KA = N*matStride;
	//			if (TM[N] != T)
	//			{
	//				memcpy(A1, A, matStride*sizeof(double));
	//				TM[N] = T;
	//				memcpy(AZ + KA, A, matStride*sizeof(double));
	//			}
	//			else
	//			{
	//				memcpy(A1, AZ + KA, matStride*sizeof(double));
	//			}
	//		}
	//	}
	//};
};