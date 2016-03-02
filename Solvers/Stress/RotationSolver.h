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
		double* _varDX;
		double* _hDDX1;
		double* _hDDX2;
		double* _hDDX3;
		double* _varX;
		double* _initX;

		double* _rframeMtx;

		size_t _vecStride; 
		size_t _vecStride2;
		size_t _matStride;

	public:

		double* GetRframeMtx(size_t elementId) const;
		double* GetAngles(size_t elementId) const;

		bool IsSingularityAngle(size_t elementId) const;

		void Update(size_t elementId, double* elementW, double* elementMtx, double timeStep, const int stageRK = 0);
		void MakeZeroVectors(size_t elementId) const;

		RotationSolver(const int nElements, int stride);
		~RotationSolver();

	protected:
		void UpdateR(const int offset, const double* elementW, const double timeStep) const;
		void UpdateR2(const int offet, const int stageRK) const;
		void UpdateMtx(const int ofset, double* elementMtx) const;


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