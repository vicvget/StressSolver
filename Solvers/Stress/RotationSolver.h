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
		double* _varR;	// ������� ����������
		double* _varDR;	// ������� ����������� �� �������
		double* _initR;	// �������� ���������� ����������� ����
		double* _hDR1;	// ��������������� ���������� (RK4)
		double* _hDR2; // ��������������� ���������� (RK4)
		double* _hDR3; // ��������������� ���������� (RK4)

// ��� ������� � ��������� ��� ���������� ����� ��������
#ifdef USE_SVML
		double* _sincache;
		double* _coscache;
#endif

		// ������� ��� ��������� ������� �������� ������� ����� ��� ����������� � �������������
		// ��� ��������� ������� ����� (x,y,z) ��� abs(y % pi) ������� � pi/2 ������� 
		// ������������� ������� ������� ��������, � ���� ����������
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
		
		// ���������� ������ ������
		bool UpdateRHS(size_t elementId) const;
		
		// ���������� ������
		void UpdateMtx(size_t elementId) const;
		void UpdateMtxs() const;

		//void UpdateR2(const int offet, const int stageRK) const;

		// ������ ��� ��������� ������������� ��������
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