#include <cstring>

#include "RotationSolver.h"
#include "Common.h"

#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"


using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#define M_PI_2 M_PI*0.5
#endif


namespace Stress
{
	void PrintVector(Vec3 vec, string comment);
	void PrintMatrix(Mat3 mat, string comment);
	void PrintMatrix(MathHelpers::Mat3x4 mat, string comment);

	double* RotationSolver::GetRframeMtx(size_t elementId) const
	{
		return _rframeMtx + _matStride * elementId;
	}

	double* RotationSolver::GetRotationMtx(size_t elementId) const
	{
		return _mtxPointer + _matStride * elementId;
	}

	double* RotationSolver::GetAngles(size_t elementId) const
	{
		return _varR + elementId * _vecStride;
	}

	double* RotationSolver::GetDerivatives(size_t elementId) const
	{
		return _varDR + elementId * _vecStride;
	}

	double* RotationSolver::GetAngularVelocity(size_t elementId) const
	{
		return _wPointer + elementId * _vecStride2;
	}

	bool RotationSolver::IsSingularityAngle(size_t elementId) const
	{
		double angle = GetAngles(elementId)[1];
		return std::abs(M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI)) < 0.1;

	}

	bool RotationSolver::IsValid() const
	{
		return _isValid;
	}

	RotationSolver::RotationSolver
		(
			int nElements,
			int stride,
			double timeStep,
			double* wPointer,
			double* mtxPointer
		) :
		_timeStep(timeStep),
		_wPointer(wPointer),
		_mtxPointer(mtxPointer),
		_nElements(nElements),
		_vecStride(stride),
		_vecStride2(stride * 2),
		_matStride(stride * 3),
		_isValid(true)
	{
		_nRVariables = nElements*_vecStride;
		
		const size_t maxRegSize = 8; // KNC - 8 doubles
		size_t matCount = ((nElements*_matStride + maxRegSize - 1) / maxRegSize) * maxRegSize;
		size_t varCount = ((nElements*_vecStride + maxRegSize - 1) / maxRegSize) * maxRegSize;
		size_t matSize = matCount*sizeof(double);
		size_t varSize = varCount*sizeof(double);

		_varDR = (double*)aligned_alloc(varSize, ALIGNMENT);
		_hDR1 = (double*)aligned_alloc(varSize, ALIGNMENT);
		_hDR2 = (double*)aligned_alloc(varSize, ALIGNMENT);
		_hDR3 = (double*)aligned_alloc(varSize, ALIGNMENT);
		_varR = (double*)aligned_alloc(varSize, ALIGNMENT);
		_initR = (double*)aligned_alloc(varSize, ALIGNMENT);

		_rframeMtx = (double*)aligned_alloc(matSize, ALIGNMENT);
		
		memset(_varR, 0, varSize);
		memset(_initR, 0, varSize);
		memset(_varDR, 0, varSize);
		memset(_hDR1, 0, varSize);
		memset(_hDR2, 0, varSize);
		memset(_hDR3, 0, varSize);

		// ���� ���������� ���� �� ������� ��������, �� ����� �� ��������� �������
		memcpy(_rframeMtx, mtxPointer, matSize);
	}

	RotationSolver::~RotationSolver()
	{
		aligned_free(_varR);
		aligned_free(_initR);
		aligned_free(_varDR);
		aligned_free(_hDR1);
		aligned_free(_hDR2);
		aligned_free(_hDR3);
		aligned_free(_rframeMtx);
	}

	void RotationSolver::MakeZeroVectors(size_t elementId) const
	{
		size_t offset = _vecStride * elementId;
		memset(_varDR + offset, 0, sizeof(double) * _vecStride);
		memset(_hDR1 + offset, 0, sizeof(double) * _vecStride);
		memset(_hDR2 + offset, 0, sizeof(double) * _vecStride);
		memset(_hDR3 + offset, 0, sizeof(double) * _vecStride);
		memset(_varR + offset, 0, sizeof(double) * _vecStride);
		memset(_initR + offset, 0, sizeof(double) * _vecStride);
	}

// corrected
	//void RotationSolver::InitIteration() const
	//{
	//	memcpy(_initR, _varR, sizeof(double)*_nRVariables);

	//	// _hDR1 = k1*h = _varDR
	//	memcpy(_hDR1, _varDR, sizeof(double)*_nRVariables);
	//}

	//void RotationSolver::Solve1()
	//{
	//	for (size_t i = 0; i < _nRVariables; i++)
	//		_varR[i] = _initR[i] + _hDR1[i] * 0.5;
	//	CalculateRHS(); // k2 = f(t+h/2,y+k1*h/2)
	//	// _hDR2 = k2*h
	//	memcpy(_hDR2, _varDR, sizeof(double)*_nRVariables);
	//	UpdateMtxs();
	//}

	//void RotationSolver::Solve2()
	//{
	//	for (size_t i = 0; i < _nRVariables; i++)
	//		_varR[i] = _initR[i] + _hDR2[i] * 0.5;
	//	CalculateRHS(); // k3 = f(t+h/2,y+k2*h/2)
	//	// _hDR3 = k3*h
	//	memcpy(_hDR3, _varDR, sizeof(double)*_nRVariables);
	//	UpdateMtxs();
	//}

	//void RotationSolver::Solve3()
	//{

	//	for (size_t i = 0; i < _nRVariables; i++)
	//		_varR[i] = _initR[i] + _hDR3[i];
	//	CalculateRHS();// k3 = f(t+h/2,y+k2*h/2)
	//	memcpy(_hDR3, _varDR, sizeof(double)*_nRVariables);
	//	UpdateMtxs();
	//}

	//void RotationSolver::Solve4()
	//{
	//	CalculateRHS();// k4 = f(t+h/2,y+k3*h)
	//	for (size_t i = 0; i < _nRVariables; i++)
	//	{
	//		//_varR[i] = _initR[i] + (_hDR1[i] + 2 * (_hDR2[i] + _hDR3[i]) + _varDR[i]) / 6.0;
	//		_varR[i] = _initR[i] + (_hDR1[i] + (_hDR2[i] + _hDR3[i])) / 6.0;
	//	}
	 
	//	// �������� �������������
	//	for (size_t _elementId = 0; _elementId < _nElements; _elementId++)
	//	{
	//		double* elementMtx = GetRotationMtx(_elementId);
	//		double* rframeMtx = GetRframeMtx(_elementId);
	//		if (IsSingularityAngle(_elementId))
	//		{
	//			MakeZeroVectors(_elementId);
	//			memcpy(rframeMtx, elementMtx, _matStride*sizeof(double));
	//		}
	//		UpdateMtx(_elementId);
	//	}
	//}

	void RotationSolver::InitialSolve()
	{
		CalculateRHS();
		UpdateMtxs();
	}


	void RotationSolver::InitIteration() const
	{
		memcpy(_initR, _varR, sizeof(double)*_nRVariables);
		// _hDR1 = k1*h = _varDR
		memcpy(_hDR1, _varDR, sizeof(double)*_nRVariables);
	}

	void RotationSolver::Solve1()
	{		
		InitIteration();
		for (size_t i = 0; i < _nRVariables; i++)
			_varR[i] = _initR[i] + _hDR1[i] * 0.5;
		CalculateRHS(); // k2 = f(t+h/2,y+k1*h/2)
		// _hDR2 = k2*h
		memcpy(_hDR2, _varDR, sizeof(double)*_nRVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve2()
	{
		for (size_t i = 0; i < _nRVariables; i++)
			_varR[i] = _initR[i] + _hDR2[i] * 0.5;
		CalculateRHS(); // k3 = f(t+h/2,y+k2*h/2)
		// _hDR3 = k3*h
		memcpy(_hDR3, _varDR, sizeof(double)*_nRVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve3()
	{
		
		for (size_t i = 0; i < _nRVariables; i++)
			_varR[i] = _initR[i] + _hDR3[i];
		CalculateRHS();// k4 = f(t+h,y+k3*h)
		//memcpy(_hDR4, _varDR, sizeof(double)*_nRVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve4()
	{
		//CalculateRHS();// k4 = f(t+h/2,y+k3*h)
		for (size_t i = 0; i < _nRVariables; i++)
		{
			_varR[i] = _initR[i] + (_hDR1[i] + 2 * (_hDR2[i] + _hDR3[i]) + _varDR[i]) / 6.0;
			//_varR[i] = _initR[i] + (_hDR1[i] + (_hDR2[i] + _hDR3[i])) / 6.0;
		}
		CalculateRHS(); // calculate  _varDR = f(t+h, y(t+h));
		// �������� �������������
		for (size_t elementId = 0; elementId < _nElements; elementId++)
		{
			double* elementMtx = GetRotationMtx(elementId);
			double* rframeMtx = GetRframeMtx(elementId);
			if (IsSingularityAngle(elementId))
			{
				MakeZeroVectors(elementId);
				memcpy(rframeMtx, elementMtx, _matStride*sizeof(double));
			}
			UpdateMtx(elementId);
		}
	}

	void RotationSolver::CalculateRHS()
	{
		for (size_t elementId = 0; elementId < _nElements; elementId++)
		{
			if (!UpdateRHS(elementId))
			{
				_isValid = false;
				break;
			}
		}
	}

	bool RotationSolver::UpdateRHS(size_t elementId) const
	{
		double* angles = GetAngles(elementId);
		double* elementW = GetAngularVelocity(elementId);

		double cosY = cos(angles[1]);
		double tanY = tan(angles[1]);
		double sinZ = sin(angles[2]);
		double cosZ = cos(angles[2]);
		double wx = elementW[0];
		double wy = elementW[1];
		double wz = elementW[2];


		// �������� ����������
		double controlAngleX = (-sinZ*wx - cosZ*wy) / cosZ;
		double controlAngleY = (wy*cosZ + wx*sinZ)*tanY;
		double maxValue = 400.0 / (_timeStep*_timeStep);
		bool result = true;
		if (std::abs(controlAngleX) > maxValue)
		{
			std::cout << "small step for 1 euler eq\n";
			result = false;
		}
		if (std::abs(controlAngleY) > maxValue)
		{
			std::cout << "small step for 3 euler eq\n";
			result = false;
		}

		double* derivatives = GetDerivatives(elementId);

		derivatives[0] = (wx*cosZ - wy*sinZ) / cosY*_timeStep;
		derivatives[1] = (wx*sinZ + wy*cosZ)*_timeStep;
		derivatives[2] = ((wy*sinZ - wx*cosZ)*tanY + wz)*_timeStep;
		return result;
	}

	void RotationSolver::UpdateMtx(size_t elementId) const
	{
		double* rotationMtx = GetRotationMtx(elementId);
		double* rframeMtx = GetRframeMtx(elementId);
		double* angles = GetAngles(elementId);
		if (_vecStride == 4)
		{
			Mat3x4 rframe(rframeMtx);
			Mat3x4 newMtx = Mat3x4::MakeXYZRotationMtx01(GetAngles(elementId));
			(rframe*newMtx).Export(rotationMtx);
		}
		else
		{
			Mat3 rframe(rframeMtx);
			Mat3 newMtx = Mat3::MakeXYZRotationMtx01(GetAngles(elementId));
			(rframe*newMtx).Export(rotationMtx);
		}
	}

	void RotationSolver::UpdateMtxs() const
	{
		for (size_t elementId = 0; elementId < _nElements; elementId++)
			UpdateMtx(elementId);
	}
}
