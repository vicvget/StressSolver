#include "RotationSolver.h"

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
		return _varX + elementId * _vecStride;
	}

	double* RotationSolver::GetAngularVelocity(size_t elementId) const
	{
		return _wPointer + elementId * _vecStride2;
	}

	bool RotationSolver::IsSingularityAngle(size_t elementId) const
	{
		double angle = GetAngles(elementId)[1];
		return std::abs(M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI) < 0.1);

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
		_nVariables = nElements*_vecStride;
		const size_t matSize = nElements*_matStride*sizeof(double);
		const size_t varSize = nElements*_vecStride*sizeof(double);

		_varDX = (double*)_aligned_malloc(varSize, alignment);
		_hDDX1 = (double*)_aligned_malloc(varSize, alignment);
		_hDDX2 = (double*)_aligned_malloc(varSize, alignment);
		_hDDX3 = (double*)_aligned_malloc(varSize, alignment);
		_varX = (double*)_aligned_malloc(varSize, alignment);
		_initX = (double*)_aligned_malloc(varSize, alignment);

		_rframeMtx = (double*)_aligned_malloc(matSize, alignment);
		
		memset(_varX, 0, varSize);
		memset(_initX, 0, varSize);
		memset(_varDX, 0, varSize);
		memset(_hDDX1, 0, varSize);
		memset(_hDDX2, 0, varSize);
		memset(_hDDX3, 0, varSize);

		// если изначально были не нулевые повороты, то будут не единичные матрицы
		memcpy(_rframeMtx, mtxPointer, matSize);
	}

	RotationSolver::~RotationSolver()
	{
		_aligned_free(_varX);
		_aligned_free(_initX);
		_aligned_free(_varDX);
		_aligned_free(_hDDX1);
		_aligned_free(_hDDX2);
		_aligned_free(_hDDX3);
		_aligned_free(_rframeMtx);
	}

	void RotationSolver::MakeZeroVectors(size_t elementId) const
	{
		size_t offset = _vecStride * elementId;
		memset(_varDX + offset, 0, sizeof(double) * _vecStride);
		memset(_hDDX1 + offset, 0, sizeof(double) * _vecStride);
		memset(_hDDX2 + offset, 0, sizeof(double) * _vecStride);
		memset(_hDDX3 + offset, 0, sizeof(double) * _vecStride);
		memset(_varX + offset, 0, sizeof(double) * _vecStride);
		memset(_initX + offset, 0, sizeof(double) * _vecStride);
	}

	void RotationSolver::InitIteration() const
	{
		memcpy(_initX, _varX, sizeof(double)*_nVariables);
	}

	void RotationSolver::Solve1()
	{
		for (size_t i = 0; i < _nVariables; i++)
			_varX[i] = _initX[i] + _hDDX1[i] * 0.5;
		CalculateRHS(); // k1 = f(t,y)
		memcpy(_hDDX1,_varDX,sizeof(double)*_nVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve2()
	{
		for (size_t i = 0; i < _nVariables; i++)
			_varX[i] = _initX[i] + _hDDX1[i] * 0.5;
		CalculateRHS(); // k2 = f(t+h/2,y+k1*h/2)
		memcpy(_hDDX2, _varDX, sizeof(double)*_nVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve3()
	{
		for (size_t i = 0; i < _nVariables; i++)
			_varX[i] = _initX[i] + _hDDX3[i];
		CalculateRHS();// k3 = f(t+h/2,y+k2*h/2)
		memcpy(_hDDX3, _varDX, sizeof(double)*_nVariables);
		UpdateMtxs();
	}

	void RotationSolver::Solve4()
	{
		CalculateRHS();// k4 = f(t+h/2,y+k3*h)
		for (size_t i = 0; i < _nVariables; i++)
			_varX[i] = _initX[i] + (_hDDX1[i] + (_hDDX2[i] + _hDDX3[i]) * 2 + _varDX[i]) / 6.0;
		
		// проверка сингулярности
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

	bool RotationSolver::UpdateRHS(int elementId) const
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


		// Проверки сходимости
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


		_varDX[elementId] =	  (wx*cosZ - wy*sinZ) / cosY*_timeStep;
		_varDX[elementId + 1] =  (wx*sinZ + wy*cosZ)*_timeStep;
		_varDX[elementId + 2] = ((wy*sinZ - wx*cosZ)*tanY + wz)*_timeStep;
		return result;
	}

	void RotationSolver::UpdateMtx(int elementId) const
	{
		if (_vecStride == 4)
		{
			Mat3x4 rframe(GetRframeMtx(elementId));
			(rframe*Mat3x4::MakeXYZRotationMtx01(_varX + elementId))
				.Export(GetRotationMtx(elementId));
		}
		else
		{
			Mat3 rframe(GetRframeMtx(elementId));
			(rframe*Mat3::MakeXYZRotationMtx01(_varX + elementId))
				.Export(GetRotationMtx(elementId));
		}
	}

	void RotationSolver::UpdateMtxs() const
	{
		for (size_t elementId = 0; elementId < _nElements; elementId++)
			UpdateMtx(elementId);
	}
}
