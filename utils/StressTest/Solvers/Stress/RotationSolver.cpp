#include "RotationSolver.h"

#include "CsrSymmetricMatrix.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"

#include <omp.h>
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

	double* RotationSolver::GetAngles(size_t elementId) const
	{
		return _varX + elementId * _vecStride;
	}

	bool RotationSolver::IsSingularityAngle(size_t elementId) const
	{
		double angle = GetAngles(elementId)[1];
		return std::abs(M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI) < 0.1);

	}

	RotationSolver::RotationSolver(const int nElements, int stride) :
		_vecStride(stride),
		_vecStride2(stride * 2),
		_matStride(stride * 3)
	{
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

		memset(_rframeMtx, 0, matSize);
		for (int i = 0; i < nElements; i++)
		{
			_rframeMtx[i*_matStride] = 1;
			_rframeMtx[i*_matStride + _vecStride + 1] = 1;
			_rframeMtx[i*_matStride + _vecStride2 + 2] = 1;
		}
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

	void RotationSolver::Update(size_t elementId, double* elementW, double* elementMtx, double timeStep, const int stageRK)
	{

		size_t offset = _vecStride * elementId;
		for (int i = offset; i < offset + _vecStride; i++)
		{
			switch (stageRK)
			{
			case 0:
				MakeZeroVectors(elementId);
				break;
			case 1:
				_varX[i] = _initX[i] + (_hDDX1[i] + (_hDDX2[i] + _hDDX3[i]) * 2 + _varDX[i]) / 6.0;
				_initX[i] = _varX[i];
				break;
			case 2:
			case 3:
				_varX[i] = _initX[i] + _hDDX1[i] * 0.5;
				break;
			case 4:
				_varX[i] = _initX[i] + _hDDX3[i];
				break;
			}
		}
		UpdateR(offset, elementW, timeStep);
		UpdateR2(offset, stageRK);
		UpdateMtx(offset, elementMtx);
	}

	void RotationSolver::UpdateR(const int id, const double* w, const double timeStep) const
	{
		double yc = cos(_varX[id + 1]);
		double yt = tan(_varX[id + 1]);
		double zs = sin(_varX[id + 2]);
		double zc = cos(_varX[id + 2]);
		double wx = w[0];
		double wy = w[1];
		double wz = w[2];


		// Проверки сходимости
		double controlAngleX = (-zs*wx - zc*wy) / zc;
		double controlAngleY = (wy*zc + wx*zs)*yt;
		double maxValue = 400.0 / (timeStep*timeStep);

		if (std::abs(controlAngleX) > maxValue)
		{
			std::cout << "small step for 1 euler eq\n";
		}
		if (std::abs(controlAngleY) > maxValue)
		{
			std::cout << "small step for 3 euler eq\n";
		}


		_varDX[id] = (wx*zc - wy*zs) / yc*timeStep;
		_varDX[id + 1] = (wx*zs + wy*zc)*timeStep;
		_varDX[id + 2] = ((wy*zs - wx*zc)*yt + wz)*timeStep;
	}

	void RotationSolver::UpdateR2(const int id, const int mets) const
	{

		for (int j = id; j < id + 3; j++)
		{
			switch (mets)
			{
			case 1:
				_hDDX1[j] = _varDX[j];
				break;
			case 2:
				_hDDX2[j] = _varDX[j];
				break;
			case 3:
				_hDDX3[j] = _varDX[j];
				break;
			}
		}
	}

	void RotationSolver::UpdateMtx(const int id, double* a) const
	{
		double xc = cos(_varX[id]);
		double yc = cos(_varX[id + 1]);
		double zc = cos(_varX[id + 2]);
		double xs = sin(_varX[id]);
		double ys = sin(_varX[id + 1]);
		double zs = sin(_varX[id + 2]);


		double* firstRow = a;
		double* secondRow = a + _vecStride;
		double* thirdRow = a + _vecStride * 2;

		firstRow[0] = yc*zc;
		firstRow[1] = -yc*zs;
		firstRow[2] = ys;

		secondRow[0] = xs*ys*zc + xc*zs;
		secondRow[1] = -xs*ys*zs + xc*zc;
		secondRow[2] = -xs*yc;

		thirdRow[0] = -xc*ys*zc + xs*zs;
		thirdRow[1] = xc*ys*zs + xs*zc;
		thirdRow[2] = xc*yc;
	}

}
