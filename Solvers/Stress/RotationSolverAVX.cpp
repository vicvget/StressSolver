#include "RotationSolverAVX.h"

namespace Stress
{

    RotationSolverAVX::RotationSolverAVX(
        int nElements,
        int stride,
        double timeStep,
        double* wPointer,
        double* mtxPointer)
        :
    RotationSolver(nElements, stride, timeStep, wPointer, mtxPointer)
    {
        cout << "AVX rotation" << std::endl;
    }


    
        RotationSolverAVX::~RotationSolverAVX()
        {
        }

        void RotationSolverAVX::CalculateRHS()
        {
            for (int elementId = 0; elementId < _nElements; elementId++)
            {
                const double* angles = _varR + elementId * _vecStride;
                const double* elementW = _wPointer + elementId * _vecStride2;

                double cosY = cos(angles[1]);
                double tanY = tan(angles[1]);
                double sinZ = sin(angles[2]);
                double cosZ = cos(angles[2]);

                double wx = elementW[0];
                double wy = elementW[1];
                double wz = elementW[2];


                // ѕроверки сходимости
                double controlAngleX = (-sinZ * wx - cosZ * wy) / cosZ;
                double controlAngleY = (wy*cosZ + wx * sinZ)*tanY;
                double maxValue = 400.0 / (_timeStep*_timeStep);

                if (std::abs(controlAngleX) > maxValue)
                {
                    std::cout << "small step for 1 euler eq\n";
                }
                if (std::abs(controlAngleY) > maxValue)
                {
                    std::cout << "small step for 3 euler eq\n";
                }

                double* derivatives = _varDR + elementId * _vecStride;

                derivatives[0] = (-wy * sinZ + wx * cosZ) / cosY * _timeStep;
                derivatives[1] = (wx*sinZ + wy * cosZ) *_timeStep;
                derivatives[2] = ((wy*sinZ - wx * cosZ) * tanY + wz) * _timeStep;
            }
        }

        void RotationSolverAVX::UpdateMtxs() const
        {
            for (size_t elementId = 0; elementId < _nElements; elementId++)
            {
                double* rotationMtx = _mtxPointer + _matStride * elementId;
                const double* rframeMtx = _rframeMtx + _matStride * elementId;
                const double* angles = _varR + elementId * _vecStride;
                double newMtx[12];
                MakeXYZRotationMtx01(4, angles, newMtx);
                matMul3x4(rframeMtx, newMtx, rotationMtx);

            }

        }

        void RotationSolverAVX::InitialSolve()
        {
            CalculateRHS();
            UpdateMtxs();
        }

        void RotationSolverAVX::InitIteration() const
        {
            memcpy(_initR, _varR, sizeof(double)*_nRVariables);
            // _hDR1 = k1*h = _varDR
            memcpy(_hDR1, _varDR, sizeof(double)*_nRVariables);
        }

        void RotationSolverAVX::Solve1()
        {
            InitIteration();
#ifdef OMP_SOLVE
#pragma omp parallel for
#endif
            for (int i = 0; i < _nRVariables; i++)
                _varR[i] = _initR[i] + _hDR1[i] * 0.5;
            CalculateRHS(); // k2 = f(t+h/2,y+k1*h/2)
            // _hDR2 = k2*h
            memcpy(_hDR2, _varDR, sizeof(double)*_nRVariables);
            UpdateMtxs();
        }

        void RotationSolverAVX::Solve2()
        {
#ifdef OMP_SOLVE
#pragma omp parallel for
#endif
            for (int i = 0; i < _nRVariables; i++)
                _varR[i] = _initR[i] + _hDR2[i] * 0.5;
            CalculateRHS(); // k3 = f(t+h/2,y+k2*h/2)
            // _hDR3 = k3*h
            memcpy(_hDR3, _varDR, sizeof(double)*_nRVariables);
            UpdateMtxs();
        }

        void RotationSolverAVX::Solve3()
        {
#ifdef OMP_SOLVE
#pragma omp parallel for
#endif
            for (int i = 0; i < _nRVariables; i++)
                _varR[i] = _initR[i] + _hDR3[i];
            CalculateRHS();// k4 = f(t+h,y+k3*h)
            //memcpy(_hDR4, _varDR, sizeof(double)*_nRVariables);
            UpdateMtxs();
        }

        void RotationSolverAVX::Solve4()
        {
            //CalculateRHS();// k4 = f(t+h/2,y+k3*h)
#ifdef OMP_SOLVE
#pragma omp parallel for
#endif
            for (int i = 0; i < _nRVariables; i++)
            {
                _varR[i] = _initR[i] + (_hDR1[i] + 2 * (_hDR2[i] + _hDR3[i]) + _varDR[i]) / 6.0;
                //_varR[i] = _initR[i] + (_hDR1[i] + (_hDR2[i] + _hDR3[i])) / 6.0;
            }
            CalculateRHS(); // calculate  _varDR = f(t+h, y(t+h));
            // проверка сингул€рности
            for (size_t elementId = 0; elementId < _nElements; elementId++)
            {
                double* elementMtx = GetRotationMtx(elementId);
                double* rframeMtx = GetRframeMtx(elementId);
                if (IsSingularityAngle(elementId))
                {
                    MakeZeroVectors(elementId);
                    memcpy(rframeMtx, elementMtx, _matStride * sizeof(double));
                }
                UpdateMtx(elementId);
            }
        }
}
