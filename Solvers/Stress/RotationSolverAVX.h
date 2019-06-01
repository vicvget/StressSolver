#pragma once
#include "RotationSolver.h"
#include "MathF.h"
namespace Stress
{
    class RotationSolverAVX :
        public RotationSolver
    {
    public:
        RotationSolverAVX(
            int nElements,
            int stride,
            double timeStep,
            double* wPointer,
            double* mtxPointer
        );
        virtual ~RotationSolverAVX();

        void RotationSolverAVX::CalculateRHS();
        void RotationSolverAVX::UpdateMtxs() const;
        void InitIteration() const;
        void InitialSolve();
        void Solve1();
        void Solve2();
        void Solve3();
        void Solve4();
    };
}

