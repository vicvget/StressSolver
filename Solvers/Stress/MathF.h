#pragma once
#include "immintrin.h"
#include "math.h"
void matMul3x4(const double * matA, const double * matB, double * res);

void vecDoubleMulSub(const double * vecA, const double * vecB, const double b, const double d, double *vecC);
void matVecMul3x4(const double * matA, const double * vec, double * res);
void matVecTMul3x4(const double * matA, const double * vec, double * res);
void vecCross(const double * vecA, const double * vecB, double * res);
void MakeXYZRotationMtx01(const size_t stride, const double* const angles, double * res);
