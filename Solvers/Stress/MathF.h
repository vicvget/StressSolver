#pragma once

void matMul3x4(const double * matA, const double * matB, double * res);

void vecDoubleMulSub(const double * vecA, const double * vecB, const double b, const double d, double *vecC);
void matVecMul3x4(const double * matA, const double * vec, double * res);
void matVecTMul3x4(const double * matA, const double * vec, double * res);
void vecCross(const double * vecA, const double * vecB, double * res);
