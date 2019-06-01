#include "MathF.h"

 void matMul3x4(const double * matA, const double * matB, double * res)
{
        __m256d matB0 = _mm256_load_pd(matB);
        __m256d matB1 = _mm256_load_pd(matB+4);
        __m256d matB2 = _mm256_load_pd(matB+8);
        __m256d res0 = _mm256_mul_pd(matB0, _mm256_broadcast_sd(matA));
        res0 = _mm256_add_pd(_mm256_mul_pd(matB1, _mm256_broadcast_sd(matA + 1)), res0);
        res0 = _mm256_add_pd(_mm256_mul_pd(matB2, _mm256_broadcast_sd(matA + 2)), res0);
        _mm256_store_pd(res, res0);

        __m256d res1 = _mm256_mul_pd(matB0, _mm256_broadcast_sd(matA + 4));
        res1 = _mm256_add_pd(_mm256_mul_pd(matB1, _mm256_broadcast_sd(matA + 1 + 4)), res1);
        res1 = _mm256_add_pd(_mm256_mul_pd(matB2, _mm256_broadcast_sd(matA + 2 + 4)), res1);
        _mm256_store_pd(res+4, res1);

        __m256d res2 = _mm256_mul_pd(matB0, _mm256_broadcast_sd(matA + 8));
        res2 = _mm256_add_pd(_mm256_mul_pd(matB1, _mm256_broadcast_sd(matA + 1 + 8)), res2);
        res2 = _mm256_add_pd(_mm256_mul_pd(matB2, _mm256_broadcast_sd(matA + 2 + 8)), res2);
        _mm256_store_pd(res+8, res2);

    //for (int i = 0; i < 3; i++)
    //{
    //    
    //    for (int j = 0; j < 3; j++)
    //    {
    //        res[i * 4 + j] = 0;
    //        for (int k = 0; k < 3; k++)
    //        {
    //            res[i * 4 + j] += matA[i * 4 + k] * matB[k * 4 + j];
    //        }
    //    }
    //}
}
/* C= -A*b - A *d */
 void vecDoubleMulSub(const double * vecA, const double * vecB, const double b, const double d, double *vecC)
{
    for (size_t i = 0; i < 3; i++)
    {
        vecC[i] = -vecA[i] * b - vecB[i] * d;
    }
}

 void matVecMul3x4(const double * matA, const double * vec, double * res)
{
     //__m256d rRes;
     //__m256d rVec = _mm256_load_pd(vec);
     //__m256i indexes = _mm256_set_epi64x(12,8,4,0);
     //__m256i indexes1 = _mm256_set_epi64x(13, 9, 5, 0);
     //__m256i indexes2 = _mm256_set_epi64x(14, 10, 6, 0);
     //rRes = _mm256_mul_pd(_mm256_i64gather_pd(matA, indexes, 4), _mm256_broadcast_sd(vec));
     //rRes = _mm256_add_pd(_mm256_mul_pd(_mm256_i64gather_pd(matA, indexes1, 8), _mm256_broadcast_sd(vec + 1)), rRes);
     //rRes = _mm256_add_pd(_mm256_mul_pd(_mm256_i64gather_pd(matA, indexes2, 8), _mm256_broadcast_sd(vec + 2)), rRes);
     //_mm256_store_pd(res, rRes);
     
    for (int i = 0; i < 3; i++)
    {
        res[i] = 0.0;
        for (int j = 0; j < 3; j++)
        {
            res[i] += matA[i * 4 + j] * vec[j];
        }
    }
}

 void matVecTMul3x4(const double * matA, const double * vec, double * res)
{
     __m256d rRes;
     __m256d rVec=_mm256_load_pd(vec);
     
     rRes = _mm256_mul_pd(_mm256_load_pd(matA), _mm256_broadcast_sd(vec));
     rRes = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(matA+4), _mm256_broadcast_sd(vec+1)), rRes);
     rRes = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(matA + 8), _mm256_broadcast_sd(vec + 2)), rRes);
     _mm256_store_pd(res, rRes);
    //for (int i = 0; i < 3; i++)
    //{
    //    res[i] = 0.0;
    //    for (int j = 0; j < 3; j++)
    //    {
    //        res[i] += matA[j * 4 + i] * vec[j];
    //    }
    //}
}

 void vecCross(const double * vecA, const double * vecB, double * res)
{
    res[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    res[1] = -vecA[0] * vecB[2] + vecA[2] * vecB[0];
    res[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
}
 void sincos(const double* angels, double *vSin, double * vCos)
 {
     for (size_t i = 0; i < 4; i++)
     {
         double x = angels[i];
         vSin[i] = x - x * x*x *(1.0/ 6.0) + x * x*x*x*x *(1.0/ 120.0);
         vCos[i] = 1 - x * x *0.5 + x * x*x*x *(1.0/ 24.0);
     }
 }
 void MakeXYZRotationMtx01(const size_t stride, const double* const angles, double * res)
 {
     double pcos[4];
     double psin[4];
     sincos(angles, psin, pcos);

     res[0 * stride + 0] = pcos[1] * pcos[2];
     res[0 * stride + 1] = -pcos[1] * psin[2];
     res[0 * stride + 2] = psin[1];
     res[1 * stride + 0] = psin[0] * psin[1]*pcos[2]+ pcos[0]* psin[2];
     res[1 * stride + 1] = -psin[0] * psin[1]*psin[2]+ pcos[0]* pcos[2];
     res[1 * stride + 2] = -psin[0] * pcos[1];
     res[2 * stride + 0] = -pcos[0] * psin[1]*pcos[2] + psin[0] * psin[2];
     res[2 * stride + 1] = psin[0] * psin[1]*psin[2]+ psin[0]* pcos[2];
     res[2 * stride + 2] = pcos[0] * pcos[1];
 }

