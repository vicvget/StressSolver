

 void matMul3x4(const double * matA, const double * matB, double * res)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res[i * 4 + j] = 0;
            for (int k = 0; k < 3; k++)
            {
                res[i * 4 + j] += matA[i * 4 + k] * matB[k * 4 + j];
            }
        }
    }
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
    for (int i = 0; i < 3; i++)
    {
        res[i] = 0.0;
        for (int j = 0; j < 3; j++)
        {
            res[i] += matA[j * 4 + i] * vec[j];
        }
    }
}

 void vecCross(const double * vecA, const double * vecB, double * res)
{
    res[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    res[1] = -vecA[0] * vecB[2] + vecA[2] * vecB[0];
    res[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
}