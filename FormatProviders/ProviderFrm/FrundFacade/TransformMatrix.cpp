#include "TransformMatrix.h"
#include "FGeometryPoint.h"

#include "../../../fcore/fcore.h"
#include "../../../fcore/Exceptions/fcExceptions.h"

#include <cmath>


TransformMatrix::TransformMatrix()
{
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			_mtx[i*4+j] = (i == j) ? 1.0 : 0.0;
		}
	}
}

TransformMatrix::TransformMatrix
	(
		const float (&mtx)[12]
	)
{
	int k;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			k = j*3+i+3;
			_mtx[i*4+j] = mtx[k];
		}
		_mtx[i*4+3] = mtx[i]*1000.;
	}
	for(int i = 0; i < 3; i++)
	{
		_mtx[i+12]=0.0;
	}
	_mtx[15]=1.0;
}

void TransformMatrix::FillRepresentation
	(
		float (&mtx)[12]
	)	const
{
	int k = 3;

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			mtx[k++] = static_cast<float>(Element(j, i));
		}
		mtx[i] = static_cast<float>(Element(i, 3) / 1000.0);
	}
}

bool TransformMatrix::operator == (const TransformMatrix& mtx) const
{
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			if(!IS_ZERO(Element(i,j)-mtx.Element(i,j)))
				return false;
		}
	return true;
}

TransformMatrix TransformMatrix::operator * (const TransformMatrix& mtx) const
{
	TransformMatrix res;
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			res[i][j]=0;
			for(int k = 0; k < 4; k++)
			{
				res[i][j]+=Element(i,k)*mtx.Element(k,j);
			}
		}
		return res;
}

void TransformMatrix::SetTranslation(const FGeometryPoint& v)
{
	for(int j = 0; j < 3; j++)
	{
		_mtx[j*4+3]=v.coords[j];
	}
}

void TransformMatrix::SetRotation(const FGeometryPoint& axis, double angle)
{
	TransformMatrix m1;
	TransformMatrix m2;
	double c0 = axis.x/axis.Magnitude();
	double c1 = axis.y/axis.Magnitude();
	double c2 = axis.z/axis.Magnitude();
	m1[0][0]=c0*c0;
	m1[1][1]=c1*c1;
	m1[2][2]=c2*c2;
	m1[0][1]=m1[1][0]=c0*c1;
	m1[0][2]=m1[2][0]=c0*c2;
	m1[1][2]=m1[2][1]=c1*c2;
	m2[0][1]=-c2;
	m2[0][2]=c1;
	m2[1][2]=-c0;
	m2[1][0]=-m2[0][1];
	m2[2][0]=-m2[0][2];
	m2[2][1]=-m2[1][2];
	m2[0][0]=m2[1][1]=m2[2][2]=0.0;

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			m1[i][j]*=(1-cos(angle));
			_mtx[i*4+j]=m1[i][j]+m2[i][j]*sin(angle);
		}
		_mtx[i*4+i]+=cos(angle);
	}
}


// надо помержить эти 2 матрицы в Mat3x3 или Mat4x3б а то уже невозможно с этим работать, каждый нагородил себе по матрице (put the brain inside)

TransformationMatrix::TransformationMatrix()
{
	for(int i=0; i < 3; i++)
		ShiftVector[i] = 0.;
	for(int i=0; i < 9; i++)
	{
		RotationMatrix[i] = 0.;
		SomeAngles[i] = 0;
	}
	RotationMatrix[0]=1.;
	RotationMatrix[4]=1.;
	RotationMatrix[8]=1.;
}

TransformationMatrix::TransformationMatrix
	(
		const float (&mtx)[12]
	)
{
	for(int i=0; i < 3; i++)
		ShiftVector[i] = mtx[i];
	for(int i=0; i < 9; i++)
	{
		RotationMatrix[i] = mtx[i+3];
		SomeAngles[i] = 0;
	}
}

void TransformationMatrix::Export
	(
		float (&mtx)[12]
	)	const
{
	for(int i=0; i < 3; i++)
		mtx[i] = ShiftVector[i];
	for(int i=0; i < 9; i++)
	{
		mtx[i+3] = RotationMatrix[i];
	}
}