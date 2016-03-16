#include "FGeometryPoint.h"

#include "TransformMatrix.h"

#include "../../../fcore/fcore.h"
#include "../../../fcore/Exceptions/fcExceptions.h"

#include <cmath>

#ifdef GNUCPP
#include "../../../fcore/wrappers/fcUnistd.h"
#else
#include <vector>
#endif


using std::vector;


FGeometryPoint::FGeometryPoint()
	:
		x(),
		y(),
		z()
{
}

FGeometryPoint::FGeometryPoint
	(
		float ix,
		float iy,
		float iz
	)
	:
		x(ix),
		y(iy),
		z(iz)
{
}

FGeometryPoint::FGeometryPoint
	(
		const FGeometryPoint& p1,
		const FGeometryPoint& p2
	)
	:
		x(p2.x - p1.x),
		y(p2.y - p1.y),
		z(p2.z - p1.z)
{
}

FGeometryPoint FGeometryPoint::Cross(const FGeometryPoint& point) const
{
	FGeometryPoint p;
	p.x = y*point.z-z*point.y;
	p.y = -x*point.z+z*point.x;
	p.z = x*point.y-y*point.x;
	return p;
}

void FGeometryPoint::Transform
	(
		const float (&transformMtx)[12],
		const FGeometryPoint& cmNode
	)
{
	float tcoords[3] = {x - cmNode.x, y - cmNode.y, z - cmNode.z};

	for (int i = 0; i < 3; i++)
	{
		coords[i] = transformMtx[i] * 1000 +
			transformMtx[i + 3] * tcoords[0] + transformMtx[i + 6] * tcoords[1] + transformMtx[i + 9] * tcoords[2];
	}
}

void FGeometryPoint::Transform
	(
		const TransformMatrix& mtx,
		const FGeometryPoint& cmNode
	)
{
	float tcoords[3] = {x - cmNode.x, y - cmNode.y, z - cmNode.z};

	for (int i = 0; i < 3; i++)
	{
		coords[i] = static_cast<float>(mtx[i][3]);
		for (int j = 0; j < 3; j++)
		{
			coords[i] += static_cast<float>(mtx[i][j] * tcoords[j]);
		}
	}
}

void FGeometryPoint::Transform
	(
		const TransformMatrix& mtx
	)
{
	Transform(mtx, FGeometryPoint());
}

void FGeometryPoint::Project
	(
		const vector<double>& projectionMtx,
		double& xx,
		double& yy
	)	const
{
	xx = coords[0]*projectionMtx[0]+coords[1]*projectionMtx[1]+coords[2]*projectionMtx[2];
	yy = coords[0]*projectionMtx[3]+coords[1]*projectionMtx[4]+coords[2]*projectionMtx[5];
}

float FGeometryPoint::Dot(const FGeometryPoint& point) const
{
	return x*point.x+y*point.y+z*point.z;
}

float FGeometryPoint::Magnitude() const
{
	return sqrt(x*x+y*y+z*z);
}

float FGeometryPoint::Angle(const FGeometryPoint& point) const
{
	float denominator = Magnitude() * point.Magnitude();
	if (IS_ZERO(denominator))
	{
		return 0.0f;
	}

	float cosAngle = Dot(point) / denominator;

	if (cosAngle < -1.0f)
	{
		return acos(-1.0f);
	}
	else
	{
		if (cosAngle > 1.0f)
		{
			return acos(1.0f);
		}
	}

	return acos(cosAngle);
}

FGeometryPoint FGeometryPoint::Normal() const
{
	FGeometryPoint p(1., 0., 0.);
	if(IS_ZERO(Cross(p).Magnitude()))
	{
		p.x = 0;
		p.y = 1;
	}
	p = Cross(p);
	return p;
}

FGeometryPoint  FGeometryPoint::operator + (const FGeometryPoint& vector) const
{
	FGeometryPoint res;
	for(int i = 0; i < 3; i++)
	{
		res.coords[i]=coords[i]+vector.coords[i];
	}
	return res;
}

FGeometryPoint  FGeometryPoint::operator - (const FGeometryPoint& vector) const
{
	FGeometryPoint res;
	for(int i = 0; i < 3; i++)
	{
		res.coords[i]=coords[i]-vector.coords[i];
	}
	return res;
}

void FGeometryPoint::operator /= (float denominator)
{
	for(int i = 0; i < 3; i++)
	{
		coords[i] /= denominator;
	}
}

bool FGeometryPoint::operator == (const FGeometryPoint& vector) const
{
	return IS_ZERO(x-vector.x) && IS_ZERO(y-vector.y) && IS_ZERO(z-vector.z);
}