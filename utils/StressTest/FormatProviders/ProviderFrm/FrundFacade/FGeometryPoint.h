#ifndef FGEOMETRY_POINT_H

#define FGEOMETRY_POINT_H


#include "TransformMatrix.h"

#ifdef GNUCPP

#include "../../../fcore/wrappers/fcUnistd.h"

#else

#include <vector>


using std::vector;

#endif


class TransformMatrix;


/** Точка в 3D декартовом пространстве
*/
union FGeometryPoint
{
public:

	struct
	{
		// координаты
		float x,y,z;
	};

	// координаты
	float coords[3];


	FGeometryPoint();

	FGeometryPoint
		(
			float ix,
			float iy,
			float iz
		);

	FGeometryPoint
		(
			const FGeometryPoint& p1,
			const FGeometryPoint& p2
		);

	void Transform
		(
			const float (&transformMtx)[12],
			const FGeometryPoint& cmNode
		);

	void Transform
		(
			const TransformMatrix& mtx,
			const FGeometryPoint& cmNode
		);

	void Transform
		(
			const TransformMatrix& mtx
		);

	void Project
		(
			const vector<double>& projectionMtx,
			double& xx,
			double& yy
		)	const;

	FGeometryPoint Cross(const FGeometryPoint& point) const;
	float Dot(const FGeometryPoint& point) const;
	float Magnitude() const;
	float Angle(const FGeometryPoint& point) const;
	FGeometryPoint Normal() const;
	FGeometryPoint  operator + (const FGeometryPoint& vector) const;
	FGeometryPoint  operator - (const FGeometryPoint& vector) const;
	bool operator == (const FGeometryPoint& vector) const;
	void operator /= (float denominator);
};


#endif // FGEOMETRY_POINT_H