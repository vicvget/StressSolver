#ifndef TRANSFORM_MATRIX_H

#define TRANSFORM_MATRIX_H


#include "FGeometryPoint.h"


union FGeometryPoint;


class TransformMatrix
{
public:

	double* operator [] (int i) {return &_mtx[i*4];}
	const double* operator [] (int i) const {return &_mtx[i*4];}
	double Element(int i, int j) const {return _mtx[i*4+j];}

	TransformMatrix();

	TransformMatrix
		(
			const float (&mtx)[12]
		);

	void FillRepresentation
		(
			float (&mtx)[12]
		)	const;

	bool operator == (const TransformMatrix& mtx) const;
	TransformMatrix operator * (const TransformMatrix& mtx) const;
	void SetTranslation(const FGeometryPoint& v);
	void SetRotation(const FGeometryPoint& axis, double angle);

private:

	double _mtx[16];

};

struct TransformationMatrix
{

	float RotationMatrix[9];
	float ShiftVector[3];
	float SomeAngles[9];

	TransformationMatrix();

	TransformationMatrix
		(
			const float (&mtx)[12]
		);

	void Export
		(
			float (&mtx)[12]
		)	const;

};


#endif // TRANSFORM_MATRIX_H