#ifndef FTRANSFORM_H

#define FTRANSFORM_H


#include "../../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../../AdditionalModules/fmath/Vector3.h"


using MathHelpers::Vec3;
using MathHelpers::Vec3CRef;
using MathHelpers::Mat3;


class FTransform
{
	Mat3 _rotation;
	Vec3 _translation;

public:

	FTransform(const double* rotation, const double* translation);

	Vec3CRef GetDirectionAxis(int direction) const;
	const Mat3& GetRotation() const;
	const Vec3& GetTranslation() const;
	
	/************************************************************************/
	/* Перевод в глобальную систему координат */
	/************************************************************************/
	Vec3 ToGCS(const Vec3& point) const;

	/************************************************************************/
	/* Перевод из глобальной системы координат */
	/************************************************************************/
	Vec3 FromGCS(const Vec3& point) const;

	/************************************************************************/
	/* Перевод в локальную систему координат lcs */
	/************************************************************************/
	Vec3 ToLCS(const Vec3& point, const FTransform& lcs) const;

	/************************************************************************/
	/* Перевод в глобальную систему координат */
	/************************************************************************/
	Vec3 ToGCSRotation(const Vec3& point) const;

	/************************************************************************/
	/* Перевод из глобальной системы координат */
	/************************************************************************/
	Vec3 FromGCSRotation(const Vec3& point) const;

	/************************************************************************/
	/* Перевод в локальную систему координат lcs */
	/************************************************************************/
	Vec3 ToLCSRotation(const Vec3& point, const FTransform& lcs) const;

};


#endif // FTRANSFORM_H
