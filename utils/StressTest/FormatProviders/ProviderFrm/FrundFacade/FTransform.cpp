#include "FTransform.h"


using MathHelpers::MakeVec3;


FTransform::FTransform(const double* rotation, const double* translation)
{
	for(int i = 0; i < 3; i++)
	{
		_translation[i] = translation[i];
		for(int j = 0; j < 3; j++)
			_rotation[i*3+j] = rotation[i*3+j];
	}
	// TODO: переписать правильно
	//_rotation = Mat3(rotation);
	//_translation.Init(translation[0], translation[1], translation[2]);
}

Vec3CRef FTransform::GetDirectionAxis(int direction) const
{
	return _rotation.Row(direction);
}

const Mat3& FTransform::GetRotation() const
{
	return _rotation;
}

const Vec3& FTransform::GetTranslation() const
{
	return _translation;
}

Vec3 FTransform::ToGCS(const Vec3& point) const
{
	return _rotation * point + _translation;
}

Vec3 FTransform::FromGCS(const Vec3& point) const
{
	return _rotation.Tmul(point - _translation);
}

Vec3 FTransform::ToLCS(const Vec3& point, const FTransform& lcs) const
{
	return lcs.FromGCS(ToGCS(point));
}

Vec3 FTransform::ToGCSRotation(const Vec3& point) const
{
	return _rotation * point;
}

Vec3 FTransform::FromGCSRotation(const Vec3& point) const
{
	return _rotation.Tmul(point);
}

Vec3 FTransform::ToLCSRotation(const Vec3& point, const FTransform& lcs) const
{
	return lcs.FromGCSRotation(ToGCSRotation(point));
}
