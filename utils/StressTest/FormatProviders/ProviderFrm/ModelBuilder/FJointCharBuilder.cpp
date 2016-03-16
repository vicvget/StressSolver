#include "FJointCharBuilder.h"

#include "FCharComponentBuilder.h"
#include "FElementBuilder.h"
#include "../FrundFacade/FElementType.h"

const FJointChar& FJointCharBuilder::Get() const
{
	return _fJointChar;
}

void FJointCharBuilder::Setup(const string& name, int id)
{
	FElementBuilder::SetupCommonProperties(_fJointChar, name, id, 
		FElementType::MakeFileId(FElementType::FETC_JointChar, id+1));
}


void FJointCylindricalCharBuilder::Setup(int directionAxis, const string& name, int id)
{
	FJointCharBuilder::Setup(name, id);
	FCharComponentTypeKinematicJointBuilder charComponentTypeKinematicJointBuilder;
	charComponentTypeKinematicJointBuilder.Setup(directionAxis, 1);
	FCharComponentTypeDampingDefaultBuilder charComponentTypeDampingDefaultBuilder;
	charComponentTypeDampingDefaultBuilder.Setup(directionAxis);
	for(int i = 1; i <= 3; i++)
	{
		if(i != directionAxis)
		{
			_fJointChar.AddCharComponent(charComponentTypeKinematicJointBuilder.Get(), charComponentTypeDampingDefaultBuilder.Get(), i);
		}
	}
}

void FJointSphericalCharBuilder::Setup(const string& name, int id)
{
	FJointCharBuilder::Setup(name, id);
	FCharComponentTypeKinematicJointBuilder charComponentTypeKinematicJointBuilder;
	charComponentTypeKinematicJointBuilder.Setup(1, 1);
	FCharComponentTypeDampingDefaultBuilder charComponentTypeDampingDefaultBuilder;
	charComponentTypeDampingDefaultBuilder.Setup(1);
	for(int i = 1; i <= 3; i++)
	{
		_fJointChar.AddCharComponent(charComponentTypeKinematicJointBuilder.Get(), charComponentTypeDampingDefaultBuilder.Get(), i);
	}
}

void FJointInPlaneCharBuilder::Setup(int direction, const string& name, int id)
{
	FJointCharBuilder::Setup(name, id);
	FCharComponentTypeKinematicJointBuilder charComponentTypeKinematicJointBuilder;
	charComponentTypeKinematicJointBuilder.Setup(direction, 2);
	FCharComponentTypeDampingDefaultBuilder charComponentTypeDampingDefaultBuilder;
	charComponentTypeDampingDefaultBuilder.Setup(direction);
	_fJointChar.AddCharComponent(charComponentTypeKinematicJointBuilder.Get(), charComponentTypeDampingDefaultBuilder.Get(), direction);
}

void FJointSpringDamperCharBuilder::Setup(const string& name, int id, double stiffness, double damping)
{
	FJointCharBuilder::Setup(name, id);
	FCharComponentTypeSpringBuilder charComponentTypeSpringBuilder;
	charComponentTypeSpringBuilder.Setup(stiffness);
	FCharComponentTypeDampingBuilder charComponentTypeDampingBuilder;
	charComponentTypeDampingBuilder.Setup(damping);
	_fJointChar.AddCharComponent(charComponentTypeSpringBuilder.Get(), charComponentTypeDampingBuilder.Get(), 1);
}

void FJointSpringType2DamperCharBuilder::Setup(const string& name, int id, double stiffness1, double stiffness2, double lInterval, double damping)
{
	FJointCharBuilder::Setup(name, id);
	FCharComponentTypeSpringType2Builder charComponentTypeSpringBuilder;
	charComponentTypeSpringBuilder.Setup(stiffness1, stiffness2, lInterval);
	FCharComponentTypeDampingBuilder charComponentTypeDampingBuilder;
	charComponentTypeDampingBuilder.Setup(damping);
	_fJointChar.AddCharComponent(charComponentTypeSpringBuilder.Get(), charComponentTypeDampingBuilder.Get(), 1);

}
