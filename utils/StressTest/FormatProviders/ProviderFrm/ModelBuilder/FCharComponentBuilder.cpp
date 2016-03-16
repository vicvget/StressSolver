#include "FCharComponentBuilder.h"
#include "../../../fcore/wrappers/StringRoutines.h"

const CharComponent& FCharComponentBuilder::Get() const
{
	return _charComponent;
}

void FCharComponentTypeGravityBuilder::Setup(int direction, string gravity)
{
	_charComponent.direction = direction;
	_charComponent.type = 3;
	_charComponent.params.clear();
	_charComponent.params.push_back(gravity);
}

void FCharComponentTypeKinematicJointBuilder::Setup(int direction, int paramsCount)
{
	_charComponent.direction = direction;
	_charComponent.type = 15;
	_charComponent.params.clear();
	for(int i = 0; i < paramsCount; i++)
		_charComponent.params.push_back("0.0");	
}

void FCharComponentTypeDampingDefaultBuilder::Setup(int direction)
{
	_charComponent.direction = direction;
	_charComponent.type = 50;
	_charComponent.params.clear();
	_charComponent.params.push_back("100.0");	
}

void FCharComponentTypeSpringBuilder::Setup(double stiffnes)
{
	_charComponent.direction = 1;
	_charComponent.type = 1;
	_charComponent.params.clear();	
	_charComponent.params.push_back(NumberToString(stiffnes));	
}

void FCharComponentTypeDampingBuilder::Setup(double damping)
{
	_charComponent.direction = 1;
	_charComponent.type = 50;
	_charComponent.params.clear();	
	_charComponent.params.push_back(NumberToString(damping));	
}

void FCharComponentTypeConstantBuilder::Setup(int direction, const string& param)
{
	_charComponent.direction = direction;
	_charComponent.type = 2;
	_charComponent.params.clear();	
	_charComponent.params.push_back(param);
}

void FCharComponentTypeConstantBuilder::Setup(int direction, double param)
{
	Setup(direction, NumberToString(param));
}

void FCharComponentTypeSpringType2Builder::Setup(double stiffnes1, double stiffnes2, double lInterval)
{
	_charComponent.direction = 1;
	_charComponent.type = 2;
	_charComponent.params.clear();	
	_charComponent.params.push_back(NumberToString(stiffnes1));	
	_charComponent.params.push_back(NumberToString(stiffnes2));	
	_charComponent.params.push_back(NumberToString(lInterval));	
	_charComponent.params.push_back("0.");	
	
}
