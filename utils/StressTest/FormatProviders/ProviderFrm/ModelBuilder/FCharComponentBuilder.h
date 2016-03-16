#ifndef FCHAR_COMPONENT_BUILDER_H

#define FCHAR_COMPONENT_BUILDER_H


#include "../../ProviderFrm/FrundFacade/BasicTypes.h"

#include <string>


using std::string;


class FCharComponentBuilder
{
protected:
	CharComponent _charComponent;
public:
	virtual const CharComponent& Get() const;
};

class FCharComponentTypeGravityBuilder: public FCharComponentBuilder
{
public:
	void Setup(int direction, string gravity = "9.81");
};

class FCharComponentTypeConstantBuilder: public FCharComponentBuilder
{
public:
	void Setup(int direction, const string& param);
	void Setup(int direction, double param);
};


class FCharComponentTypeKinematicJointBuilder: public FCharComponentBuilder
{
public:
	void Setup(int direction, int paramsCount);
};

class FCharComponentTypeDampingDefaultBuilder: public FCharComponentBuilder
{
public:
	void Setup(int direction);
};

class FCharComponentTypeSpringBuilder: public FCharComponentBuilder
{
public:
	void Setup(double stiffnes);
};

class FCharComponentTypeSpringType2Builder: public FCharComponentBuilder
{
public:
	void Setup(double stiffnes1, double stiffnes2, double lInterval);
};

class FCharComponentTypeDampingBuilder: public FCharComponentBuilder
{
public:
	void Setup(double damping);
};


#endif // FCHAR_COMPONENT_BUILDER_H