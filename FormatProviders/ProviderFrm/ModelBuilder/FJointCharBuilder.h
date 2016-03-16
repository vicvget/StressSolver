#ifndef FJOINT_CHAR_BUILDER_H

#define FJOINT_CHAR_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FJointChar.h"

#include <string>


using std::string;


class FJointCharBuilder
{
protected:
	FJointChar _fJointChar;
	void Setup(const string& name, int id);
public:
	virtual const FJointChar& Get() const;
};

class FJointCylindricalCharBuilder: public FJointCharBuilder
{
public:
	void Setup(int directionAxis, const string& name, int id);
};

class FJointSphericalCharBuilder: public FJointCharBuilder
{
public:
	void Setup(const string& name, int id);
};

class FJointInPlaneCharBuilder: public FJointCharBuilder
{
public:
	void Setup(int directionAxis, const string& name, int id);
};

class FJointSpringDamperCharBuilder: public FJointCharBuilder
{
public:
	void Setup(const string& name, int id, double stiffness, double damping);
};

class FJointSpringType2DamperCharBuilder: public FJointCharBuilder
{
public:
	void Setup(const string& name, int id, double stiffness1, double stiffness2, double lInterval, double damping);
};


#endif // FJOINT_CHAR_BUILDER_H