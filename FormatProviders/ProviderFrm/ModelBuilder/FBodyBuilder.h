#ifndef FBODY_BUILDER_H

#define FBODY_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FBody.h"
#include "../../../AdditionalModules/fmath/Vector3.h"

#include <string>


using std::string;
using MathHelpers::Vec3;


class FBodyBuilder
{
protected:
	FBody _fBody;

public:
	virtual const FBody& Get() const;
};

class FFreeBodyBuilder: public FBodyBuilder
{
public:
	void Setup(const string& name, int id, double mass, Vec3 inertia, int geometryId);
};


class FFixedBodyBuilder: public FBodyBuilder
{
public:
	void Setup(const string& name, int id, double mass, Vec3 inertia, int geometryId);
};


#endif // FBODY_BUILDER_H