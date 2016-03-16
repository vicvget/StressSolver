#ifndef FFORCE_BUILDER_H

#define FFORCE_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FForce.h"

#include <string>


using std::string;


class FForceBuilder
{
protected:
	FForce _force;
public:
	virtual const FForce& Get() const;
	virtual void Setup(const string& name, int id, int bodyNumber, int nodeNumber, int charNumber);
	virtual void Setup(const string& name, int id, const BodyNodeNumber& bodyNodeNumber, int charNumber);
};


#endif // FFORCE_BUILDER_H