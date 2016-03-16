#ifndef FJOINT_BUILDER_H

#define FJOINT_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FJoint.h"

#include <string>


using std::string;


class FJointBuilder
{
protected:
	FJoint _fJoint;
public:
	virtual const FJoint& Get() const;
	virtual void SetSpringType(int springType);
	virtual void Setup(const string& name, int id, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, int charNumber);
	virtual void Setup(const string& name, int id, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, int charNumber);
	void SetLength(double length);
};


#endif // FJOINT_BUILDER_H