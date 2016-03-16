#include "FJointBuilder.h"
#include "FElementBuilder.h"
#include "../FrundFacade/FElementType.h"
#include "../../../fcore/wrappers/StringRoutines.h"
const FJoint& FJointBuilder::Get() const
{
	return _fJoint;
}

void FJointBuilder::Setup(const string& name, int id, 
						  int bodyNumber1, int nodeNumber1, 
						  int bodyNumber2, int nodeNumber2, int charNumber)
{
	FElementBuilder::SetupCommonProperties(_fJoint, name, id, 
		FElementType::MakeFileId(FElementType::FETC_Joint, id+1));

	_fJoint.BodyNodeNumber1(BodyNodeNumber(bodyNumber1,nodeNumber1));
	_fJoint.BodyNodeNumber2(BodyNodeNumber(bodyNumber2,nodeNumber2));
	_fJoint.CharNumber(charNumber);
}

void FJointBuilder::Setup(const string& name, int id, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, int charNumber)
{
	Setup(name, id, bodyNode1.bodyNumber, bodyNode1.nodeNumber, bodyNode2.bodyNumber, bodyNode2.nodeNumber, charNumber);
}

void FJointBuilder::SetSpringType(int springType)
{
	_fJoint.SpringType(springType);
}

void FJointBuilder::SetLength(double length)
{
	_fJoint.SpringLength(length);
	_fJoint.SspringLength(NumberToString(length));
}
