#include "FElementBuilder.h"
#include "FForceBuilder.h"
#include "../FrundFacade/FElementType.h"

const FForce& FForceBuilder::Get() const
{
	return _force;
}

void FForceBuilder::Setup(const string& name, int id, const BodyNodeNumber& bodyNodeNumber, int charNumber)
{
	Setup(name, id, bodyNodeNumber.bodyNumber, bodyNodeNumber.nodeNumber, charNumber);
}

void FForceBuilder::Setup(const string& name, int id, int bodyNumber, int nodeNumber, int charNumber)
{
	FElementBuilder::SetupCommonProperties(_force, name, id, 
		FElementType::MakeFileId(FElementType::FETC_Force, id+1));

	_force.BodyNumber(bodyNumber);
	_force.NodeNumber(nodeNumber);
	_force.CharNumber(charNumber);
}
