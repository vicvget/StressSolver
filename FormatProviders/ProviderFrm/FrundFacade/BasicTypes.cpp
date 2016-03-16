#include "BasicTypes.h"

#include "../../../AdditionalModules/RealsComparing/RealComparer.h"


CharComponent::CharComponent():
	stiffness(0.), 
	type(1), 
	direction(1)
{
}

BodyNodeNumber::BodyNodeNumber()
	:
		bodyNumber(),
		nodeNumber()
{
}

BodyNodeNumber::BodyNodeNumber
	(
		int bodyNumber,
		int nodeNumber
	)
	:
		bodyNumber(bodyNumber),
		nodeNumber(nodeNumber)
{
}


// struct FDebugParams

bool FDebugParams::Load()
{
	isLoaded = false;
	if (fs::FileExist(FILE_DEBUG_PARAMS))
	{
		FDebugParams tmpDebugParams;
		ifstream ifs(FILE_DEBUG_PARAMS);

		if (!ifs.is_open())
		{
			return false;
		}

		ifs >> tmpDebugParams.minMass >> tmpDebugParams.minInertia >> tmpDebugParams.maxStiffness;
		ifs.close();

		DoubleComparer comparer;

		if (comparer.GreaterThanZero(tmpDebugParams.minMass))
		{
			minMass = tmpDebugParams.minMass;
			isLoaded = true;
		}
		if (comparer.GreaterThanZero(tmpDebugParams.minInertia))
		{
			minInertia = tmpDebugParams.minInertia;
			isLoaded = true;
		}
		if (comparer.GreaterThanZero(tmpDebugParams.maxStiffness))
		{
			maxStiffness = tmpDebugParams.maxStiffness;
			isLoaded = true;
		}
	}

	return true;
}

bool FDebugParams::Save() const
{
	if (isLoaded)
	{
		ofstream ofs(FILE_DEBUG_PARAMS);

		if (!ofs.is_open())
		{
			return false;
		}
		ofs
			<< minMass << std::endl
			<< minInertia << std::endl
			<< (maxStiffness >= MaxStiffness ? 0.0 : maxStiffness) << std::endl;
		ofs.close();
	}

	return true;
}


FGeometryLink::FGeometryLink()
	:
		node1(),
		node2()
{
}

FGeometryLink::FGeometryLink
	(
		int n1,
		int n2
	)
	:
		node1(n1),
		node2(n2)
{
}

bool operator==(const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2)
{
	return (bodyNodeNumber1.bodyNumber == bodyNodeNumber2.bodyNumber) 
		&& (bodyNodeNumber1.nodeNumber == bodyNodeNumber2.nodeNumber);
}
