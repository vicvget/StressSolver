#include "FForce.h"

#include "../../../fcore/fcore.h"


using std::setw;
using std::setprecision;
using std::endl;


FForce::FForce()
	:
		_charNumber()
{
}

FForce::FForce
	(
		int orderIndex
	)
	:
		FElement(orderIndex),
		_charNumber()
{
}

FForce::FForce
	(
		const FForce& src
	)
	:
		FElement(src),
		_bodyNodeNumber(src._bodyNodeNumber),
		_node(src._node),
		_charNumber(src._charNumber)
{
}

FForce::FForce
	(
		const char* path
	)
	:
		_charNumber()
{
	FElement::Load(path);
}

FForce::FForce
	(
		ifstream& ifs
	)
	:
		_charNumber()
{
	Load(ifs);
}

bool FForce::Load(ifstream& stream)
{
	int bodyIndex;

	fs::ReadLineString(stream, _name);
	stream >> _number;
	_charNumber = _number;
	stream >> _bodyNodeNumber.bodyNumber >> bodyIndex >> _bodyNodeNumber.nodeNumber;
	stream >> _node.x >> _node.y >> _node.z;

	return true;
}

void FForce::Save(ofstream& ofs) const
{
	int bodyIndex = 0;
	ofs << _name << endl
		<< _number << endl
		<< _bodyNodeNumber.bodyNumber << ' '
		<< bodyIndex << ' '
		<< _bodyNodeNumber.nodeNumber << ' '
		<< setiosflags(ios_base::fixed)
		<< setprecision(6)
		<< _node.x << ' ' 
		<< _node.y << ' '
		<< _node.z << std::endl;
}

void FForce::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	if (!isRType)
	{
		// BL19 NSS    NT   NUZ
		int l = 4;
		stream << "@ \'" << _name << "\'\n";
		stream << setw(9) << _number << ' ' <<
				setw(l) << _bodyNodeNumber.bodyNumber << ' ' <<
				setw(l) << _bodyNodeNumber.nodeNumber << '\n';
	}
}
