#include "FJoint.h"

#include "../../../fcore/fcore.h"
#include "../../../fcore/wrappers/FileRoutines.h"

#include <sstream>
#include <string>
#include <iomanip>


using std::setw;
using std::setprecision;



FJoint::FJoint()
	:
		charNumber(),
		_springType(),
		_springLength(),
		_sortId()
{
	_excludedDirections.resize(6, false);
}

FJoint::FJoint
	(
		int orderIndex
	)
	:
		FElement(orderIndex),
		charNumber(),
		_springType(),
		_springLength(),
		_sortId()
{
	_excludedDirections.resize(6, false);
}

FJoint::FJoint
	(
		const FJoint& src
	)
	:
		FElement(src),
		_bodyNodeNumber1(src._bodyNodeNumber1),
		_bodyNodeNumber2(src._bodyNodeNumber2),
		node1(src.node1),
		node2(src.node2),
		_excludedDirections(src._excludedDirections),
		charNumber(src.charNumber),
		_springType(src._springType),
		_sspringLength(src._sspringLength),
		_springLength(src._springLength),
		_tag(src._tag),
		_sortId(src._sortId)
{
}

bool FJoint::Load(ifstream& stream)
{
	_springType = 0;

	int bodyIndex;

	fs::ReadLineString(stream, _name);
	stream >> _number;
	fs::ReadLineString(stream, _tag); // остаток строки пихаем в Tag

// stream >> number;
// stream.getline(buf, STRLIMIT); // 3ds
	int bodyNumber1, nodeNumber1;
	int bodyNumber2, nodeNumber2;

	stream >> bodyNumber1 >> bodyIndex >> nodeNumber1
		>> node1.x >> node1.y >> node1.z
		>> bodyNumber2 >> bodyIndex >> nodeNumber2
		>> node2.x >> node2.y >> node2.z
		>> charNumber >> _springType >> _springLength >> _sspringLength;
	BodyNodeNumber1(BodyNodeNumber(bodyNumber1, nodeNumber1));
	BodyNodeNumber2(BodyNodeNumber(bodyNumber2, nodeNumber2));

	int tmp;

	_excludedDirections.clear();
	for (int i = 0; i < 6; i++)
	{
		stream >> tmp;
		_excludedDirections.push_back(tmp);
	}
	// TODO: еще есть информация о длине в свободном состоянии и куча каких-то чисел

	return true;
}

void FJoint::Save(ofstream& ofs) const
{
	int bodyIndex = 0; // TODO: что это?
	ofs << _name << std::endl;
	ofs << _number << _tag << std::endl; // TODO: &tag
	ofs << BodyNumber1() << ' '
		<< bodyIndex << ' ' 
		<< NodeNumber1() << ' '
		<< setiosflags(ios_base::fixed) 
		<< setprecision(6)
		<< node1.x << ' ' 
		<< node1.y << ' '
		<< node1.z << std::endl;
	ofs << BodyNumber2() << ' '
		<< bodyIndex << ' ' 
		<< NodeNumber2() << ' '
		<< setiosflags(ios_base::fixed) 
		<< setprecision(6)
		<< node2.x << ' ' 
		<< node2.y << ' '
		<< node2.z << std::endl;
	ofs << charNumber << std::endl;	
	ofs << _springType << ' ';
	if(_springType == 1)
		ofs << setiosflags(ios_base::fixed) 
		<< setprecision(6)
		<< _springLength << std::endl
		<< _sspringLength << std::endl;
	else
		ofs << setiosflags(ios_base::fixed) 
		<< setprecision(6)
		<< 0.0f << std::endl
		<< _sspringLength << std::endl;
	for(int i = 0; i < 6; i++)
	{
		int tmp = _excludedDirections[i] ? 1:0;
		ofs << tmp << std::endl;	
	}
	// TODO: missing clinks logic
	ofs << 0 << std::endl;		
}

void FJoint::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	if (!isRType)
	{
		// BL17 NSE   NT1  NUZ1   NT2  NUZ2   NMX
		int l = 4;
		stream << "@ \'" << _name << "\'\n";
		stream << setw(9) << _number << ' ' <<
				setw(l) << BodyNumber1() << ' ' <<
				setw(l) << NodeNumber1() << ' ' <<
				setw(l) << BodyNumber2() << ' ' <<
				setw(l) << NodeNumber2() << ' ' <<
				setw(l) << charNumber << '\n';
    
		if(_springType != 0)
		{
			stringstream ss;
			ss << "2 " 
				<< setw(4) << _number << ' ' 
				<< setw(4) <<  _springType << ' '
				<< setw(4) << 0 << ' '
				<< setw(4) << 0 << ' '
				<< setw(4) << 0 << ' '
				<< std::endl;
			fs::AppendLineToFile(FILE_MHL_SPRINGS, ss);
		}
	}
	else
	{
		if(_springType != 0)
		{
			stringstream ss;
			ss << "3 " 
				<< setw(4) << _number << ' ' 
				<< setw(4) << 0 << ' '
				<< setw(4) << 0 << ' '
				<< setw(4) << 0 << ' '
				<< setw(4) << _springLength << ' '
				<< setw(4) << 0 << ' '
				<< setw(4) << 0 << ' '
				<< std::endl;
			fs::AppendLineToFile(FILE_RHL_SPRINGS, ss);
		}
	}
}

void FJoint::Eval(Calc::Calculator &calc) const
{
	_springLength = calc.Eval(_sspringLength.c_str());
}

