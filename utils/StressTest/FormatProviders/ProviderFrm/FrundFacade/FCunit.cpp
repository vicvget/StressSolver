#include "FCunit.h"

#include "../../../fcore/fcore.h"

#include <string>


using std::setw;


FCunit::FCunit()
	:
		_x(),
		_y()
{
}

FCunit::FCunit
	(
		int orderIndex
	)
	:
		FElement(orderIndex),
		_x(),
		_y()
{
}

bool FCunit::Load(ifstream& stream)
{
	fs::ReadLineString(stream, _name);
	fs::ReadLineValue(stream, _number);
	fs::ReadLineValue(stream, _id);
	stream >> _x >> _y;
	return true;
}

void FCunit::Save(ofstream& ofs) const
{
	ofs << _name << std::endl
		<< _number << std::endl
		<< _id << std::endl;
	ofs << _x << ' ' << _y << std::endl;
}

void FCunit::Output(ofstream& stream, int isRType) const
{
	if(isRType) 
	{
		// BL31 NPDS    MASSA      JX      JY      JZ	

		//int l = 8;
		//stream << "@ \'" << name << "\'\n";
		//stream << setw(9) << number << ' ' <<
		//	setw(l) << strm << ' ' <<
		//	setw(l) << strjx << ' ' << 
		//	setw(l) << strjy << ' ' << 
		//	setw(l) << strjz << '\n';
	}
	else
	{
		// BL14 NPDS NTIP   LX   LY   LZ   UX   UY   UZ
		int l = 4;
		stream << "@ \'" << _name << "\'\n";
		stream << setw(9) << _number << ' ' <<
			setw(l) << _id << ' ' <<
			setw(l) << 0 << ' ' << 
			setw(l) << 0 << ' ' << 
			setw(l) << 0 << ' ' <<
			setw(l) << 0 << ' ' << 
			setw(l) << 0 << ' ' << 
			setw(l) << 0 << '\n';
	}
}