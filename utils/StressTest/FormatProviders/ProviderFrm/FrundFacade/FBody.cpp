#include "FBody.h"

#include "../../../fcore/fcore.h"
#include "../../../fcore/Exceptions/fcExceptions.h"

#include <string>
#include <sstream>
#include <cstring>


using std::setw;
using std::endl;


FBodyDof::FBodyDof()
	:
		x(3),
		y(3),
		z(3),
		rx(1),
		ry(1),
		rz(1)
{
}

void FBodyDof::Init(int x, int y, int z, int rx, int ry, int rz)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->rx = rx;
	this->ry = ry;
	this->rz = rz;
}

bool FBodyDof::IsSmallMovementsMode() const
{
	return
		!((rx == 3) || (ry == 3) || (rz == 3)) &&
		((rx == 1) || (ry == 1) || (rz == 1));
}

bool FBodyDof::IsFixed() const
{
	return x == 0 && y == 0 && z == 0 && rx == 0 && ry == 0 && rz == 0;
}


FBodyInertia::FBodyInertia()
	:
		m(1.0),
		jx(1.0),
		jy(1.0),
		jz(1.0),
		strm("1"),
		strjx("1"),
		strjy("1"),
		strjz("1"),
		minMass(0.0),
		minInertia(0.0)
{
}

void FBodyInertia::Eval(Calc::Calculator &calc) const
{
	m = calc.Eval(strm.c_str());
	jx = calc.Eval(strjx.c_str());
	jy = calc.Eval(strjy.c_str());
	jz = calc.Eval(strjz.c_str());
}

void FBodyInertia::ToStrings(string& strm_out, string& strjx_out, string& strjy_out, string& strjz_out) const
{
	strm_out = ToStringByMinValue(strm, m, minMass);
	strjx_out = ToStringByMinValue(strjx, jx, minInertia);
	strjy_out = ToStringByMinValue(strjy, jy, minInertia);
	strjz_out = ToStringByMinValue(strjz, jz, minInertia);
}

// static
string FBodyInertia::ToStringByMinValue(const string& str, double value, double minValue)
{
	string res = str;

	if (value < minValue)
	{
		std::stringstream sstr;

		sstr << minValue;
		sstr >> res;
	}

	return res;
}


FFlexibleBodyParams::FFlexibleBodyParams()
	:
		isFlexible(false),
		flexFormNumber(0),
		flexFormPath(),
		node1(0),
		node2(0),
		node3(0),
		freeze()
{
}


FBody::FBody()
	:
		_geometryId(),
		_mph(nullptr)
{
}

FBody::FBody
	(
		const FBody& src
	)
	:
		FElement(src),
		_dof(src._dof),
		_inertia(src._inertia),
		_flexible(src._flexible),
		_geometryId(src._geometryId),
		_mph(src._mph)
{
}

void FBody::Load(const char* path)
{
	DefaultLoad(path);

	string strPath(path);

	int pos = strPath.find_first_of(EXT_MODEL_ELEMENT);

	if (pos == string::npos)
	{
		exceptions::ThrowFileInvalidExtension(strPath);
	}
	_fileId = strPath.substr(0, pos);

	//this->_geometry = nullptr;
	_geometryId = _id; // 1-based, 0 = not initialized
	_mph = nullptr;

	string fileMph = _fileId + EXT_MULTIPHYSICS;
	ifstream streamMph(fileMph);

	if (streamMph.is_open())
	{
		//_multiphysics = new FMultiphysics(streamMph);
		_mph = new FMph(streamMph);
		streamMph.close();
	}
	MinMass(0);
	MinInertia(0);
}

bool FBody::Load(ifstream& ifs)
{
	//_geometry = nullptr;
	fs::ReadLineString(ifs, _name);
	fs::ReadLineValue(ifs, _flexible.flexFormNumber);
	_flexible.isFlexible = _flexible.flexFormNumber ? 1 : 0;

	int x, y, z, rx, ry, rz;

	fs::ReadLineValue(ifs, x);
	fs::ReadLineValue(ifs, y);
	fs::ReadLineValue(ifs, z);
	fs::ReadLineValue(ifs, rx);
	fs::ReadLineValue(ifs, ry);
	fs::ReadLineValue(ifs, rz);
	_dof.Init(x,y,z,rx,ry,rz);

	fs::ReadLineValue(ifs, _inertia.m);
	fs::ReadLineValue(ifs, _inertia.jx);
	fs::ReadLineValue(ifs, _inertia.jy);
	fs::ReadLineValue(ifs, _inertia.jz);
	fs::ReadLineString(ifs, _inertia.strm);
	fs::ReadLineString(ifs, _inertia.strjx);
	fs::ReadLineString(ifs, _inertia.strjy);
	fs::ReadLineString(ifs, _inertia.strjz);
	fs::ReadLineString(ifs, _flexible.flexFormPath);

	if (ifs.eof())
	{
		return false;
	}
	ifs >> _number >> _flexible.freeze;

	return true;
}

void FBody::SaveByIndex() const
{
	FElement::SaveByIndex();
	
	//if(this->_geometry != nullptr)
	//{
	//	string path = this->_fileId + EXT_BODY_GEOMETRY;
	//	ofstream ofs(path);
	//	if(ofs.is_open())
	//	{
	//		this->_geometry->Save(ofs);
	//		ofs.close();
	//	}
	//	else
	//	{
	//		exceptions::ThrowFileNotOpened(path);
	//	}
	//}

//	if(this->_multiphysics != nullptr)
	
	if (_mph != nullptr)
	{
		string path = _fileId + EXT_MULTIPHYSICS;
		ofstream ofs(path);

		if (ofs.is_open())
		{
			//this->_multiphysics->Save(ofs);
			this->_mph->Save(ofs);
			ofs.close();
		}
		else
		{
			exceptions::ThrowFileNotOpened(path);
		}
	}
	else
	{
		string path = _fileId + EXT_MULTIPHYSICS;

		remove(path.c_str());
	}
}

void FBody::Save(ofstream& ofs) const
{
	ofs << _name << endl;
	ofs << IsFlexible() << endl;
	ofs << X() << endl 
		<< Y() << endl 
		<< Z() << endl 
		<< Rx() << endl 
		<< Ry() << endl
		<< Rz() << endl;
	ofs.setf(std::ios_base::fixed);
	ofs.precision(6);
	ofs << M() << endl 
		<< Jx() << endl 
		<< Jy() << endl 
		<< Jz() << endl
		<< Strm() << endl 
		<< Strjx() << endl 
		<< Strjy() << endl 
		<< Strjz() << endl
		<< FlexFormPath() << endl
		<< _number << ' ' << Freeze() << endl;
}

void FBody::OutputLba() const
{
	if (!_fileId.empty())
	{
		string fileName = _fileId + ".lba";
		ofstream stream(fileName);

		if (stream.is_open())
		{
			OutputLba(stream);
			stream.close();
		}
	}
}

void FBody::OutputLba(ofstream& stream) const
{
	if (_flexible.isFlexible)
	{
		exceptions::ThrowMessage("Flexible is Not Supported");
	}
	else
	{
		stream.setf(ios_base::fixed | ios_base::showpoint);
		stream.precision(6);
		stream << 2 << '\n' << 
			GetIdToOutput() << '\n' << 
			_inertia.m << '\n' << 
			_name << '\n' << 
			_inertia.jx << ' ' << 
			_inertia.jy << ' ' << 
			_inertia.jz << ' ' << 
			"0.0" << ' ' << "0.0" << ' ' << "0.0" << '\n';
	}
}

void FBody::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	if (isRType)
	{
		// BL31 NPDS    MASSA      JX      JY      JZ
		string strm_out, strjx_out, strjy_out, strjz_out;

		_inertia.ToStrings(strm_out, strjx_out, strjy_out, strjz_out);

		int l = 8;

		stream << "@ \'" << _name << "\'\n";
		stream << setw(9) << _number << ' ' <<
			setw(l) << strm_out << ' ' <<
			setw(l) << strjx_out << ' ' <<
			setw(l) << strjy_out << ' ' <<
			setw(l) << strjz_out << '\n';
	}
	else
	{
		// BL14 NPDS NTIP   LX   LY   LZ   UX   UY   UZ
		int l = 4;

		stream << "@ \'" << _name << "\'\n";
		stream << setw(9) << _number << ' ' <<
			setw(l) << GetIdToOutput() << ' ' <<
			setw(l) << _dof.x << ' ' << 
			setw(l) << _dof.y << ' ' << 
			setw(l) << _dof.z << ' ' <<
			setw(l) << _dof.rx  << ' ' << 
			setw(l) << _dof.ry << ' ' << 
			setw(l) << _dof.rz << '\n';
	}
}

void FBody::Eval(Calc::Calculator &calc) const
{
	_inertia.Eval(calc);
}

SolverParamsBase* FBody::AddNewSolver(SolverTypes solverType)
{
	if (GetMph() == nullptr)
	{
		_mph = new FMph(_fileId);
	}

	return _mph->AddSolver(solverType);
}

/**
* Признак малых движений
*/
bool FBody::IsSmallMovementsMode() const
{
	return _dof.IsSmallMovementsMode();
}

int FBody::GetIdToOutput() const
{
	return _flexible.isFlexible ? _id + 2001 : _id + 101;
}
