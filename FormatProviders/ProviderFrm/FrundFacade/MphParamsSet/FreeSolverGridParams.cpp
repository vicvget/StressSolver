#include "FreeSolverGridParams.h"
#include "../../../GridProvider/MeshDataSolverFormatter.h"
#include "../../../../Fcore/wrappers/FileRoutines.h"
#include <stdexcept>

FreeSolverGridParams::FreeSolverGridParams(): ParamsSetBase()
{
	Init();
}

FreeSolverGridParams::FreeSolverGridParams(const string &name)
	: ParamsSetBase(name)
{
	Init();
}

FreeSolverGridParams::FreeSolverGridParams(const FreeSolverGridParams &freeSolverGridParams)
	: ParamsSetBase(freeSolverGridParams.GetName()), _gridStep(freeSolverGridParams.GridStep())
{
	_bcMapper = new BcMapper(freeSolverGridParams.GetBcMapper());
}

FreeSolverGridParams::FreeSolverGridParams(SolverGridParams &solverGridParams) : ParamsSetBase(solverGridParams.GetName())
{	
	_gridStep = solverGridParams.GridStep();
	//delete _bcMapper;
	_bcMapper = new BcMapper(solverGridParams.GetCurrentMapper().GetName());
	_bcMapper->SetRestBcId(solverGridParams.GetCurrentMapper().GetRestBcId());
	_bcMapper->SetIsRestBcEnabled(solverGridParams.GetCurrentMapper().IsRestBcEnabled());
	map<size_t, size_t> tmp_map;
	for (const auto &j : solverGridParams.GetCurrentMapper().GetBcParams())
	{
		tmp_map[j.first] = j.second;
	}
	_bcMapper->SetBcParams(tmp_map);
	vector<BcSurface> tmp_bcsurface_vector;
	for (auto &l : solverGridParams.GetCurrentMapper().GetBcSurfaces())
	{
		BcSurface tmp_surface(l.GetName());
		if (l.IsEnabled())
		{
			tmp_surface.Enable();
		}
		else
		{
			tmp_surface.Disable();
		}
		for (auto &k : l.GetBcSurfacesIds())
		{
			tmp_surface.AddSurface(k);
		}
		tmp_bcsurface_vector.push_back(tmp_surface);
	}
	_bcMapper->SetBcSurfaces(tmp_bcsurface_vector);
}

void FreeSolverGridParams::Init()
{
	_gridStep = 1;
	_bcMapper = new BcMapper();
}

void FreeSolverGridParams::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	ofs.precision(10);
	ofs << std::fixed << _gridStep << std::endl; // для вывода шага сетки нужна большая точность!!!	
	/*текст из _bcMapper->Save(ofs), только без мапы*/
	/*_bcMapper->Save(ofs);*/
	ofs << _bcMapper->GetName() << std::endl;
	fs::WriteCommentedLine(ofs, _bcMapper->GetRestBcId(), "rest bc params id");
	int isEnabled = _bcMapper->IsRestBcEnabled() ? 1 : 0;
	fs::WriteCommentedLine(ofs, isEnabled, "flag rest bc enabled");
	fs::WriteCommentedLine(ofs, _bcMapper->CountBcs(), "count bcs");
	for (size_t i = 0; i < _bcMapper->CountBcs(); i++)
	{
		//fs::WriteCommentedLine(ofs, _bcParams.find(i)->second, "bc params id");
		_bcMapper->GetBcSurfaces()[i].Save(ofs);
	}
}

void FreeSolverGridParams::SaveBinary(ofstream &ofs)
{
	// saves name
	int tmp_int = _name.length();
	ofs.write((char *)&tmp_int, sizeof(tmp_int));
	ofs.write(_name.c_str(), tmp_int);
	// saves step	
	ofs.write(reinterpret_cast<char*>(&_gridStep), sizeof(_gridStep));	
	// saves bcMapper
	// saves name of the mapper
	tmp_int = _bcMapper->GetName().length();
	ofs.write((char *)&tmp_int, sizeof(tmp_int));
	ofs.write(_bcMapper->GetName().c_str(), tmp_int);
	// saves _isRestBcEnabled of the mapper	
	bool tmp_bool = (bool)_bcMapper->IsRestBcEnabled();
	ofs.write(reinterpret_cast<char*>(&tmp_bool), sizeof(tmp_bool));
	// saves _restBcId
	unsigned int tmp_size_t = static_cast<unsigned int>(_bcMapper->GetRestBcId());
	ofs.write(reinterpret_cast<char*>(&tmp_size_t), sizeof(tmp_size_t));
	// save map or not??? (map<size_t, size_t> _bcParams)
	// saves bcSurfaces
	tmp_int = _bcMapper->GetBcSurfaces().size();
	ofs.write((char *)&tmp_int, sizeof(tmp_int));
	for (BcSurface &j : _bcMapper->GetBcSurfaces())
	{
		// saves name of the surface set
		tmp_int = j.GetName().length();
		ofs.write((char *)&tmp_int, sizeof(tmp_int));
		ofs.write(j.GetName().c_str(), tmp_int);
		// saves _isEnabled of the surface	
		tmp_bool = j.IsEnabled();
		ofs.write(reinterpret_cast<char*>(&tmp_bool), sizeof(tmp_bool));
		// saves surface ids (names)
		tmp_int = j.CountSurfaces();
		ofs.write((char *)&tmp_int, sizeof(tmp_int));
		for (string &k : j.GetBcSurfacesIds())
		{
			tmp_int = k.length();
			ofs.write((char *)&tmp_int, sizeof(tmp_int));
			ofs.write(k.c_str(), tmp_int);
		}
	}

}

void FreeSolverGridParams::LoadBinary(ifstream &ifs)
{
	try
	{
		// проверяем, не пустой ли файл
		if (ifs.peek() == ifstream::traits_type::eof())
		{
			return;
		}

		int tmp_int = 0;

		// reads name
		ifs.read(reinterpret_cast<char*>(&tmp_int), sizeof(tmp_int));
		char *buf = new char[tmp_int];
		ifs.read(buf, tmp_int);
		_name.assign(buf, tmp_int);
		delete[] buf;
		// reads grid step
		ifs.read(reinterpret_cast<char*>(&_gridStep), sizeof(_gridStep));
		// reads bcMapper	
		// reads name of the mapper
		ifs.read(reinterpret_cast<char*>(&tmp_int), sizeof(tmp_int));
		buf = new char[tmp_int];
		ifs.read(buf, tmp_int);
		string tmp_string;
		tmp_string.assign(buf, tmp_int);
		_bcMapper->SetName(tmp_string);
		delete[] buf;
		// reads _isRestBcEnabled of the mapper	
		bool tmp_bool;
		ifs.read(reinterpret_cast<char*>(&tmp_bool), sizeof(tmp_bool));
		_bcMapper->SetIsRestBcEnabled(tmp_bool);
		// reads _restBcId
		unsigned int tmp_size_t;
		ifs.read(reinterpret_cast<char*>(&tmp_size_t), sizeof(tmp_size_t));
		_bcMapper->SetRestBcId(static_cast<size_t>(tmp_size_t));
		// reads all bcSurfaces
		int bcSurfaceCounter;
		ifs.read(reinterpret_cast<char*>(&bcSurfaceCounter), sizeof(bcSurfaceCounter));
		vector<BcSurface> tmp_bc_surface_vector;
		for (int j = 0; j < bcSurfaceCounter; j++)
		{
			BcSurface tmp_surface;
			// reads name of the surface set
			ifs.read(reinterpret_cast<char*>(&tmp_int), sizeof(tmp_int));
			buf = new char[tmp_int];
			ifs.read(buf, tmp_int);
			tmp_string.assign(buf, tmp_int);
			tmp_surface.SetName(tmp_string);
			delete[] buf;
			// reads _isEnabled of the surface set
			ifs.read(reinterpret_cast<char*>(&tmp_bool), sizeof(tmp_bool));
			tmp_bool ? tmp_surface.Enable() : tmp_surface.Disable();
			// reads all bcSurface ids (names)
			int bcSurfaceIdsCounter;
			ifs.read(reinterpret_cast<char*>(&bcSurfaceIdsCounter), sizeof(bcSurfaceIdsCounter));
			for (int k = 0; k < bcSurfaceIdsCounter; k++)
			{
				ifs.read(reinterpret_cast<char*>(&tmp_int), sizeof(tmp_int));
				buf = new char[tmp_int];
				ifs.read(buf, tmp_int);
				tmp_string.assign(buf, tmp_int);
				tmp_surface.AddSurface(tmp_string);
				delete[] buf;
			}
			tmp_bc_surface_vector.push_back(tmp_surface);
		}
		_bcMapper->SetBcSurfaces(tmp_bc_surface_vector);
	}
	catch (const std::exception&)
	{
		delete _bcMapper;
		Init();
	}
}

void FreeSolverGridParams::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	fs::ReadLineValue(ifs, _gridStep);
	/*текст из _bcMapper->Load(ifs), только без мапы и слегка переделанный*/
	//_bcMapper->Load(ifs);
	string tmp_string;
	fs::ReadLineString(ifs, tmp_string);
	_bcMapper->SetName(tmp_string);
	size_t tmp_size_t;
	fs::ReadLineValue(ifs, tmp_size_t);
	_bcMapper->SetRestBcId(tmp_size_t);
	bool tmp_bool;
	fs::ReadLineValue(ifs, tmp_bool);
	_bcMapper->SetIsRestBcEnabled(tmp_bool);
	size_t tmpCount;
	fs::ReadLineValue(ifs, tmpCount);
	vector<BcSurface> tmp_vect;
	tmp_vect.resize(tmpCount);
	for (size_t i = 0; i < tmpCount; i++)
	{
		tmp_vect[i].Load(ifs);
	}
	_bcMapper->SetBcSurfaces(tmp_vect);
}