#include "BcMapper.h"
#include "../../../../fcore/wrappers/FileRoutines.h"
#include "../../../../fcore/Exceptions/fcExceptions.h"

BcMapper::BcMapper()
{
	Init();
}

BcMapper::BcMapper
	(
		const std::string& name
	)
	:
		ParamsSetBase(name)
{
	Init();
}

void BcMapper::Init()
{
	_restBcId = 0;
	_isRestBcEnabled = false;
}

void BcMapper::Clear()
{
	this->_bcSurfaces.clear();
	this->_bcParams.clear();
}

// TODO: разобраться, что это за параметр bcType
void BcMapper::AddBc(int /* bcType */)
{
	BcSurface bcSurface;
	_bcSurfaces.push_back(bcSurface);
	_bcParams[_bcSurfaces.size()-1] = 0;
}

// TODO: разобраться, что это за параметр bcType
BcSurface& BcMapper::AddNewBc(int /* bcType */)
{
	BcSurface bcSurface;
	_bcSurfaces.push_back(bcSurface);
	_bcParams[_bcSurfaces.size()-1] = 0;
	return _bcSurfaces[_bcSurfaces.size()-1]; 
}


void BcMapper::RemoveBc(size_t bcId)
{
	if(bcId < _bcSurfaces.size())
	{
		_bcSurfaces.erase(_bcSurfaces.begin() + bcId);
		if(_bcParams.find(bcId) != _bcParams.end())
			_bcParams.erase(_bcParams.find(bcId));
	}
}

void BcMapper::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	//ofs << GetRestBcId() << std::endl;
	fs::WriteCommentedLine(ofs, GetRestBcId(), "rest bc params id");
	int isEnabled = IsRestBcEnabled() ? 1 : 0;
	fs::WriteCommentedLine(ofs, isEnabled, "flag rest bc enabled");
	fs::WriteCommentedLine(ofs, CountBcs(), "count bcs");
	
	//ofs << isEnabled << std::endl;		
	//ofs << CountBcs() << std::endl;
	for(size_t i = 0; i < CountBcs(); i++)
	{
		//ofs << _bcParams[i] << std::endl;
		fs::WriteCommentedLine(ofs, _bcParams.find(i)->second, "bc params id");
		_bcSurfaces[i].Save(ofs);
	}
}

void BcMapper::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	fs::ReadLineValue(ifs, _restBcId);
	fs::ReadLineValue(ifs, _isRestBcEnabled);	

	size_t tmpCount, tmpId;
	fs::ReadLineValue(ifs, tmpCount);
	_bcSurfaces.resize(tmpCount);
	for(size_t i = 0; i < tmpCount; i++)
	{
		fs::ReadLineValue(ifs, tmpId);
		_bcSurfaces[i].Load(ifs);
		_bcParams[i] = tmpId;
	}
}

BcSurface& BcMapper::GetBc(size_t bcId)
{
	if(bcId >= _bcSurfaces.size())
		exceptions::ThrowCollectionOutOfBounds("_bcSurfaces");
	return _bcSurfaces[bcId];
}

size_t BcMapper::GetBcParamsId(size_t bcId)
{
	if(_bcParams.find(bcId)==_bcParams.end())
		exceptions::ThrowCollectionOutOfBounds("_bcParams");
	return _bcParams[bcId];
}

