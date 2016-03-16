#include "BcSurface.h"

#include "../../../../fcore/wrappers/FileRoutines.h"


// очистка списка поверностей
void BcSurface::Clear()
{
	_bcSurfacesIds.clear();
}
// добавить
void BcSurface::AddSurface(string surfaceId)
{
	_bcSurfacesIds.push_back(surfaceId);
}

// получить поверхность по id
string BcSurface::GetSurface(size_t id) const
{
	return _bcSurfacesIds[id];
}

// размер
size_t BcSurface::CountSurfaces() const
{
	return _bcSurfacesIds.size();
}

void BcSurface::Save(ofstream& ofs) const
{	
	DefaultSave(ofs);
	int isEnabled = _isEnabled ? 1 : 0;
	fs::WriteCommentedLine(ofs, isEnabled, "flag enabled");
	//ofs << isEnabled << std::endl;
	ofs << CountSurfaces() << std::endl;
	for(int i = 0; i < CountSurfaces(); i++)
	{
		ofs << _bcSurfacesIds[i] << std::endl;
	}
}

void BcSurface::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	size_t tmpCount;
	fs::ReadLineValue(ifs, _isEnabled);
	fs::ReadLineValue(ifs, tmpCount);
	_bcSurfacesIds.resize(tmpCount);
	for(size_t i = 0; i < tmpCount; i++)
	{
		fs::ReadLineString(ifs, _bcSurfacesIds[i]);
	}
}
