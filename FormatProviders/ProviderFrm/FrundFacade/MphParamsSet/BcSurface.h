#ifndef BC_SURFACE_H

#define BC_SURFACE_H


#include "ParamsSetBase.h"

#include <vector>
#include <string>


using std::vector;
using std::string;


/** Граничная поверхность, которой соответствует
* набор геометрических поверхностей в CAD геометрии
* и один конкретный тип граничных условий, 
* привязываемый через BcMapper
*
* @author Getmanskiy Victor
*/
class BcSurface : public ParamsSetBase
{
	// номера граничных поверхностей
	vector<string> _bcSurfacesIds;
	// активация
	bool _isEnabled;

public:

	BcSurface()
		:
			ParamsSetBase(),
			_isEnabled()
	{
	}

	BcSurface
		(
			const string& name
		)
		:
			ParamsSetBase(name),
			_isEnabled()
	{
	}

	void Enable() {_isEnabled = true;}
	void Disable() {_isEnabled = false;}
	bool IsEnabled() const {return _isEnabled;}

	// очистка списка поверностей
	void Clear();
	// добавить
	void AddSurface(string surfaceId);

	/** получить id поверхности по ее номеру в векторе
	* @param id - индекс в векторе
	* @return id поверхности
	*/
	string GetSurface(size_t id) const;
	// размер
	size_t CountSurfaces() const;
	
	vector<string>& GetSurfaceIds() {return _bcSurfacesIds;}
#pragma region accessors
	vector<string>& GetBcSurfacesIds() {return _bcSurfacesIds;}
#pragma endregion

#pragma region overriden	
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);
#pragma endregion
};


#endif // BC_SURFACE_H