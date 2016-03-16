#ifndef BC_MAPPER_H

#define BC_MAPPER_H


#include "ParamsSetBase.h"
#include "BcSurface.h"

#include <map>


using std::map;


/** Параметры граничных условий
* класс для связывания наборов параметров для разных типов границ
* с объектами BcSurface
*
* @author Getmanskiy Victor
*/
class BcMapper: public ParamsSetBase
{
	// граничные наборы поверхностей
	vector<BcSurface> _bcSurfaces;
	// соответствие BcGeometryId - BcParamsId
	map<size_t, size_t> _bcParams;
	// набор параметров для свободных граничных условий
	size_t _restBcId;
	// включение оставшихся точек в граничные условия
	bool _isRestBcEnabled;
public:

	BcMapper();

	BcMapper(const string &name);


	/** очистка списка поверностей
	*/
	void Clear();
	
	/** Добавить пустую граничную поверхность (не содержащую геометрических поверхностей)
	* и прявязать к параметрам по умолчанию (то есть к 0)
	* @param bcType тип граничных условий
	*/
	void AddBc(int bcType);
	
	/** Добавить пустую граничную поверхность (не содержащую геометрических поверхностей)
	* и прявязать к параметрам по умолчанию (то есть к 0)
	* @param bcType - тип граничных условий
	* @return граничная поверхность для добавления геомтерических поверхностей
	*/
	BcSurface& AddNewBc(int bcType);
	
	/** Получить граничную поверхность
	* @param bcId - номер граничной поверхности
	* @return граничная поверхность
	*/	
	BcSurface& GetBc(size_t bcId);

	/** получить Id параметров BC в заголовке
	* @param bcId - номер граничной поверхности
	* @return Id параметров BC в заголовке
	*/
	size_t GetBcParamsId(size_t bcId);

	/** Ассоциация параметров с набором граничных поверхностей
	*/
	void SetBcParamId(size_t surfaceId, size_t paramsId) 
	{
		_bcParams[surfaceId]=paramsId;
	}
	
	void EnableBc(int surfaceId = -1)
	{
		if(surfaceId == -1)
		{
			_isRestBcEnabled = true;
		}
		else if(surfaceId < static_cast<int>(_bcSurfaces.size()))
		{
			_bcSurfaces[surfaceId].Enable();
		}
	}

	void DisableBc(int surfaceId = -1)
	{
		if(surfaceId == -1)
		{
			_isRestBcEnabled = false;
		}
		else if(surfaceId < static_cast<int>(_bcSurfaces.size()))
		{
			_bcSurfaces[surfaceId].Disable();
		}
	}

	size_t GetRestBcId() const {return _restBcId;}
	void SetRestBcId(size_t id) {_restBcId = id;}

	bool IsRestBcEnabled() const { return _isRestBcEnabled; }
	void SetIsRestBcEnabled(bool value) {_isRestBcEnabled = value;}
//	void EnableRestBc {_isRestBcEnabled = true;}
//	void DisablRestBc {_isRestBcEnabled = false;}

	/** Удалить Bc
	*/
	void RemoveBc(size_t bcId);
	
	/** Количество поверхностей с параметрами
	*/
	size_t CountBcs() const {return _bcSurfaces.size();}

#pragma region accessors
	vector<BcSurface> GetBcSurfaces()const {return _bcSurfaces;}
	void SetBcSurfaces(vector<BcSurface> val) {_bcSurfaces = val;}
	map<size_t, size_t> GetBcParams()const {return _bcParams;}
	void SetBcParams(map<size_t, size_t> val){_bcParams = val;}
#pragma endregion

#pragma region overriden	
	void Init();
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);

#pragma endregion

};


#endif // BC_MAPPER_H
