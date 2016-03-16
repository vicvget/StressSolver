#ifndef SOLVER_GRID_PARAMS_H

#define SOLVER_GRID_PARAMS_H


#include "ParamsSetBase.h"
#include "BcMapper.h"
#include "../BasicTypes.h"
#include "../FGeometryPoint.h"

#include <utility>
#include <vector>


using std::pair;
using std::vector;


typedef pair<size_t, FGeometryPoint> IndexPointPair;
typedef vector <IndexPointPair> PointsMapper;


/** Набор параметров сетки,
* общих для всех междисциплинарных решателей
*
* @author Getmanskiy Victor
*/
class SolverGridParams: public ParamsSetBase
{
protected:
	// имя CAD файла
	string _pathCAD;
	// имя Grid файл
	string _pathGrid;
	// шаг сетки
	double _gridStep;
	// номер набора параметров граничных условий из заголовка (фильтруются по типу решателя)
	size_t _bcMapperId;
	// набор граничных условий и параметров (соответствие)
	vector<BcMapper> _bcMappers;

	// совмещение по 3 точкам с CAD-геометрией
	// вектор состоит из N (0-3) точек, выбранных для совмещения
	// если точек 0, то совмещение не требуется
	// map устроен так, что ключом является индекс точки в массиве transformedNodes в FGeometry
	// key0 - point0 - первая пара точек, которая сопоставляется (параллельный перенос)
	// key1 - point1 - вторая пара точек, дающая сопоставление по линии (поворот в плоскости)
	// key2 - point2 - третья пара точек, которая используется для произвольной 3D трансформации
	//map<size_t, FGeometryPoint> _transformPointsMapper;
	PointsMapper _transformPointsMapper;
	bool _isMirrored;

//	vector<FGeometryPoint> _transformPoints;
//	vector<int> _transformPointIds;
public:
	SolverGridParams(const string &name);
	SolverGridParams();
	SolverGridParams(const SolverGridParams& gridParams):
		ParamsSetBase(gridParams),
		_pathCAD(gridParams._pathCAD),
		_pathGrid(gridParams._pathGrid),
		_gridStep(gridParams._gridStep),
		_bcMapperId(gridParams._bcMapperId),
		_bcMappers(gridParams._bcMappers),
		_transformPointsMapper(gridParams._transformPointsMapper),
		_isMirrored(gridParams._isMirrored)
	{
	};
	
	~SolverGridParams() {};
	
	/** Получить вектор размерностей границ
	*/
	void GetNodesCountsInBcs(vector<size_t> &fullBoundaryNodesCounts);

	/** Создает маппер и делает его текущим
	* @return маппер
	*/
	BcMapper& CreateNewBcMapper();
	
	size_t CountTransformPoints() const {return _transformPointsMapper.size();}
	
	/** Получить точку для совмещения геометрии
	* @param id - id соответствующей точки из FGeometry.transformedNodes
	* @return точка
	*/
	//FGeometryPoint& GetPoint(size_t id) {return _transformPointsMapper[id];}
	
	/** Добавить точку для совмещения CAD геометрии
	* @param id - id соответствующей точки из FGeometry.transformedNodes
	* @param point - точка на CAD геометрии
	*/
	void AddPoint(int id, FGeometryPoint& point)
	{
		_transformPointsMapper.push_back(IndexPointPair(id, point));
	}

	/** Получить текущий маппер граничных условий
	* return текущий маппер граничных условий
	*/
	BcMapper& GetCurrentMapper() {return _bcMappers[_bcMapperId];};
	
	/** Установить текущий маппер
	* @param mapperId - номер маппера в _bcMappers
	*/
	void SetCurrentMapper(size_t mapperId) {_bcMapperId = mapperId;};

	size_t CountBcMappers() const {return _bcMappers.size();};
#pragma region accessors
	const string& PathCAD() const {return _pathCAD;}
	void PathCAD(string val){_pathCAD = val;}
	const string& PathGrid() const {return _pathGrid;}
	void PathGrid(string val){_pathGrid = val;}
	double GridStep() const {return _gridStep;}
	void GridStep(double val){_gridStep = val;}
	size_t GetBcMapperId() const {return _bcMapperId;}
	void SetBcMapperId(int bcMapperId) {_bcMapperId = bcMapperId;}
	vector<BcMapper> GetBcMapper() const {return _bcMappers;}
	void SetBcMapper(vector<BcMapper> val) { _bcMappers = val; }

	bool IsMirrored() const {return _isMirrored;}
	void SetIsMirrored(bool isMirrored) {_isMirrored = isMirrored;}


	PointsMapper TransformPointsMapper() const { return _transformPointsMapper; }
	void TransformPointsMapper(PointsMapper val) { _transformPointsMapper = val; }
#pragma endregion

#pragma region overriden	
	void Init();
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);
#pragma endregion
};


#endif // SOLVER_GRID_PARAMS_H