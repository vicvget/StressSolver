#ifndef FGEOMETRY_H

#define FGEOMETRY_H


#include "BasicTypes.h"
#include "FElement.h"
#include "FGeometryPoint.h"
#include "TransformMatrix.h"
#include "FTransform.h"

//#include "../FortranDllsHandling/AvmodelInterface.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания геометрии тел модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FGeometry : public FElement
{
	string _code60;
	// Номер узла с центром масс (номер в массиве nodes 0-based);
	int _cmNodeNumber;

	vector<FGeometryPoint>	_nodes;
	vector<FGeometryLink>	_links;
	vector<FGeometryCircle>	_circles;
	vector<FGeometryPoint>	_transformedNodes;
	vector<FGeometryLink>	_transformedLinks;

	float _transformMtx[12];
	bool _isTransformed;

public:
	FGeometry();
	FGeometry(ifstream& stream);
	FGeometry(ifstream& stream, float scale);
	FGeometry(const char* path){Load(path);};

	// См. комментарии в базовом классе
#pragma region overriden

	void Load(const char* path);
	//int  Load(ifstream& stream);
	//void Save(ofstream& ofs) const;

	//void Output(ofstream& stream, int stage);
	//void SaveByIndex() const;

	// Загрузка геометрии из .elg формата
	void Load(ifstream& stream, float scale);

	// Загрузка геометрии из .elg формата
	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	// Сохранение геометрии в .elg формат
	void Save(ofstream& ofs) const;
	// Сохранение геометрии в .lba формат
	void OutputLba(ofstream& stream) const;

#pragma endregion


	void AddNode(const FGeometryPoint& node);
	void AddLink(const FGeometryLink& link);
//	void SetCmNode(int id) {cmNodeNumber = id;}

	const string& Code60() const
	{
		return _code60;
	}

	void SetCode60(const string& val)
	{
		_code60 = val;
	}
	
	// Преобразование поворота и сдвига (первые 3 компонента - смещение, остальные 9 - матрица поворота)
	void GlobalTransform
		(
			const float (&transform)[12]
		);

	const vector<FGeometryPoint>& Nodes() const
	{
		return _nodes;
	}

	const vector<FGeometryCircle>& Circles() const
	{
		return _circles;
	}
	
	// Генерация окружностеи (в массивы points и links)	
	void GenerateAllNodes();

	/** Получение узла
	* @param id - номер узла
	* @return узел
	*/
	const FGeometryPoint& GetNode(size_t id) const
	{
		return _nodes[id];
	}

	/** Получение стержня
	* @param id - номер стержня
	* @return стержень
	*/
	const FGeometryLink& GetLink(size_t id) const
	{
		return _links[id];
	}

	/** Получение окружности
	* @param id - номер окружности
	* @return окружность
	*/
	const FGeometryCircle& GetCircle(size_t id) const
	{
		return _circles[id];
	}

	const vector<FGeometryLink>& Links() const
	{
		return _links;
	}

	const vector<FGeometryPoint>& TransformedNodes() const
	{
		return _transformedNodes;
	}

	const vector<FGeometryLink>& TransformedLinks() const
	{
		return _transformedLinks;
	}

	float (&GetTransformArray())[12]
	{
		return _transformMtx;
	}

	const float (&GetTransformArray() const)[12]
	{
		return _transformMtx;
	}

	void GetTransformArray
		(
			float (&mtx)[12]
		)	const
	{
		std::copy(_transformMtx, _transformMtx + 12, mtx);
	}

	void SetTransformArray
		(
			const float (&mtx)[12]
		)
	{
		memcpy(_transformMtx, mtx, 12 * sizeof(float));
		_isTransformed = true;
	}

	TransformMatrix GetTransformMatrix() const
	{
		// TODO: фигурные скобки нужны, чтобы код компилировался в g++
		return {_transformMtx};
	}

	TransformationMatrix GetTransformationMatrix() const
	{
		// TODO: фигурные скобки нужны, чтобы код компилировался в g++
		return {_transformMtx};
	}

	void SetTransformationMatrix
		(
			const TransformationMatrix& mtx
		)
	{
		mtx.Export(_transformMtx);
		_isTransformed = true;
	}

	FTransform GetTransform() const;


	void GetBound(MathHelpers::Vec3& minPoint, MathHelpers::Vec3& maxPoint) const;
	void GetPoints(std::vector<MathHelpers::Vec3>& points, bool toGCS=false) const;

	bool IsTransformed() const
	{
		return _isTransformed;
	}

	void SetTransformed(bool val){ _isTransformed = val; }

	int CmNodeNumber() const {return _cmNodeNumber;}
	void SetCmNodeNumber(int val){_cmNodeNumber = val;}

	const FGeometryPoint& CmNode() const {return _nodes[_cmNodeNumber];}

	const FGeometryPoint& GetTransformedNode(size_t id) const
	{
		return _transformedNodes[id];
	}
	const FGeometryPoint& TransformedCmNode() const {return _transformedNodes[_cmNodeNumber];}
	const FGeometryLink& GetTransformedLink(size_t id) const {return _transformedLinks[id];}

	/** Возвращает узел с учетом локальной трансформации (для начальных условий)
	*/
	FGeometryPoint GetLocalTransformedNode(size_t id) const;


	void SetNodes(vector<FGeometryPoint> val){_nodes = val;}
	void SetLinks(vector<FGeometryLink> val){_links = val;}
	void SetCircles(vector<FGeometryCircle> val){_circles = val;}
	void SetTransformedNodes(vector<FGeometryPoint> val){_transformedNodes = val;}
	void SetTransformedLinks(vector<FGeometryLink> val){_transformedLinks = val;}

	size_t CountLinks() const {return _links.size();};
	size_t CountTransformedLinks() const {return _transformedLinks.size();};
	size_t CountNodes() const {return _nodes.size();};
	size_t CountTransformedNodes() const {return _transformedNodes.size();};
		
	bool Combine(int* ids, FGeometryPoint* points, int count);

	Vec3 GetGeometryPoint(int pointId) const;

	// TODO: вынести в Mat3
	static
	void AnglesToMtx
		(
			double a,
			double b,
			double c,
			double* mtx
		);
};


#endif // FGEOMETRY_H
