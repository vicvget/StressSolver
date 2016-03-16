#ifndef FGEOMETRY_BUILDER_H

#define FGEOMETRY_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FGeometry.h"
#include "../../../AdditionalModules/fmath/Matrix3x3.h"

#include <memory>


using MathHelpers::Vec3;
using MathHelpers::Mat3;


class FGeometryBuilder
{
protected:	
	FGeometry _fGeometry;

public:
	FGeometryBuilder(const FGeometry& fGeometry): _fGeometry(fGeometry){}
	FGeometryBuilder(){}

	virtual const FGeometry& Get() const;
	
	/**
	* ƒобавл€ет точку и возвращает ее номер (1-based)
	*/
	int AddPoint(const Vec3& point);
	int AddCmPoint(const Vec3& point);

	void AddLink(int pointId1, int pointId2);
	void SetTransform(const Mat3& mtxTransform, const Vec3& cmNode);
	void Setup(const string& name, int id, int number);
};


class FGeometryBoundBuilder: public FGeometryBuilder
{
	void BuildParallelepiped(const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Vec3* points);
	void BuildMarkersBound(int countPoints);
public:
	void Setup(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint);
	
	// «адает ограничивающий параллелепипед точками, координаты которых уже преобразованы
	// так как вращение всегда происходит относительно ц.м., а транспонированна€ матрица поворота соответствует обратной,
	// преобразование координат P перед добавлением выполн€етс€ как mtxTransform^T * (P - cmNode);
	void SetupTransformed(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, const Mat3& mtxTransform, const Vec3& cmNode);

	// «адает ограничивающий параллелепипед точками, координаты которых не преобразованы
	// так как вращение всегда происходит относительно ц.м.
	void SetupUntransformed(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Mat3 mtxTransform, Vec3 cmNode);

	// «адает ограничивающий параллелепипед точками с предварительным обратным преобразованием
	// так как вращение всегда происходит относительно ц.м.
	void SetupTransformedBound(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Mat3 mtxTransform, Vec3 cmNode);

	//ƒобавление Linkов дл€ Markerов
	void SetupMarkersBound(const string& name, int id, int number, vector<Vec3> points, Mat3 mtxTransform, Vec3 cmNode);
};

typedef std::unique_ptr<FGeometryBuilder> PGeometryBuilder;


#endif // FGEOMETRY_BUILDER_H