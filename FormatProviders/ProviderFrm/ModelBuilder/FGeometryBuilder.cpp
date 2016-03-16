#include "FGeometryBuilder.h"

#include "FElementBuilder.h"
#include "../FrundFacade/FElementType.h"


using MathHelpers::MakeVec3;


const FGeometry& FGeometryBuilder::Get() const
{
	return _fGeometry;
}

int FGeometryBuilder::AddPoint(const Vec3& point)
{
	const double multiplier = 1e3;
	_fGeometry.AddNode
		(
			FGeometryPoint
				(
					static_cast<float>(point.X() * multiplier),
					static_cast<float>(point.Y() * multiplier),
					static_cast<float>(point.Z() * multiplier)
				)
		); // stored in mm

	return _fGeometry.CountNodes();
}

int FGeometryBuilder::AddCmPoint(const Vec3& point)
{
	_fGeometry.SetCmNodeNumber(_fGeometry.CountNodes()); // 1-based
	return AddPoint(point);
}

void FGeometryBuilder::AddLink(int pointId1, int pointId2)
{
	_fGeometry.AddLink(FGeometryLink(pointId1+1, pointId2+1)); // 1-based
}

void FGeometryBuilder::SetTransform(const Mat3& mtxTransform, const Vec3& cmNode)
{
	float mtx[12];
	cmNode.Export(mtx);
	mtxTransform.Export(mtx+3);
	_fGeometry.SetTransformArray(mtx);
}

void FGeometryBuilder::Setup(const string& name, int id, int number)
{
	FElementBuilder::SetupCommonBodyProperties(_fGeometry, name, id, number,
		FElementType::MakeFileId(FElementType::FETC_Body, id+1));
	
}

void FGeometryBoundBuilder::BuildParallelepiped(const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Vec3* points)
{
	Vec3 midPoint = (maxCoordsPoint + minCoordsPoint) * 0.5;
	Vec3 halfLength = (maxCoordsPoint - minCoordsPoint).Abs() * 0.5;
	points[0].Init(midPoint.X() - halfLength.X(), midPoint.Y() - halfLength.Y(), midPoint.Z() - halfLength.Z());
	points[1].Init(midPoint.X() + halfLength.X(), midPoint.Y() - halfLength.Y(), midPoint.Z() - halfLength.Z());
	points[2].Init(midPoint.X() - halfLength.X(), midPoint.Y() + halfLength.Y(), midPoint.Z() - halfLength.Z());
	points[3].Init(midPoint.X() + halfLength.X(), midPoint.Y() + halfLength.Y(), midPoint.Z() - halfLength.Z());
	points[4].Init(midPoint.X() - halfLength.X(), midPoint.Y() - halfLength.Y(), midPoint.Z() + halfLength.Z());
	points[5].Init(midPoint.X() + halfLength.X(), midPoint.Y() - halfLength.Y(), midPoint.Z() + halfLength.Z());
	points[6].Init(midPoint.X() - halfLength.X(), midPoint.Y() + halfLength.Y(), midPoint.Z() + halfLength.Z());
	points[7].Init(midPoint.X() + halfLength.X(), midPoint.Y() + halfLength.Y(), midPoint.Z() + halfLength.Z());
	
	AddLink(0,1);
	AddLink(2,3);
	AddLink(4,5);
	AddLink(6,7);

	AddLink(0,2);
	AddLink(1,3);
	AddLink(4,6);
	AddLink(5,7);

	AddLink(0,4);
	AddLink(1,5);
	AddLink(2,6);
	AddLink(3,7);
}

void FGeometryBoundBuilder::BuildMarkersBound(int countPoints)
{
	for(int i = 0 ; i < countPoints ; i++)
	{
		AddLink(i, i+1);
	}
}

void FGeometryBoundBuilder::Setup(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint)
{
	FGeometryBuilder::Setup(name, id, number);
	
	Vec3 points[8];		
	BuildParallelepiped(minCoordsPoint, maxCoordsPoint, points);
	Vec3 midPoint = (maxCoordsPoint + minCoordsPoint) * 0.5;
	AddCmPoint(midPoint);
	
	for(int i = 0; i < 8; i++)
	{
		AddPoint(points[i]);
	}	
}

void FGeometryBoundBuilder::SetupTransformed(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, const Mat3& mtxTransform, const Vec3& cmNode)
{
	FGeometryBuilder::Setup(name, id, number);

	Vec3 points[8];		
	BuildParallelepiped(minCoordsPoint, maxCoordsPoint, points);
	for(int i = 0; i < 8; i++)
	{
		points[i] = mtxTransform.Tmul(points[i]-cmNode);
		AddPoint(points[i]);
	}
	AddCmPoint(cmNode);
	SetTransform(mtxTransform, cmNode);	
}

void FGeometryBoundBuilder::SetupUntransformed(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Mat3 mtxTransform, Vec3 cmNode)
{
	FElementBuilder::SetupCommonBodyProperties(_fGeometry, name, id, number,
		FElementType::MakeFileId(FElementType::FETC_Body, id+1));

	Vec3 points[8];
	BuildParallelepiped(minCoordsPoint, maxCoordsPoint, points);
	for(int i = 0; i < 8; i++)
	{
		points[i] = mtxTransform*points[i] + cmNode;
		AddPoint(points[i]);
	}
	AddCmPoint(cmNode);
	SetTransform(mtxTransform, cmNode);	
}

void FGeometryBoundBuilder::SetupTransformedBound(const string& name, int id, int number, const Vec3& minCoordsPoint, const Vec3& maxCoordsPoint, Mat3 mtxTransform, Vec3 cmNode)
{
	FElementBuilder::SetupCommonBodyProperties(_fGeometry, name, id, number,
		FElementType::MakeFileId(FElementType::FETC_Body, id+1));

	Vec3 points[8];		
	Vec3 minCoordsPointUntransformed = mtxTransform.Tmul(minCoordsPoint-cmNode);
	Vec3 maxCoordsPointUntransformed = mtxTransform.Tmul(maxCoordsPoint-cmNode);
	
	BuildParallelepiped(minCoordsPointUntransformed, maxCoordsPointUntransformed, points);
	for(int i = 0; i < 8; i++)
	{
		AddPoint(points[i]);
	}
	AddCmPoint(MakeVec3(0.,0.,0.));
	SetTransform(mtxTransform, cmNode);	
}

void FGeometryBoundBuilder::SetupMarkersBound(const string& name, int id, int number, vector<Vec3> points, Mat3 mtxTransform, Vec3 cmNode)
{
	FElementBuilder::SetupCommonBodyProperties(_fGeometry, name, id, number,
		FElementType::MakeFileId(FElementType::FETC_Body, id+1));

	for(size_t i = 0 ; i < points.size() ; i++)
		points.at(i) = mtxTransform.Tmul(points.at(i) - cmNode);

	BuildMarkersBound(points.size());

	for(size_t i = 0; i < points.size(); i++)
	{
		AddPoint(points.at(i));
	}
	AddCmPoint(MakeVec3(0.,0.,0.));
	SetTransform(mtxTransform, cmNode);	
}