#include "FElementBuilder.h"
#include "FBodyBuilder.h"
#include "../FrundFacade/FElementType.h"
#include "../../../fcore/wrappers/StringRoutines.h"

const FBody& FBodyBuilder::Get() const
{
	return _fBody;
}

void FFreeBodyBuilder::Setup(const string& name, int id, double mass, Vec3 inertia, int geometryId)
{
	int number = id < 48 ? id + 1 : id + 2;
	FElementBuilder::SetupCommonBodyProperties(_fBody, name, id, number, FElementType::MakeFileId(FElementType::FETC_Body, id+1));

	_fBody.M(mass);
	_fBody.Jx(inertia.X());
	_fBody.Jy(inertia.Y());
	_fBody.Jz(inertia.Z());

	_fBody.Strm(NumberToString(mass));
	_fBody.Strjx(NumberToString(inertia.X()));
	_fBody.Strjy(NumberToString(inertia.Y()));
	_fBody.Strjz(NumberToString(inertia.Z()));

	_fBody.X(1);
	_fBody.Y(1);
	_fBody.Z(1);
	_fBody.Rx(3);
	_fBody.Ry(3);
	_fBody.Rz(3);

	//FGeometry* pGeometry = new FGeometry(geometry); // TODO: утечка памяти, переделат на индексы
	_fBody.GeometryId(geometryId);
}

void FFixedBodyBuilder::Setup(const string& name, int id, double mass, Vec3 inertia, int geometryId)
{
	int number = id < 48 ? id + 1 : id + 2;
	FElementBuilder::SetupCommonBodyProperties(_fBody, name, id, number,
		FElementType::MakeFileId(FElementType::FETC_Body, id+1));

	_fBody.M(mass);
	_fBody.Jx(inertia.X());
	_fBody.Jy(inertia.Y());
	_fBody.Jz(inertia.Z());

	_fBody.Strm(NumberToString(mass));
	_fBody.Strjx(NumberToString(inertia.X()));
	_fBody.Strjy(NumberToString(inertia.Y()));
	_fBody.Strjz(NumberToString(inertia.Z()));

	_fBody.X(0);
	_fBody.Y(0);
	_fBody.Z(0);
	_fBody.Rx(0);
	_fBody.Ry(0);
	_fBody.Rz(0);

	//FGeometry* pGeometry = new FGeometry(geometry); // TODO: утечка памяти, переделат на индексы
	_fBody.GeometryId(geometryId  );
}