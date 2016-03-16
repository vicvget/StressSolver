#include "SurfaceParam.h"


SurfaceParam::SurfaceParam
	(
		const string& n,
		const string& v
	)
	:
		name(n),
		var(v)
{
}

void SurfaceParam::ClearPoints()
{
	for (Point3D* child : SurfaceChildren)
	{
		delete child;
	}
	SurfaceChildren.clear();
}
//добавление граничной точки
void SurfaceParam::AddBoundaryPoints(int pntx,int pnty,int pntz)
{
	Point3D *pnt = new Point3D();
	pnt->x=pntx;
	pnt->y=pnty;
	pnt->z=pntz;
	SurfaceChildren.push_back(pnt);
}

//вернуть граничныю точку под индексом npnt
int SurfaceParam::GetBoundaryPoints(int npnt, int *x, int *y, int *z)
{
	if (npnt < static_cast<int>(SurfaceChildren.size()))
	{
		*x = SurfaceChildren.at(npnt)->x;
		*y = SurfaceChildren.at(npnt)->y;
		*z = SurfaceChildren.at(npnt)->z;
		return 1;
	}
	return -1;
}

//вернуть количество граничных точек
size_t SurfaceParam::GetNumBP() const
{
	return SurfaceChildren.size();
}

void SurfaceParam::SetIndex
	(
		int ind
	)
{
	index = ind;
}

int SurfaceParam::GetIndex() const
{
	return index;
}