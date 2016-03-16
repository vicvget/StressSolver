#include "ArrayVertexToGridFormatter.h"

#include <cmath>


ArrayVertexToGridFormatter::ArrayVertexToGridFormatter
	(
		const vector<Vertex*>& myArray
	)
	:
		_array(myArray)
{
	
}
ArrayVertexToGridFormatter::ArrayVertexToGridFormatter
	(
		const vector<Vertex*>& myArray,
		const vector<BoundaryCondition*>& bc
	)
	:
		_array(myArray),
		_bc(bc)
{
}

/**
* Функция преобразовывает массив указателей объектов типа Vertex в сетку
* @return возвращет 1 если преобразование успешно
*/
int ArrayVertexToGridFormatter::Transform()
{
	//вычисление размерности сетки, размера элементов и координат начального элемента
	GetSize();

	_grid=new OccRectilinearGrid(nx,ny,nz);
	_grid->SetElementSize(_dx,_dy,_dz);
	_grid->SetStartPoint(startPointX,startPointY,startPointZ);

	for (Vertex* tmpPoint : _array)
	{
		int i = (int)((tmpPoint->pt.x-startPointX)/_dx);
		int j = (int)((tmpPoint->pt.y-startPointY)/_dy);
		int k = (int)((tmpPoint->pt.z-startPointZ)/_dz);

		if((i >= nx) || (j >= ny) || (k >= nz))
		{
			std::cout << "ERROR: Node outside grid: (" << i << ' ' << j << ' ' << k << ')' << std::endl;
		}
		else
		{
			(*_grid)[i][j][k] = 1;
		}
	}
	//TODO: bcs

	return 1;
}

/**
* Функция возвращает размер элемента сетки
* @return возвращет 1 если размер элемента найден
*/
int ArrayVertexToGridFormatter::GetSize()
{
	double spacingsx = 1E10;
	vector<Vertex *>::iterator it = _array.begin();
	Vertex *tmpPoint = *it;

	++it;
	const double epsilon = 1e-12;
	while(it != _array.end())
	{
		double ddx = fabs((*it)->pt.x - tmpPoint->pt.x);
		double ddy = fabs((*it)->pt.y - tmpPoint->pt.y);
		double ddz = fabs((*it)->pt.z - tmpPoint->pt.z);
		if((ddx > epsilon)&&(ddx < spacingsx))
				spacingsx = ddx;
		if((ddy > epsilon)&&(ddy < spacingsx))
				spacingsx = ddy;
		if((ddz > epsilon)&&(ddz < spacingsx))
				spacingsx = ddz;
		tmpPoint = *it; 
		++it;
	}

	_dx = _dy = _dz = spacingsx;
	double sx,sy,sz;
	GetBBox(&sx,&sy,&sz);
	nx = (int)(sx / _dx + 1);
	ny = (int)(sy / _dy + 1);
	nz = (int)(sz / _dz + 1);
	return 1;
}
/**
* Функция возвращает размер сетки
* @return возвращет 1 если размер найден
*/
int ArrayVertexToGridFormatter::GetBBox(double *sx, double *sy, double *sz)
{
	size_t i = 0;
	size_t nPoints = _array.size();
	// TODO: переписать код так, чтобы не выдавалось предупреждений
	// о возможности использования непроинициализированных min_/max_
	double maxx{}, maxy{}, maxz{};
	double minx{}, miny{}, minz{};
	Vertex *tmpPoint;

	while (i < nPoints)
	{
		tmpPoint = _array.at(i);
		if (i == 0)
		{
			minx = maxx = tmpPoint->pt.x;
			miny = maxy = tmpPoint->pt.y;
			minz = maxz = tmpPoint->pt.z;
		}
		if (tmpPoint->pt.x < minx) minx = tmpPoint->pt.x;
		if (tmpPoint->pt.y < miny) miny = tmpPoint->pt.y;
		if (tmpPoint->pt.z < minz) minz = tmpPoint->pt.z;
		if (tmpPoint->pt.x > maxx) maxx = tmpPoint->pt.x;
		if (tmpPoint->pt.y > maxy) maxy = tmpPoint->pt.y;
		if (tmpPoint->pt.z > maxz) maxz = tmpPoint->pt.z;
		i++;
	}
	*sx = maxx - minx;
	*sy = maxy - miny;
	*sz = maxz - minz;
	startPointX = minx;
	startPointY = miny;
	startPointZ = minz;
	return 1;
}
