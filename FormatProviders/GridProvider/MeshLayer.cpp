#include "MeshLayer.h"

MeshLayer::MeshLayer()
{
}
MeshLayer::MeshLayer(int nl,int nr)
{
	for(int i=0;i<nl;i++)
	{
		Layer.push_back(new MeshRow(nr));
	}
}

MeshLayer::~MeshLayer(void)
{
	for(MeshRow* row : Layer)
	{
		delete row;
	}
}

//Добавить строку
void MeshLayer::AddRow(MeshRow *row)
{
	Layer.push_back(row);
}

const MeshRow& MeshLayer::operator[]
	(
		int index
	)	const
{
	return *(Layer[index]);
}

MeshRow& MeshLayer::operator[]
	(
		int index
	)
{
	return *(Layer[index]);
}