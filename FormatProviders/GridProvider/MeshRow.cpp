
#include "MeshRow.h"

MeshRow::MeshRow(int nr)
{
	n=nr;
	Row=new int[n];
	for(int i=0;i<nr;i++)
		Row[i]=0;
}

MeshRow::~MeshRow()
{
	if(n>0)
		delete[] Row;
}

// ���������� ��������� [].
const int& MeshRow::operator[]
	(
		int index
	)	const
{
	return Row[index];
}

// ���������� ��������� [].
int& MeshRow::operator[]
	(
		int index
	)
{
	return Row[index];
}