#ifndef MESHLAYER_H
#define MESHLAYER_H

//#pragma GCC visibility push(hidden)
#include "MeshRow.h"

#include <vector>
#include <cstddef>


using std::vector;


/*Класс служит для хранения и работы со слоем сетки */
class MeshLayer
{
protected:
	//fields
	//строки
	vector<MeshRow *> Layer;
	//methods
public:
	//
	MeshLayer();
	MeshLayer(int nl,int nr);
	~MeshLayer(void);
	//methods
	//вернуть количество строк
	size_t GetNLayer() const {return Layer.size();}
	//Добавить строку
	void AddRow(MeshRow *row);

	// Перегрузка оператора []
	const MeshRow& operator[]
		(
			int index
		)	const;

	// Перегрузка оператора []
	MeshRow& operator[]
		(
			int index
		);

	//fields
};
//#pragma GCC visibility pop

#endif
