#ifndef MESHROW_H
#define MESHROW_H
//#pragma GCC visibility push(hidden)
/*Класс служит для хранения и работы с одной строчкой элементов*/
class MeshRow
{
protected:
	//fields
	//количество элементов в строке
	int n;
	//элементы строки
	int *Row;
	//methods

public:

	//
	MeshRow
		(
			int nr
		);

	MeshRow(const MeshRow&) = delete;

	~MeshRow();

	//methods
	//вернуть количество элементов
	int GetNRow() const {return n;}
	//установить 1
	void Set1To(int i){Row[i]=1;}
	//установить 0
	void Set0To(int i){Row[i]=0;}
	//вернуть значение элемента i
	int GetVRow(int i) const {return Row[i];}

	// Перегрузка оператора [].
	const int& operator[]
		(
			int index
		)	const;

	// Перегрузка оператора [].
	int& operator[]
		(
			int index
		);

	//fields
};
//#pragma GCC visibility pop
#endif
