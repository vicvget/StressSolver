#ifndef FFORCE_CHAR_H

#define FFORCE_CHAR_H


#include "FElement.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания характеристик нагрузок модели ФРУНД
* характеристика нагрузки для каждой нагрузки своя, оптимальнее 
* было бы реализовать принцип ссылок, как в шарнирах
*
* @author Getmanskiy Victor
*/
class FForceChar : public FElement
{
private:
	vector<CharComponent> _params;
public:
	FForceChar();
	FForceChar(int orderIndex);
	FForceChar(const char* path){FElement::Load(path);};
	FForceChar(ifstream& ifs){Load(ifs);};
	void AddCharComponent(const CharComponent& charComponent)
	{
		_params.push_back(charComponent);
	}
// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	void Save(ofstream& ofs) const;
	
	/** Возвращвет вектор из 6 компонент, true, если есть данная степень свободы
	* первые 3 поступательные, последние три - вращательные.
	* поступательные изображаются как стрелки, вращательные - как моменты (дуга вокруг оси со стрелкой)
	*/
	void GetComponentsVector(vector<bool>& components) const;

	virtual
	void Output
		(
			ofstream& stream,
			int stage
		)	const override;
	
	vector<CharComponent> Params() const { return _params; }
	const vector<CharComponent>& GetParams() const { return _params; }
	void Params(vector<CharComponent> val) { _params = val; }
#pragma endregion

};


#endif // FFORCE_CHAR_H