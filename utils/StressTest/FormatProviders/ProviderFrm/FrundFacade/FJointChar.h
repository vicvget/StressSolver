#ifndef FJOINT_CHAR_H

#define FJOINT_CHAR_H


#include "FElement.h"

#include "../../../fcore/Calculator/Calculator.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания характеристики шарнира модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FJointChar : public FElement
{
private:
	vector<CharComponent> _springParams;
	// параметры диссипативной характеристики
	vector<CharComponent> _dampingParams;
	// ограничение по мах жесткости
	mutable double _maxStiffness;

public:

	FJointChar();

	FJointChar(int orderIndex);

	FJointChar
		(
			const char* path
		)
		:
			_maxStiffness()
	{
		FElement::Load(path);
	}

	FJointChar
		(
			ifstream& ifs
		)
		:
			_maxStiffness()
	{
		Load(ifs);
	}

	void AddCharComponent(const CharComponent& charSpringComponent, const CharComponent& charDampingComponent);

	void AddCharComponent(const CharComponent& charSpringComponent, const CharComponent& charDampingComponent, int direction)
	{
		CharComponent spring = charSpringComponent;
		CharComponent damping = charDampingComponent;
		spring.direction = direction;
		damping.direction = direction;
		_springParams.push_back(spring);
		_dampingParams.push_back(damping);
	}



	const vector<CharComponent>& GetSpringParams() const
	{
		return _springParams;
	}

	void GetSpringParamsType(int dof, int& type, size_t& subtype) const
	{
		type = _springParams[dof].type;
		subtype = _springParams[dof].params.size();
	}

	string GetSpringParam(int dof, int paramId) const
	{
		return _springParams[dof].params[paramId];
	}

	const vector<CharComponent>& GetDampingParams() const
	{
		return _dampingParams;
	}
	

	void SetSpringParams(vector<CharComponent> val) { _springParams = val; }
	void SetDampingParams(vector<CharComponent> val) { _dampingParams = val; }
	void SetMaxStiffness(double val) const { _maxStiffness = val; }

// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	void Save(ofstream& ofs) const;

	virtual
	void Output
		(
			ofstream& stream,
			int stage
		)	const override;

#pragma endregion

	/**
	* Заменить значения выражений для жесткости их численными значениями
	* @param calc - калькулятор
	*/
	void Eval(Calc::Calculator &calc) const;

};


#endif // FJOINT_CHAR_H