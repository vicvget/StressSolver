#ifndef MODEL_APPEARANCE_H

#define MODEL_APPEARANCE_H


#include "../../../ExternalModules/TinyXml/tinyxml.h"

#include <map>
#include <string>


/*
* Структура, содержащая компоненты цвета
*/
class Color final
{
public:

	// Константы

	// количество компонент
	static
	const size_t ComponentsCount = 4;


	// Типы

	// компонента цвета
	using Component = unsigned char;

	// массив компонент цвета
	using Components = Component[ComponentsCount];

	// вещественное число в диапазоне [0, 1], используемое в VTK в качестве компоненты цвета
	using VTKComponent = double;

	// массив значений, используемых в VTK в качестве компонент цвета
	using VTKComponents = VTKComponent[ComponentsCount];


	// Константы

	// максимальное значение компоненты
	static
	const Component MaxComponentValue;


	// Конструкторы и деструкторы

	/**
	* Конструктор по умолчанию
	*/
	Color();

	/**
	* Конструктор
	* @param alpha - прозрачность
	* @param red - красная компонента
	* @param green - зеленая компонента
	* @param blue - синяя компонента
	*/
	Color
		(
			Component alpha,
			Component red,
			Component green,
			Component blue
		);

	/**
	* Конструктор
	* @param components - массив компонент цвета
	*/
	Color
		(
			const Components& components
		);

	/**
	* Конструктор
	* @param componentsVTK - массив вещественных значений в диапазоне [0, 1],
	* используемых в VTK в качестве компонент цвета
	*/
	Color
		(
			const VTKComponents& componentsVTK
		);


	// Селекторы

	/**
	* Получить прозрачность
	* @return прозрачность
	*/
	Component GetAlpha() const;

	/**
	* Получить непрозрачность цвета
	* @return непрозрачность цвета
	*/
	VTKComponent GetOpacity() const;

	/**
	* Получить массив компонент цвета
	* @return массив компонент цвета
	*/
	const Components& GetComponents() const;

	/**
	* Получить массив значений в диапазоне [0, 1], используемых в VTK в качестве компонент цвета
	* @param componentsVTK - массив значений в диапазоне [0, 1], используемых в VTK в качестве компонент цвета
	*/
	void GetVTKComponents
		(
			VTKComponents& componentsVTK
		)	const;


	// Чтение из файла

	/**
	* Чтение структуры с компонентами цвета из xml-документа
	* @param colorElement - элемент xml-документа, содержащий структуру с компонентами цвета
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlElement
		(
			const TiXmlElement* colorElement
		);


private:

	// компоненты цвета
	union
	{

		// набор компонент
		struct
		{

			// прозрачность
			Component _alpha;

			// красная компонента
			Component _red;

			// зеленая компонента
			Component _green;

			// синяя компонента
			Component _blue;

		};

		// массив компонент
		Components _components;

	};


	// Вспомогательные функции

	/*
	* Сформировать цвет из строки, хранящей представление цвета в HTML-формате
	* @param htmlColor - строка, хранящая представление цвета в HTML-формате
	* @return признак успешного (true) или неуспешного (false) формирования
	*/
	bool FromHtmlColor
		(
			const std::string& htmlColor
		);

	/**
	* Преобразовать значение компоненты цвета в вещественное число в диапазоне [0, 1] (используется в VTK)
	* @param component - значение компоненты цвета
	* @return соответствующее вещественное число в диапазоне [0, 1]
	*/
	static
	VTKComponent FromColorComponentToVTK
		(
			Component component
		);

	/**
	* Преобразовать вещественное число в диапазоне [0, 1] (используется в VTK) в значение компоненты цвета
	* @param componentVTK - вещественное число в диапазоне [0, 1]
	* @return соответствующее значение компоненты цвета
	*/
	static
	Component FromVTKToColorComponent
		(
			VTKComponent componentVTK
		);

	/**
	* Преобразовать массив компонент цвета в массив вещественных чисел в диапазоне [0, 1] (используются в VTK)
	* @param components - массив компонент цвета
	* @param componentsVTK - соответствующий массив вещественных чисел в диапазоне [0, 1]
	*/
	static
	void FromColorComponentsToVTK
		(
			const Components& components,
			VTKComponents& componentsVTK
		);

	/**
	* Преобразовать массив вещественных чисел в диапазоне [0, 1] (используются в VTK) в массив компонент цвета
	* @param componentsVTK - массив вещественных чисел в диапазоне [0, 1]
	* @param components - соответствующий массив компонент цвета
	*/
	static
	void FromVTKToColorComponents
		(
			const VTKComponents& componentsVTK,
			Components& components
		);

};


/**
* Внешний вид тела
*/
class BodyAppearance final
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	BodyAppearance() = default;

	/**
	* Конструктор
	* @param bodyNumber - номер тела
	* @param color - цвет тела
	*/
	BodyAppearance
		(
			size_t bodyNumber,
			const Color& color
		);


	// Селекторы

	/**
	* Получить номер тела
	* @return номер тела
	*/
	size_t GetBodyNumber() const;

	/**
	* Получить цвет тела
	* @return цвет тела
	*/
	const Color& GetColor() const;


	// Чтение из файла

	/**
	* Чтение внешнего вида тела из xml-документа
	* @param bodyAppearanceElement - элемент xml-документа, содержащий внешний вид тела
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlElement
		(
			const TiXmlElement* bodyAppearanceElement
		);


private:

	// номер тела
	size_t _bodyNumber{};

	// цвет тела
	Color _color;

};

// список объектов, содержащих внешний вид тел
using BodyAppearances = std::map<size_t, BodyAppearance>;


// Чтение из файла

/**
* Чтение списка объектов, содержащих внешний вид тел, из xml-документа
* @param bodyAppearancesElement - элемент xml-документа, содержащий список параметров макроса
* @param bodyAppearances - считанный список объектов, содержащих внешний вид тел
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlElement
	(
		const TiXmlElement* bodyAppearancesElement,
		BodyAppearances& bodyAppearances
	);


/**
* Внешний вид модели
*/
class ModelAppearance final
{
public:

	// Селекторы

	/**
	* Получить список объектов, содержащих внешний вид тел
	* @return список объектов, содержащих внешний вид тел
	*/
	const BodyAppearances& GetBodyAppearances() const;


	// Модификаторы

	/**
	* Добавить объект, содержащий внешний вид тела
	* @param bodyAppearance - добавляемый объект, содеражщий внешний вид тела
	*/
	void AddBodyAppearance
		(
			const BodyAppearance& bodyAppearance
		);


	// Чтение из файла

	/**
	* Чтение внешнего вида модели из xml-документа
	* @param xmlFileName - наименование xml-документа, содержащиго внешний вид модели
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlDocument
		(
			const std::string& xmlFileName
		);

	/**
	* Чтение внешнего вида модели из файла с внешним видом модели ФРУНД
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromFile();


private:

	// наименование файла с внешним видом модели
	static
	const std::string _modelAppearanceFileName;

	// список объектов, содержащих внешний вид тел
	BodyAppearances _bodyAppearances;

};


#endif // MODEL_APPEARANCE_H