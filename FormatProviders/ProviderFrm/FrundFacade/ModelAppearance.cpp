#include "ModelAppearance.h"

#include "../../../fcore/exceptions/fcExceptions.h"
#include "../../../fcore/wrappers/FileRoutines.h"

#include <algorithm>


using std::string;


// class Color

// максимальное значение компоненты
// static
const Color::Component Color::MaxComponentValue = std::numeric_limits<Component>::max();

/**
* Конструктор по умолчанию
*/
Color::Color()
	:
		_alpha(),
		_red(),
		_green(),
		_blue()
{
}

/**
* Конструктор
* @param alpha - прозрачность
* @param red - красная компонента
* @param green - зеленая компонента
* @param blue - синяя компонента
*/
Color::Color
	(
		Component alpha,
		Component red,
		Component green,
		Component blue
	)
	:
		_alpha(alpha),
		_red(red),
		_green(green),
		_blue(blue)
{
}

/**
* Конструктор
* @param components - массив компонент цвета
*/
Color::Color
	(
		const Components& components
	)
{
	std::copy_n(components, ComponentsCount, _components);
}

/**
* Конструктор
* @param componentsVTK - массив вещественных значений в диапазоне [0, 1],
* используемых в VTK в качестве компонент цвета
*/
Color::Color
	(
		const VTKComponents& componentsVTK
	)
{
	FromVTKToColorComponents(componentsVTK, _components);
}

/**
* Получить прозрачность
* @return прозрачность
*/
Color::Component Color::GetAlpha() const
{
	return _alpha;
}

/**
* Получить непрозрачность цвета
* @return непрозрачность цвета
*/
Color::VTKComponent Color::GetOpacity() const
{
	return FromColorComponentToVTK(_alpha);
}

/**
* Получить массив компонент цвета
* @return массив компонент цвета
*/
const Color::Components& Color::GetComponents() const
{
	return _components;
}

/**
* Получить массив значений в диапазоне [0, 1], используемых в VTK в качестве компонент цвета
* @param componentsVTK - массив значений в диапазоне [0, 1], используемых в VTK в качестве компонент цвета
*/
void Color::GetVTKComponents
	(
		VTKComponents& componentsVTK
	)	const
{
	FromColorComponentsToVTK(_components, componentsVTK);
}

/**
* Чтение структуры с компонентами цвета из xml-документа
* @param colorElement - элемент xml-документа, содержащий структуру с компонентами цвета
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool Color::FromXmlElement
	(
		const TiXmlElement* colorElement
	)
{
	if (string(colorElement->Value()) != "Color")
	{
		return false;
	}

	const char* colorText;

	colorText = colorElement->GetText();
	if (colorText == nullptr)
	{
		return false;
	}
	if (!FromHtmlColor(colorText))
	{
		return false;
	}

	return true;
}

/**
* Сформировать цвет из строки, хранящей представление цвета в HTML-формате
* @param htmlColor - строка, хранящая представление цвета в HTML-формате
* @return признак успешного (true) или неуспешного (false) формирования
*/
bool Color::FromHtmlColor
	(
		const std::string& htmlColor
	)
{
	static const size_t digitsCount = 2;
	static const size_t htmlColorLength = 1 + digitsCount * ComponentsCount;

	if (htmlColor.size() != htmlColorLength)
	{
		return false;
	}
	if (htmlColor[0] != '#')
	{
		return false;
	}

	const char* htmlComponent = htmlColor.c_str() + 1;
	Component components[ComponentsCount];

	try
	{
		for
			(
				size_t componentIndex = 0;
				componentIndex < ComponentsCount;
				++componentIndex,
				htmlComponent += digitsCount
			)
		{
			string componentString(htmlComponent, digitsCount);
			size_t component = static_cast<size_t>(std::stoull(componentString, nullptr, 16));

			if (component > static_cast<size_t>(MaxComponentValue))
			{
				return false;
			}
			components[componentIndex] = static_cast<Component>(component);
		}
	}
	catch (const std::invalid_argument&)
	{
		return false;
	}
	catch (const std::out_of_range&)
	{
		return false;
	}
	std::copy(components, components + ComponentsCount, _components);

	return true;
}

/**
* Преобразовать значение компоненты цвета в вещественное число в диапазоне [0, 1] (используется в VTK)
* @param component - значение компоненты цвета
* @return соответствующее вещественное число в диапазоне [0, 1]
*/
// static
Color::VTKComponent Color::FromColorComponentToVTK
	(
		Component component
	)
{
	return component * 1.0 / Color::MaxComponentValue;
}

/**
* Преобразовать вещественное число в диапазоне [0, 1] (используется в VTK) в значение компоненты цвета
* @param componentVTK - вещественное число в диапазоне [0, 1]
* @return соответствующее значение компоненты цвета
*/
// static
Color::Component Color::FromVTKToColorComponent
	(
		VTKComponent componentVTK
	)
{
	return static_cast<Color::Component>(componentVTK * Color::MaxComponentValue);
}

/**
* Преобразовать массив компонент цвета в массив вещественных чисел в диапазоне [0, 1] (используются в VTK)
* @param components - массив компонент цвета
* @param componentsVTK - соответствующий массив вещественных чисел в диапазоне [0, 1]
*/
// static
void Color::FromColorComponentsToVTK
	(
		const Components& components,
		VTKComponents& componentsVTK
	)
{
	std::transform
		(
			components,
			components + Color::ComponentsCount,
			componentsVTK,
			FromColorComponentToVTK
		);
}

/**
* Преобразовать массив вещественных чисел в диапазоне [0, 1] (используются в VTK) в массив компонент цвета
* @param componentsVTK - массив вещественных чисел в диапазоне [0, 1]
* @param components - соответствующий массив компонент цвета
*/
// static
void Color::FromVTKToColorComponents
	(
		const VTKComponents& componentsVTK,
		Components& components
	)
{
	std::transform
		(
			componentsVTK,
			componentsVTK + Color::ComponentsCount,
			components,
			FromVTKToColorComponent
		);
}


// class BodyAppearance

/**
* Конструктор
* @param bodyNumber - номер тела
* @param color - цвет тела
*/
BodyAppearance::BodyAppearance
	(
		size_t bodyNumber,
		const Color& color
	)
	:
		_bodyNumber(bodyNumber),
		_color(color)
{
}

/**
* Получить номер тела
* @return номер тела
*/
size_t BodyAppearance::GetBodyNumber() const
{
	return _bodyNumber;
}

/**
* Получить цвет тела
* @return цвет тела
*/
const Color& BodyAppearance::GetColor() const
{
	return _color;
}

/**
* Чтение внешнего вида тела из xml-документа
* @param bodyAppearanceElement - элемент xml-документа, содержащий внешний вид тела
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool BodyAppearance::FromXmlElement
	(
		const TiXmlElement* bodyAppearanceElement
	)
{
	if (string(bodyAppearanceElement->Value()) != "BodyAppearance")
	{
		return false;
	}

	const char* bodyNumberAsString;

	bodyNumberAsString = bodyAppearanceElement->Attribute("BodyNumber");
	if (bodyNumberAsString == nullptr)
	{
		return false;
	}

	size_t bodyNumber;

	try
	{
		bodyNumber = static_cast<size_t>(std::stoull(bodyNumberAsString));
	}
	catch (const std::invalid_argument&)
	{
		return false;
	}
	catch (const std::out_of_range&)
	{
		return false;
	}

	const TiXmlElement* colorElement;

	colorElement = bodyAppearanceElement->FirstChildElement("Color");
	if (colorElement == nullptr)
	{
		return false;
	}

	Color color;

	if (!color.FromXmlElement(colorElement))
	{
		return false;
	}
	_bodyNumber = bodyNumber;
	_color = color;

	return true;
}

/**
* Чтение списка объектов, содержащих внешний вид тел, из xml-документа
* @param bodyAppearancesElement - элемент xml-документа, содержащий список параметров макроса
* @param bodyAppearances - считанный список объектов, содержащих внешний вид тел
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* bodyAppearancesElement,
		BodyAppearances& bodyAppearances
	)
{
	if (string(bodyAppearancesElement->Value()) != "BodyAppearances")
	{
		return false;
	}

	BodyAppearances bodyAppearancesTmp;
	const TiXmlElement* bodyAppearanceElement;

	bodyAppearanceElement = bodyAppearancesElement->FirstChildElement();
	while (bodyAppearanceElement != nullptr)
	{
		BodyAppearance bodyAppearance;

		if (!bodyAppearance.FromXmlElement(bodyAppearanceElement))
		{
			return false;
		}

		size_t bodyNumber = bodyAppearance.GetBodyNumber();

		if (bodyAppearancesTmp.find(bodyNumber) != bodyAppearancesTmp.end())
		{
			return false;
		}
		bodyAppearancesTmp.emplace(bodyNumber, bodyAppearance);
		bodyAppearanceElement = bodyAppearanceElement->NextSiblingElement();
	}
	bodyAppearances = bodyAppearancesTmp;

	return true;
}

/**
* Получить список объектов, содержащих внешний вид тел
* @return список объектов, содержащих внешний вид тел
*/
const BodyAppearances& ModelAppearance::GetBodyAppearances() const
{
	return _bodyAppearances;
}

/**
* Добавить объект, содержащий внешний вид тела
* @param bodyAppearance - добавляемый объект, содеражщий внешний вид тела
*/
void ModelAppearance::AddBodyAppearance
	(
		const BodyAppearance& bodyAppearance
	)
{
	_bodyAppearances.emplace(bodyAppearance.GetBodyNumber(), bodyAppearance);
}

/**
* Чтение внешнего вида модели из xml-документа
* @param xmlFileName - наименование xml-документа, содержащиго внешний вид модели
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool ModelAppearance::FromXmlDocument
	(
		const std::string& xmlFileName
	)
{
	TiXmlDocument xmlFile(xmlFileName.c_str());

	if (!xmlFile.LoadFile())
	{
		exceptions::ThrowFileInvalidFormat(xmlFileName);
	}

	const TiXmlElement* modelAppearanceElement;

	modelAppearanceElement = xmlFile.RootElement();
	if (modelAppearanceElement == nullptr)
	{
		return false;
	}
	if (string(modelAppearanceElement->Value()) != "ModelAppearance")
	{
		return false;
	}

	const TiXmlElement* bodyAppearancesElement;

	bodyAppearancesElement = modelAppearanceElement->FirstChildElement("BodyAppearances");
	if (bodyAppearancesElement == nullptr)
	{
		return false;
	}

	return ::FromXmlDocument(bodyAppearancesElement, _bodyAppearances);
}

/**
* Чтение внешнего вида модели из файла с внешним видом модели ФРУНД
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool ModelAppearance::FromFile()
{
	try
	{
		if (!fs::FileExist(_modelAppearanceFileName))
		{
			exceptions::ThrowFileNotFound(_modelAppearanceFileName);
		}
		if (!FromXmlDocument(_modelAppearanceFileName))
		{
			exceptions::ThrowFileNotOpened(_modelAppearanceFileName);
		}
	}
	catch (const exceptions::CoreException& exception)
	{
		exceptions::Output(exception);

		return false;
	}

	return true;
}

// наименование файла с внешним видом модели
// static
const std::string ModelAppearance::_modelAppearanceFileName = "ModelAppearance.xml";