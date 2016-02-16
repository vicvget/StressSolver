#ifndef COMPARISON_TAG_H

#define COMPARISON_TAG_H


#include <cstddef>


// тег сравнения в качестве элемента перечисления тегов сравнения
#define COMPARISON_TAG_ENUM_ELEMENT(tagName) \
	tagName,

// определелить макрос, отвечающий за развертывания тега сравнения
#define COMPARISON_TAG COMPARISON_TAG_ENUM_ELEMENT

// перечисление тегов сравнения,
// используемых в обертке над функциями сравнения для определения типа сравнения
enum class ComparisonTag : size_t
{
#	include "ComparisonTagList.h"
};

#undef COMPARISON_TAG

#undef COMPARISON_TAG_ENUM_ELEMENT


// макрос для определения константы тега сравнения
#define DEFINE_COMPARISON_TAG(tagName) \
	const ComparisonTag tagName##Tag = ComparisonTag::tagName;

// определелить макрос, отвечающий за развертывания тега сравнения
#define COMPARISON_TAG DEFINE_COMPARISON_TAG

#include "ComparisonTagList.h"

#undef COMPARISON_TAG

#undef DEFINE_COMPARISON_TAG


#endif // COMPARISON_TAG_H