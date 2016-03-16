#include "FForceCharBuilder.h"

#include "../FrundFacade/FElementType.h"
#include "FCharComponentBuilder.h"
#include "FElementBuilder.h"
#include "FModelBuilder.h"
#include "../../../AdditionalModules/RealsComparing/RealComparer.h"


const FForceChar& FForceCharBuilder::Get() const
{
	return _forceChar;
}

void FForceGravityCharBuilder::Setup(int direction, const string& name, int id)
{
	FElementBuilder::SetupCommonProperties(_forceChar, name, id, 
		FElementType::MakeFileId(FElementType::FETC_ForceChar, id+1));

	FCharComponentTypeGravityBuilder charComponentTypeGravityBuilder;
	charComponentTypeGravityBuilder.Setup(direction);
	_forceChar.AddCharComponent(charComponentTypeGravityBuilder.Get());
}


void FForceConstantCharBuilder::Setup(int direction, const string& name, int id, double param)
{
	FElementBuilder::SetupCommonProperties(_forceChar, name, id, 
		FElementType::MakeFileId(FElementType::FETC_ForceChar, id+1));

	FCharComponentTypeConstantBuilder charComponentTypeBuilder;
	charComponentTypeBuilder.Setup(direction, param);
	_forceChar.AddCharComponent(charComponentTypeBuilder.Get());

}

void FForceConstant3DCharBuilder::Setup(const string& name, int id, const MathHelpers::Vec3& force)
{
	FElementBuilder::SetupCommonProperties
		(
			_forceChar,
			name,
			id, 
			FElementType::MakeFileId(FElementType::FETC_ForceChar, id + 1)
		);

	// TODO: если будут проблемы с компаратором, задать нужное значение абсолютной погрешности вручную
	// с помощью метода RealComparer::SetAbsoluteError
	DoubleComparer comparer;

	for (int i = 0; i < 3; i++)
	{
		if (comparer.IsNotZero(force[i]))
		{
			FCharComponentTypeConstantBuilder charComponentTypeBuilder;

			charComponentTypeBuilder.Setup(i + 1, force[i]);
			_forceChar.AddCharComponent(charComponentTypeBuilder.Get());
		}
	}
}