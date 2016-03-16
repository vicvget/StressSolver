#ifndef FFORCE_CHAR_BUILDER_H

#define FFORCE_CHAR_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FForceChar.h"
#include "../../../AdditionalModules/fmath/Vector3.h"

#include <string>
#include "FModelBuilder.h"


using std::string;


class FForceCharBuilder
{
protected:
	FForceChar _forceChar;

public:
	virtual const FForceChar& Get() const;
};

class FForceGravityCharBuilder: public FForceCharBuilder
{
public:
	void Setup(int direction, const string& name, int id);
};

class FForceConstantCharBuilder: public FForceCharBuilder
{
public:
	void Setup(int direction, const string& name, int id, double param);
};

class FForceConstant3DCharBuilder: public FForceCharBuilder
{
public:
	void Setup(const string& name, int id, const MathHelpers::Vec3& force);
};

#endif // FFORCE_CHAR_BUILDER_H