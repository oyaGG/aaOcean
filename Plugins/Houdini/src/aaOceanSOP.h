
#ifndef __aaOceanSOP_h__
#define __aaOceanSOP_h__

#include <SOP/SOP_Node.h>
#include "aaOceanClass.cpp"

class aaOceanSOP : public SOP_Node
{
public:
	aaOceanSOP(OP_Network *net, const char *name, OP_Operator *op);
	virtual ~aaOceanSOP();

	static PRM_Template myTemplateList[];
	static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

protected:
	virtual const char	*inputLabel(unsigned idx) const;

	/// Method to cook geometry for the SOP
	virtual OP_ERROR cookMySop(OP_Context &context);

private:
	//void	getGroups(UT_String &str)	{ evalString(str, "group", 0, 0); }
	int		RESOLUTION()			{ return evalInt("resolution", 0, 0); }
	int		SEED()					{ return evalInt("seed", 0, 0); }
	fpreal	OCEANSCALE(fpreal t)	{ return evalFloat("oceanScale", 0, t); }

	fpreal	VELOCITY(fpreal t)		{ return evalFloat("velocity", 0, t); }
	fpreal	CUTOFF(fpreal t)		{ return evalFloat("cutoff", 0, t); }
	fpreal	WINDDIR(fpreal t)		{ return evalFloat("windDir", 0, t); }
	int		WINDALIGN()				{ return evalInt("windAlign", 0, 0); }

	fpreal	DAMP(fpreal t)			{ return evalFloat("damp", 0, t); }
	fpreal	WAVESPEED(fpreal t)		{ return evalFloat("waveSpeed", 0, t); }
	fpreal	WAVEHEIGHT(fpreal t)	{ return evalFloat("waveHeight", 0, t); }
	fpreal	CHOP(fpreal t)			{ return evalFloat("chop", 0, t); }
	int		ENABLEFOAM()			{ return evalInt("enableFoam", 0, 0); }
	fpreal	TIMEOFFSET(fpreal t)	{ return evalFloat("timeOffset", 0, t); }
	
	/// This variable is used together with the call to the "checkInputChanged"
	/// routine to notify the handles (if any) if the input has changed.
	GU_DetailGroupPair	 myDetailGroupPair;

	/// This is the group of geometry to be manipulated by this SOP and cooked
	/// by the method "cookInputGroups".
	const GA_PointGroup	*myGroup;

	aaOcean *pOcean;
};

#endif