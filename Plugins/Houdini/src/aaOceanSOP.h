// aaOcean v2.6 Houdini SOP Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

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
	// working variables
	bool enableEigens;
	char eVecPlusName[10];
	char eVecMinusName[10];
	char eValuesName[10];
	UT_String	UvAttribute;
	

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
    int		ENABLEEIGENS()			{ return evalInt("enableEigens", 0, 0); }
    fpreal	TIMEOFFSET(fpreal t)	{ return evalFloat("timeOffset", 0, t); }

    void	getUVAttributeName(UT_String &str){ evalString(str, "uvAttribute", 0, 0); }
    void	getEigenPlusAttribute(UT_String &str){ evalString(str, "eigenPlusAttribute", 0, 0); }
    void	getEigenMinusAttribute(UT_String &str){ evalString(str, "eigenMinusAttribute", 0, 0); }
    void	getEigenValueAttribute(UT_String &str){ evalString(str, "eigenValueAttribute", 0, 0); }
    
    /// This variable is used together with the call to the "checkInputChanged"
    /// routine to notify the handles (if any) if the input has changed.
    GU_DetailGroupPair	 myDetailGroupPair;

    // pointer to our aaOcean class
    aaOcean *pOcean;
};

#endif