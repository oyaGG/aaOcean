// aaOcean v2.6 Houdini SOP Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifdef MSVC_DEFINES // defined in MSVC project settings, preprocessor section
#define VERSION "12.1.125"
#define I386 
#define WIN32 
#define SWAP_BITFIELDS 
#define _WIN32_WINNT 0x0501
#define WINVER 0x0501
#define NOMINMAX 
#define SESI_LITTLE_ENDIAN 
#define NEED_SPECIALIZATION_STORAGE 
#define AMD64 
#define SIZEOF_VOID_P 8
#define MAKING_DSO
#define SESI_TAGINFO "Produced by: Amaan Akram"
#endif

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Math.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include <UT/UT_Options.h>
#include "aaOceanSOP.h"


void newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator("aaOceanSOP",
        "aaOceanSOP",
        aaOceanSOP::myConstructor,
        aaOceanSOP::myTemplateList,
        1,
        1,
        0));
}

static PRM_Name names[] = 
{
    PRM_Name("resolution",		"Resolution"),
    PRM_Name("seed",			"Seed"),
    PRM_Name("oceanScale",		"Ocean Scale"),

    PRM_Name("velocity",		"Wave Size"),
    PRM_Name("cutoff",			"Wave Smooth"),
    PRM_Name("windDir",			"Wind Dir"),
    PRM_Name("windAlign",		"Wind Align"),

    PRM_Name("damp",			"Reflected Waves"),
    PRM_Name("waveSpeed",		"Wave Speed"),
    PRM_Name("waveHeight",		"Wave Height"),
    PRM_Name("chop",			"Chop Amount"),
    PRM_Name("enableEigens",	"Output Eigens Attributes"),
    PRM_Name("timeOffset",		"Time Offset"),
    PRM_Name("uvAttribute",		"UV Attribute"),

};

PRM_Template aaOceanSOP::myTemplateList[] = 
{	
    PRM_Template(PRM_INT_E,	1, &names[0],  PRMfourDefaults, 0, &PRMfrequencyRange),		// resolution
    PRM_Template(PRM_FLT_J,	1, &names[2],  PRM100Defaults,  0, &PRMpositiveRange),		// oceanScale
    PRM_Template(PRM_INT_E,	1, &names[1],  PRMoneDefaults,  0, &PRMfrequency10Range),	// seed
    PRM_Template(PRM_FLT_J,	1, &names[12], PRMzeroDefaults, 0, &PRMscaleRange),			// timeOffset

    PRM_Template(PRM_FLT_J,	1, &names[9],  PRMoneDefaults,  0, &PRMdivision0Range),		// waveHeight
    PRM_Template(PRM_FLT_J,	1, &names[3],  PRMfourDefaults, 0, &PRMdivision0Range),		// velocity (Wave Size)
    PRM_Template(PRM_FLT_J,	1, &names[8],  PRMoneDefaults,  0, &PRMdivision0Range),		// waveSpeed
    PRM_Template(PRM_FLT_J,	1, &names[10], PRMzeroDefaults, 0, &PRMrolloffRange),		// chop
    PRM_Template(PRM_FLT_J,	1, &names[4],  PRMzeroDefaults, 0, &PRMdivision0Range),		// cutoff (Wave Smooth)

    PRM_Template(PRM_FLT_J,	1, &names[5],  PRMzeroDefaults, 0, &PRMangleRange),			// windDir
    PRM_Template(PRM_FLT_J,	1, &names[7],  PRMzeroDefaults, 0, &PRMunitRange),			// damp
    PRM_Template(PRM_INT_E,	1, &names[6],  PRMzeroDefaults, 0, &PRMdivision0Range),		// windAlign

    PRM_Template(PRM_TOGGLE,1, &names[11]),												// enable Foam
    PRM_Template(PRM_STRING,1, &names[13], 0),											// UV Attribute

    PRM_Template(),
};


OP_Node *aaOceanSOP::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new aaOceanSOP(net, name, op);
}

aaOceanSOP::aaOceanSOP(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
{
	sprintf(eVecPlusName,	"eVecPlus");
	sprintf(eVecMinusName,	"eVecMinus");
	sprintf(eValuesName,	"eValues");
	enableEigens = FALSE;

    pOcean = new aaOcean;
}

aaOceanSOP::~aaOceanSOP() 
{
    if(pOcean)
        delete pOcean;
}

OP_ERROR aaOceanSOP::cookMySop(OP_Context &context)
{
    fpreal		now		= context.getTime();
    int			npts	= gdp->points().entries();
    float u, v;
	UT_Vector3 eigenVectorPlusValue, eigenVectorMinusValue, eigenValuesValue;
	UT_Vector4	PtValue;
	GEO_AttributeHandle	UvHandle, PtHandle, eVecPlusHandle, eVecMinusHandle, eValuesHandle;
	
    if (lockInputs(context) >= UT_ERROR_ABORT)
        return error();

    duplicateSource(0, context);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();

    // Flag the SOP as being time dependent (i.e. cook on time changes)
    flags().timeDep = 1;

    // start pulling in SOP inputs and send to aaOcean 
	enableEigens = (ENABLEEIGENS() != 0);
	if(pOcean->isChoppy() && enableEigens)
		enableEigens = TRUE;
    now = now + TIMEOFFSET(now);
    pOcean->input(	RESOLUTION(), 
                    SEED(),
                    OCEANSCALE(now), 
                    VELOCITY(now), 
                    CUTOFF(now), 
                    WINDDIR(now), 
                    WINDALIGN(), 
                    DAMP(now), 
                    WAVESPEED(now), 
                    WAVEHEIGHT(now),
                    CHOP(now), 
                    now,
                    enableEigens);

    if(pOcean->isValid() == FALSE)
    {
		char msg[256];
		sprintf(msg, "[aaOcean] Failed to allocate Ocean. Bad input", msg);
		addError(SOP_MESSAGE, ""); 
		cout<<msg;
		cout.flush();
        unlockInputs();
        return error();
    }

    // get the user-specified attribute that holds uv-data
    getUVAttributeName(UvAttribute);
    if(UvAttribute.length() == 0)
        UvAttribute = "uv";
    const char* attribName = (const char *)UvAttribute;
    UvHandle = gdp->getAttribute(GEO_POINT_DICT, attribName);

    if(UvHandle.isAttributeValid() == FALSE)
    {
        // uv attribute not found
		char msg[256];
		sprintf(msg, "[aaOcean] Specified UV attribute \'%s\' not found on geometry.\
					 \nUV's are required for aaOcean to cook", attribName);
		cout<<msg;
		cout.flush();
		addError(SOP_MESSAGE, msg); 
        unlockInputs();
        return error();
    }
	
	// setup local variables to output Eigens
	if(enableEigens)
	{
		gdp->addFloatTuple(GA_ATTRIB_POINT, eVecPlusName,	3);
		gdp->addFloatTuple(GA_ATTRIB_POINT, eVecMinusName,	3);
		gdp->addFloatTuple(GA_ATTRIB_POINT, eValuesName,	3);
	}
    
    // inputs validated. Begin writing ocean data to output handles
    #pragma omp parallel for private(PtHandle, PtValue, u, v, UvHandle, eVecPlusHandle, eVecMinusHandle, eValuesHandle, eigenVectorPlusValue, eigenVectorMinusValue, eigenValuesValue)
    for (int i = 0; i < npts; ++i)
    {
        UvHandle = gdp->getAttribute(GEO_POINT_DICT, attribName);
        PtHandle  = gdp->getAttribute(GEO_POINT_DICT, "P");

        UvHandle.setElement(gdp->points()(i));
        PtHandle.setElement(gdp->points()(i));

        u = UvHandle.getV3().x();
        // Houdini V coord runs in opposite direction compared to Softimage/Maya
        // Conforming with other apps to make ocean shape consistent across apps
        v = 1.0 - (fmod(UvHandle.getV3().y(), 1.0f));

        PtValue = PtHandle.getV3();

        PtValue.y() += pOcean->getOceanData(u, v, aaOcean::eHEIGHTFIELD);
        if(pOcean->isChoppy())
        {
            PtValue.x() += pOcean->getOceanData(u, v, aaOcean::eCHOPX);
            PtValue.z() += pOcean->getOceanData(u, v, aaOcean::eCHOPZ);
        }
        PtHandle.setV3(PtValue);

		if(enableEigens)
		{
			eVecPlusHandle.setElement(gdp->points()(i));
			eVecMinusHandle.setElement(gdp->points()(i));
			eValuesHandle.setElement(gdp->points()(i));

			eVecPlusHandle	= gdp->getAttribute(GEO_POINT_DICT, eVecPlusName);
			eVecMinusHandle = gdp->getAttribute(GEO_POINT_DICT, eVecMinusName);
			eValuesHandle	= gdp->getAttribute(GEO_POINT_DICT, eValuesName);

			eigenVectorPlusValue.x() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSX);
			eigenVectorPlusValue.y() =  0.0f;
			eigenVectorPlusValue.z() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSZ);

			eigenVectorMinusValue.x() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSX);
			eigenVectorMinusValue.y() = 0.0f;
			eigenVectorMinusValue.z() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSZ);
 
			eigenValuesValue.x() = pOcean->getOceanData(u, v, aaOcean::eFOAM);
			eigenValuesValue.z() = eigenValuesValue.y() = eigenValuesValue.x();

			eVecPlusHandle.setV3(eigenVectorPlusValue);
			eVecMinusHandle.setV3(eigenVectorMinusValue);
			eValuesHandle.setV3(eigenValuesValue);
		}
    }
	
    // Notify the display cache
    gdp->notifyCache(GU_CACHE_ALL);
    unlockInputs();
    return error();
}


const char *aaOceanSOP::inputLabel(unsigned) const
{
    return "UV'ed Geometry to simulate ocean on";
}
