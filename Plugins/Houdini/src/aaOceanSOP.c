// aaOcean v2.6 Houdini SOP Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef __GNUC__ // defined in MSVC project settings, preprocessor section
#define VERSION "12.5.376"
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
#endif

#include "aaOceanSOP.h"
#include "timer/Timer.cpp"

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Math.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include <UT/UT_Options.h>
#include <GEO/GEO_AttributeHandle.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>

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
	PRM_Name("oceanDepth",		"Ocean Depth"),
	PRM_Name("surfaceTension",	"Surface Tension"),

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
	PRM_Name("loopTime",		"Loop Time"),
    PRM_Name("uvAttribute",		"UV Attribute"),

};

// defining some custom ranges and defaults

static PRM_Range		resolutionRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_RESTRICTED, 6);
static PRM_Default		resolutionDefault(4);

static PRM_Range		oceanScaleRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_UI, 200.0);
static PRM_Default		oceanScaleDefault(100.0);

static PRM_Range		oceanDepthRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_RESTRICTED, 10000.0);
static PRM_Default		oceanDepthDefault(10000.0);

static PRM_Range		seedRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_RESTRICTED, 15);
static PRM_Default		seedDefault(1);

static PRM_Range		velocityRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_RESTRICTED, 30.0);
static PRM_Default		velocityDefault(4.0);

static PRM_Range		loopTimeRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_UI, 1000.0);
static PRM_Default		loopTimeDefault(1000.0);

PRM_Template aaOceanSOP::myTemplateList[] = 
{	
    PRM_Template(PRM_INT_E,	1, &names[0],  &resolutionDefault,	0, &resolutionRange),		// resolution	// 0
    PRM_Template(PRM_FLT_J,	1, &names[2],  &oceanScaleDefault,  0, &oceanScaleRange),		// oceanScale	// 2
	PRM_Template(PRM_FLT_J,	1, &names[3],  &oceanDepthDefault,	0, &oceanDepthRange),		// oceanDepth	// 3
    PRM_Template(PRM_INT_E,	1, &names[1],  &seedDefault,		0, &seedRange),				// seed			// 1
    PRM_Template(PRM_FLT_J,	1, &names[14], PRMzeroDefaults,		0, &PRMscaleRange),			// timeOffset	// 14
	PRM_Template(PRM_FLT_J,	1, &names[15], &loopTimeDefault,	0, &loopTimeRange),			// loop time	// 15

    PRM_Template(PRM_FLT_J,	1, &names[11], PRMoneDefaults,		0, &PRMdivision0Range),		// waveHeight	// 11
    PRM_Template(PRM_FLT_J,	1, &names[5],  &velocityDefault,	0, &velocityRange),			// velocity (Wave Size) //5
    PRM_Template(PRM_FLT_J,	1, &names[10], PRMoneDefaults,		0, &PRMdivision0Range),		// waveSpeed	// 10
    PRM_Template(PRM_FLT_J,	1, &names[12], PRMzeroDefaults,		0, &PRMrolloffRange),		// chop			// 12
    PRM_Template(PRM_FLT_J,	1, &names[6],  PRMzeroDefaults,		0, &PRMdivision0Range),		// cutoff (Wave Smooth) // 6

    PRM_Template(PRM_FLT_J,	1, &names[7],  PRMzeroDefaults,		0, &PRMangleRange),			// windDir		// 7
    PRM_Template(PRM_FLT_J,	1, &names[9],  PRMzeroDefaults,		0, &PRMunitRange),			// damp			// 9
    PRM_Template(PRM_INT_E,	1, &names[8],  PRMzeroDefaults,		0, &PRMdivision0Range),		// windAlign	// 8

    PRM_Template(PRM_TOGGLE,1, &names[13]),													// enable Foam  // 13
    PRM_Template(PRM_STRING,1, &names[16], 0),												// UV Attribute	// 16

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
	if (lockInputs(context) >= UT_ERROR_ABORT)
		return error();

    duplicateSource(0, context);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();

	// variable declarations
	float now  = context.getTime();

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
					OCEANDEPTH(now),
					SURFACETENSION(now),
					VELOCITY(now), 
					CUTOFF(now), 
					WINDDIR(now), 
					WINDALIGN(), 
					DAMP(now), 
					WAVESPEED(now), 
					WAVEHEIGHT(now),
					CHOP(now), 
					now,
					LOOPTIME(now),
					enableEigens,
					FALSE);

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
    const char* UVAttribName = (const char *)UvAttribute;
	uvRef = gdp->findFloatTuple(GA_ATTRIB_POINT, UVAttribName, 3);

    if(uvRef.isValid() == TRUE)
    {
		uvAttribute = uvRef.getAttribute();
		uvTuple = uvRef.getAIFTuple(); 
    }
	else
	{
		// uv attribute not found
        char msg[256];
        sprintf(msg, "[aaOcean] Specified UV attribute \'%s\' not found on geometry.\
                     \nUV's are required for aaOcean to cook", UVAttribName);
        cout<<msg;
        cout.flush();
        addError(SOP_MESSAGE, msg); 
        unlockInputs();
        return error();
	}

    // setup local variables to output Eigens
    if(enableEigens)
    {
        eVecPlusRef  = gdp->addFloatTuple(GA_ATTRIB_POINT, eVecPlusName,	3);
        eVecMinusRef = gdp->addFloatTuple(GA_ATTRIB_POINT, eVecMinusName,	3);
        eValuesRef   = gdp->addFloatTuple(GA_ATTRIB_POINT, eValuesName,		1);

		eVecPlusHandle	= GA_RWHandleV3::GA_RWHandleT(eVecPlusRef.getAttribute());
		eVecMinusHandle = GA_RWHandleV3::GA_RWHandleT(eVecMinusRef.getAttribute());
		eValuesHandle	= GA_RWHandleF::GA_RWHandleT(eValuesRef.getAttribute());
    }
	
    // inputs validated. Begin writing ocean data to output handles
	int npts = gdp->getNumPoints();
	#pragma omp parallel for 
	for (int pt_offset = 0; pt_offset < npts; ++pt_offset)
    {
		UT_Vector3F pos = gdp->getPos3(pt_offset);
		UT_Vector3F UV;
		
        uvTuple->get(uvAttribute, pt_offset, UV.data(), 3);
        // Houdini V coord runs in opposite direction compared to Softimage/Maya
        // Conforming with other apps to make ocean shape consistent across apps
		float u = UV.x();
		float v = 1.0f - (fmod(UV.y(), 1.0f));

		pos.y() += pOcean->getOceanData(u, v, aaOcean::eHEIGHTFIELD);
        if(pOcean->isChoppy())
        {
            pos.x() += pOcean->getOceanData(u, v, aaOcean::eCHOPX);
            pos.z() += pOcean->getOceanData(u, v, aaOcean::eCHOPZ);
        }
		gdp->setPos3(pt_offset, pos);

       if(enableEigens)
        {
			UT_Vector3F eigenVectorPlusValue;
			UT_Vector3F eigenVectorMinusValue;
			float eigenValue;

            eigenVectorPlusValue.x() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSX);
            eigenVectorPlusValue.y() =  0.0f;
            eigenVectorPlusValue.z() =  pOcean->getOceanData(u, v, aaOcean::eEIGENPLUSZ);

            eigenVectorMinusValue.x() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSX);
            eigenVectorMinusValue.y() = 0.0f;
            eigenVectorMinusValue.z() = pOcean->getOceanData(u, v, aaOcean::eEIGENMINUSZ);

            eigenValue = pOcean->getOceanData(u, v, aaOcean::eFOAM);

			eVecPlusHandle.set(pt_offset,eigenVectorPlusValue);
			eVecMinusHandle.set(pt_offset,eigenVectorMinusValue);
			eValuesHandle.set(pt_offset,eigenValue);
        }
    }
    unlockInputs();

    return error();
}


const char *aaOceanSOP::inputLabel(unsigned) const
{
    return "UV'ed Geometry to simulate ocean on";
}
