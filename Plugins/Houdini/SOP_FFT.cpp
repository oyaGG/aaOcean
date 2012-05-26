// Houdini FFT SOP
// FORWARD FFT NOT TESTED
// Author: Amaan Akram 
// www.amaanakram.com
// IFFT and FFT routines based on Paul Bourke
// http://paulbourke.net/miscellaneous/dft/


#define DLLEXPORT __declspec(dllexport)
#define SESI_LITTLE_ENDIAN
#define VERSION "10.0.295"
#define I386
#define SWAP_BITFIELDS
#define MAKING_DSO

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Math.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Guide.h>
#include <GEO/GEO_Detail.h>

#include "SOP_FFT.h"
#include "fft_funcs.h"
#include <cmath>


#define PRM_MENU_CHOICES	(PRM_ChoiceListType)(PRM_CHOICELIST_EXCLUSIVE |\
						     PRM_CHOICELIST_REPLACE)
static PRM_Name names[] = {
	PRM_Name("dimension",		"FFT Dimension"),
	PRM_Name("direction",		"FFT Direction"),
	PRM_Name("real_attrib",		"Real attribute"),
	PRM_Name("imag_attrib",		"Imaginary attribute"),
	PRM_Name("invert_alternate","Invert Alternate Points"),
	PRM_Name("tile_points",		"Tile Points"),
	PRM_Name("real_output",		"Real output to"),
	PRM_Name(0)
};

static PRM_Name	FFTDimMenuNames[] = {
	PRM_Name("1d",	"1D"),
	PRM_Name("2d",	"2D"),
	PRM_Name(0)
};

static PRM_Name	FFTDirMenuNames[] = {
	PRM_Name("forward",	"Forward"),
	PRM_Name("reverse",	"Reverse"),
	PRM_Name(0)
};

static PRM_Name	AxesMenuNames[] = {
	PRM_Name("none",	"none"),
	PRM_Name("x-axis",	"X"),
	PRM_Name("y-axis",	"Y"),
	PRM_Name("z-axis",	"Z"),
	PRM_Name(0)
};


PRM_ChoiceList	SOP_FFT::FFTDirMenu(PRM_MENU_CHOICES, FFTDirMenuNames); // direction menu
PRM_ChoiceList	SOP_FFT::FFTDimMenu(PRM_MENU_CHOICES, FFTDimMenuNames); // dimension menu
PRM_ChoiceList	SOP_FFT::AxisMenu(PRM_MENU_CHOICES, AxesMenuNames); // dimension menu

PRM_Template
SOP_FFT::myTemplateList[] = {
    PRM_Template(PRM_ORD, 1,	&names[0], PRMzeroDefaults, &FFTDimMenu, 0,0,0,1,"\nSelect whether you want a one or two dimensional FFT",0),
	PRM_Template(PRM_ORD, 1,	&names[1], PRMzeroDefaults, &FFTDirMenu, 0,0,0,1,"\nForward is Real-to-Complex\nReverse is Complex-to-Real",0),
	PRM_Template(PRM_STRING, 1, &names[2], PRMzeroDefaults,	0,0,0,0,1,"\nType the name of the point-class attribute that holds the REAL component data\nFFT Output is returned to this same attribute",0),
	PRM_Template(PRM_STRING, 1, &names[3], PRMzeroDefaults,	0,0,0,0,1,"\nType the name of the point-class attribute that holds the IMAGINARY component data\nFFT Output is returned to this same attribute",0),
	PRM_Template(PRM_TOGGLE, 1, &names[4], PRMzeroDefaults,	0,0,0,0,1,"\nMultiply alternate points with -1",0),
	PRM_Template(PRM_TOGGLE, 1, &names[5], PRMzeroDefaults,	0,0,0,0,1,"\nTile boundary points since our grid does not strictly contain 2^N points",0),
    PRM_Template(PRM_ORD, 1,	&names[6], PRMzeroDefaults, &AxisMenu, 0,0,0,1,"\nUse Real component of the output and apply to a selected direction for all points",0),
	PRM_Template(),
};

void newSopOperator(OP_OperatorTable *table)
{
     table->addOperator(new OP_Operator("fft",
					"FFT",
					 SOP_FFT::myConstructor,
					 SOP_FFT::myTemplateList,
					 1,
					 1,
					 0));
}	

OP_Node *SOP_FFT::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_FFT(net, name, op);
}

SOP_FFT::SOP_FFT(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op), myGroup(0)
{    
	myCurrPoint = -1;			// To prevent garbage values from being returned
    myTotalPoints = 0;			// Set the NPT local variable value
}

SOP_FFT::~SOP_FFT() {}


OP_ERROR SOP_FFT::cookMySop(OP_Context &context)
{
    GEO_Point	 *ppt;
    float		 now;
    UT_Vector3   normal, p;

	//my vars
	UT_String	 fft_dir;
	UT_String	 fft_dim;
	UT_String	 fft_real_name;
	UT_String	 fft_imag_name;
	UT_String	 axis_name;
	bool _2D_flag = false;
	COMPLEX **comp2D=NULL;
	double *real_1D,*imag_1D;	
	float *real_val=NULL;
	float *imag_val=NULL;
	int dir = 1; //default to forward FFT dir
	int axis, N;
	float sqroot;

    if (lockInputs(context) >= UT_ERROR_ABORT)
		return error();
    now = context.myTime;
    duplicateSource(0, context);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();
	
		
	//start set up vars
	myTotalPoints = gdp->points().entries();
	getFFTDir(fft_dir);
	getFFTDim(fft_dim);
	getRealAttribute(fft_real_name);
	getImagAttribute(fft_imag_name);
	getAxisName(axis_name);
	if(fft_dir=="reverse")
		dir = -1;
	if(axis_name == "x-axis")
		axis = 0;
	else if (axis_name == "y-axis")
		axis = 1;
	else if (axis_name == "z-axis")
		axis = 2;
	else
		axis = -1;
	
	if(fft_dim == "2d")
		_2D_flag = true;
	if(_2D_flag) 
		sqroot = sqrt((float)myTotalPoints) -1;
	else
		sqroot = myTotalPoints;
	//end setting up of vars
	
	//do some error checking to see if we can run an FFT -- need 2^n points	
	if (std::floor(sqroot) == sqroot)
	{
		if(PowerOfTwo(sqroot))
		{
			//printf("2D FFT ready with pt count of: %d", myTotalPoints);
			N = (int)sqroot + 1;
			
			if(_2D_flag)
				create_2Dcomplex_array(comp2D,N);
			else
			{
				real_1D = new double[myTotalPoints];
				imag_1D = new double[myTotalPoints];
			}
		}
		else
			return error();
	}
	else
		return error();


	//find the attributes that the user has typed in the parameter interface for the FFT SOP node
	int real_index = gdp->findPointAttrib((const char *)fft_real_name, sizeof(float), GB_ATTRIB_FLOAT);
	int imag_index = gdp->findPointAttrib((const char *)fft_imag_name, sizeof(float), GB_ATTRIB_FLOAT);
	
    
	//collect data for FFT if the attributes are found
	if (error() < UT_ERROR_ABORT && cookInputGroups(context) < UT_ERROR_ABORT && real_index!=-1 && imag_index!=-1)
    {
		int i=0; int j=0; 
		for (int k = 0; k < myTotalPoints; k++)
		{
			ppt = gdp->points()(k);
			real_val = (float *)ppt->getAttribData(real_index);
			imag_val = (float *)ppt->getAttribData(imag_index);
			
			if (i==N)
			{
				i=0;
				j++;
			}

			if(_2D_flag){
				comp2D[i][j].real = *real_val;
				comp2D[i][j].imag = *imag_val;
			}
			else
			{
				real_1D[k] = *real_val;
				imag_1D[k] = *imag_val;
			}
			i++;
		}
		
		//do the FFT
		if(_2D_flag)
			fft_2D(comp2D, N-1, N-1, dir);
		else
			fft_1D(myTotalPoints,real_1D,imag_1D);

		
		//send data back to Houdini
		i=0; j=0;
		for (int k = 0; k < myTotalPoints; k++)
		{
			if (i==N)
			{
				i=0;
				j++;
			}
			ppt = gdp->points()(k);
			real_val = (float *)ppt->getAttribData(real_index);
			imag_val = (float *)ppt->getAttribData(imag_index);

			if (alternate())
			{
				if(_2D_flag)
				{
					comp2D[i][j].real *= pow(-1.0f,k);
					comp2D[i][j].imag *= pow(-1.0f,k);
				}
				else
				{
					real_1D[k]*= pow(-1.0f,k);
					imag_1D[k]*= pow(-1.0f,k);
				}
			}

			if (tilepoints() && _2D_flag)
			{
				if( i<N-1 && j<N-1)	{	    
					*real_val = comp2D[i][j].real;
					*imag_val = comp2D[i][j].imag;
					if(axis!=-1) ppt->getPos()(axis) += comp2D[i][j].real;
				}
				if(i==N-1){			
					*real_val = comp2D[0][j].real;
					*imag_val = comp2D[0][j].imag;
					if(axis!=-1) ppt->getPos()(axis) += comp2D[0][j].real;
				}
				if(j==N-1){			
					*real_val = comp2D[i][0].real;
					*imag_val = comp2D[i][0].imag;
					if(axis!=-1) ppt->getPos()(axis) += comp2D[i][0].real;
				}
				if(j==N-1 && i==N-1){	
					*real_val = comp2D[0][0].real;
					*imag_val = comp2D[0][0].imag; 
					if(axis!=-1) ppt->getPos()(axis) += comp2D[0][0].real;
				}
			}
			else 
			{
				if(_2D_flag)
				{
					*real_val = comp2D[i][j].real;
					*imag_val = comp2D[i][j].imag;
					if(axis!=-1) ppt->getPos()(axis) += comp2D[i][j].real;
				}
				else
				{
					*real_val = real_1D[k];
					*imag_val = imag_1D[k];
					if(axis!=-1) ppt->getPos()(axis) += real_1D[k];
				}
			}
			
			i++;
		}
    }

	//clean up
    unlockInputs();
    resetLocalVarRefs();
	if(_2D_flag)
		delete_2Dcomplex_array(comp2D,N);
	else{
		delete [] real_1D;
		delete [] imag_1D;
	}
   
    return error();
}

