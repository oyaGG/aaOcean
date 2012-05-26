// Houdini FFT SOP
// FORWARD FFT NOT TESTED
// Author: Amaan Akram 
// www.amaanakram.com
// IFFT and FFT routines based on Paul Bourke
// http://paulbourke.net/miscellaneous/dft/

#ifndef __SOP_FFT_h__
#define __SOP_FFT_h__

#include <SOP/SOP_Node.h>

class SOP_FFT : public SOP_Node
{
public:
	     
	SOP_FFT(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_FFT();

    static PRM_Template		 myTemplateList[];
    static OP_Node		*myConstructor(OP_Network*, const char *,  OP_Operator *);
	

protected:
    // Method to cook geometry for the SOP
    virtual OP_ERROR		 cookMySop(OP_Context &context);

	static PRM_ChoiceList	 FFTDirMenu;
    static PRM_ChoiceList	 FFTDimMenu;
	static PRM_ChoiceList	 AxisMenu;

private:
    void	getFFTDim(UT_String &str)			{ evalString(str, 0, 0, 0); }
	void	getFFTDir(UT_String &str)			{ evalString(str, 1, 0, 0); }
	void	getRealAttribute(UT_String &str)	{ evalString(str, 2, 0, 0); }
	void	getImagAttribute(UT_String &str)	{ evalString(str, 3, 0, 0); }
	int		alternate(void)						{ return  evalInt(4, 0, 0); }
	int		tilepoints(void)					{ return  evalInt(5, 0, 0); }
	void	getAxisName(UT_String &str)			{ evalString(str, 6, 0, 0); }

    // This variable is used together with the call to the "checkInputChanged"
    // routine to notify the handles (if any) if the input has changed.
    GU_DetailGroupPair	 myDetailGroupPair;

    // This is the group of geometry to be manipulated by this SOP and cooked
    // by the method "cookInputGroups".
    const GB_PointGroup	*myGroup;

	// "Local" vars
    int		myCurrPoint;
    int		myTotalPoints;
};

#endif
