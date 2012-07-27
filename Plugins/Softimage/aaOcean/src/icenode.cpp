// aaOcean v2.5 Softimage ICE Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#include "xsi_includes.h"
#include "aaOceanClass.cpp"
#include "icenode_portdefs.h"
#include "icenode_registration.h"
#include "functions.h"


SICALLBACK aaOcean_Init( CRef& in_ctxt )
{
	fftwf_init_threads();
	aaOcean *ICEocean	= new aaOcean; 
	Context ctxt(in_ctxt);
	ctxt.PutUserData( (CValue::siPtrType)ICEocean);
	Application().LogMessage("[aaOcean ICE] : Successfully allocated memory");
	return CStatus::OK;
}

SICALLBACK aaOcean_Term( CRef& in_ctxt )
{
	Context ctxt( in_ctxt );
	CValue userData = ctxt.GetUserData( );
	aaOcean *ICEocean = (aaOcean *)(CValue::siPtrType)ctxt.GetUserData();
	if(ICEocean)
		delete ICEocean;	
	ICEocean = 0;
	ctxt.PutUserData( CValue() );
	Application().LogMessage(L"[aaOcean ICE] : Successfully cleaned up memory");
	return CStatus::OK;
}

SICALLBACK aaOcean_BeginEvaluate( ICENodeContext& in_ctxt )
{
	aaOcean *ICEocean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

	CDataArrayLong PointID( in_ctxt, ID_IN_PointID);
	int resolution = int_sqrt(PointID.GetCount()) - 1;		

	ICEocean->reInit(resolution);				
	fetchICEUserInput(in_ctxt, ICEocean);

	// not using aaOcean->Input. Using ICE dirty states, via fetchICEUserInput, since they are faster
	/*CDataArrayLong seed( in_ctxt, ID_IN_SEED);
	CDataArrayFloat	waveHeight( in_ctxt, ID_IN_WAVE_HEIGHT);
	CDataArrayFloat	waveSpeed( in_ctxt, ID_IN_WAVESPEED);
	CDataArrayFloat	chop( in_ctxt, ID_IN_CHOP);
	CDataArrayFloat	oceanScale( in_ctxt, ID_IN_OCEAN_SCALE );
	CDataArrayFloat windDir( in_ctxt, ID_IN_WINDDIR	);
	CDataArrayFloat cutoff( in_ctxt, ID_IN_CUTOFF);
	CDataArrayFloat velocity( in_ctxt, ID_IN_WINDVELOCITY);
	CDataArrayLong	windAlign( in_ctxt, ID_IN_WINDALIGN );
	CDataArrayFloat damp( in_ctxt, ID_IN_DAMP);*/
	//ICEocean->input(resolution, seed[0],oceanScale[0], velocity[0], cutoff[0], windDir[0], 
	//			windAlign[0], damp[0], waveSpeed[0], waveHeight[0],chop[0], (float)in_ctxt.GetTime().GetTime(CTime::Seconds));
	return CStatus::OK;
}

SICALLBACK aaOcean_Evaluate( ICENodeContext& in_ctxt )
{
	aaOcean *ICEocean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

	if(!ICEocean->m_isValid )
	{
		ICEocean = NULL;
		return CStatus::OK;
	}
	
	CIndexSet indexSet(in_ctxt, ID_IN_PointID );
	CDataArrayLong PointID(in_ctxt, ID_IN_PointID );
	CDataArrayBool bEnable( in_ctxt, ID_IN_ENABLE);
	CDataArrayFloat gridU( in_ctxt, ID_IN_GRID_LENGTH_U);
	CDataArrayFloat gridV( in_ctxt, ID_IN_GRID_LENGTH_V);
	int count = PointID.GetCount();
	ULONG out_portID = in_ctxt.GetEvaluatedOutputPortID( );	
	switch( out_portID )
	{
		case ID_OUT_OCEAN:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(bEnable[0])
				ICEocean->prepareOcean(1,1,0,0,0);			
			displayHeightField(outData,ICEocean, bEnable[0], gridU[0], gridV[0]);
		}
		break;

		case ID_OUT_FOAM:
		{
			if(ICEocean->m_isAllocated)
			{				
				ICEocean->prepareOcean(0,0,1,0,0);
				CDataArrayFloat outData( in_ctxt );
				displayFoam(outData,ICEocean);
			}
		}
		break;

		case ID_OUT_EIGEN_MINUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if( ICEocean->m_isAllocated && ICEocean->m_chopAmount > 0.0 && bEnable[0])
			{
				ICEocean->prepareOcean(0,0,1,0,0);				
				displayEigenMinus(outData,ICEocean);
			}
			else
			{
				for(int i = 0; i<count; ++i){
					outData[i].PutX(0);outData[i].PutY(0);outData[i].PutZ(0);}
			}
		}
		break;

		case ID_OUT_EIGEN_PLUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(ICEocean->m_isAllocated && ICEocean->m_chopAmount > 0.0 && bEnable[0])
			{
				ICEocean->prepareOcean(0,0,1,0,0);				
				displayEigenPlus(outData,ICEocean);
			}
			else
			{
				for(int i = 0; i<count; ++i)
				{
					outData[i].PutX(0);
					outData[i].PutY(0);
					outData[i].PutZ(0);
				}
			}
		}
		break;
	}

	return CStatus::OK;
}

