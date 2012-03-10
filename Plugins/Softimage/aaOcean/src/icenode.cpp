#include "xsi_includes.h"
#include "aaOceanClass.cpp"
#include "icenode_portdefs.h"
#include "icenode_registration.h"
#include "functions.h"


SICALLBACK aaOcean_Init( CRef& in_ctxt )
{
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

	CDataArrayBool bRender( in_ctxt, ID_IN_RENDER_READY);
	CDataArrayLong pointID( in_ctxt, ID_IN_PointID);

	ICEocean->m_pointCount = pointID.GetCount();
	ICEocean->m_renderReady = bRender[0];
	
	int gridRes = int_sqrt(pointID.GetCount()) - 1;		
	ICEocean->reInit(gridRes);				
	fetchICEUserInput(in_ctxt, ICEocean);
	
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

	ULONG out_portID = in_ctxt.GetEvaluatedOutputPortID( );	
	switch( out_portID )
	{
		case ID_OUT_OCEAN:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(!ICEocean->m_renderReady)
				ICEocean->prepareOcean(1,1,0,0);			
				displayHeightField(outData,ICEocean);
		}
		break;

		case ID_OUT_FOAM:
		{
			if(ICEocean->m_isAllocated)
			{				
				ICEocean->prepareOcean(0,0,1,0);
				CDataArrayFloat outData( in_ctxt );
				displayFoam(outData,ICEocean);
			}
		}
		break;

		case ID_OUT_EIGEN_MINUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if( ICEocean->m_isAllocated && ICEocean->m_chopAmount > 0.0 && !ICEocean->m_renderReady)
			{
				ICEocean->prepareOcean(0,0,1,0);				
				displayEigenMinus(outData,ICEocean);
			}
			else
			{
				for(int i = 0; i<ICEocean->m_pointCount; i++){
					outData[i].PutX(0);outData[i].PutY(0);outData[i].PutZ(0);}
			}
		}
		break;

		case ID_OUT_EIGEN_PLUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(ICEocean->m_isAllocated && ICEocean->m_chopAmount > 0.0 && !ICEocean->m_renderReady)
			{
				ICEocean->prepareOcean(0,0,1,0);				
				displayEigenPlus(outData,ICEocean);
			}
			else
			{
				for(int i = 0; i<ICEocean->m_pointCount; i++)
				{
					outData[i].PutX(0);
					outData[i].PutY(0);
					outData[i].PutZ(0);
				}
			}
		}
		break;

		case ID_OUT_NORMALS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(ICEocean->m_isAllocated)
			{
				ICEocean->prepareOcean(1,1,0,1);
				displayNormals(outData,ICEocean);
			}
		}
		break;
	}

	return CStatus::OK;
}

