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

SICALLBACK aaOcean_Init( CRef& in_ctxt )
{
	// initialize fftw threads routines
	fftwf_init_threads();

	// create new aaOceanClass instance
	aaOcean *pOcean	= new aaOcean;

	// store ocean pointer in user-data
	Context ctxt(in_ctxt);
	ctxt.PutUserData( (CValue::siPtrType)pOcean);

	Application().LogMessage("[aaOcean ICE] : Successfully allocated memory");
	return CStatus::OK;
}

SICALLBACK aaOcean_Term( CRef& in_ctxt )
{
	Context ctxt( in_ctxt );

	// get ocean pointer from user-data and do clean up
	CValue userData = ctxt.GetUserData( );
	aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)ctxt.GetUserData();
	if(pOcean)
		delete pOcean;	
	pOcean = 0;
	ctxt.PutUserData( CValue() );
	Application().LogMessage(L"[aaOcean ICE] : Successfully cleaned up memory");

	// call fftw cleanup routines
	fftwf_cleanup_threads();
	fftwf_cleanup();

	return CStatus::OK;
}

SICALLBACK aaOcean_BeginEvaluate( ICENodeContext& in_ctxt )
{
	// get ocean pointer from user-data
	aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

	// get ICE node input port arrays
	CDataArrayLong PointID( in_ctxt, ID_IN_PointID);
	CDataArrayLong resolution( in_ctxt, ID_IN_RESOLUTION);		
	CDataArrayLong seed( in_ctxt, ID_IN_SEED);
	CDataArrayFloat	waveHeight( in_ctxt, ID_IN_WAVE_HEIGHT);
	CDataArrayFloat	waveSpeed( in_ctxt, ID_IN_WAVESPEED);
	CDataArrayFloat	chop( in_ctxt, ID_IN_CHOP);
	CDataArrayFloat	oceanScale( in_ctxt, ID_IN_OCEAN_SCALE );
	CDataArrayFloat windDir( in_ctxt, ID_IN_WINDDIR	);
	CDataArrayFloat cutoff( in_ctxt, ID_IN_CUTOFF);
	CDataArrayFloat velocity( in_ctxt, ID_IN_WINDVELOCITY);
	CDataArrayLong	windAlign( in_ctxt, ID_IN_WINDALIGN );
	CDataArrayFloat damp( in_ctxt, ID_IN_DAMP);
	CDataArrayBool enableFoam( in_ctxt, ID_IN_ENABLEFOAM);
	CDataArrayFloat time( in_ctxt, ID_IN_TIME);

	pOcean->input(resolution[0], 
		seed[0],
		oceanScale[0], 
		velocity[0], 
		cutoff[0], 
		windDir[0], 
		windAlign[0], 
		damp[0], 
		waveSpeed[0], 
		waveHeight[0],
		chop[0], 
		time[0],
		enableFoam[0]);

	return CStatus::OK;
}

SICALLBACK aaOcean_Evaluate( ICENodeContext& in_ctxt )
{
	aaOcean *pOcean = (aaOcean *)(CValue::siPtrType)in_ctxt.GetUserData();

	if(!pOcean->isValid() )
	{
		pOcean = NULL;
		Application().LogMessage(CString("Invalid Ocean Input"), siErrorMsg);
		return CStatus::OK;
	}
	
	CIndexSet indexSet(in_ctxt, ID_IN_PointID );
	CDataArrayLong PointID(in_ctxt, ID_IN_PointID );
	CDataArrayBool bEnable(in_ctxt, ID_IN_ENABLE);
	CDataArrayFloat uCoord(in_ctxt, ID_IN_U);
	CDataArrayFloat vCoord(in_ctxt, ID_IN_V);
	CDataArrayMatrix4f transform(in_ctxt, ID_IN_TRANSFORM);

	float worldSpaceDisplacementVec[3] = {0.0f, 0.0f, 0.0f};
	float oceanLocalSpace[3];
	const int count = PointID.GetCount();
	int transformArraySize = 0;
	bool transformSingleton = TRUE;

	if(transform.GetCount() > 1)
		transformSingleton = FALSE;

	ULONG out_portID = in_ctxt.GetEvaluatedOutputPortID( );	
	switch( out_portID )
	{
		case ID_OUT_OCEAN:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(bEnable[0])
			{
				#pragma omp parallel for private(worldSpaceDisplacementVec, oceanLocalSpace)
				for(int i = 0; i<count; ++i)
				{
					// get ocean displacement vector
					worldSpaceDisplacementVec[1] = pOcean->getOceanData(uCoord[i], vCoord[i], HEIGHTFIELD);
					if(pOcean->isChoppy())
					{
						worldSpaceDisplacementVec[0] = pOcean->getOceanData(uCoord[i], vCoord[i], CHOPX);
						worldSpaceDisplacementVec[2] = pOcean->getOceanData(uCoord[i], vCoord[i], CHOPZ);
					}

					// multiply displacement vector by input transform matrix
					if(!transformSingleton)
						transformArraySize = i;

					oceanLocalSpace[0] =	worldSpaceDisplacementVec[0] * transform[transformArraySize].GetValue(0,0) + 
											worldSpaceDisplacementVec[1] * transform[transformArraySize].GetValue(1,0) + 
											worldSpaceDisplacementVec[2] * transform[transformArraySize].GetValue(2,0);

					oceanLocalSpace[1] =	worldSpaceDisplacementVec[0] * transform[transformArraySize].GetValue(0,1) + 
											worldSpaceDisplacementVec[1] * transform[transformArraySize].GetValue(1,1) + 
											worldSpaceDisplacementVec[2] * transform[transformArraySize].GetValue(2,1);

					oceanLocalSpace[2] =	worldSpaceDisplacementVec[0] * transform[transformArraySize].GetValue(0,2) + 
											worldSpaceDisplacementVec[1] * transform[transformArraySize].GetValue(1,2) + 
											worldSpaceDisplacementVec[2] * transform[transformArraySize].GetValue(2,2);

					outData[i].PutX(oceanLocalSpace[0]);
					outData[i].PutY(oceanLocalSpace[1]);
					outData[i].PutZ(oceanLocalSpace[2]);
				}
			}
		}
		break;

		case ID_OUT_FOAM:
		{
			if(pOcean->isChoppy() && bEnable[0])
			{				
				CDataArrayFloat outData( in_ctxt );

				// output raw (unscaled) foam in ICE deformer
				#pragma omp parallel for
				for(int i = 0; i<count; ++i)
					outData[i] = pOcean->getOceanData(uCoord[i], vCoord[i], FOAM);
			}
		}
		break;

		case ID_OUT_EIGEN_MINUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(pOcean->isChoppy() && bEnable[0])
			{
				#pragma omp parallel for
				for(int i = 0; i<count; ++i)
				{
					outData[i].PutX(pOcean->getOceanData(uCoord[i], vCoord[i], EIGENPLUSX));
					outData[i].PutY(0.0f);
					outData[i].PutZ(pOcean->getOceanData(uCoord[i], vCoord[i], EIGENPLUSZ));
				}
			}
			else
			{
				#pragma omp parallel for
				for(int i = 0; i<count; ++i){
					outData[i].PutX(0);outData[i].PutY(0);outData[i].PutZ(0);}
			}
		}
		break;

		case ID_OUT_EIGEN_PLUS:
		{
			CDataArrayVector3f outData( in_ctxt );
			if(pOcean->isChoppy() && bEnable[0])
			{
				#pragma omp parallel for
				for(int i = 0; i<count; ++i)
				{
					outData[i].PutX(pOcean->getOceanData(uCoord[i], vCoord[i], EIGENMINUSX));
					outData[i].PutY(0.0f);
					outData[i].PutZ(pOcean->getOceanData(uCoord[i], vCoord[i], EIGENMINUSZ));
				}
			}
			else
			{
				#pragma omp parallel for
				for(int i = 0; i<count; ++i){
					outData[i].PutX(0); outData[i].PutY(0); outData[i].PutZ(0);}
			}
		}
		break;
	}

	return CStatus::OK;
}

