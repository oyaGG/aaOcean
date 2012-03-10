#include "ai.h"
#include "aaOcean.h"
#include "utility_funcs.h"
#include <limits>

AI_SHADER_NODE_EXPORT_METHODS(aaOceanArnFoamMethods);

enum aaOceanArnFoamParams
{
	p_renderResolution,
	p_oceanSize,
	p_seed,
	p_waveHeight,
	p_waveSize,
	p_waveSpeed,
	p_waveChop,
	p_waveSmooth,
	p_waveDirection,
	p_waveReflection,
	p_waveAlign,
	p_uv_coords,
	p_time,
	p_gamma,
	p_brightness,
	p_fade,
	p_raw,
	p_foamMin,
	p_foamMax,
};

#include "shader_funcs.h"

node_parameters
{
	AiParameterINT ( "renderResolution"	, 4 );
	AiParameterFLT ( "oceanSize"		, 100.0f );
	AiParameterINT ( "seed"				, 1 );
	AiParameterFLT ( "waveHeight"		, 5.0f);
	AiParameterFLT ( "waveSize"			, 4.5f);
	AiParameterFLT ( "waveSpeed"		, 1.0f);
	AiParameterFLT ( "waveChop"			, 1.0f);
	AiParameterFLT ( "waveSmooth"		, 0.0f);
	AiParameterFLT ( "waveDirection"	, 45.0f);
	AiParameterFLT ( "waveReflection"	, 0.985f);
	AiParameterINT ( "waveAlign"		, 0);
	AiParameterVEC ( "uv_coords",		1.0f, 1.0f, 1.0f);
	AiParameterFLT ( "time"				, 0.1f);
	AiParameterFLT ( "gamma"			, 1.0f);
	AiParameterFLT ( "brightness"		, 1.0f);
	AiParameterFLT ( "fade"				, 0.0f);
	AiParameterBOOL( "raw"				, 0);
	AiParameterFLT ( "foamMin"			, 0.0f);
	AiParameterFLT ( "foamMax"			, 0.0f);
}

node_update
{
	aaOceanClass *pOcean= (aaOceanClass *)AiNodeGetLocalData(node);
	//get user input
	AtBoolean rawOutput	= params[p_raw].BOOL;
	bool isChanged = fetchInput(params,pOcean);
	int renderRes = params[p_renderResolution].INT  + 4;
	renderRes = (int)pow(2.0f,renderRes);

	//initialize ocean with user input
	if(pOcean->reInit(renderRes))
	{
		if(isChanged)
		{
			pOcean->isChoppy = 0;
			if((float)pOcean->chopAmount > 0.0f)
				pOcean->isChoppy = 1;
			
			pOcean->prepareOcean(pOcean->redoHoK,1,1,0);
			
			//prepare ocean arrays to display in arnold
			//we copy, tile and rotate to make the array compatible with arnold storage format
			copy_and_tile(pOcean->fft_htField,pOcean);
			if(pOcean->isChoppy) 
				copy_and_tile(pOcean->fft_chopX,pOcean);
			if(pOcean->isChoppy) 
				copy_and_tile(pOcean->fft_chopZ,pOcean);
			
			#pragma omp sections
			{
				#pragma omp section
				{
					rotateArray(pOcean->fft_htField,pOcean);
				}
				#pragma omp section
				{
					if(pOcean->isChoppy) rotateArray(pOcean->fft_chopX,pOcean);
				}
				#pragma omp section
				{
					if(pOcean->isChoppy) rotateArray(pOcean->fft_chopZ,pOcean);
				}
			}
			copy_and_tile(pOcean->fft_jxz,pOcean);		
			rotateArray(pOcean->fft_jxz,pOcean);
			getArrayBounds(pOcean->fft_jxz,1,pOcean->resolution, pOcean->foamMin , pOcean->foamMax);
			AiMsgWarning("[aaOcean Foam] New foam bounds: min = %f, max = %f", pOcean->foamMin, pOcean->foamMax);
		}

		if(!rawOutput)
		{
			 if( params[p_foamMin].FLT != 0.0f || params[p_foamMax].FLT != 0.0f )
			 {	
				 pOcean->foamMin = params[p_foamMin].FLT;
				 pOcean->foamMax = params[p_foamMax].FLT;
			 }
			 else 
				getArrayBounds(pOcean->fft_jxz,1,pOcean->resolution, pOcean->foamMin , pOcean->foamMax);
		}
	}

}

node_initialize
{
	aaOceanClass *pOcean;
	pOcean = new aaOceanClass;
	sprintf(pOcean->name,"Foam");
	AiMsgInfo("[aaOcean Foam] Created new aaOcean Foam data");
	AiNodeSetLocalData(node, pOcean);
}

node_finish
{
	delete (aaOceanClass*) AiNodeGetLocalData(node);
	AiMsgInfo("[aaOcean Foam] Deleted aaOcean data");		
}

shader_evaluate
{
	aaOceanClass *pOcean = (aaOceanClass *)AiNodeGetLocalData(node);
	
	AtFloat	gamma		= AiShaderEvalParamFlt(p_gamma);
	AtFloat	brightness  = AiShaderEvalParamFlt(p_brightness);
	AtFloat	fade		= AiShaderEvalParamFlt(p_fade);
	AtBoolean rawOutput	= AiShaderEvalParamBool(p_raw);
	AtVector myuv		= AiShaderEvalParamVec(p_uv_coords);

	if(pOcean->isValid)
	{
		AtPoint uvPt;
		AtFloat foam;
		uvPt.x = sg->u;
		uvPt.y = sg->v;
		foam = catromPrep( pOcean,  pOcean->fft_jxz, uvPt);

		if(!rawOutput)
		{
			// apply colour correction
			foam  = ( foam + fabs( pOcean->foamMin ) + 1.0f ) / ( fabs( pOcean->foamMin ) + ( pOcean->foamMax ) + 1.1f) ; 
			foam  = max(foam,0.0f);
			foam  = min(foam,1.0f);
			foam  = pow(foam, gamma);
			foam *= brightness;	
			foam *= (1.0f - fade);
		}
		
		sg->out.FLT = foam;
	}
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaOceanArnFoamMethods;
   node->output_type  = AI_TYPE_FLOAT;
   node->name         = "aaOceanArnFoam";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}