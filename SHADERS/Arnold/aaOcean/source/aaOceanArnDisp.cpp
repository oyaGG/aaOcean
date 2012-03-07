#include "ai.h"
#include "aaOcean.h"
#include "utility_funcs.h"

AI_SHADER_NODE_EXPORT_METHODS(aaOceanArnDispMethods);

enum aaOceanArnDispParams
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
}

node_update
{
	aaOceanClass *pOcean = (aaOceanClass*)AiNodeGetLocalData(node);
	//get user input
	bool isChanged = fetchInput(params,pOcean);
	int renderRes = params[p_renderResolution].INT  + 4;
	renderRes = (int)pow(2.0f,renderRes);

	//initialize ocean with user input
	if(pOcean->reInit(renderRes))
	{
		if(isChanged)
		{
			AiMsgWarning("displacement changed!");
			pOcean->isChoppy = 0;
			if((float)pOcean->chopAmount > 0.0f)
				pOcean->isChoppy = 1;
			
			pOcean->prepareOcean(pOcean->redoHoK,1,0,0);
			
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
		}
	}
}

node_initialize
{
	aaOceanClass *pOcean;
	pOcean = new aaOceanClass;
	sprintf(pOcean->name, "Displacement");
	AiMsgInfo("[aaOcean Displacement] Created new aaOcean data");
	AiNodeSetLocalData(node, pOcean);
}

node_finish
{
	delete (aaOceanClass*)AiNodeGetLocalData(node);
	AiMsgInfo("[aaOcean Displacement] Deleted aaOcean data");
}

shader_evaluate
{
	AtVector myuv = AiShaderEvalParamVec(p_uv_coords);

	aaOceanClass *pOcean = (aaOceanClass *)AiNodeGetLocalData(node);
	if(pOcean->isValid)
	{
		AtPoint uvPt;
		uvPt.x = sg->u;
		uvPt.y = sg->v;
		sg->out.VEC.y = catromPrep( pOcean,  pOcean->fft_htField, uvPt);
		if((float)pOcean->chopAmount > 0.0)
		{
			sg->out.VEC.x = -catromPrep( pOcean, pOcean->fft_chopX,  uvPt);
			sg->out.VEC.z = -catromPrep( pOcean, pOcean->fft_chopZ,  uvPt);
		}
	}
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaOceanArnDispMethods;
   node->output_type  = AI_TYPE_VECTOR;
   node->name         = "aaOceanArnDisp";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}