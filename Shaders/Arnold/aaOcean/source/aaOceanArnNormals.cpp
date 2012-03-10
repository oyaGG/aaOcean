#include "ai.h"
#include "aaOcean.h"
#include "utility_funcs.h"
#include <cmath>
#include <limits>

AI_SHADER_NODE_EXPORT_METHODS(aaOceanArnNormalsMethods);

enum aaOceanArnNormalsParams
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

node_initialize
{
	aaOceanClass *pOcean = NULL;
	//AiMsgWarning("[aaOceanFoam] In node-init");
	if (node->local_data == NULL) 
	{
		pOcean = new aaOceanClass;
		AiMsgWarning("[aaOceanNormals] Created new aaOcean Normals data");
		node->local_data = pOcean;
	}
	else
		pOcean = (aaOceanClass *)node->local_data;

	// Initializing Critical Section
	AiCritSecInit(&pOcean->critical_section);

	fetchInput(params,pOcean);
	int renderRes = params[p_renderResolution].INT  + 4;
	renderRes = (int)pow(2.0f,renderRes);

	if(pOcean->reInit(renderRes,true,false))
	{
		pOcean->normalsResolution = renderRes;
		pOcean->prepareOcean(pOcean->redoHoK,1,0,1);
		copy_and_tile(pOcean->fft_normX, pOcean);
		rotateArray(pOcean->fft_normX,	pOcean);
		copy_and_tile(pOcean->fft_normY, pOcean);
		rotateArray(pOcean->fft_normY,	pOcean);
		copy_and_tile(pOcean->fft_normZ, pOcean);
		rotateArray(pOcean->fft_normZ,   pOcean);
	}
	else
		pOcean->isValid = FALSE;
}

node_finish
{
	aaOceanClass *pOcean = (aaOceanClass*) node->local_data;
	if (!destroy_node) 
	{
		node->initialized =1;
	}
	else 
	{
		// Destroying Critical Section
		AiCritSecClose(&pOcean->critical_section);
		delete pOcean;
		AiMsgInfo("[aaOceanFoam] Deleted aaOcean data");		
	}
}

shader_evaluate
{
	AtVector myuv = AiShaderEvalParamVec(p_uv_coords);
	aaOceanClass *pOcean= (aaOceanClass *)node->local_data;

	AiCritSecEnter(&pOcean->critical_section);
		if(fetchInputEval(sg, node,pOcean)==1 && pOcean->isArnoldInit == 0)
		{
			int renderRes = AiShaderEvalParamInt(p_renderResolution)  + 4;
			renderRes = (int)pow(2.0f,renderRes);

			if(pOcean->reInit(renderRes,true,false))
			{
				pOcean->normalsResolution = renderRes;
				pOcean->prepareOcean(pOcean->redoHoK,1,0,1);
				
				copy_and_tile(pOcean->fft_normX, pOcean);
				rotateArray(pOcean->fft_normX,	pOcean);
				copy_and_tile(pOcean->fft_normY, pOcean);
				rotateArray(pOcean->fft_normY,	pOcean);
				copy_and_tile(pOcean->fft_normZ, pOcean);
				rotateArray(pOcean->fft_normZ,   pOcean);
			}
		}
		pOcean->isArnoldInit = 0;
	AiCritSecLeave(&pOcean->critical_section);	

	if(pOcean->isValid)
	{
		AtPoint uvPt;
		uvPt.x = sg->u;
		uvPt.y = sg->v;
		sg->Ns.x = sg->N.x = sg->out.VEC.x = catromPrep( pOcean, pOcean->fft_normX, uvPt);
		sg->Ns.y = sg->N.y = sg->out.VEC.y = catromPrep( pOcean, pOcean->fft_normY, uvPt);
		sg->Ns.z = sg->N.z = sg->out.VEC.z = catromPrep( pOcean, pOcean->fft_normZ, uvPt);		
	}
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaOceanArnNormalsMethods;
   node->output_type  = AI_TYPE_VECTOR;
   node->name         = "aaOceanArnNormals";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}