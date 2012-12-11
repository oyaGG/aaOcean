// aaOcean Arnold Shader
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
// Run RGB through a Color to Vector and feed into sta_displacement
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#include <limits>
#include <string>

AI_SHADER_NODE_EXPORT_METHODS(aaOceanArnoldMethods);

enum aaOceanArnoldParams
{
	p_uv_coords,
	p_useUVInput,
	p_fade,
	p_resolution,
	p_oceanScale,
	p_seed,
	p_waveHeight,
	p_velocity,
	p_waveSpeed,
	p_chopAmount,
	p_cutoff,
	p_windDir,
	p_damp,
	p_windAlign,
	p_time,
	p_gamma,
	p_brightness,
	p_raw,
	p_fMin,
	p_fMax,
	p_writeFile,
	p_outputFolder,
	p_postfix,
	p_currentFrame,
	p_rotateUV,
	p_transform
};

//#include "shader_funcs.h"


node_parameters
{
	AtMatrix matrix44;
	AiM4Identity(matrix44);
	AiParameterVEC ( "uv_coords"		, 1.0f, 1.0f, 1.0f);
	AiParameterBOOL( "use_uv_input"		, 0);
	AiParameterFLT ( "fade"				, 0.0f);
	AiParameterINT ( "resolution"		, 4);
	AiParameterFLT ( "oceanScale"		, 100.0f);
	AiParameterINT ( "seed"				, 1);
	AiParameterFLT ( "waveHeight"		, 5.0f);
	AiParameterFLT ( "velocity"			, 4.5f);
	AiParameterFLT ( "waveSpeed"		, 1.0f);
	AiParameterFLT ( "chopAmount"		, 1.0f);
	AiParameterFLT ( "cutoff"			, 0.0f);
	AiParameterFLT ( "windDir"			, 45.0f);
	AiParameterFLT ( "damp"				, 0.985f);
	AiParameterINT ( "windAlign"		, 0);
	AiParameterFLT ( "time"				, 0.1f);
	AiParameterFLT ( "gamma"			, 1.0f);
	AiParameterFLT ( "brightness"		, 1.0f);
	AiParameterBOOL( "raw"				, 0);
	AiParameterFLT ( "fMin"				, 0.0f);
	AiParameterFLT ( "fMax"				, 0.0f);
	AiParameterBOOL( "writeFile"		, 0);
	AiParameterStr ( "outputFolder"     , "");
	AiParameterStr ( "postfix"			, "");
	AiParameterINT ( "currentFrame"		, 1);
	AiParameterBOOL( "rotateUV"			, 0);
	AiParameterMTX ( "transform"		, matrix44);
}

node_update
{
	aaOcean *ocean = (aaOcean *)AiNodeGetLocalData(node);

	// main input function
	ocean->input(params[p_resolution].INT,
		params[p_seed].INT,
		params[p_oceanScale].FLT,
		params[p_velocity].FLT,
		params[p_cutoff].FLT,
		params[p_windDir].FLT,
		params[p_windAlign].INT,
		params[p_damp].FLT,
		params[p_waveSpeed].FLT,
		params[p_waveHeight].FLT,
		params[p_chopAmount].FLT,
		params[p_time].FLT,
		TRUE);

	if(!ocean->isValid())
	{
		AiMsgWarning("%s", ocean->getState());	
		return;
	}
	
	// see if user has requested normalized or raw foam values
	float rawOutput	= params[p_raw].BOOL;
	if(ocean->isChoppy() && !rawOutput)
	{
		float outMin, outMax;
		float	fmin		= params[p_fMin].FLT;
		float	fmax		= params[p_fMax].FLT;
		ocean->getFoamBounds(fmin, fmax, outMin, outMax);

		float epsilon = 1e-3f;
		if(!isfEqual(fmin, outMin, epsilon) || !isfEqual(fmax, outMax, epsilon) )
			AiMsgWarning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
						outMin, outMax);
	}
	//clear arrays that are not required
	ocean->clearResidualArrays();
}

shader_evaluate
{
	aaOcean *ocean = (aaOcean*)AiNodeGetLocalData(node);

	if(!ocean->isValid())
		return;

	// get our UV's
	AtPoint uvPt;
	if(AiShaderEvalParamBool(p_useUVInput))
		uvPt = AiShaderEvalParamVec(p_uv_coords);
	else
	{
		uvPt.x = sg->u;
		uvPt.y = sg->v;
	}

	// retrieve heightfield in world space
	AtPoint oceanWorldSpace;
	oceanWorldSpace.y = ocean->getOceanData(uvPt.x, uvPt.y, HEIGHTFIELD);
	oceanWorldSpace.x = oceanWorldSpace.z = 0.0f;

	if(ocean->isChoppy())
	{
		// retrieve chop displacement
		oceanWorldSpace.x = ocean->getOceanData(uvPt.x, uvPt.y, CHOPX);
		oceanWorldSpace.z = ocean->getOceanData(uvPt.x, uvPt.y, CHOPZ);

		// retrieve foam
		float foam = ocean->getOceanData(uvPt.x, uvPt.y, FOAM);

		// see if user has requested normalized or raw foam values
		if(AiShaderEvalParamBool(p_raw))
			sg->out.RGBA.a = foam;
		else
		{
			// get normalization weights
			float	fmin = AiShaderEvalParamFlt(p_fMin);
			float	fmax = AiShaderEvalParamFlt(p_fMax);
			
			// fitting to 0 - 1 range using rescale
			foam  = rescale(foam, fmin, fmax, 0.0f, 1.0f);
			// removing negative leftovers
			foam  = maximum<float>(foam, 0.0f);
			// gamma
			foam  = pow(foam, AiShaderEvalParamFlt(p_gamma));
			foam *= AiShaderEvalParamFlt(p_brightness);
		}
	}
	else
		sg->out.RGBA.a = 0.f;

	// get space-transform matrix
	AtMatrix matrix;
	AiM4Identity(matrix);
	AiM4Copy(matrix, *AiShaderEvalParamMtx(p_transform));

	// local space point
	AtPoint oceanLocalSpace;
	oceanLocalSpace.x = oceanLocalSpace.z = 0.f;

	// convert to local space
	AiM4PointByMatrixMult(&oceanLocalSpace, matrix, &oceanWorldSpace);

	// store result in output
	sg->out.RGBA.r = oceanLocalSpace.x;
	sg->out.RGBA.g = oceanLocalSpace.y;
	sg->out.RGBA.b = oceanLocalSpace.z;
}

node_initialize
{
	fftwf_init_threads();

	aaOcean *ocean;
	ocean = new aaOcean;
	AiMsgInfo("[aaOcean Arnold] Created new aaOcean data");
	AiNodeSetLocalData(node, ocean);
}

node_finish
{
	aaOcean *ocean = (aaOcean *)AiNodeGetLocalData(node);

	// cleanup
	delete ocean;
	AiMsgInfo("[aaOcean Arnold] Deleted aaOcean data");

	fftwf_cleanup_threads();
	fftwf_cleanup();
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaOceanArnoldMethods;
   node->output_type  = AI_TYPE_RGBA;
   node->name         = "aaOceanArnold";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}
