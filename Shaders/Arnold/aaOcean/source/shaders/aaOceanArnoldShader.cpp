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
	p_outputFileName
};

#include "shader_funcs.h"

node_parameters
{
	AiParameterVEC ( "uv_coords"		, 1.0f, 1.0f, 1.0f);
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
	AiParameterStr ( "outputFileName"   , "");
}

node_update
{
	aaOcean *ocean = (aaOcean *)AiNodeGetLocalData(node);
	if(!fetchInput(params, ocean))
		return;

	int renderRes = (int)pow(2.0f, params[p_resolution].INT + 4);
	if(ocean->reInit(renderRes))
	{
		ocean->prepareOcean(TRUE, TRUE, TRUE, FALSE, TRUE);

		if(ocean->m_chopAmount > 0.0f)
		{
			// get foam array bounds to generate normalization information
			getArrayBounds(ocean->m_fft_jxz, 1, ocean->m_resolution, ocean->m_fmin, ocean->m_fmax);

			float	fmin		= params[p_fMin].FLT;
			float	fmax		= params[p_fMax].FLT;
			float	rawOutput	= params[p_raw].BOOL;

			float epsilon = 1e-3f;
			if(( !isfEqual(fmin, ocean->m_fmin, epsilon) || !isfEqual(fmax, ocean->m_fmax, epsilon) ) && !rawOutput)
				AiMsgWarning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
							ocean->m_fmin, ocean->m_fmax);
		}
	}
	else // bad input to ocean class
		AiMsgInfo("%s", ocean->m_state);	
}

shader_evaluate
{
	aaOcean *ocean = (aaOcean*)AiNodeGetLocalData(node);

	if(ocean->m_isValid)
	{
		AtPoint uvPt;
		uvPt.x = sg->v * -1.0f;
		uvPt.y = sg->u;

		sg->out.RGBA.g = catromPrep( ocean,  ocean->m_fft_htField, uvPt);

		if((float)ocean->m_chopAmount > 0.0f)
		{
			sg->out.RGBA.r = -catromPrep( ocean, ocean->m_fft_chopX, uvPt);
			sg->out.RGBA.b = -catromPrep( ocean, ocean->m_fft_chopZ, uvPt);

			float	gamma		= AiShaderEvalParamFlt(p_gamma);
			float	brightness  = AiShaderEvalParamFlt(p_brightness);
			bool	rawOutput	= AiShaderEvalParamBool(p_raw);
			float	fmin		= AiShaderEvalParamFlt(p_fMin);
			float	fmax		= AiShaderEvalParamFlt(p_fMax);

			float foam = catromPrep(ocean, ocean->m_fft_jxz, uvPt);

			if(!rawOutput)
			{
				foam  = rescale(foam, fmin, fmax, 0.0f, 1.0f);
				foam  = maximum<float>(foam, FLT_MIN);
				foam  = pow(foam, gamma);
				foam *= brightness;
				foam *= (1.0f - AiShaderEvalParamFlt(p_fade));
			}
			sg->out.RGBA.a = foam;
		}
		else
			sg->out.RGBA.a = 0.0f;
	}
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

	const char *outputFileName = AiNodeGetStr(node, "outputFileName");
	if(strcmp(outputFileName, "") != 0)
	{
		// dump ocean data to open-exr file
		AiMsgInfo("[aaOcean Arnold] Written file to disk"); // log file path and name to console
	}

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
