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
	p_rotateUV
};

#include "shader_funcs.h"

node_parameters
{
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
}

node_update
{
	aaOcean *ocean = (aaOcean *)AiNodeGetLocalData(node);

	int renderRes = (int)pow(2.0f, params[p_resolution].INT + 4);

	ocean->input(renderRes,
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
		params[p_time].FLT);

	if(ocean->m_isValid)
	{
		bool chop = false;
		if((float)ocean->m_chopAmount > 0.0f)
			chop = true;

		ocean->prepareOcean(1,chop,chop,1,1);

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
	else // bad input to ocean class
		AiMsgWarning("%s", ocean->m_state);	
}

shader_evaluate
{
	aaOcean *ocean = (aaOcean*)AiNodeGetLocalData(node);

	if(ocean->m_isValid)
	{
		AtPoint uvPt;
		bool rotate = AiShaderEvalParamBool(p_rotateUV);

		if(AiShaderEvalParamBool(p_useUVInput))
		{
			uvPt = AiShaderEvalParamVec(p_uv_coords);
			// Rotating UV space to be in line with deformer
			if(rotate)
			{
				uvPt.z = uvPt.y * -1.0f;
				uvPt.y = uvPt.x;
				uvPt.x = uvPt.z;
			}
		}
		else
		{
			// Rotating UV space to be in line with deformer
			if(rotate)
			{
				uvPt.x = sg->v * -1.0f;
				uvPt.y = sg->u;
			}
			else
			{
				uvPt.x = sg->u;
				uvPt.y = sg->v;
			}
		}

#ifdef VECTORSHADER
		sg->out.VEC.y = catromPrep( ocean,  ocean->m_fft_htField, uvPt);
#else
		sg->out.RGBA.g = catromPrep( ocean,  ocean->m_fft_htField, uvPt);
#endif

		if((float)ocean->m_chopAmount > 0.0f)
		{
			#ifdef VECTORSHADER
			sg->out.VEC.x = -catromPrep( ocean, ocean->m_fft_chopX, uvPt);
			sg->out.VEC.z = -catromPrep( ocean, ocean->m_fft_chopZ, uvPt);
			#else
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
			#endif
		}
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

	if(AiNodeGetBool(node, "writeFile"))
	{
		const char *outputFolder = AiNodeGetStr(node, "outputFolder");
		if(!dirExists(outputFolder))
			AiMsgWarning("[aaOcean Arnold] Invalid folder path: %s", outputFolder);
		else
		{
			const char* postfix = AiNodeGetStr(node, "postfix");
			int frame = AiNodeGetInt(node, "currentFrame");
			char outputFileName[255];
			genFullFilePath(&outputFileName[0], &outputFolder[0], &postfix[0], frame);
			writeExr(&outputFileName[0], ocean);
			AiMsgInfo("[aaOcean Arnold] Image written to %s", outputFileName);
		}
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
   #ifdef VECTORSHADER
   node->output_type  = AI_TYPE_VECTOR;
   #else
   node->output_type  = AI_TYPE_RGBA;
   #endif
   node->name         = "aaOceanArnold";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}
