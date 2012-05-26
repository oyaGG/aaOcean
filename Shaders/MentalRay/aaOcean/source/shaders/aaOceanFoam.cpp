// aaOcean Mental Ray Foam Shader
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#include "aaOceanFoam.h"

extern "C" DLLEXPORT miBoolean
aaOceanFoamShader( miColor *result, miState *state, aaOceanFoamShader_t *params)
{
	oceanStore** os;
	mi_query(miQ_FUNC_USERPTR, state, 0, (void *)&os);
	aaOcean *ocean = (*os)->ocean;
	
	if(ocean->m_chopAmount > 0.0f)
	{
		miVector	*coord		=  mi_eval_vector(&params->uv_coords);
		miScalar	gamma		= *mi_eval_scalar(&params->gamma);
		miScalar	brightness  = *mi_eval_scalar(&params->brightness);
		miScalar	fade		= *mi_eval_scalar(&params->fade);
		miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);
		miScalar	fmin		= *mi_eval_scalar(&params->fmin);
		miScalar	fmax		= *mi_eval_scalar(&params->fmax);

		miScalar foam = catromPrep(ocean, ocean->m_fft_jxz,state,coord);

		if(!rawOutput)
		{
			foam  = rescale(foam, fmin, fmax, 0.0f, 1.0f);
			foam  = maxim(foam, 0.0f);
			foam  = pow((float)foam, (float)gamma);
			foam *= brightness;
		}
		result->b = result->g = result->r = foam * (1.0f - fade);
	}
	else
		result->b = result->g = result->r = 0;

	return( miTRUE );
}

extern "C" DLLEXPORT void
aaOceanFoamShader_init(miState *state, aaOceanFoamShader_t *params,	miBoolean *inst_init_req)
{
	if( !params )
	{
		*inst_init_req = miTRUE;
		fftwf_init_threads();
		return;
	}
	Timer timer;
	timer.start();

	oceanStore** os = NULL;
	mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
	*os = (oceanStore *)mi_mem_allocate(sizeof(oceanStore));	
	(*os)->ocean = new aaOcean;
	aaOcean *ocean = (*os)->ocean;

	float velocity			= *mi_eval_scalar(&params->velocity);
	int resolution			= (int)pow(2.0f,(4 + (*mi_eval_integer(&params->resolution))));
	ocean->m_oceanScale		= maximum<float>(*mi_eval_scalar(&params->oceanScale), 1e-5f);
	ocean->m_seed			= *mi_eval_integer(&params->seed);
	ocean->m_waveHeight		= *mi_eval_scalar(&params->waveHeight) * 0.01f;
	ocean->m_velocity		= maximum<float>(((velocity * velocity) / (aa_GRAVITY)), 1e-5f);
	ocean->m_waveSpeed		= *mi_eval_scalar(&params->waveSpeed);
	ocean->m_chopAmount		= *mi_eval_scalar(&params->chopAmount) * 0.01f;
	ocean->m_cutoff			= fabs(*mi_eval_scalar(&params->cutoff) * 0.01f);
	ocean->m_windDir		= DegsToRads(*mi_eval_scalar(&params->windDir));
	ocean->m_windAlign		= maximum<int>(((*mi_eval_integer(&params->windAlign) + 1) * 2), 2); 
	ocean->m_damp			= minimum<float>(*mi_eval_scalar(&params->damp),1);
	ocean->m_time			= *mi_eval_scalar(&params->time);
	
	miTag shaderInst; // tag of the currently running shader
	mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);

	//start main shader work
	if((float)ocean->m_chopAmount > 0.0f)
	{		
		if(!ocean->reInit(resolution))
		{
			mi_error("%s. Shader ID: %d", ocean->m_state, shaderInst);	
			return;
		}
		mi_info("%s. Shader ID: %d", ocean->m_state, shaderInst);	

		ocean->prepareOcean(TRUE, TRUE, TRUE, FALSE);

		copy_and_tile(ocean->m_fft_jxz, ocean);
		rotateArray(ocean->m_fft_jxz, ocean);
		getArrayBounds(ocean->m_fft_jxz, 1, ocean->m_resolution, ocean->m_fmin, ocean->m_fmax);

		miScalar	fmin		= *mi_eval_scalar(&params->fmin);
		miScalar	fmax		= *mi_eval_scalar(&params->fmax);
		miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);

		float epsilon = 1e-3f;
		if(( !isfEqual(fmin, ocean->m_fmin, epsilon) || !isfEqual(fmax, ocean->m_fmax, epsilon) ) && !rawOutput)
			mi_warning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
						ocean->m_fmin, ocean->m_fmax);

		ocean->clearResidualArrays();
		timer.stop();

		mi_info("[aaOcean Shader] Foam shader initiated in %f seconds. Shader ID: %d, location: %p", 
				timer.getElapsedTimeInSec(), shaderInst, ocean);		
	}
	else
		mi_error("[aaOcean Shader] Choppiness is 0. Increase chopiness to see foam");
}

extern "C" DLLEXPORT void 
aaOceanFoamShader_exit(	miState	*state,	aaOceanFoamShader_t *params)
{
	if( params )
	{
		oceanStore** os = NULL;
		if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os ))
		{
			aaOcean *ocean = (*os)->ocean;
			if((*os)->ocean)
				delete ((*os)->ocean);
			mi_mem_release( (void *)(*os) );
			
			miTag shaderInst; // tag of the currently running shader
			mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
			mi_info("[aaOcean Shader] Displacement Shader ID %d, terminated at location: %p",shaderInst, ocean);
		}
	}
	else
	{
		//fftwf_cleanup_threads();
		//fftwf_cleanup();
	}
}

extern "C" DLLEXPORT int
aaOceanFoamShader_version( )
{
	return( 1 );
}
