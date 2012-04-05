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
			if(fabs(fmin) != 0.0f && fabs(fmax) != 0.0f)
			{
				foam  = rescale(foam, fmin, fmax, 0.0f, 1.0f);
				foam  = maxim(foam, 0.0f);
				foam  = pow((float)foam, (float)gamma);
				foam *= brightness;
			}
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

	int resolution			= (int)pow(2.0f,(4 + (*mi_eval_integer(&params->resolution))));
	ocean->m_oceanScale		= maximum<float>(*mi_eval_scalar(&params->oceanScale),0.00001f);
	ocean->m_seed			= *mi_eval_integer(&params->seed);
	ocean->m_waveHeight		= *mi_eval_scalar(&params->waveHeight) * .01f;
	ocean->m_velocity		= maximum<float>((((*mi_eval_scalar(&params->velocity)) * (*mi_eval_scalar(&params->velocity))) / (aa_GRAVITY)),0.00001f);
	ocean->m_waveSpeed		= *mi_eval_scalar(&params->waveSpeed);
	ocean->m_chopAmount		= *mi_eval_scalar(&params->chopAmount) * 0.01f;
	ocean->m_cutoff			= fabs(*mi_eval_scalar(&params->cutoff) * 0.01f);
	ocean->m_windDir		= DegsToRads(*mi_eval_scalar(&params->windDir));
	ocean->m_windAlign		= maximum<int>(((*mi_eval_integer(&params->windAlign) + 1) * 2), 2); 
	ocean->m_damp			= minimum<float>(*mi_eval_scalar(&params->damp),1);
	ocean->m_time			= *mi_eval_scalar(&params->time);
	
	//start main shader work
	if((float)ocean->m_chopAmount > 0.0f)
	{		
		if(!ocean->reInit(resolution))
		{
			ocean->m_isValid = FALSE;
			mi_error("[aaOcean] Foam shader failed to initialize.");	
			return;
		}

		ocean->m_redoHoK = TRUE;		
		ocean->prepareOcean(TRUE, TRUE, TRUE, FALSE);
		copy_and_tile(ocean->m_fft_jxz, ocean);
		rotateArray(ocean->m_fft_jxz, ocean);
		
		getArrayBounds(ocean->m_fft_jxz,1,ocean->m_resolution, ocean->m_fmin, ocean->m_fmax);
		mi_warning("[aaOcean] Foam array bounds for aaOceanFoam shader, min: %f, max: %f", ocean->m_fmin, ocean->m_fmax);

		ocean->clearResidualArrays();
		timer.stop();

		miTag shaderInst; // tag of the currently running shader
		mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
		mi_info("[aaOcean] Foam shader initiated in %f seconds, at %dx%d, Shader Instance ID: %d, location: %p", 
				timer.getElapsedTimeInSec(), resolution, resolution, shaderInst, ocean);		
	}
	else
		mi_error("[aaOcean] Choppiness is 0. Increase chopiness to see foam");
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
			mi_info("[aaOcean] Displacement Shader Instance ID %d, terminated at location: %p",shaderInst, ocean);
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
