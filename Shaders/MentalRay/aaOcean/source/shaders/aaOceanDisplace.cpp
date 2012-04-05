#include "aaOceanDisplace.h"

extern "C" DLLEXPORT 
miBoolean aaOceanDisplaceShader(miScalar *result, miState *state, aaOceanDisplaceShader_t *params)
{
	oceanStore** os;
	mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
	aaOcean *ocean = (*os)->ocean;
	
	miScalar	layerOcean	= *mi_eval_scalar(&params->layerOcean);
	miVector	*coord		=  mi_eval_vector(&params->uv_coords);
	miScalar	fade		= *mi_eval_scalar(&params->fade);	
	
	miVector point;	
	point.y = catromPrep(ocean, ocean->m_fft_htField,  state,  coord) * (1.0f - fade);
	if(ocean->m_chopAmount > 0.0f)
	{
		point.x = catromPrep( ocean, ocean->m_fft_chopX,  state,  coord)  *  (1.0f - fade);
		point.z = catromPrep( ocean, ocean->m_fft_chopZ,  state,  coord)  *  (1.0f - fade);
	}

	state->point.x -= point.x;
	state->point.y += point.y;
	state->point.z -= point.z;

	return( miTRUE );
}

extern "C" DLLEXPORT 
void aaOceanDisplaceShader_init(miState *state, aaOceanDisplaceShader_t *params, miBoolean *inst_init_req)
{
	if( !params )
	{
		*inst_init_req = miTRUE;
		fftwf_init_threads();
	}
	else
	{
		Timer timer;
		timer.start();

		miScalar layerOcean	= *mi_eval_scalar(&params->layerOcean);

		oceanStore** os = NULL;
		mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
		*os = (oceanStore *)mi_mem_allocate( sizeof(oceanStore));	
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

		if(!ocean->reInit(resolution))
		{
			ocean->m_isValid = FALSE;
			mi_error("[aaOcean] Displacement shader failed to initialize");	
			return;
		}
		ocean->m_redoHoK = TRUE;	
		ocean->prepareOcean(TRUE, TRUE, FALSE, FALSE);

		copy_and_tile(ocean->m_fft_htField, ocean);
		rotateArray(ocean->m_fft_htField, ocean);
		if(ocean->m_chopAmount > 0.0f)
		{
			copy_and_tile(ocean->m_fft_chopX, ocean);
			copy_and_tile(ocean->m_fft_chopZ, ocean);
			rotateArray(ocean->m_fft_chopX,   ocean);
			rotateArray(ocean->m_fft_chopZ,   ocean);
		}

		//clear arrays that are not required
		ocean->clearResidualArrays();
		timer.stop();

		miTag shaderInst; // tag of the currently running shader
		mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
		mi_info("[aaOcean] Displacement shader initiated in %f seconds, at %dx%d, Shader Instance ID: %d, location: %p", 
				timer.getElapsedTimeInSec(), resolution, resolution, shaderInst, ocean);	
	}
}

extern "C" DLLEXPORT 
void aaOceanDisplaceShader_exit(miState	*state,	aaOceanDisplaceShader_t *params)
{
	if( params )
	{
		oceanStore** os = NULL;
		if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os))
		{
			aaOcean *ocean = (*os)->ocean;
			if(ocean)
				delete ocean;
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

extern "C" DLLEXPORT int aaOceanDisplaceShader_version( )
{
	return( 1 );
}
