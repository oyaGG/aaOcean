#include "aaOceanNormals.h"

extern "C" DLLEXPORT miBoolean aaOceanNormalsShader( miVector *result, miState *state, aaOceanNormalsShader_t *params)
{
	oceanStore** os = NULL;
	if(!mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os ))
		return( miTRUE);

	aaOceanClass *ocean = (*os)->MRocean;

	miVector *coord	 =  mi_eval_vector(&params->uv_coords);
	result->x = catromPrep(ocean,ocean->fft_normX,state,coord);
	result->y = catromPrep(ocean,ocean->fft_normY,state,coord);
	result->z = catromPrep(ocean,ocean->fft_normZ,state,coord);

	miVector *layerNormals	 =  mi_eval_vector(&params->layerNormals);
	result->x = (result->x + layerNormals->x) / 2;
	result->y = (result->y + layerNormals->y) / 2;
	result->z = (result->z + layerNormals->z) / 2;

	miScalar fade = *mi_eval_scalar(&params->fade);
	result->x = result->x * (1.0 - fade) + state->normal.x * fade;
	result->y = result->y * (1.0 - fade) + state->normal.y * fade;
	result->z = result->z * (1.0 - fade) + state->normal.z * fade;

	mi_vector_normalize(result);
	
	return( miTRUE );
}


extern "C" DLLEXPORT void aaOceanNormalsShader_init(miState *state, aaOceanNormalsShader_t *params,	miBoolean *inst_init_req)
{
	if( params == NULL )
	{
		// TODO: Shader global initialization code goes here (if needed)
		
		// Request a per-instance shader initialization as well (set to miFALSE if not needed)
		*inst_init_req = miTRUE;
	}
	else
	{
		miLock *lock;
		mi_query(miQ_FUNC_LOCK, state, miNULLTAG, &lock);
		mi_lock(*lock);

		miVector *layerNormals	 =  mi_eval_vector(&params->layerNormals);
 		oceanStore** os = NULL;
		if(!shader_init(os, state))
			return;	
		
		(*os)->MRocean = new aaOceanClass;
		aaOceanClass *ocean = (*os)->MRocean;	

		int resolution		= pow(2.0f,(4+ (*mi_eval_integer(&params->resolution))));
		ocean->oceanScale	= maximum<float>(*mi_eval_scalar(&params->oceanScale),0.00001);
		ocean->seed			= *mi_eval_integer(&params->seed);
		ocean->waveHeight	= *mi_eval_scalar(&params->waveHeight) * .01;
		ocean->velocity		= maximum<float>((((*mi_eval_scalar(&params->velocity))  * (*mi_eval_scalar(&params->velocity))) / (9.81)),0.00001);
		ocean->waveSpeed	= *mi_eval_scalar(&params->waveSpeed);
		ocean->chopAmount	= *mi_eval_scalar(&params->chopAmount) * .01;
		ocean->cutoff		= fabs(*mi_eval_scalar(&params->cutoff) / 100);
		ocean->windDir		= ((*mi_eval_scalar(&params->windDir))/180) * 3.141592653589793238462643383;
		ocean->windAlign	= maximum<int>	 (((*mi_eval_integer(&params->windAlign) + 1) * 2),2); 
		ocean->damp			= minimum<float>(*mi_eval_scalar(&params->damp),1);
		ocean->time			= *mi_eval_scalar(&params->time);
			
		if(!ocean->reInit( resolution ,true, 1))
		{
			(*os)->isValid = false;
			mi_error("aaOcean: Normals shader failed to initialize.");	
			return;
		}
		
		ocean->normalsResolution = resolution;	
		if(ocean->MRredoHoK)
			ocean->redoHoK = true;
		ocean->prepareOcean(0,1,1,0,1);
		ocean->MRredoHoK = false;

		copy_and_tile(ocean->fft_normX, ocean);
		rotateArray(ocean->fft_normX,	ocean);
		copy_and_tile(ocean->fft_normY, ocean);
		rotateArray(ocean->fft_normY,	ocean);
		copy_and_tile(ocean->fft_normZ, ocean);
		rotateArray(ocean->fft_normZ,   ocean);

		//clear arrays that are not required
		shaderCleanup(ocean);
		mi_warning("aaOcean: Normals shader initiated successfully at %dx%d", resolution, resolution);	

		mi_unlock(*lock);
	}
}

extern "C" DLLEXPORT void 
aaOceanNormalsShader_exit(	miState	*state,	aaOceanNormalsShader_t *params)
{
	if( params == NULL )
	{
		// TODO: Shader global cleanup code goes here (if needed)
	}
	else
	{
		// TODO: Shader instance-specific cleanup code goes here (if needed)
		oceanStore** os = NULL;
		if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os ))
		{
			if((*os)->MRocean)
				delete ((*os)->MRocean);
			mi_mem_release( (void *)(*os) );
			mi_warning("aaOcean: Normals shader terminated successfully");
		}
	}
}

extern "C" DLLEXPORT int
aaOceanNormalsShader_version( )
{
	return( 1 );
}
