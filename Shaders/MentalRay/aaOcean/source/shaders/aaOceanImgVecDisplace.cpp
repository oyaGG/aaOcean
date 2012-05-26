// aaOcean Vector Displacement shader
// Now redundant because Mental Ray ships with one of its own
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html			

#include "aaOceanImgVecDisplace.h"

extern "C" DLLEXPORT miBoolean aaOceanImgVecDisplace(miScalar *result, miState *state, aaOceanImgVecDisplace_t *params )
{
	// TODO: Shader main code goes here
	miTag	 tex		  = *mi_eval_tag(&params->texture);
	miVector *coord		  =  mi_eval_vector(&params->uv_coords);
	miScalar scale_height = *mi_eval_scalar(&params->scale_height);
	miScalar scale_chop	  = *mi_eval_scalar(&params->scale_chop);
	miScalar fade		  = *mi_eval_scalar(&params->fade);
	
	miColor color;
	mi_lookup_color_texture(&color, state, tex, coord);

	miVector vec;
	vec.x = color.r;
	vec.y = color.g;
	vec.z = color.b;

	mi_point_from_object(state,&vec,&vec);  // convert object-space vector from texture to internal space

	state->point.x += vec.x * scale_chop	*	(1.0f - fade); 
	state->point.y += vec.y * scale_height	*	(1.0f - fade); 
	state->point.z += vec.z * scale_chop	*	(1.0f - fade); 

	return( miTRUE );
}

extern "C" DLLEXPORT void aaOceanImgVecDisplace_init( miState *state, aaOceanImgVecDisplace_t *params, miBoolean *inst_init_req)
{
	if( params == NULL )
		*inst_init_req = miTRUE;
	else
	{
		// TODO: Shader instance-specific initialization code goes here (if needed)
	}
}

extern "C" DLLEXPORT void aaOceanImgVecDisplace_exit(miState *state, aaOceanImgVecDisplace_t *params)
{
	if( params == NULL )
	{
		// TODO: Shader global cleanup code goes here (if needed)
	}
	else
	{
		// TODO: Shader instance-specific cleanup code goes here (if needed)
	}
}

extern "C" DLLEXPORT int aaOceanImgVecDisplace_version( )
{
	return( 1 );
}
