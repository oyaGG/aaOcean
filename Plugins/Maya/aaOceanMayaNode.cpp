// aaOcean v2.5 Maya Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#include <string.h>
#include <maya/MIOStream.h>
#include <math.h>
#include <maya/MGlobal.h>

#include <maya/MPxDeformerNode.h> 
#include <maya/MItGeometry.h>
#include <maya/MPxLocatorNode.h> 

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMatrixData.h>

#include <maya/MFnPlugin.h>
#include <maya/MFnDependencyNode.h>

#include <maya/MTypeId.h> 
#include <maya/MPlug.h>

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>

#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>

#include <maya/MAnimControl.h>
#include <maya/MPointArray.h>

#include <maya/MThreadUtils.h>
#include <maya/MDagModifier.h>

#include <maya/MFloatMatrix.h>
#include <maya/MFloatArray.h>
#include <maya/MGeometryManager.h>
#include <maya/MGeometry.h>
#include <maya/MGeometryData.h>
#include <maya/MGeometryPrimitive.h>
#include <maya/MGeometry.h>
#include <maya/MGeometryData.h>
#include <maya/MFnMesh.h>

#include "aaOceanMayaNode.h"

MStatus aaOceanMaya::deform( MDataBlock& block,	MItGeometry& iter,	const MMatrix& mat, unsigned int multiIndex)
{
	// get maya points array
	MPointArray verts;
	int numVerts = iter.count();
	iter.allPositions(verts);

	// The following block of code is from Chad Vernon's site
	// It provides access to the geometry the deformer is applied to
	// so that we can query data on it, such as UVs
	MStatus status;
    MArrayDataHandle hInput = block.outputArrayValue( input, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status )
	if (status != MS::kSuccess) 		// Make sure we didn't fail.
	    return status ;
    status = hInput.jumpToElement( multiIndex );
    CHECK_MSTATUS_AND_RETURN_IT( status )
	if (status != MS::kSuccess) 		// Make sure we didn't fail.
	    return status ;
    MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
    MFnMesh mesh( oInputGeom );
	
	MFloatArray uCoord(numVerts,0.0);
	MFloatArray vCoord(numVerts,0.0);

	// only getting one UV for now
	// change this later to user-specified
	MString uvSetName;
	mesh.getCurrentUVSetName( uvSetName );
	mesh.getUVs(uCoord, vCoord, &uvSetName);

	// main input function
	pOcean->input(	block.inputValue(resolution).asInt(),
					block.inputValue(seed).asInt(),
					block.inputValue(oceanSize).asFloat(),
					block.inputValue(waveSize).asFloat(),
					block.inputValue(waveSmooth).asFloat(),
					block.inputValue(waveDirection).asFloat(),
					block.inputValue(waveAlign).asInt(),
					block.inputValue(waveReflection).asFloat(),
					block.inputValue(waveSpeed).asFloat(),
					block.inputValue(waveHeight).asFloat(),
					block.inputValue(waveChop).asFloat(),
					block.inputValue(currTime).asFloat(),
					block.inputValue(doFoam).asInt(),
					1);

	if(pOcean->isValid())
	{
		MPoint worldSpaceDisplacementVec;
		MPoint oceanLocalSpace;

		// the following matrix holds junk values
		MDataHandle matData = block.inputValue(inTransform);
		MMatrix transform = matData.asMatrix();

		#pragma omp parallel for private(worldSpaceDisplacementVec, oceanLocalSpace)
		for(int i = 0; i < verts.length(); i++)
		{
			// get height field
			worldSpaceDisplacementVec[1] = pOcean->getOceanData(uCoord[i], vCoord[i], HEIGHTFIELD);
			if(pOcean->isChoppy())
			{
				// get x and z displacement
				worldSpaceDisplacementVec[0] = pOcean->getOceanData(uCoord[i], vCoord[i], CHOPX);
				worldSpaceDisplacementVec[2] = pOcean->getOceanData(uCoord[i], vCoord[i], CHOPZ);
			}

			// oceanLocalSpace = worldSpaceDisplacementVec * transform;
			verts[i] += worldSpaceDisplacementVec;
		}
		iter.setAllPositions(verts);
		return MS::kSuccess;
	}
	else
		return MS::kFailure;
}


