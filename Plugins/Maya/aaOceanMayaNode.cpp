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

#include "aaOceanMayaNode.h"

MStatus aaOceanMaya::deform( MDataBlock& block,	MItGeometry& iter,	const MMatrix& mat, unsigned int multiIndex)
{
	// get maya points array
	MPointArray verts;
	iter.allPositions(verts);

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
		float worldSpaceDisplacementVec[3] = {0.0f, 0.0f, 0.0f};
		float transformedVector[3];

		#pragma omp parallel for private(worldSpaceDisplacementVec, transformedVector)
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
		}
		iter.setAllPositions(verts);
		return MS::kSuccess;
	}
	else
		return MS::kFailure;
}


