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
	MPointArray verts;
	iter.allPositions(verts);

	int gridRes = int_sqrt(verts.length()) - 1;	
	fetchInput(block);
	if(pOcean->reInit(gridRes))
	{
		bool chop = false;
		if((float)pOcean->m_chopAmount > 0.0001f)
			chop = true;

		pOcean->prepareOcean(1,chop,0,0,1,1);

		#pragma omp parallel for
		for(int i = 0; i<verts.length(); i++)
		{
			/*
			// this block of code for vec(x, y) -> vec(-y, x)
			// fix this later
			verts[i].y -= pOcean->m_fft_htField[i][1]; // why flip y?
			verts[i].x += pOcean->m_fft_chopZ[i][1];
			verts[i].z -= pOcean->m_fft_chopX[i][1];*/

			verts[i].y += pOcean->m_fft_htField[i][1];
			verts[i].x -= pOcean->m_fft_chopZ[i][1];
			verts[i].z -= pOcean->m_fft_chopX[i][1];
		}

		iter.setAllPositions(verts);
	}

	return MS::kSuccess;
}


