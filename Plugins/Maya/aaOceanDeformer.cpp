// aaOcean v2.5 Maya Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html


#include "aaOceanDeformer.h"

bool aaOceanDeformer::getUVs(MFnMesh &mesh, MDataBlock &block)
{
    MString userUV = block.inputValue(uvMap).asString();
    MStringArray uvSetNames;
    mesh.getUVSetNames(uvSetNames);

    int uvIndex = 0;
    bool uvFound = FALSE;
    int numUVs = uvSetNames.length();

    for(; uvIndex < numUVs; ++uvIndex)
    {
        if(userUV == uvSetNames[uvIndex])
        {
            uvFound = TRUE;
            break;
        }
    }
    if(!uvFound)
    {
        char errMsg[64];
        sprintf(errMsg, "[aaOcean] uvset \'%s\' not found on mesh", userUV.asChar());
        MGlobal::displayError( MString(errMsg));
        return FALSE;
    }
    else
    {
        mesh.getUVs(u, v, &uvSetNames[uvIndex]);
        return TRUE;
    }
}

void aaOceanDeformer::getColorSets(MFnMesh &mesh, MDataBlock &block)
{
    // query Color maps
    MString userEigenVec = block.inputValue(eigenVectorMap).asString();
    MString userEigenVal = block.inputValue(eigenValueMap).asString();
    MStringArray colSetNames;
    foundEigenVector = foundEigenValue = FALSE;
    MStatus status;

    int numcolors = mesh.numColorSets();
    mesh.getColorSetNames(colSetNames);
    if(numcolors == 0)
        MGlobal::displayInfo( "[aaOcean] No Color-At-Vertex Maps found. Skipping EigenVectors/Values");
    else
    {
        for(unsigned int i = 0; i < colSetNames.length(); ++i)
        {
            if(userEigenVec == colSetNames[i])
                foundEigenVector = TRUE;
            if(userEigenVal == colSetNames[i])
                foundEigenValue = TRUE;
        }
        if(foundEigenVector)
        {
            status = mesh.createColorSetDataMesh(userEigenVec);
            if(status != MStatus::kSuccess)
                MGlobal::displayInfo( status.errorString());
            colArrayEigenVector.setLength(mesh.numVertices());
        }
        else
            MGlobal::displayError( "[aaOcean] Specified Color Set for Eigen Vectors not found on mesh. Skipping Eigen Vectors");

        if(foundEigenValue)
        {
            status = mesh.createColorSetDataMesh(userEigenVal);
            if(status != MStatus::kSuccess)
                MGlobal::displayInfo( status.errorString());
            colArrayEigenValue.setLength(mesh.numVertices());
        }
        else
            MGlobal::displayError( "[aaOcean] Specified Color Set for Eigen Values not found on mesh. Skipping Eigen Values");
    }
}

void aaOceanDeformer::setColorSets(MFnMesh &mesh, MDataBlock &block)
{
    if(foundEigenVector || foundEigenValue)
    {
        int numPolys = mesh.numPolygons();
        if(numPolygons != numPolys)
            faceColorID.setLength(mesh.numFaceVertices());

        int index;
        MIntArray vertexList;
        
        #pragma omp parallel for private(vertexList, index)
        for(int i = 0; i < numPolys; ++i)
        {
            vertexList.clear();                 
            mesh.getPolygonVertices(i,vertexList);
            
            for (unsigned int j = 0; j < vertexList.length(); ++j)
            {
                mesh.getFaceVertexColorIndex(i, j, index);
                faceColorID[index] = vertexList[j];                         
            }
        }

        MStatus status;
        if(foundEigenVector)
        {
            MString colSet = block.inputValue(eigenVectorMap).asString();
            status = mesh.setColors(colArrayEigenVector, &colSet);
            status = mesh.assignColors(faceColorID, &colSet);
        }

        if(foundEigenValue)
        {
            MString colSet = block.inputValue(eigenValueMap).asString();
            status = mesh.setColors(colArrayEigenValue, &colSet);
            status = mesh.assignColors(faceColorID, &colSet);
        }
    }
}

MDataHandle aaOceanDeformer::getMeshHandle(const MPlug& plug, MDataBlock &block)
{
    // using compute() in deformer as described in
    // http://download.autodesk.com/us/maya/2011help/API/class_m_px_deformer_node.html
    // get the input corresponding to this output
    unsigned int index = plug.logicalIndex();
    MObject thisNode = this->thisMObject();
    MPlug inPlug(thisNode, input);
    inPlug.selectAncestorLogicalIndex(index, input);
    MDataHandle hInput = block.inputValue(inPlug);

    // get the input geometry and input groupId
    MDataHandle hGeom = hInput.child(inputGeom);
    MDataHandle hGroup = hInput.child(groupId);
    MDataHandle hOutput = block.outputValue(plug);
    hOutput.copy(hGeom);
    return hOutput;
}

MStatus aaOceanDeformer::compute(const MPlug& plug, MDataBlock& block)
{
    MStatus status = MStatus::kUnknownParameter;

    if (plug.attribute() == outputGeom) 
    {
        MDataHandle hOutput = getMeshHandle(plug, block);
        MFnMesh mesh(hOutput.asMesh(), &status);

        // if no UVs on mesh, return
        if(getUVs(mesh, block) == FALSE)
            return MStatus::kNotFound;
        
        // pull in some attribute values for convenience
        float timeOffsetValue = block.inputValue(timeOffset).asFloat();
        MTime appTime = block.inputValue(time).asTime();
        float currentTime = (float)appTime.as(MTime::kSeconds) + timeOffsetValue;
        bool foam = block.inputValue(doFoam).asBool();
        bool invert = block.inputValue(invertFoam).asBool();
        MMatrix transform = block.inputValue(inTransform).asMatrix();

        // main ocean input function
        pOcean->input(  block.inputValue(resolution).asInt(),
                        block.inputValue(seed).asInt(),
                        block.inputValue(oceanSize).asFloat(),
                        block.inputValue(oceanDepth).asFloat(),
                        block.inputValue(surfaceTension).asFloat(),
                        block.inputValue(waveSize).asFloat(),
                        block.inputValue(waveSmooth).asFloat(),
                        block.inputValue(waveDirection).asFloat(),
                        block.inputValue(waveAlign).asInt(),
                        block.inputValue(waveReflection).asFloat(),
                        block.inputValue(waveSpeed).asFloat(),
                        block.inputValue(waveHeight).asFloat(),
                        block.inputValue(waveChop).asFloat(),
                        currentTime,
                        block.inputValue(repeatTime).asFloat(),
                        foam,
                        FALSE);

        if(foam)
            getColorSets(mesh, block);

        MPoint worldSpaceVec;
        MPoint localSpaceVec;
        MPointArray verts;
        mesh.getPoints(verts);

        int numVertices = mesh.numVertices();
        float r, g, b, a = 0.f;
        #pragma omp parallel for private(worldSpaceVec, localSpaceVec, r, g, b, a)
        for(int i = 0; i < numVertices; ++i) 
        {
            // get height field
            worldSpaceVec[1] = pOcean->getOceanData(u[i], v[i], aaOcean::eHEIGHTFIELD);
            if(pOcean->isChoppy())
            {
                // get x and z displacement
                worldSpaceVec[0] = pOcean->getOceanData(u[i], v[i], aaOcean::eCHOPX);
                worldSpaceVec[2] = pOcean->getOceanData(u[i], v[i], aaOcean::eCHOPZ);

                if(foam)
                {
                    if(foundEigenVector)
                    {
                        r = pOcean->getOceanData(u[i], v[i], aaOcean::eEIGENMINUSX);
                        g = pOcean->getOceanData(u[i], v[i], aaOcean::eEIGENMINUSZ);
                        b = pOcean->getOceanData(u[i], v[i], aaOcean::eEIGENPLUSX);
                        a = pOcean->getOceanData(u[i], v[i], aaOcean::eEIGENPLUSZ);
                        colArrayEigenVector.set(i, r, g, b, a);
                    }
                    if(foundEigenValue)
                    {
                        if(invert)
                            r = 1.0f - pOcean->getOceanData(u[i], v[i], aaOcean::eFOAM);
                        else
                            r = pOcean->getOceanData(u[i], v[i], aaOcean::eFOAM);

                        colArrayEigenValue.set(i, r, r, r);
                    }
                }
            }

            localSpaceVec = worldSpaceVec * transform;
            verts[i] += localSpaceVec;
        }

        if(foam)
            setColorSets(mesh, block);
        mesh.setPoints(verts);
    }
    block.setClean(plug);
    return status;
}
