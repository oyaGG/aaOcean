// aaOcean v2.5 Maya Deformer
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#include "aaOceanClass.cpp"

//Maya Node ID 0x20B6EF34  (548859700)  -- Randomly generated. Change if this conflicts

class aaOceanMaya : public MPxDeformerNode
{
public:
						aaOceanMaya();
	virtual				~aaOceanMaya();
	static  void*		creator();
	static  MStatus		initialize();

    virtual MStatus  deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex);

	static  MObject  viewRes;	
	static  MObject  oceanSize;	
	static  MObject  seed;	

	static  MObject  waveHeight;
	static  MObject  waveSize;
	static  MObject  waveSpeed;
	static  MObject  waveChop;
	static  MObject  waveSmooth;

	static  MObject  waveDirection;
	static  MObject  waveReflection;
	static  MObject  waveAlign;
	static  MObject  currTime;
	static  MTypeId	 id;

	static aaOcean* pOcean;

	void fetchInput(MDataBlock& block);

};

MObject		aaOceanMaya::viewRes;	
MObject		aaOceanMaya::oceanSize;	
MObject		aaOceanMaya::seed;	
MObject		aaOceanMaya::waveHeight;
MObject		aaOceanMaya::waveSize;
MObject		aaOceanMaya::waveSpeed;
MObject		aaOceanMaya::waveChop;
MObject		aaOceanMaya::waveSmooth;
MObject		aaOceanMaya::waveDirection;
MObject		aaOceanMaya::waveReflection;
MObject		aaOceanMaya::waveAlign;
MObject		aaOceanMaya::currTime;
MTypeId     aaOceanMaya::id( 0x20B6EF34 ); //Maya Node ID 548859700

aaOcean* aaOceanMaya::pOcean;


MStatus aaOceanMaya::initialize()
{
	MFnNumericAttribute nAttrViewRes;
	viewRes = nAttrViewRes.create( "Resolution", "viewRes", MFnNumericData::kInt, 2 );
    nAttrViewRes.setKeyable (true);	
	nAttrViewRes.setWritable(true);
	nAttrViewRes.setSoftMin(1);	
	nAttrViewRes.setSoftMax(6);	
	nAttrViewRes.setMin(1);	
	nAttrViewRes.setMax(6);	
    addAttribute( viewRes );
    attributeAffects( aaOceanMaya::viewRes, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrOceanSize;
	oceanSize= nAttrOceanSize.create( "oceanSize", "oceanSize", MFnNumericData::kFloat, 100.f );
    nAttrOceanSize.setKeyable(  true );	
	nAttrOceanSize.setWritable(true);
	nAttrOceanSize.setSoftMin(1.f);	
	nAttrOceanSize.setSoftMax(1000.f);		
    addAttribute( oceanSize );
    attributeAffects( aaOceanMaya::oceanSize, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrSeed;
	seed = nAttrSeed.create( "Seed", "seed", MFnNumericData::kInt, 1 );
    nAttrSeed.setKeyable(  true );	
	nAttrSeed.setWritable(true);
	nAttrSeed.setSoftMin(1);	
	nAttrSeed.setSoftMax(32);	
	nAttrSeed.setMin(1);	
	nAttrSeed.setMax(64);	
    addAttribute( seed );
    attributeAffects( aaOceanMaya::seed, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveHeight;
	waveHeight = nAttrWaveHeight.create( "waveHeight", "waveHeight", MFnNumericData::kFloat, 2.0f );
    nAttrWaveHeight.setKeyable(  true );	
	nAttrWaveHeight.setWritable(true);
	nAttrWaveHeight.setSoftMin(0.001f);	
	nAttrWaveHeight.setSoftMax(100.f);	
	nAttrWaveHeight.setMin(0.001f);		
    addAttribute( waveHeight );
    attributeAffects( aaOceanMaya::waveHeight, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveSize;
	waveSize = nAttrWaveSize.create( "waveSize", "waveSize", MFnNumericData::kFloat, 4.0f );
    nAttrWaveSize.setKeyable(  true );	
	nAttrWaveSize.setWritable(true);
	nAttrWaveSize.setSoftMin(1.f);	
	nAttrWaveSize.setSoftMax(30.f);		
    addAttribute( waveSize );
    attributeAffects( aaOceanMaya::waveSize, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveSpeed;
	waveSpeed = nAttrWaveSpeed.create( "waveSpeed", "waveSpeed", MFnNumericData::kFloat, 1.0f );
    nAttrWaveSpeed.setKeyable(  true );	
	nAttrWaveSpeed.setWritable(true);
	nAttrWaveSpeed.setSoftMin(1.f);	
	nAttrWaveSpeed.setSoftMax(10.f);	
    addAttribute( waveSpeed );
    attributeAffects( aaOceanMaya::waveSpeed, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveChop;
	waveChop = nAttrWaveChop.create( "waveChop", "waveChop", MFnNumericData::kFloat, 2.0f );
    nAttrWaveChop.setKeyable(  true );	
	nAttrWaveChop.setWritable(true);
	nAttrWaveChop.setSoftMin(0.0f);	
	nAttrWaveChop.setSoftMax(6.0f);
    addAttribute( waveChop );
    attributeAffects( aaOceanMaya::waveChop, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveSmooth;
	waveSmooth = nAttrWaveSmooth.create( "waveSmooth", "waveSmooth", MFnNumericData::kFloat, 0.0f);
    nAttrWaveSmooth.setKeyable(  true );	
	nAttrWaveSmooth.setWritable(true);
	nAttrWaveSmooth.setSoftMin(0.0f);	
	nAttrWaveSmooth.setSoftMax(20.f);
    addAttribute( waveSmooth );
    attributeAffects( aaOceanMaya::waveSmooth, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveDirection;
	waveDirection = nAttrWaveDirection.create( "waveDirection", "waveDirection", MFnNumericData::kFloat, 45.0f );
    nAttrWaveDirection.setKeyable(  true );	
	nAttrWaveDirection.setWritable(true);
	nAttrWaveDirection.setSoftMin(0.0f);	
	nAttrWaveDirection.setSoftMax(360.0f);	
	nAttrWaveDirection.setMin(0.0f);	
	nAttrWaveDirection.setMax(360.0f);	
    addAttribute( waveDirection );
    attributeAffects( aaOceanMaya::waveDirection, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveReflection;
	waveReflection = nAttrWaveReflection.create( "waveReflection", "waveReflection", MFnNumericData::kFloat, 0.0f );
    nAttrWaveReflection.setKeyable(  true );	
	nAttrWaveReflection.setWritable(true);
	nAttrWaveReflection.setSoftMin(0.0f);	
	nAttrWaveReflection.setSoftMax(1.0f);	
	nAttrWaveReflection.setMin(0.0f);	
	nAttrWaveReflection.setMax(1.0f);	
    addAttribute( waveReflection );
    attributeAffects( aaOceanMaya::waveReflection, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrWaveAlign;
	waveAlign = nAttrWaveAlign.create( "waveAlign", "waveAlign", MFnNumericData::kInt, 1 );
    nAttrWaveAlign.setKeyable(  true );	
	nAttrWaveAlign.setWritable(true);
	nAttrWaveAlign.setSoftMin(0);	
	nAttrWaveAlign.setSoftMax(10);	
	nAttrWaveAlign.setMin(0);	
	nAttrWaveAlign.setMax(10);	
    addAttribute( waveAlign );
    attributeAffects( aaOceanMaya::waveAlign, aaOceanMaya::outputGeom);

	MFnNumericAttribute nAttrCurrTime;
	currTime = nAttrCurrTime.create( "currTime", "currTime", MFnNumericData::kFloat, 0.042f );
    nAttrCurrTime.setKeyable(  true );
	nAttrWaveAlign.setWritable(true);
    addAttribute( currTime );
    attributeAffects( aaOceanMaya::currTime, aaOceanMaya::outputGeom);

	return MStatus::kSuccess;
}

aaOceanMaya::aaOceanMaya() 
{
	pOcean = new aaOcean;
	MGlobal::displayInfo( "[aaOcean Maya] Created a new ocean patch" );
}
aaOceanMaya::~aaOceanMaya() 
{
	if(pOcean)
	{
		delete pOcean;
		MGlobal::displayInfo( "[aaOcean Maya] Deleted ocean patch" );
	}
}
void* aaOceanMaya::creator()
{
	return new aaOceanMaya();
}


MStatus initializePlugin( MObject obj )
{
	MStatus result;
	MFnPlugin plugin( obj, "Amaan Akram", "2.0", "Any");
	result = plugin.registerNode( "aaOceanMaya", aaOceanMaya::id, aaOceanMaya::creator, 
								  aaOceanMaya::initialize, MPxNode::kDeformerNode );

	return result;
}

MStatus uninitializePlugin( MObject obj)
{
	MStatus result;
	MFnPlugin plugin( obj );
	result = plugin.deregisterNode( aaOceanMaya::id );
	return result;
}

void aaOceanMaya::fetchInput(MDataBlock& block)
{	
	float temp;
	int temp1;

	temp =  DegsToRads(block.inputValue(waveDirection).asFloat());
	if(pOcean->m_windDir != temp)
	{
		pOcean->m_windDir = temp;
		pOcean->m_redoHoK = true;
	}
	temp = (block.inputValue(waveSmooth).asFloat() * 0.01f);
	if(pOcean->m_cutoff != temp)
	{
		pOcean->m_cutoff = temp;
		pOcean->m_redoHoK = true;
	}
	temp = maximum<float>(block.inputValue(oceanSize).asFloat(),0.00001f);
	if(pOcean->m_oceanScale	!= temp)
	{
		pOcean->m_oceanScale = temp;
		pOcean->m_redoHoK = true;
	}
	temp = maximum<float>(((block.inputValue(waveSize).asFloat()  * block.inputValue(waveSize).asFloat()) / (aa_GRAVITY)),0.00001f);
	if(pOcean->m_velocity !=  temp)
	{
		pOcean->m_velocity = temp;
		pOcean->m_redoHoK = true;
	}
	temp1 = maximum<int>(((block.inputValue(waveAlign).asInt() + 1) * 2),2);
	if(pOcean->m_windAlign != temp1)
	{
		pOcean->m_windAlign = temp1;
		pOcean->m_redoHoK = true;
	}
	temp = block.inputValue(waveReflection).asFloat();
	if(pOcean->m_damp != temp)
	{
		pOcean->m_damp = temp;
		pOcean->m_redoHoK = true;
	}

	temp1 =  block.inputValue(seed).asInt();
	if(pOcean->m_seed	!= temp1)
	{
		pOcean->m_seed	= temp1;
		pOcean->m_redoHoK = true;
		pOcean->setup_grid(); 
	}
	pOcean->m_chopAmount	= block.inputValue(waveChop).asFloat()  * .01f;		//divided by scale for better ocean control;
	pOcean->m_waveHeight	= block.inputValue(waveHeight).asFloat()* .01f;		//divided by scale for better ocean control;
	pOcean->m_waveSpeed		= block.inputValue(waveSpeed).asFloat();
	pOcean->m_time			= block.inputValue(currTime).asFloat(); // time in seconds
}