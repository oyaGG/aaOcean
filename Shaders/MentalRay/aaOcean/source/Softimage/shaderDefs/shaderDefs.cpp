SICALLBACK aaOceanMR_oceanDataShader_1_8_DefineInfo(const CRef& in_ctxt)
{
	Context ctxt(in_ctxt);

	ctxt.PutAttribute(L"Category", CValue(L"aaOceanMR"));
	ctxt.PutAttribute(L"DisplayName", CValue(L"aaOcean Data Shader"));

	return CStatus::OK;
}

SICALLBACK aaOceanMR_oceanDataShader_1_8_Define(const CRef& in_ctxt)
{
	Application app;
	Context ctxt(in_ctxt);

	ShaderDef sdef(ctxt.GetAttribute(L"Definition"));
	sdef.AddShaderFamily(siShaderFamilyTexture, true);

	// PARAMETERS
	Factory fact = app.GetFactory();
	ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
	ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

	// Output
	CRef poutopts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
	outpdefs.AddParamDef(L"out", siShaderDataTypeColor4, poutopts);

	#include "shaderDefs_common.h"

	CRef gamma_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions gamma_opts = ShaderParamDefOptions(gamma_popts);
	gamma_opts.SetTexturable(true);
	gamma_opts.SetInspectable(true);
	gamma_opts.SetHardLimit(0, 1000);
	gamma_opts.SetSoftLimit(0, 10);
	gamma_opts.SetDefaultValue(1);
	gamma_opts.SetLongName("Gamma");
	inpdefs.AddParamDef(L"gamma", siShaderDataTypeScalar, gamma_popts);

	CRef brightness_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions brightness_opts = ShaderParamDefOptions(brightness_popts);
	brightness_opts.SetTexturable(true);
	brightness_opts.SetInspectable(true);
	brightness_opts.SetHardLimit(0, 1000);
	brightness_opts.SetSoftLimit(0, 10);
	brightness_opts.SetDefaultValue(1);
	brightness_opts.SetLongName("Brightness");
	inpdefs.AddParamDef(L"brightness", siShaderDataTypeScalar, brightness_popts);

	CRef rawOutput_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions rawOutput_opts = ShaderParamDefOptions(rawOutput_popts);
	rawOutput_opts.SetTexturable(true);
	rawOutput_opts.SetInspectable(true);
	rawOutput_opts.SetDefaultValue(false);
	rawOutput_opts.SetLongName("Raw Output");
	inpdefs.AddParamDef(L"rawOutput", siShaderDataTypeBoolean, rawOutput_popts);

	CRef fmin_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fmin_opts = ShaderParamDefOptions(fmin_popts);
	fmin_opts.SetTexturable(true);
	fmin_opts.SetInspectable(true);
	fmin_opts.SetDefaultValue(0.0);
	fmin_opts.SetLongName("Min");
	inpdefs.AddParamDef(L"fmin", siShaderDataTypeScalar, fmin_popts);

	CRef fmax_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fmax_opts = ShaderParamDefOptions(fmax_popts);
	fmax_opts.SetTexturable(true);
	fmax_opts.SetInspectable(true);
	fmax_opts.SetDefaultValue(0.0f);
	fmax_opts.SetLongName("Max");
	inpdefs.AddParamDef(L"fmax", siShaderDataTypeScalar, fmax_popts);

	CRef layerOcean_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions layerOcean_opts = ShaderParamDefOptions(layerOcean_popts);
	layerOcean_opts.SetTexturable(true);
	layerOcean_opts.SetInspectable(false);
	layerOcean_opts.SetDefaultValue(0);
	layerOcean_opts.SetLongName("Layer Ocean");
	inpdefs.AddParamDef(L"layerOcean", siShaderDataTypeScalar, layerOcean_popts);

	PPGLayout oPPGLayout = sdef.GetPPGLayout();

	oPPGLayout.AddTab("Ocean Parameters");
	oPPGLayout.AddGroup("Ocean Parameters");
		oPPGLayout.AddItem("resolution");
		oPPGLayout.AddItem("oceanScale");
		oPPGLayout.AddItem("seed");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wave Parameters");
		oPPGLayout.AddItem("waveHeight");
		oPPGLayout.AddItem("velocity");
		oPPGLayout.AddItem("waveSpeed");
		oPPGLayout.AddItem("chopAmount");
		oPPGLayout.AddItem("cutoff");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wind Parameters");
		oPPGLayout.AddItem("windDir");
		oPPGLayout.AddItem("damp");
		oPPGLayout.AddItem("windAlign");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Misc");
		oPPGLayout.AddItem("fade");
		oPPGLayout.AddItem("time");
	oPPGLayout.EndGroup();

	oPPGLayout.AddTab("Foam Parameters");
	oPPGLayout.AddItem("rawOutput");
	oPPGLayout.AddGroup("Non raw options");
	oPPGLayout.AddItem("gamma");
	oPPGLayout.AddItem("brightness");
	oPPGLayout.AddItem("fmin");
	oPPGLayout.AddItem("fmax");

	// RENDERERS
	MetaShaderRendererDef rend = sdef.AddRendererDef(L"mental ray");
	rend.PutSymbolName(L"aaOceanDataShader");
	rend.PutCodePath(L"aaOceanMentalRay");
	ValueMap ropts = rend.GetRendererOptions();
	ropts.Set(L"version", CValue(1));
	return CStatus::OK;

	return CStatus::OK;
}

SICALLBACK aaOceanMR_displacementShader_1_8_DefineInfo( const CRef& in_ctxt )
{
	Context ctxt(in_ctxt);

	ctxt.PutAttribute(L"Category", CValue(L"aaOceanMR"));
	ctxt.PutAttribute(L"DisplayName", CValue(L"aaOcean Displacement Shader"));

	return CStatus::OK;
}

SICALLBACK aaOceanMR_displacementShader_1_8_Define( const CRef& in_ctxt )
{
	Application app;
	Context ctxt(in_ctxt);

	ShaderDef sdef(ctxt.GetAttribute(L"Definition"));
	sdef.AddShaderFamily(siShaderFamilyTexture);

	// PARAMETERS
	Factory fact = app.GetFactory();
	ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
	ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

	// Output
	CRef poutopts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
	outpdefs.AddParamDef(L"out", siShaderDataTypeScalar, poutopts);

	#include "shaderDefs_common.h"

	CRef layerOcean_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions layerOcean_opts = ShaderParamDefOptions(layerOcean_popts);
	layerOcean_opts.SetTexturable(true);
	layerOcean_opts.SetInspectable(false);
	layerOcean_opts.SetDefaultValue(0);
	layerOcean_opts.SetLongName("Layer Ocean");
	inpdefs.AddParamDef(L"layerOcean", siShaderDataTypeScalar, layerOcean_popts);

	//Layout Defs
	PPGLayout oPPGLayout = sdef.GetPPGLayout();
	oPPGLayout.AddGroup("Ocean Parameters");
		oPPGLayout.AddItem("resolution");
		oPPGLayout.AddItem("oceanScale");
		oPPGLayout.AddItem("seed");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wave Parameters");
		oPPGLayout.AddItem("waveHeight");
		oPPGLayout.AddItem("velocity");
		oPPGLayout.AddItem("waveSpeed");
		oPPGLayout.AddItem("chopAmount");
		oPPGLayout.AddItem("cutoff");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wind Parameters");
		oPPGLayout.AddItem("windDir");
		oPPGLayout.AddItem("damp");
		oPPGLayout.AddItem("windAlign");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Misc");
		oPPGLayout.AddItem("fade");
		oPPGLayout.AddItem("time");
	oPPGLayout.EndGroup();
	
	// RENDERERS
	MetaShaderRendererDef rend = sdef.AddRendererDef(L"mental ray");
	rend.PutSymbolName(L"aaOceanDisplaceShader");
	rend.PutCodePath(L"aaOceanMentalRay");
	ValueMap ropts = rend.GetRendererOptions();
	ropts.Set(L"version", CValue(1));

	return CStatus::OK;
}

SICALLBACK aaOceanMR_foamShader_1_8_DefineInfo( const CRef& in_ctxt )
{
	Context ctxt(in_ctxt);
	ctxt.PutAttribute(L"Category", CValue(L"aaOceanMR"));
	ctxt.PutAttribute(L"DisplayName", CValue(L"aaOcean Foam Shader"));
	return CStatus::OK;
}

SICALLBACK aaOceanMR_foamShader_1_8_Define( const CRef& in_ctxt )
{
	Application app;
	Context ctxt(in_ctxt);

	ShaderDef sdef(ctxt.GetAttribute(L"Definition"));
	sdef.AddShaderFamily(siShaderFamilyTexture, true);

	// PARAMETERS
	Factory fact = app.GetFactory();
	ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
	ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

	// Output
	CRef poutopts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
	outpdefs.AddParamDef(L"out", siShaderDataTypeColor4, poutopts);

	// Input
	#include "shaderDefs_common.h"

	CRef gamma_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions gamma_opts = ShaderParamDefOptions(gamma_popts);
	gamma_opts.SetTexturable(true);
	gamma_opts.SetInspectable(true);
	gamma_opts.SetHardLimit(0, 1000);
	gamma_opts.SetSoftLimit(0, 10);
	gamma_opts.SetDefaultValue(1);
	gamma_opts.SetLongName("Gamma");
	inpdefs.AddParamDef(L"gamma", siShaderDataTypeScalar, gamma_popts);

	CRef brightness_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions brightness_opts = ShaderParamDefOptions(brightness_popts);
	brightness_opts.SetTexturable(true);
	brightness_opts.SetInspectable(true);
	brightness_opts.SetHardLimit(0, 1000);
	brightness_opts.SetSoftLimit(0, 10);
	brightness_opts.SetDefaultValue(1);
	brightness_opts.SetLongName("Brightness");
	inpdefs.AddParamDef(L"brightness", siShaderDataTypeScalar, brightness_popts);

	CRef rawOutput_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions rawOutput_opts = ShaderParamDefOptions(rawOutput_popts);
	rawOutput_opts.SetTexturable(true);
	rawOutput_opts.SetInspectable(true);
	rawOutput_opts.SetDefaultValue(false);
	rawOutput_opts.SetLongName("Raw Output");
	inpdefs.AddParamDef(L"rawOutput", siShaderDataTypeBoolean, rawOutput_popts);

	CRef fmin_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fmin_opts = ShaderParamDefOptions(fmin_popts);
	fmin_opts.SetTexturable(true);
	fmin_opts.SetInspectable(true);
	fmin_opts.SetDefaultValue(0.0);
	fmin_opts.SetLongName("Min");
	inpdefs.AddParamDef(L"fmin", siShaderDataTypeScalar, fmin_popts);

	CRef fmax_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fmax_opts = ShaderParamDefOptions(fmax_popts);
	fmax_opts.SetTexturable(true);
	fmax_opts.SetInspectable(true);
	fmax_opts.SetDefaultValue(0.0);
	fmax_opts.SetLongName("Max");
	inpdefs.AddParamDef(L"fmax", siShaderDataTypeScalar, fmax_popts);

	PPGLayout oPPGLayout = sdef.GetPPGLayout();
	oPPGLayout.AddGroup("Ocean Parameters");
		oPPGLayout.AddItem("resolution");
		oPPGLayout.AddItem("oceanScale");
		oPPGLayout.AddItem("seed");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wave Parameters");
		oPPGLayout.AddItem("waveHeight");
		oPPGLayout.AddItem("velocity");
		oPPGLayout.AddItem("waveSpeed");
		oPPGLayout.AddItem("chopAmount");
		oPPGLayout.AddItem("cutoff");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wind Parameters");
		oPPGLayout.AddItem("windDir");
		oPPGLayout.AddItem("damp");
		oPPGLayout.AddItem("windAlign");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Post Processing");
	oPPGLayout.AddItem("gamma");
	oPPGLayout.AddItem("brightness");
	oPPGLayout.AddItem("fade");
	oPPGLayout.AddItem("fmin");
	oPPGLayout.AddItem("fmax");
	oPPGLayout.AddItem("rawOutput");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Misc");
		oPPGLayout.AddItem("time");
	oPPGLayout.EndGroup();

	// RENDERERS
	MetaShaderRendererDef rend = sdef.AddRendererDef(L"mental ray");
	rend.PutSymbolName(L"aaOceanFoamShader");
	rend.PutCodePath(L"aaOceanMentalRay");
	ValueMap ropts = rend.GetRendererOptions();
	ropts.Set(L"version", CValue(1));
	return CStatus::OK;
}

SICALLBACK aaOceanMR_normalsShader_1_8_DefineInfo( const CRef& in_ctxt )
{
	Context ctxt(in_ctxt);
	ctxt.PutAttribute(L"Category", CValue(L"aaOceanMR"));
	ctxt.PutAttribute(L"DisplayName", CValue(L"aaOcean Normals Shader"));
	return CStatus::OK;
}

SICALLBACK aaOceanMR_normalsShader_1_8_Define( const CRef& in_ctxt )
{
	Application app;
	Context ctxt(in_ctxt);

	ShaderDef sdef(ctxt.GetAttribute(L"Definition"));
	sdef.AddShaderFamily(siShaderFamilyTexture, true);

	// PARAMETERS
	Factory fact = app.GetFactory();
	ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
	ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

	// Output
	CRef poutopts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
	outpdefs.AddParamDef(L"out", siShaderDataTypeVector3, poutopts);

	// Input
	#include "shaderDefs_common.h"

	CRef layerNormals_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions layerNormals_opts = ShaderParamDefOptions(layerNormals_popts);
	layerNormals_opts.SetTexturable(true);
	layerNormals_opts.SetInspectable(false);
	layerNormals_opts.SetDefaultValue(0);
	layerNormals_opts.SetLongName("Layer Normals");
	inpdefs.AddParamDef(L"layerNormals", siShaderDataTypeVector3, layerNormals_popts);

	PPGLayout oPPGLayout = sdef.GetPPGLayout();
	oPPGLayout.AddGroup("Ocean Parameters");
		oPPGLayout.AddItem("resolution");
		oPPGLayout.AddItem("oceanScale");
		oPPGLayout.AddItem("seed");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wave Parameters");
		oPPGLayout.AddItem("waveHeight");
		oPPGLayout.AddItem("velocity");
		oPPGLayout.AddItem("waveSpeed");
		oPPGLayout.AddItem("chopAmount");
		oPPGLayout.AddItem("cutoff");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Wind Parameters");
		oPPGLayout.AddItem("windDir");
		oPPGLayout.AddItem("damp");
		oPPGLayout.AddItem("windAlign");
	oPPGLayout.EndGroup();
	oPPGLayout.AddGroup("Misc");
		oPPGLayout.AddItem("time");
	oPPGLayout.EndGroup();

	// RENDERERS
	MetaShaderRendererDef rend = sdef.AddRendererDef(L"mental ray");
	rend.PutSymbolName(L"aaOceanNormalsShader");
	rend.PutCodePath(L"aaOceanMentalRay");
	ValueMap ropts = rend.GetRendererOptions();
	ropts.Set(L"version", CValue(1));
	return CStatus::OK;
}