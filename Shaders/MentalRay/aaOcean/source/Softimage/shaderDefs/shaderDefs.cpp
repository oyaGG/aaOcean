// Includes specific to this Shader Definition

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
	sdef.AddShaderFamily(siDisplacementShaderFamily,true);

	// PARAMETERS
	Factory fact = app.GetFactory();
	ShaderParamDefContainer inpdefs = sdef.GetInputParamDefs();
	ShaderParamDefContainer outpdefs = sdef.GetOutputParamDefs();

	// Output
	CRef poutopts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions outopts = ShaderParamDefOptions(poutopts);
	outpdefs.AddParamDef(L"out", siShaderDataTypeScalar, poutopts);

	// Input
	CRef uv_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions uv_opts = ShaderParamDefOptions(uv_popts);
	uv_opts.SetTexturable(true);
	uv_opts.SetInspectable(false);
	uv_opts.SetDefaultValue(0);
	uv_opts.SetLongName("UV Coords");
	inpdefs.AddParamDef(L"uv_coords", siShaderDataTypeVector3, uv_popts);

	CRef fade_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fade_opts = ShaderParamDefOptions(fade_popts);
	fade_opts.SetTexturable(true);
	fade_opts.SetInspectable(true);
	fade_opts.SetHardLimit(0, 1);
	fade_opts.SetDefaultValue(0);
	fade_opts.SetLongName("Fade");
	inpdefs.AddParamDef(L"fade", siShaderDataTypeScalar, fade_popts);

	CRef resolution_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions resolution_opts = ShaderParamDefOptions(resolution_popts);
	resolution_opts.SetTexturable(true);
	resolution_opts.SetInspectable(true);
	resolution_opts.SetHardLimit(3, 8);
	resolution_opts.SetDefaultValue(4);
	resolution_opts.SetLongName("Resolution");
	inpdefs.AddParamDef(L"resolution", siShaderDataTypeInteger, resolution_popts);

	CRef oceanScale_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions oceanScale_opts = ShaderParamDefOptions(oceanScale_popts);
	oceanScale_opts.SetTexturable(true);
	oceanScale_opts.SetInspectable(true);
	oceanScale_opts.SetDefaultValue(100);
	oceanScale_opts.SetHardLimit(0.0, 200000);
	oceanScale_opts.SetSoftLimit(0.0001,15.0);
	oceanScale_opts.SetLongName("Ocean Size");
	inpdefs.AddParamDef(L"oceanScale", siShaderDataTypeScalar, oceanScale_popts);

	CRef seed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions seed_opts = ShaderParamDefOptions(seed_popts);
	seed_opts.SetTexturable(true);
	seed_opts.SetInspectable(true);
	seed_opts.SetHardLimit(1, 20);
	seed_opts.SetDefaultValue(1);
	seed_opts.SetLongName("Seed");
	inpdefs.AddParamDef(L"seed", siShaderDataTypeInteger, seed_popts);

	CRef waveHeight_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveHeight_opts = ShaderParamDefOptions(waveHeight_popts);
	waveHeight_opts.SetTexturable(true);
	waveHeight_opts.SetInspectable(true);
	waveHeight_opts.SetHardLimit(0.0, 10000);
	waveHeight_opts.SetSoftLimit(0.001,15.0);
	waveHeight_opts.SetDefaultValue(1);
	waveHeight_opts.SetLongName("Wave Height");
	inpdefs.AddParamDef(L"waveHeight", siShaderDataTypeScalar, waveHeight_popts);

	CRef velocity_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions velocity_opts = ShaderParamDefOptions(velocity_popts);
	velocity_opts.SetTexturable(true);
	velocity_opts.SetInspectable(true);
	velocity_opts.SetHardLimit(0.0, 30);
	velocity_opts.SetDefaultValue(5);
	velocity_opts.SetLongName("Wave Size");
	inpdefs.AddParamDef(L"velocity", siShaderDataTypeScalar, velocity_popts);

	CRef waveSpeed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveSpeed_opts = ShaderParamDefOptions(waveSpeed_popts);
	waveSpeed_opts.SetTexturable(true);
	waveSpeed_opts.SetInspectable(true);
	waveSpeed_opts.SetHardLimit(0.00001, 1000);
	waveSpeed_opts.SetSoftLimit(0.001,15.0);
	waveSpeed_opts.SetDefaultValue(1);
	waveSpeed_opts.SetLongName("Wave Speed");
	inpdefs.AddParamDef(L"waveSpeed", siShaderDataTypeScalar, waveSpeed_popts);

	CRef chopAmount_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions chopAmount_opts = ShaderParamDefOptions(chopAmount_popts);
	chopAmount_opts.SetTexturable(true);
	chopAmount_opts.SetInspectable(true);
	chopAmount_opts.SetHardLimit(0.0, 100000);
	chopAmount_opts.SetDefaultValue(1);
	chopAmount_opts.SetSoftLimit(0.001,15.0);
	chopAmount_opts.SetLongName("Chop Amount");
	inpdefs.AddParamDef(L"chopAmount", siShaderDataTypeScalar, chopAmount_popts);

	CRef cutoff_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions cutoff_opts = ShaderParamDefOptions(cutoff_popts);
	cutoff_opts.SetTexturable(true);
	cutoff_opts.SetInspectable(true);
	cutoff_opts.SetHardLimit(0.0, 1000);
	cutoff_opts.SetSoftLimit(0.0,15.0);
	cutoff_opts.SetDefaultValue(0);
	cutoff_opts.SetLongName("Smooth Amount");
	inpdefs.AddParamDef(L"cutoff", siShaderDataTypeScalar, cutoff_popts);

	CRef windDir_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windDir_opts = ShaderParamDefOptions(windDir_popts);
	windDir_opts.SetTexturable(true);
	windDir_opts.SetInspectable(true);
	windDir_opts.SetHardLimit(0.0, 360.0);
	windDir_opts.SetDefaultValue(45);
	windDir_opts.SetLongName("Wind Direction");
	inpdefs.AddParamDef(L"windDir", siShaderDataTypeScalar, windDir_popts);

	CRef damp_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions damp_opts = ShaderParamDefOptions(damp_popts);
	damp_opts.SetTexturable(true);
	damp_opts.SetInspectable(true);
	damp_opts.SetHardLimit(0.0, 1.0);
	damp_opts.SetDefaultValue(0.985);
	damp_opts.SetLongName("Reflected Waves");
	inpdefs.AddParamDef(L"damp", siShaderDataTypeScalar, damp_popts);

	CRef windAlign_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windAlign_opts = ShaderParamDefOptions(windAlign_popts);
	windAlign_opts.SetTexturable(true);
	windAlign_opts.SetInspectable(true);
	windAlign_opts.SetHardLimit(0, 200);
	windAlign_opts.SetSoftLimit(0,10);
	windAlign_opts.SetDefaultValue(1);
	windAlign_opts.SetLongName("Wind Align");
	inpdefs.AddParamDef(L"windAlign", siShaderDataTypeInteger, windAlign_popts);

	CRef time_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions time_opts = ShaderParamDefOptions(time_popts);
	time_opts.SetTexturable(true);
	time_opts.SetInspectable(true);
	time_opts.SetDefaultValue(0.042);
	time_opts.SetLongName("Current Time");
	inpdefs.AddParamDef(L"time", siShaderDataTypeScalar, time_popts);

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
	CRef uv_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions uv_opts = ShaderParamDefOptions(uv_popts);
	uv_opts.SetTexturable(true);
	uv_opts.SetInspectable(true);
	uv_opts.SetDefaultValue(0);
	uv_opts.SetLongName("UV Coords");
	inpdefs.AddParamDef(L"uv_coords", siShaderDataTypeVector3, uv_popts);

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

	CRef fade_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fade_opts = ShaderParamDefOptions(fade_popts);
	fade_opts.SetTexturable(true);
	fade_opts.SetInspectable(true);
	fade_opts.SetHardLimit(0,1);
	fade_opts.SetDefaultValue(0);
	fade_opts.SetLongName("Fade");
	inpdefs.AddParamDef(L"fade", siShaderDataTypeScalar, fade_popts);

	CRef rawOutput_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions rawOutput_opts = ShaderParamDefOptions(rawOutput_popts);
	rawOutput_opts.SetTexturable(true);
	rawOutput_opts.SetInspectable(true);
	rawOutput_opts.SetDefaultValue(false);
	rawOutput_opts.SetLongName("Raw Output");
	inpdefs.AddParamDef(L"rawOutput", siShaderDataTypeBoolean, rawOutput_popts);

	CRef resolution_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions resolution_opts = ShaderParamDefOptions(resolution_popts);
	resolution_opts.SetTexturable(true);
	resolution_opts.SetInspectable(true);
	resolution_opts.SetHardLimit(3, 8);
	resolution_opts.SetDefaultValue(4);
	resolution_opts.SetLongName("Resolution");
	inpdefs.AddParamDef(L"resolution", siShaderDataTypeInteger, resolution_popts);

	CRef oceanScale_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions oceanScale_opts = ShaderParamDefOptions(oceanScale_popts);
	oceanScale_opts.SetTexturable(true);
	oceanScale_opts.SetInspectable(true);
	oceanScale_opts.SetDefaultValue(100);
	oceanScale_opts.SetHardLimit(0.0, 200000);
	oceanScale_opts.SetSoftLimit(0.0001,15.0);
	oceanScale_opts.SetLongName("Ocean Size");
	inpdefs.AddParamDef(L"oceanScale", siShaderDataTypeScalar, oceanScale_popts);

	CRef seed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions seed_opts = ShaderParamDefOptions(seed_popts);
	seed_opts.SetTexturable(true);
	seed_opts.SetInspectable(true);
	seed_opts.SetHardLimit(1, 20);
	seed_opts.SetDefaultValue(1);
	seed_opts.SetLongName("Seed");
	inpdefs.AddParamDef(L"seed", siShaderDataTypeInteger, seed_popts);

	CRef waveHeight_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveHeight_opts = ShaderParamDefOptions(waveHeight_popts);
	waveHeight_opts.SetTexturable(true);
	waveHeight_opts.SetInspectable(true);
	waveHeight_opts.SetHardLimit(0.0, 10000);
	waveHeight_opts.SetSoftLimit(0.001,15.0);
	waveHeight_opts.SetDefaultValue(1);
	waveHeight_opts.SetLongName("Wave Height");
	inpdefs.AddParamDef(L"waveHeight", siShaderDataTypeScalar, waveHeight_popts);

	CRef velocity_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions velocity_opts = ShaderParamDefOptions(velocity_popts);
	velocity_opts.SetTexturable(true);
	velocity_opts.SetInspectable(true);
	velocity_opts.SetHardLimit(0.0, 30);
	velocity_opts.SetDefaultValue(5);
	velocity_opts.SetLongName("Wave Size");
	inpdefs.AddParamDef(L"velocity", siShaderDataTypeScalar, velocity_popts);

	CRef waveSpeed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveSpeed_opts = ShaderParamDefOptions(waveSpeed_popts);
	waveSpeed_opts.SetTexturable(true);
	waveSpeed_opts.SetInspectable(true);
	waveSpeed_opts.SetHardLimit(0.00001, 1000);
	waveSpeed_opts.SetSoftLimit(0.001,15.0);
	waveSpeed_opts.SetDefaultValue(1);
	waveSpeed_opts.SetLongName("Wave Speed");
	inpdefs.AddParamDef(L"waveSpeed", siShaderDataTypeScalar, waveSpeed_popts);

	CRef chopAmount_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions chopAmount_opts = ShaderParamDefOptions(chopAmount_popts);
	chopAmount_opts.SetTexturable(true);
	chopAmount_opts.SetInspectable(true);
	chopAmount_opts.SetHardLimit(0.0, 100000);
	chopAmount_opts.SetDefaultValue(1);
	chopAmount_opts.SetSoftLimit(0.001,15.0);
	chopAmount_opts.SetLongName("Chop Amount");
	inpdefs.AddParamDef(L"chopAmount", siShaderDataTypeScalar, chopAmount_popts);

	CRef cutoff_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions cutoff_opts = ShaderParamDefOptions(cutoff_popts);
	cutoff_opts.SetTexturable(true);
	cutoff_opts.SetInspectable(true);
	cutoff_opts.SetHardLimit(0.0, 1000);
	cutoff_opts.SetSoftLimit(0.0,15.0);
	cutoff_opts.SetDefaultValue(0);
	cutoff_opts.SetLongName("Smooth Amount");
	inpdefs.AddParamDef(L"cutoff", siShaderDataTypeScalar, cutoff_popts);

	CRef windDir_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windDir_opts = ShaderParamDefOptions(windDir_popts);
	windDir_opts.SetTexturable(true);
	windDir_opts.SetInspectable(true);
	windDir_opts.SetHardLimit(0.0, 360.0);
	windDir_opts.SetDefaultValue(45);
	windDir_opts.SetLongName("Wind Direction");
	inpdefs.AddParamDef(L"windDir", siShaderDataTypeScalar, windDir_popts);

	CRef damp_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions damp_opts = ShaderParamDefOptions(damp_popts);
	damp_opts.SetTexturable(true);
	damp_opts.SetInspectable(true);
	damp_opts.SetHardLimit(0.0, 1.0);
	damp_opts.SetDefaultValue(0.985);
	damp_opts.SetLongName("Reflected Waves");
	inpdefs.AddParamDef(L"damp", siShaderDataTypeScalar, damp_popts);

	CRef windAlign_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windAlign_opts = ShaderParamDefOptions(windAlign_popts);
	windAlign_opts.SetTexturable(true);
	windAlign_opts.SetInspectable(true);
	windAlign_opts.SetHardLimit(0, 200);
	windAlign_opts.SetSoftLimit(0,10);
	windAlign_opts.SetDefaultValue(1);
	windAlign_opts.SetLongName("Wind Align");
	inpdefs.AddParamDef(L"windAlign", siShaderDataTypeInteger, windAlign_popts);

	CRef time_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions time_opts = ShaderParamDefOptions(time_popts);
	time_opts.SetTexturable(true);
	time_opts.SetInspectable(true);
	time_opts.SetDefaultValue(0.042);
	time_opts.SetLongName("Current Time");
	inpdefs.AddParamDef(L"time", siShaderDataTypeScalar, time_popts);

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
	CRef uv_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions uv_opts = ShaderParamDefOptions(uv_popts);
	uv_opts.SetTexturable(true);
	uv_opts.SetInspectable(true);
	uv_opts.SetDefaultValue(0);
	uv_opts.SetLongName("UV Coords");
	inpdefs.AddParamDef(L"uv_coords", siShaderDataTypeVector3, uv_popts);

	CRef fade_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions fade_opts = ShaderParamDefOptions(fade_popts);
	fade_opts.SetTexturable(true);
	fade_opts.SetInspectable(true);
	fade_opts.SetHardLimit(0, 1);
	fade_opts.SetDefaultValue(0);
	fade_opts.SetLongName("Fade");
	inpdefs.AddParamDef(L"fade", siShaderDataTypeScalar, fade_popts);

	CRef resolution_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions resolution_opts = ShaderParamDefOptions(resolution_popts);
	resolution_opts.SetTexturable(true);
	resolution_opts.SetInspectable(true);
	resolution_opts.SetHardLimit(3, 8);
	resolution_opts.SetDefaultValue(4);
	resolution_opts.SetLongName("Resolution");
	inpdefs.AddParamDef(L"resolution", siShaderDataTypeInteger, resolution_popts);

	CRef oceanScale_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions oceanScale_opts = ShaderParamDefOptions(oceanScale_popts);
	oceanScale_opts.SetTexturable(true);
	oceanScale_opts.SetInspectable(true);
	oceanScale_opts.SetDefaultValue(100);
	oceanScale_opts.SetHardLimit(0.0, 200000);
	oceanScale_opts.SetSoftLimit(0.0001,15.0);
	oceanScale_opts.SetLongName("Ocean Size");
	inpdefs.AddParamDef(L"oceanScale", siShaderDataTypeScalar, oceanScale_popts);

	CRef seed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions seed_opts = ShaderParamDefOptions(seed_popts);
	seed_opts.SetTexturable(true);
	seed_opts.SetInspectable(true);
	seed_opts.SetHardLimit(1, 20);
	seed_opts.SetDefaultValue(1);
	seed_opts.SetLongName("Seed");
	inpdefs.AddParamDef(L"seed", siShaderDataTypeInteger, seed_popts);

	CRef waveHeight_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveHeight_opts = ShaderParamDefOptions(waveHeight_popts);
	waveHeight_opts.SetTexturable(true);
	waveHeight_opts.SetInspectable(true);
	waveHeight_opts.SetHardLimit(0.0, 10000);
	waveHeight_opts.SetSoftLimit(0.001,15.0);
	waveHeight_opts.SetDefaultValue(1);
	waveHeight_opts.SetLongName("Wave Height");
	inpdefs.AddParamDef(L"waveHeight", siShaderDataTypeScalar, waveHeight_popts);

	CRef velocity_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions velocity_opts = ShaderParamDefOptions(velocity_popts);
	velocity_opts.SetTexturable(true);
	velocity_opts.SetInspectable(true);
	velocity_opts.SetHardLimit(0.0, 30);
	velocity_opts.SetDefaultValue(5);
	velocity_opts.SetLongName("Wave Size");
	inpdefs.AddParamDef(L"velocity", siShaderDataTypeScalar, velocity_popts);

	CRef waveSpeed_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions waveSpeed_opts = ShaderParamDefOptions(waveSpeed_popts);
	waveSpeed_opts.SetTexturable(true);
	waveSpeed_opts.SetInspectable(true);
	waveSpeed_opts.SetHardLimit(0.00001, 1000);
	waveSpeed_opts.SetSoftLimit(0.001,15.0);
	waveSpeed_opts.SetDefaultValue(1);
	waveSpeed_opts.SetLongName("Wave Speed");
	inpdefs.AddParamDef(L"waveSpeed", siShaderDataTypeScalar, waveSpeed_popts);

	CRef chopAmount_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions chopAmount_opts = ShaderParamDefOptions(chopAmount_popts);
	chopAmount_opts.SetTexturable(true);
	chopAmount_opts.SetInspectable(true);
	chopAmount_opts.SetHardLimit(0.0, 100000);
	chopAmount_opts.SetDefaultValue(1);
	chopAmount_opts.SetSoftLimit(0.001,15.0);
	chopAmount_opts.SetLongName("Chop Amount");
	inpdefs.AddParamDef(L"chopAmount", siShaderDataTypeScalar, chopAmount_popts);

	CRef cutoff_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions cutoff_opts = ShaderParamDefOptions(cutoff_popts);
	cutoff_opts.SetTexturable(true);
	cutoff_opts.SetInspectable(true);
	cutoff_opts.SetHardLimit(0.0, 1000);
	cutoff_opts.SetSoftLimit(0.0,15.0);
	cutoff_opts.SetDefaultValue(0);
	cutoff_opts.SetLongName("Smooth Amount");
	inpdefs.AddParamDef(L"cutoff", siShaderDataTypeScalar, cutoff_popts);

	CRef windDir_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windDir_opts = ShaderParamDefOptions(windDir_popts);
	windDir_opts.SetTexturable(true);
	windDir_opts.SetInspectable(true);
	windDir_opts.SetHardLimit(0.0, 360.0);
	windDir_opts.SetDefaultValue(45);
	windDir_opts.SetLongName("Wind Direction");
	inpdefs.AddParamDef(L"windDir", siShaderDataTypeScalar, windDir_popts);

	CRef damp_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions damp_opts = ShaderParamDefOptions(damp_popts);
	damp_opts.SetTexturable(true);
	damp_opts.SetInspectable(true);
	damp_opts.SetHardLimit(0.0, 1.0);
	damp_opts.SetDefaultValue(0.985);
	damp_opts.SetLongName("Reflected Waves");
	inpdefs.AddParamDef(L"damp", siShaderDataTypeScalar, damp_popts);

	CRef windAlign_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions windAlign_opts = ShaderParamDefOptions(windAlign_popts);
	windAlign_opts.SetTexturable(true);
	windAlign_opts.SetInspectable(true);
	windAlign_opts.SetHardLimit(0, 200);
	windAlign_opts.SetSoftLimit(0,10);
	windAlign_opts.SetDefaultValue(1);
	windAlign_opts.SetLongName("Wind Align");
	inpdefs.AddParamDef(L"windAlign", siShaderDataTypeInteger, windAlign_popts);

	CRef time_popts = fact.CreateShaderParamDefOptions();
	ShaderParamDefOptions time_opts = ShaderParamDefOptions(time_popts);
	time_opts.SetTexturable(true);
	time_opts.SetInspectable(true);
	time_opts.SetDefaultValue(0.042);
	time_opts.SetLongName("Current Time");
	inpdefs.AddParamDef(L"time", siShaderDataTypeScalar, time_popts);

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