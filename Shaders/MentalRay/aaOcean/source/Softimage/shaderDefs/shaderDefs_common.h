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
	resolution_opts.SetHardLimit(3, 10);
	resolution_opts.SetSoftLimit(4, 8);
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

	