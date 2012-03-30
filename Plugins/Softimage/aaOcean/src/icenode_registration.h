XSI::CStatus RegisteraaOcean( XSI::PluginRegistrar& in_reg );

SICALLBACK XSILoadPlugin( PluginRegistrar& in_reg )
{
	in_reg.PutAuthor(L"Amaan Akram");
	in_reg.PutName(L"aaOceanICENODE");
	in_reg.PutEmail(L"amaan@amaanakram.com");
	in_reg.PutURL(L"http://www.amaanakram.com");
	in_reg.PutVersion(2,0);

	//in_reg.RegisterProperty(L"aaOcean_Menu");
	//in_reg.RegisterMenu(siMenuTbGetPropertyID,L"aaOcean_Menu_Menu",false,false);

	RegisteraaOcean( in_reg );

	return CStatus::OK;
}

SICALLBACK XSIUnloadPlugin( const PluginRegistrar& in_reg )
{
	CString strPluginName;
	strPluginName = in_reg.GetName();
	Application().LogMessage(strPluginName + L" has been unloaded.");
	return CStatus::OK;
}