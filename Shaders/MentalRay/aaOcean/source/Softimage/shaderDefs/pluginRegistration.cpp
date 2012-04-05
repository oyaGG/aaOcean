#include <xsi_application.h>
#include <xsi_context.h>
#include <xsi_pluginregistrar.h>
#include <xsi_status.h>
#include <xsi_factory.h>
#include <xsi_shaderdef.h>
#include <xsi_shaderparamdefoptions.h>
#include <xsi_shaderparamdefcontainer.h>
#include <xsi_metashaderrendererdef.h>
#include <xsi_valuemap.h>
using namespace XSI; 

XSI::CStatus RegisteraaOcean( XSI::PluginRegistrar& in_reg );

SICALLBACK XSILoadPlugin( PluginRegistrar& in_reg )
{
	in_reg.PutAuthor(L"Amaan Akram");
	in_reg.PutName(L"aaOceanMR");
	in_reg.PutEmail(L"amaan@amaanakram.com");
	in_reg.PutURL(L"http://www.amaanakram.com");
	in_reg.PutVersion(1,8);

	in_reg.RegisterShader(L"displacementShader", 1, 8);
	in_reg.RegisterShader(L"foamShader", 1, 8);
	in_reg.RegisterShader(L"normalsShader", 1, 8);
	//RegisteraaOcean( in_reg );

	return CStatus::OK;
}

SICALLBACK XSIUnloadPlugin( const PluginRegistrar& in_reg )
{
	CString strPluginName;
	strPluginName = in_reg.GetName();
	Application().LogMessage(strPluginName + L" has been unloaded.",siVerboseMsg);
	return CStatus::OK;
}