// aaOcean Mental Ray Shaders, Softimage Shader Definitions
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

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
	in_reg.PutName(L"aaOceanMentalRay");
	in_reg.PutEmail(L"amaan@amaanakram.com");
	in_reg.PutURL(L"http://www.amaanakram.com");
	in_reg.PutVersion(2,5);

	//in_reg.RegisterShader(L"displacementShader", 1, 8);
	//in_reg.RegisterShader(L"foamShader", 1, 8);
	//in_reg.RegisterShader(L"normalsShader", 1, 8);
	in_reg.RegisterShader(L"oceanDataShader", 2, 5);
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