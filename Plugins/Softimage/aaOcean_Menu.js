// aaOcean v2.5 Softimage aaOcean Menu
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

function XSILoadPlugin( in_reg )
{
	in_reg.Author = "Amaan Akram";
	in_reg.Name = "aaOcean_Menu";
	in_reg.Email = "amaan@amaanakram.com";
	in_reg.URL = "http://www.amaanakram.com";
	in_reg.Major = 2;
	in_reg.Minor = 5;

	in_reg.RegisterProperty("aaOcean_Menu");
	in_reg.RegisterMenu(siMenuTbGetPropertyID,"aaOcean_Menu_Menu",false,false);

	return true;
}


function aaOcean_Menu_Define( in_ctxt )
{
	var path = XSIUtils.BuildPath(Application.InstallationPath(siProjectPath), "Pictures");
	
	var oCustomProperty;
	oCustomProperty = in_ctxt.Source;
	oCustomProperty.AddParameter2("Enable",siBool,false,null,null,null,null,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Enable");
	oCustomProperty.AddParameter2("ViewportResolution",siInt4,2,1,7,1,7,siClassifUnknown,siPersistable | siKeyable, "Resolution");
	oCustomProperty.AddParameter2("RenderResolution",siInt4,4,3,8,3,8,siClassifUnknown,siPersistable | siKeyable, "Render Res.");
	oCustomProperty.AddParameter2("ViewSubdivs",siInt4,32,4,1024,32,1024,siClassifUnknown,siPersistable | siKeyable | siReadOnly, "View SubDs");
	oCustomProperty.AddParameter2("RenderSubdivs",siInt4,256,4,4096,256,4096,siClassifUnknown,siPersistable | siKeyable | siReadOnly, "Render SubDs");
	oCustomProperty.AddParameter2("Ocean_Size",siFloat,100,1,100000,16,1000,siClassifUnknown,siPersistable | siKeyable , "Ocean Size");
	oCustomProperty.AddParameter2("Wind_Speed",siFloat,4.5,1,50,1,30,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Wave Size");
	oCustomProperty.AddParameter2("Wave_height",siFloat,5,0.0001,200,0.001,30,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Height");
	oCustomProperty.AddParameter2("Wave_Speed",siFloat,1,0.001,100,0.1,50,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Speed");
	oCustomProperty.AddParameter2("Wind_dir",siFloat,45,0,360,0,360,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Direction");
	oCustomProperty.AddParameter2("Wind_Align",siInt4,0,0,10,1,10,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Align");
	oCustomProperty.AddParameter2("Reflected_Waves",siFloat,0.985,0,1,0,1,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Reduce Reflected Waves");
	oCustomProperty.AddParameter2("Filter_small_waves",siFloat,0,0,2000,0,10,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Smooth");
	oCustomProperty.AddParameter2("Chop_Amount",siFloat,1,0,1000,0,5,siClassifUnknown,siPersistable | siKeyable | siAnimatable, "Chop Amount");
	oCustomProperty.AddParameter2("Seed",siInt4,1,1,20,1,20,siClassifUnknown,siPersistable | siKeyable , "Seed");
		
	return true;
}

function aaOcean_Menu_DefineLayout( in_ctxt )
{
	var oLayout;
	var oItem;
	oLayout = in_ctxt.Source;
	oLayout.Clear();
	oLayout.Clear();

	oLayout.AddTab("Ocean Controls");
	oLayout.AddSpacer(00 , 10);
	
		oLayout.AddItem("Enable");
		oLayout.AddGroup("Ocean Parameters");
			oLayout.AddRow();
				oLayout.AddItem("ViewportResolution");
				oLayout.AddItem("ViewSubdivs");
			oLayout.EndRow();
			oLayout.AddItem("Ocean_Size");
			oLayout.AddItem("Seed");
		oLayout.EndGroup();
		
		oLayout.AddGroup("Wave Parameters");
			oLayout.AddItem("Wave_height");
			oLayout.AddItem("Wind_Speed");
			oLayout.AddItem("Wave_speed");
			oLayout.AddItem("Chop_Amount"); 
			oLayout.AddItem("Filter_small_waves");
		oLayout.EndGroup();
		
		oLayout.AddGroup("Wind Parameters");
			oLayout.AddItem("Wind_dir");
			oLayout.AddItem("Reflected_Waves");
			oLayout.AddItem("Wind_Align");
		oLayout.EndGroup();

	return true;
}

function aaOcean_Menu_OnInit( )
{
	Application.LogMessage("aaOcean_Menu_OnInit called",siVerbose); //

	var Param = PPG.ViewSubdivs;
	var oParam = PPG.ViewportResolution;
	Param.Value = Math.pow(2,(oParam.Value+4));
	oParam.Parent3DObject.Parameters( "subdivu" ).Value = PPG.ViewSubdivs.Value;
	oParam.Parent3DObject.Parameters( "subdivv" ).Value = PPG.ViewSubdivs.Value;
	
	var rParam = PPG.RenderSubdivs;
	var rParam1 = PPG.RenderResolution;
	rParam.Value = Math.pow(2,(rParam1.Value+4));
	
	oParam.Parent3DObject.Parameters( "ulength" ).Value = PPG.Ocean_Size.Value;
	oParam.Parent3DObject.Parameters( "vlength" ).Value = PPG.Ocean_Size.Value;
	
}

function aaOcean_Menu_Ocean_Size_OnChanged( )
{
	var oParam = PPG.Ocean_Size;
	oParam.Parent3DObject.Parameters( "ulength" ).Value = PPG.Ocean_Size.Value;
	oParam.Parent3DObject.Parameters( "vlength" ).Value = PPG.Ocean_Size.Value;
}

function aaOcean_Menu_ViewportResolution_OnChanged( )
{
	var oParam = PPG.ViewportResolution;
	var paramVal = oParam.Value;
	var subdivParam = PPG.ViewSubdivs
	subdivParam.Value = Math.pow(2,(oParam.Value+4));
	
	oParam.Parent3DObject.Parameters( "subdivu" ).Value = PPG.ViewSubdivs.Value;
	oParam.Parent3DObject.Parameters( "subdivv" ).Value = PPG.ViewSubdivs.Value;
}

function aaOcean_Menu_RenderResolution_OnChanged( )
{
	var rParam = PPG.RenderSubdivs;
	var rParam1 = PPG.RenderResolution;
	rParam.Value = Math.pow(2,(rParam1.Value+4));
}

function aaOcean_Menu_Menu_Init( in_ctxt )
{
	var oMenu;
	oMenu = in_ctxt.Source;
	oMenu.AddCallbackItem("aaOcean_Menu","OnaaOcean_MenuMenuClicked");
	return true;
}

function OnaaOcean_MenuMenuClicked( in_ctxt )
{
	var oProp;
	var oParent;
	if (Selection.Count == 0)
	{
		XSIUIToolkit.MsgBox( "Please select a Primitive Grid object and try again" ) ; 
		return;
	}
	oParent = Selection.Item(0);
	oProp = oParent.AddProperty("aaOcean_Menu");
	InspectObj(oProp);
	return 1;
}


