
// Defines port, group and map identifiers used for registering the ICENode
enum IDs
{
	ID_IN_PointID = 0,
	ID_IN_OCEAN_SCALE,
	ID_IN_WINDDIR,
	ID_IN_CUTOFF,
	ID_IN_WINDVELOCITY,
	ID_IN_WINDALIGN,
	ID_IN_DAMP,
	ID_IN_WAVESPEED = 7,
	ID_IN_CHOP = 9,
	ID_IN_GRID_SCALE = 600,
	ID_IN_WAVE_HEIGHT = 25,
	ID_IN_SEED = 29,
	ID_IN_RENDER_READY = 30,
	ID_IN_GRID_LENGTH_U = 31,
	ID_IN_GRID_LENGTH_V = 32,

	ID_IN_RANDOM_TYPE = 40,
	
	ID_G_100 = 100,
	ID_G_101,
	ID_G_102,

	ID_G_200 = 200,
	ID_G_201 = 201,
	ID_G_202 = 202,
	ID_G_203 = 203,
	ID_G_204 = 204,

	ID_G_300 = 300,
	ID_OUT_OCEAN,
	ID_OUT_FOAM,
	ID_OUT_EIGEN_MINUS,
	ID_OUT_EIGEN_PLUS,
	ID_OUT_NORMALS,

	ID_UNDEF = ULONG_MAX
};

CStatus RegisteraaOcean( PluginRegistrar& in_reg )
{
	ICENodeDef nodeDef;
	nodeDef = Application().GetFactory().CreateICENodeDef(L"aaOcean");

	CStatus st;
	st = nodeDef.PutColor(133,164,192);
	st.AssertSucceeded();

	st = nodeDef.PutThreadingModel(XSI::siICENodeSingleThreading);
	st.AssertSucceeded();

	// Add input ports and groups.
	st = nodeDef.AddPortGroup(ID_G_100,1,1,L"Ocean Params");
	st.AssertSucceeded();
	st = nodeDef.AddPortGroup(ID_G_101,1,1,L"Export Params");
	st.AssertSucceeded();
	st = nodeDef.AddPortGroup(ID_G_102,1,1,L"Misc"); 
	st.AssertSucceeded();

	st = nodeDef.AddInputPort(	ID_IN_OCEAN_SCALE,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Ocean Size",L"Ocean_Size",
								100.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_WAVE_HEIGHT,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Wave Height",L"Wave_Height",
								1.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_WINDVELOCITY,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Wave Size",L"Wave_Size",
								5.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_WAVESPEED,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Wave Speed",L"Wave_Speed",
								1.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_CHOP,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Chopiness",L"Chopiness",
								1.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_CUTOFF,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Smooth",L"Smooth",
								0.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_WINDDIR,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Wind Direction",L"Wind_Direction",
								45.0f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_DAMP,
								ID_G_100,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Reflected Waves",L"Reflected_Waves",
								0.985f);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_WINDALIGN,
								ID_G_100,
								siICENodeDataLong,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Wind Align",L"Wind_Align",
								0);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_SEED,
								ID_G_100,
								siICENodeDataLong,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Seed",L"Seed",
								1);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_GRID_LENGTH_U,
								ID_G_101,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Grid Length U",L"Grid_Length_U",
								10);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_GRID_LENGTH_V,
								ID_G_101,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Grid Length V",L"Grid_Length_V",
								10.f);
	st.AssertSucceeded( ) ;
	

	st = nodeDef.AddInputPort(	ID_IN_RANDOM_TYPE,
								ID_G_101,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Random Type",L"Random_Type",
								0.f,0.0f, 1.0f, UINT_MAX,UINT_MAX,UINT_MAX);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_RENDER_READY,
								ID_G_102,
								siICENodeDataBool,siICENodeStructureSingle,siICENodeContextSingleton,
								L"Ready to Render",L"Ready_to_Render",
								false);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddInputPort(	ID_IN_PointID,
								ID_G_102,
								siICENodeDataLong,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"PointID",L"PointID",
								1);
	st.AssertSucceeded( ) ;

//-----------------------------------START OUTPUT PORTS----------------

	st = nodeDef.AddPortGroup(ID_G_200);
	st.AssertSucceeded( ) ;
	st = nodeDef.AddOutputPort(	ID_OUT_OCEAN,
								ID_G_200,
								siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"Ocean",L"Ocean");
	st.AssertSucceeded( ) ;

	st = nodeDef.AddPortGroup(ID_G_201);
	st.AssertSucceeded( ) ;
	st = nodeDef.AddOutputPort(	ID_OUT_FOAM,
								ID_G_201,
								siICENodeDataFloat,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"Foam",L"Foam");
	st.AssertSucceeded( ) ;

	st = nodeDef.AddPortGroup(ID_G_202);
	st.AssertSucceeded( ) ;
	st = nodeDef.AddOutputPort(	ID_OUT_EIGEN_MINUS,
								ID_G_202,
								siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"egnVector Minus",L"egnVector_Minus");
	st.AssertSucceeded( ) ;

	st = nodeDef.AddPortGroup(ID_G_203);
	st.AssertSucceeded( ) ;

	st = nodeDef.AddOutputPort(	ID_OUT_EIGEN_PLUS,
								ID_G_203,
								siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"egnVector Plus",L"egnVector_Plus");
	st.AssertSucceeded( ) ;

	st = nodeDef.AddPortGroup(ID_G_204);
	st.AssertSucceeded( ) ;
	st = nodeDef.AddOutputPort(	ID_OUT_NORMALS,
								ID_G_204,
								siICENodeDataVector3,siICENodeStructureSingle,siICENodeContextComponent0D,
								L"Normals",L"Normals");
	st.AssertSucceeded( ) ;

	PluginItem nodeItem = in_reg.RegisterICENode(nodeDef);
	nodeItem.PutCategories(L"aaOcean Suite");

	return CStatus::OK;
}
