#ifndef ICE_FUNCS_H
#define ICE_FUNCS_H

inline bool isPortDirty( CICEPortState& in_port)
{
	if( in_port.IsDirty(CICEPortState::siAnyDirtyState) )
	{
		in_port.ClearState();		
		return true;
	}
	else
		return false;
}

inline void fetchICEUserInput(ICENodeContext& in_ctxt, aaOcean *&ICEocean)
{
	//Check port states that require a re-evaluation of HoK
	CICEPortState pState_oceanScale( in_ctxt, ID_IN_OCEAN_SCALE);
	CICEPortState pState_windDir( in_ctxt, ID_IN_WINDDIR);
	CICEPortState pState_cutoff( in_ctxt, ID_IN_CUTOFF);
	CICEPortState pState_velocity( in_ctxt, ID_IN_WINDVELOCITY);
	CICEPortState pState_windAlign( in_ctxt, ID_IN_WINDALIGN );
	CICEPortState pState_damp( in_ctxt, ID_IN_DAMP);

	if( isPortDirty(pState_oceanScale) || isPortDirty(pState_windDir)   || isPortDirty(pState_cutoff) ||
		isPortDirty(pState_velocity)   || isPortDirty(pState_windAlign) || isPortDirty(pState_damp) )
	{
		CDataArrayFloat	oceanScale( in_ctxt, ID_IN_OCEAN_SCALE );
		CDataArrayFloat windDir( in_ctxt, ID_IN_WINDDIR	);
		CDataArrayFloat cutoff( in_ctxt, ID_IN_CUTOFF);
		CDataArrayFloat velocity( in_ctxt, ID_IN_WINDVELOCITY);
		CDataArrayLong	windAlign( in_ctxt, ID_IN_WINDALIGN );
		CDataArrayFloat damp( in_ctxt, ID_IN_DAMP);

		ICEocean->m_windDir		= DegsToRads(windDir[0]);
		ICEocean->m_cutoff		= fabs(cutoff[0] * 0.01f);
		ICEocean->m_oceanScale	= maximum<float>(oceanScale[0],0.00001f);
		ICEocean->m_velocity	= maximum<float>(((velocity[0]  * velocity[0]) / (aa_GRAVITY)),0.00001f);
		ICEocean->m_windAlign	= maximum<int>(((windAlign[0] + 1) * 2),2); 
		ICEocean->m_damp		= minimum<float>(damp[0],1.f);
		
		//only re-evaluate the expensive HoK Tessendorf function 
		//if user input changes for these ports
		ICEocean->m_redoHoK = true;
	}

	CICEPortState pState_waveHeight( in_ctxt, ID_IN_WAVE_HEIGHT);
	CICEPortState pState_waveSpeed( in_ctxt, ID_IN_WAVESPEED);
	CICEPortState pState_chop( in_ctxt, ID_IN_CHOP);
	if( isPortDirty(pState_waveHeight) || isPortDirty(pState_waveSpeed) || isPortDirty(pState_chop) )
	{
		CDataArrayFloat	waveHeight( in_ctxt, ID_IN_WAVE_HEIGHT);
		CDataArrayFloat	waveSpeed( in_ctxt, ID_IN_WAVESPEED);
		CDataArrayFloat	chop( in_ctxt, ID_IN_CHOP);
		ICEocean->m_chopAmount =  chop[0] * .01f;		//divided by scale for better ocean control
		ICEocean->m_waveHeight =  waveHeight[0] * .01f; //divided by scale for better ocean control
		ICEocean->m_waveSpeed	 =  waveSpeed[0];
	}
	
	CICEPortState pState_seed( in_ctxt, ID_IN_SEED);
	CICEPortState pState_ranType( in_ctxt, ID_IN_RANDOM_TYPE);
	if( isPortDirty(pState_seed) || isPortDirty(pState_ranType))
	{
		CDataArrayLong seed( in_ctxt, ID_IN_SEED);
		CDataArrayFloat ranType( in_ctxt, ID_IN_RANDOM_TYPE);
		ICEocean->m_randomType = ranType[0];
		ICEocean->m_seed	= seed[0];
		ICEocean->m_redoHoK = TRUE;
		ICEocean->setup_grid(); 	
	}
	//get grid dimensions 
	CDataArrayFloat gridU( in_ctxt, ID_IN_GRID_LENGTH_U);
	CDataArrayFloat gridV( in_ctxt, ID_IN_GRID_LENGTH_V);
	ICEocean->m_3DGridULength = gridU[0];
	ICEocean->m_3DGridVLength = gridV[0];

	//set current time and frame values
	ICEocean->m_time  = (float)in_ctxt.GetTime().GetTime(CTime::Seconds) * -1.0f;
}

void displayHeightField(CDataArrayVector3f &outData, aaOcean *&ICEocean)
{
	bool chop = FALSE;
	int index, i, j;	
	float x,z,mult1, mult2, halfRes;
	int gridRes = int_sqrt(ICEocean->m_pointCount) - 1;
	halfRes = (float)gridRes/2;
	mult1 = float(ICEocean->m_3DGridULength) / float(gridRes);
	mult2 = float(ICEocean->m_3DGridVLength) / float(gridRes);

	if(ICEocean->m_chopAmount > 0.0f)
		chop = TRUE;

	if(!ICEocean->m_renderReady)
	{
		int n = gridRes;
		if(chop)
		{
			#pragma omp parallel for private(i,j,index, x, z) 
			for(i = 0; i< n+1; i++)
			{					
				for(j = 0; j< n+1; j++)
				{
					index = i*(n+1)+j;
					x = (-halfRes + i) * mult1 ;
					z =	(-halfRes + j) * mult2 ;
					if(i==n)  //copy left-most col to right-most col
					{
						outData[index].PutX(x - ICEocean->m_fft_chopX[j][0]);
						outData[index].PutZ(z - ICEocean->m_fft_chopZ[j][0]);
						outData[index].PutY(ICEocean->m_fft_htField[j][0]);
					}
					if(j==n) // copy top row to bottom row
					{
						outData[index].PutX(x - ICEocean->m_fft_chopX[i*n][0]);
						outData[index].PutZ(z - ICEocean->m_fft_chopZ[i*n][0]);
						outData[index].PutY(ICEocean->m_fft_htField[i*n][0]);
					}
					if(i==n && j==n) // copy top-left corner to bottom-left
					{
						outData[index].PutX(x - ICEocean->m_fft_chopX[0][0]);
						outData[index].PutZ(z - ICEocean->m_fft_chopZ[0][0]);
						outData[index].PutY(ICEocean->m_fft_htField[0][0]);
					}
					if( i<n && j<n) // regular array copy
					{
						outData[index].PutX(x - ICEocean->m_fft_chopX[i*n+j][0]);
						outData[index].PutZ(z - ICEocean->m_fft_chopZ[i*n+j][0]);
						outData[index].PutY(ICEocean->m_fft_htField[i*n+j][0]);
					}
				}
			}
		}
		else
		{
			//NO CHOPINESS
			#pragma omp parallel for private(i,j,index, x, z) 
			for(i = 0; i< n+1; i++)
			{					
				for(j = 0; j< n+1; j++)
				{
					index = i*(n+1)+j;
					x = (-halfRes + i) * mult1 ;
					z =	(-halfRes + j) * mult2 ;
					outData[index].PutX(x);
					outData[index].PutZ(z);
					if(i==n)  //copy left-most col to right-most col
						outData[index].PutY(ICEocean->m_fft_htField[j][0]);
					if(j==n) // copy top row to bottom row
						outData[index].PutY(ICEocean->m_fft_htField[i*n][0]);
					if(i==n && j==n) // copy top-left corner to bottom-left
						outData[index].PutY(ICEocean->m_fft_htField[0][0]);
					if( i<n && j<n) // regular array copy
						outData[index].PutY(ICEocean->m_fft_htField[i*n+j][0]);
				}
			}
		}
	}
	else
	{
		//Rendering, display flat mesh
		int n = int_sqrt(ICEocean->m_pointCount);
		#pragma omp parallel for private(i,j,index, x, z) 
		for(i = 0; i< n; i++)
		{					
			for(j = 0; j< n; j++)
			{
				index = (i*n)+j;
				x = (-halfRes + i) * mult1 ;
				z =	(-halfRes + j) * mult2 ;
				outData[index].PutX(x);
				outData[index].PutZ(z);
				outData[index].PutY(0);				
			}
		}
	}
}

void displayFoam(CDataArrayFloat &outData, aaOcean *&ICEocean)
{
	int index;
	const int n = ICEocean->m_resolution;

	#pragma omp parallel for private(index)
	for(int i = 0; i< (n+1); i++)
	{					
		for(int j = 0; j< (n+1); j++)
		{
			index = i*(n+1)+j;
			if(i==n)  //copy left col to right col
				outData[index] = ICEocean->m_fft_jxz[j][0];
			if(j==n) // copy top row to bottom row
				outData[index] = ICEocean->m_fft_jxz[i*n][0];
			if(i==n && j==n)
				outData[index] = ICEocean->m_fft_jxz[0][0];
			if( i<n && j<n)
				outData[index] = ICEocean->m_fft_jxz[i*n+j][0];
		}
	}
}

void displayEigenMinus(CDataArrayVector3f &outData, aaOcean *&ICEocean)
{
	int index;
	const int n = ICEocean->m_resolution;

	#pragma omp parallel for private(index)
	for(int i = 0; i< (n+1); i++)
	{					
		for(int j = 0; j< (n+1); j++)
		{
			index = i*(n+1)+j;
			outData[index].PutY(0.0f);
			if(i==n)  //copy left col to right col
			{
				outData[index].PutX(ICEocean->m_eigenMinusX[j]);
				outData[index].PutZ(ICEocean->m_eigenMinusZ[j]);
			}
			if(j==n) // copy top row to bottom row
			{
				outData[index].PutX(ICEocean->m_eigenMinusX[i*n]);
				outData[index].PutZ(ICEocean->m_eigenMinusZ[i*n]);
			}
			if(i==n && j==n)
			{
				outData[index].PutX(ICEocean->m_eigenMinusX[0]);
				outData[index].PutZ(ICEocean->m_eigenMinusZ[0]);
			}
			if( i<n && j<n)
			{
				outData[index].PutX(ICEocean->m_eigenMinusX[i*n+j]);
				outData[index].PutZ(ICEocean->m_eigenMinusZ[i*n+j]);
			}			
		}
	}
}

void displayEigenPlus(CDataArrayVector3f &outData, aaOcean *&ICEocean)
{
	int index;
	const int n = ICEocean->m_resolution;

	#pragma omp parallel for private(index)
	for(int i = 0; i< (n+1); i++)
	{					
		for(int j = 0; j< (n+1); j++)
		{
			index = i*(n+1)+j;
			outData[index].PutY(0.0f);				
			if(i==n)  //copy left col to right col
			{
				outData[index].PutX(ICEocean->m_eigenPlusX[j]);
				outData[index].PutZ(ICEocean->m_eigenPlusZ[j]);
			}
			if(j==n) // copy top row to bottom row
			{
				outData[index].PutX(ICEocean->m_eigenPlusX[i*n]);
				outData[index].PutZ(ICEocean->m_eigenPlusZ[i*n]);
			}
			if(i==n && j==n)
			{
				outData[index].PutX(ICEocean->m_eigenPlusX[0]);
				outData[index].PutZ(ICEocean->m_eigenPlusZ[0]);
			}
			if( i<n && j<n)
			{
				outData[index].PutX(ICEocean->m_eigenPlusX[i*n+j]);
				outData[index].PutZ(ICEocean->m_eigenPlusZ[i*n+j]);
			}		
		}
	}
}

void displayNormals(CDataArrayVector3f &outData, aaOcean *&ICEocean)
{
	int index;
	const int n = ICEocean->m_resolution;

	#pragma omp parallel for private(index)
	for(int i = 0; i< (n+1); i++)
	{					
		for(int j = 0; j< (n+1); j++)
		{
			index = i*(n+1)+j;
			
			if(i==n)  //copy left col to right col
			{
				outData[index].PutY(ICEocean->m_fft_normY[j][0]);
				outData[index].PutX(ICEocean->m_fft_normX[j][0]);
				outData[index].PutZ(ICEocean->m_fft_normZ[j][0]);
			}
			if(j==n) // copy top row to bottom row
			{
				outData[index].PutY(ICEocean->m_fft_normY[i*n][0]);
				outData[index].PutX(ICEocean->m_fft_normX[i*n][0]);
				outData[index].PutZ(ICEocean->m_fft_normZ[i*n][0]);
			}
			if(i==n && j==n)
			{
				outData[index].PutY(ICEocean->m_fft_normY[0][0]);
				outData[index].PutX(ICEocean->m_fft_normX[0][0]);
				outData[index].PutZ(ICEocean->m_fft_normZ[0][0]);
			}
			if( i<n && j<n)
			{
				outData[index].PutY(ICEocean->m_fft_normY[i*n+j][0]);
				outData[index].PutX(ICEocean->m_fft_normX[i*n+j][0]);
				outData[index].PutZ(ICEocean->m_fft_normZ[i*n+j][0]);
			}			
		}
	}
}

#endif /*ICE_FUNCS*/
