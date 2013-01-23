#ifndef OPENEXR_OUTPUT_H
#define OPENEXR_OUTPUT_H

#include <ImfArray.h>
#include <half.h>
#include <ImfRgbaFile.h>

void genFullFilePath(char* dest, const char* outputFolder, const char* postfix, const int frame)
{
	#ifdef _MSC_VER
	char* slash = "\\";
	#else
	char* slash = "/";
	#endif

	char cFrame[10];
	char prefixFrame[10];
	if ( frame < 10)
		strcpy(prefixFrame, "000");
	else if ( frame < 100) 
		strcpy(prefixFrame, "00");
	else if ( frame < 1000) 
		strcpy(prefixFrame, "0");
	sprintf(cFrame, "%s%d", prefixFrame, frame);
					
	strcpy(dest, outputFolder);
	strcat(dest, slash);

	if(strcmp(postfix,"") == 0)
		strcat(dest, "aaOceanData");
	else
		strcat(dest, "aaOceanData_");

	strcat(dest, postfix);
	strcat(dest, ".");
	strcat(dest, cFrame);
	strcat(dest, ".exr");
}

void writeExr(const char* outputFileName, int dimension, float *&red, float *&green, float *&blue, float *&alpha )
{
	Imf::Array2D<Imf::Rgba> pixels(dimension, dimension);

	#pragma omp parallel for
	for(int i = 0; i < dimension; i++)
	{
		for(int j = 0; j < dimension; j++)
		{
			pixels[j][i].g = green[i*dimension+j];
			if(red)
			{
				pixels[j][i].r = red[i*dimension+j];
				pixels[j][i].b = blue[i*dimension+j];
				pixels[j][i].a = alpha[i*dimension+j];
			}
			else
				pixels[j][i].r = pixels[j][i].b = pixels[j][i].a = 0.0f;
		}
	}
	Imf::RgbaOutputFile outA(&outputFileName[0], Imf::Header(dimension, dimension), Imf::WRITE_RGBA);
	outA.setFrameBuffer (&pixels[0][0], 1, dimension);
	outA.writePixels (dimension);
}

#endif /* OPENEXR_OUTPUT_H */