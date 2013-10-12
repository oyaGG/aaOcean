#ifndef OPENEXR_OUTPUT_H
#define OPENEXR_OUTPUT_H

#include <ImfNamespace.h>
#include <ImfOutputFile.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <ImfArray.h>
#include <half.h>
#include <ImfRgbaFile.h>

namespace CustomImf = OPENEXR_IMF_NAMESPACE;
using namespace CustomImf;
using namespace IMATH_NAMESPACE;

#if defined(_MSC_VER)
	#include <io.h> // For access().
#else
	#include <sys/io.h> // For access().
#endif

inline bool dirExists(const char* path)
{
	#ifdef _MSC_VER
		if ( _access(path, 0) == 0 )
		{
			return 1;
		}
		else
    		return 0;
	#else
		if ( access(path, 0) == 0 )
		{
			return 1;
		}
		else
    		return 0;
	#endif
}

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

void writeFullFloatExr (const char fileName[], const float *rPixels,  const float *gPixels, const float *bPixels, const float *aPixels,
	  int width,
	  int height)
{
 
    Header header (width, height);
    header.channels().insert ("R", Channel (FLOAT));
	header.channels().insert ("G", Channel (FLOAT));
	header.channels().insert ("B", Channel (FLOAT));
	header.channels().insert ("A", Channel (FLOAT));

    OutputFile file (fileName, header);

    FrameBuffer frameBuffer;

    frameBuffer.insert ("R",					// name
		        Slice (FLOAT,					// type
			       (char *) rPixels,			// base
			       sizeof (*rPixels) * 1,		// xStride
			       sizeof (*rPixels) * width));	// yStride
	
	frameBuffer.insert ("G",					// name
		        Slice (FLOAT,					// type
			       (char *) gPixels,			// base
			       sizeof (*gPixels) * 1,		// xStride
			       sizeof (*gPixels) * width));	// yStride

	frameBuffer.insert ("B",					// name
		        Slice (FLOAT,					// type
			       (char *) bPixels,			// base
			       sizeof (*bPixels) * 1,		// xStride
			       sizeof (*bPixels) * width));	// yStride

	frameBuffer.insert ("A",					// name
		        Slice (FLOAT,					// type
			       (char *) aPixels,			// base
			       sizeof (*aPixels) * 1,		// xStride
			       sizeof (*aPixels) * width));	// yStride


    file.setFrameBuffer (frameBuffer);
    file.writePixels (height);
}

void writeExr(const char* outputFileName, int dimension, float *&red, float *&green, float *&blue, float *&alpha )
{
	Array2D<float> rPixels (dimension, dimension);
	Array2D<float> gPixels (dimension, dimension);
	Array2D<float> bPixels (dimension, dimension);
	Array2D<float> aPixels (dimension, dimension);
	
	#pragma omp parallel for
	for(int i = 0; i < dimension; i++)
	{
		for(int j = 0; j < dimension; j++)
		{
			gPixels[j][i] = green[i*dimension+j];
			if(red)
			{
				rPixels[j][i] = red[i*dimension+j];
				bPixels[j][i] = blue[i*dimension+j];
				aPixels[j][i] = alpha[i*dimension+j];
			}
			else
				rPixels[j][i] = bPixels[j][i] = aPixels[j][i] = 0.0f;
		}
	}

	writeFullFloatExr(&outputFileName[0], &rPixels[0][0], &gPixels[0][0], &bPixels[0][0], &aPixels[0][0], dimension, dimension);
}
#endif /* OPENEXR_OUTPUT_H */