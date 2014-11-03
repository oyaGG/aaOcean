#include "RixPattern.h"
#include "RixShadingUtils.h"
#include <string.h>

#include "aaOceanClass.cpp"
#ifdef WRITE_OPENEXR
#include "openEXROutput.h"
#endif /* WRITE_EXR */

// getUniformParam(param) (sctx->EvalParam(k_param, -1, &param, &m_param))
#define getUniformParam(param) (sctx->EvalParam(k_##param, -1, &param))

class aaOceanShader : public RixPattern
{
public:

    aaOceanShader();
    virtual ~aaOceanShader();

    virtual int Init(RixContext &, char const *pluginpath);
    virtual RixSCParamInfo const *GetParamTable();
    virtual void Finalize(RixContext &);

    virtual int ComputeOutputParams(RixShadingContext const *,
                                    RtInt *noutputs, 
                                    OutputSpec **outputs,
                                    RtConstPointer instanceData,
                                    RixSCParamInfo const *);
private:
    RtInt       m_writeFile;
    RtInt       m_currentFrame;
    RtInt       m_rotateUV;
    char        m_outputFolder[256];
    char        m_postfix[256];

    RixMessages *m_msg;
    RixMutex *m_mutex;

    bool isOceanReady;
    aaOcean *m_ocean;
};


aaOceanShader::aaOceanShader() :
    m_writeFile(0),
    m_currentFrame(1),
    m_msg(NULL),
    m_mutex(0L),
    isOceanReady(0)
{
    strcpy (m_outputFolder, "");
    strcpy (m_postfix, "");
}

aaOceanShader::~aaOceanShader()
{
    
}

enum paramId
{
    k_resultRGB=0,      // vector displacement output
    k_resultEigenvalue, // eigenvalue (foam/spray) output
    k_uv_coords,
    k_useUVInput,
    k_fade,
    k_resolution,
    k_oceanScale,
    k_oceanDepth,
    k_surfaceTension,
    k_seed,
    k_waveHeight,
    k_velocity,
    k_waveSpeed,
    k_chopAmount,
    k_cutoff,
    k_windDir,
    k_damp,
    k_windAlign,
    k_time,
    k_repeatTime,
    k_gamma,
    k_brightness,
    k_raw,
    k_invertFoam,
    k_fMin,
    k_fMax,
    k_writeFile,
    k_outputFolder,
    k_postfix,
    k_currentFrame,
    k_rotateUV,
    k_transform,
    k_numParams
};

RixSCParamInfo const * aaOceanShader::GetParamTable()
{
    static RixSCParamInfo s_ptable[] =
    {
        //outputs
        RixSCParamInfo("resultRGB",         k_RixSCColor, k_RixSCOutput),
        RixSCParamInfo("resultEigenvalue",  k_RixSCFloat, k_RixSCOutput),

        // inputs - constant
        RixSCParamInfo("uv_coords",      k_RixSCVector),
        RixSCParamInfo("useUVInput",     k_RixSCInteger),
        RixSCParamInfo("fade",           k_RixSCFloat),
        RixSCParamInfo("resolution",     k_RixSCInteger),
        RixSCParamInfo("oceanScale",     k_RixSCFloat),
        RixSCParamInfo("oceanDepth",     k_RixSCFloat),
        RixSCParamInfo("surfaceTension", k_RixSCFloat),
        RixSCParamInfo("seed",           k_RixSCInteger),
        RixSCParamInfo("waveHeight",     k_RixSCFloat),
        RixSCParamInfo("velocity",       k_RixSCFloat),
        RixSCParamInfo("waveSpeed",      k_RixSCFloat),
        RixSCParamInfo("chopAmount",     k_RixSCFloat),
        RixSCParamInfo("cutoff",         k_RixSCFloat),
        RixSCParamInfo("windDir",        k_RixSCFloat),
        RixSCParamInfo("damp",           k_RixSCFloat),
        RixSCParamInfo("windAlign",      k_RixSCInteger),
        RixSCParamInfo("time",           k_RixSCFloat),
        RixSCParamInfo("repeatTime",     k_RixSCFloat),
        RixSCParamInfo("gamma",          k_RixSCFloat),
        RixSCParamInfo("brightness",     k_RixSCFloat),
        RixSCParamInfo("raw",            k_RixSCInteger),
        RixSCParamInfo("invertFoam",     k_RixSCInteger),
        RixSCParamInfo("fMin",           k_RixSCFloat),
        RixSCParamInfo("fMax",           k_RixSCFloat),
        RixSCParamInfo("writeFile",      k_RixSCInteger),
        RixSCParamInfo("outputFolder",   k_RixSCString),
        RixSCParamInfo("postfix",        k_RixSCString),
        RixSCParamInfo("currentFrame",   k_RixSCInteger),
        RixSCParamInfo("rotateUV",       k_RixSCInteger),
        RixSCParamInfo("transform",      k_RixSCMatrix),

        RixSCParamInfo(), // empty, end of table
    };
    return &s_ptable[0];
}

// Init: Called when the plugin is first loaded by the renderer. 
// The plugin will remain loaded for the lifetime of the render. 
// Any global work that would be shared by all instances of a plugin 
// should be done here. Init returns 0 if there was no error initializing the plugin.
int aaOceanShader::Init(RixContext &ctx, char const *pluginpath)
{
    m_msg = (RixMessages*)ctx.GetRixInterface(k_RixMessages);
    if (!m_msg)
        return 1;

    RixThreadUtils* tu = (RixThreadUtils*)ctx.GetRixInterface(k_RixThreadUtils);
    m_mutex = tu->NewMutex();

    m_ocean = new aaOcean;

    return 0;
}

// Finalize: Called when the plugin is unloaded from memory by the renderer.
void aaOceanShader::Finalize(RixContext &ctx)
{
 #ifdef WRITE_OPENEXR
        if(m_writeFile)
        {
            if(!dirExists(&m_outputFolder[0]))
                m_msg->WarningAlways("[aaOcean] Invalid folder path: %s", m_outputFolder);
            else
            {
                char outputFileName[255];
                sprintf(outputFileName, "none");
                const char* postfix = &m_postfix[0];
                oceanDataToEXR(m_ocean, &m_outputFolder[0], &postfix[0], m_currentFrame, &outputFileName[0]);
                m_msg->WarningAlways("[aaOcean] Image written to %s", outputFileName);
            }
        }
        else
            m_msg->WarningAlways("[aaOcean] Write File false" );
#endif
    delete m_mutex;
    delete m_ocean;
}

// ComputeOutputParams: This  reads the input \
// parameters and computes the output parameters. It is called once per graph execution. \
// All outputs must be computed during this one call. The renderer provides a list of the \
// outputs it expects the plugin to compute. Most often, this is exactly the same as the \
// outputs declared in the parameter table.
int aaOceanShader::ComputeOutputParams(RixShadingContext const *sctx,
                                RtInt *noutputs, 
                                OutputSpec **outputs,
                                RtConstPointer instanceData,
                                RixSCParamInfo const *ignored)
{
    RixSCType type;
    RixSCConnectionInfo cInfo;

    // evaluate input parameters
    RtInt   const *resolution;
    RtFloat const *oceanScale;
    RtFloat const *oceanDepth;
    RtFloat const *surfaceTension;
    RtInt   const *seed;
    RtFloat const *waveHeight;
    RtFloat const *velocity;
    RtFloat const *waveSpeed;
    RtFloat const *chopAmount;
    RtFloat const *cutoff;
    RtFloat const *windDir;
    RtFloat const *damp;
    RtInt   const *windAlign;
    RtFloat const *time;
    RtFloat const *repeatTime;
    RtFloat const *gamma;
    RtFloat const *brightness;
    RtInt     const *raw;
    RtInt     const *invertFoam;
    RtFloat   const *fMin;
    RtFloat   const *fMax;
    RtInt     const *writeFile;
    RtConstString* outputFolder = NULL;
    RtConstString* postfix = NULL;
    RtInt     const *currentFrame;
    RtInt     const *rotateUV;
    RtMatrix4x4  const *transform;

    // retrieving critical ocean inputs first
    getUniformParam(resolution);
    getUniformParam(oceanScale);
    getUniformParam(oceanDepth);
    getUniformParam(surfaceTension);
    getUniformParam(seed);
    getUniformParam(waveHeight);
    getUniformParam(velocity);
    getUniformParam(waveSpeed);
    getUniformParam(chopAmount);
    getUniformParam(cutoff);
    getUniformParam(windDir);
    getUniformParam(damp);
    getUniformParam(windAlign);
    getUniformParam(time);
    getUniformParam(repeatTime);

    // retrieving non-critical ocean inputs
    getUniformParam(gamma);
    getUniformParam(brightness);
    getUniformParam(raw);
    getUniformParam(invertFoam);
    getUniformParam(fMin);
    getUniformParam(fMax);

    // openexr output params
    getUniformParam(writeFile);
    getUniformParam(currentFrame);
    getUniformParam(outputFolder);
    getUniformParam(postfix);
    m_writeFile = *writeFile;
    m_currentFrame = *currentFrame;
    strcpy(m_outputFolder, *outputFolder);
    strcpy(m_postfix, *postfix);

    // tranform parameters
    getUniformParam(rotateUV);
    sctx->EvalParam(k_transform, -1, &transform);
    //getUniformParam(transform);

    // check if ocean needs re-initialization
    // TODO: WAIT STATE! find a better way to do this
    m_mutex->Lock();
    if (isOceanReady == 0)
    {
        // initialize_ocean
        m_ocean->input(resolution[0], 
        seed[0],
        oceanScale[0], 
        oceanDepth[0],
        surfaceTension[0],
        velocity[0], 
        cutoff[0], 
        windDir[0], 
        windAlign[0], 
        damp[0], 
        waveSpeed[0], 
        waveHeight[0],
        chopAmount[0], 
        time[0],
        repeatTime[0],
        TRUE,
        FALSE);

        isOceanReady = 1;
        m_msg->WarningAlways("[aaOcean] Created new aaOcean %dx%d tile", (int)pow(2.0f, (4 + abs(resolution[0]))));

    }
    m_mutex->Unlock();

    // Find the number of outputs
    RixSCParamInfo const* paramTable = GetParamTable();
    int numOutputs = -1;
    while (paramTable[++numOutputs].access == k_RixSCOutput) {}

    // Allocate and bind our outputs
    RixShadingContext::Allocator pool(sctx);
    OutputSpec* out = pool.AllocForPattern<OutputSpec>(numOutputs);
    *outputs = out;
    *noutputs = numOutputs;

    // looping through the different output ids
    for (int i = 0; i < numOutputs; ++i)
    {
        out[i].paramId = i;
        out[i].detail = k_RixSCInvalidDetail;
        out[i].value = NULL;

        type = paramTable[i].type;    // we know this

        sctx->GetParamInfo(i, &type, &cInfo);
        if(cInfo == k_RixSCNetworkValue)
        {
            if( type == k_RixSCColor )
            {
                out[i].detail = k_RixSCVarying;
                out[i].value = pool.AllocForPattern<RtColorRGB>(sctx->numPts);
            }
            else if( type == k_RixSCFloat )
            {
                out[i].detail = k_RixSCVarying;
                out[i].value = pool.AllocForPattern<RtFloat>(sctx->numPts);
            }
        }
    }

    RtColorRGB* resultRGB = (RtColorRGB*) out[k_resultRGB].value;
    if(!resultRGB)
    {
        // make sure the resultRGB space is allocated because it
        // will store the composite color results.
        resultRGB = pool.AllocForPattern<RtColorRGB>(sctx->numPts);
    }

    RtFloat* resultEigenvalue = (RtFloat*) out[k_resultEigenvalue].value;
    //if(!resultEigenvalue)
    //    resultEigenvalue = pool.AllocForPattern<RtFloat>(sctx->numPts);

    // pulling UVs from primvars
    RtFloat2 const *uv2;
    RtFloat const *uv2Width;
    RtFloat2 fill (0.f, 0.f);
    sctx->GetPrimVar("st", fill, &uv2, &uv2Width);

    // loop over our points and shade them
    for (int i = 0; i < sctx->numPts; i++)
    {
        if(resultRGB)
        {
            resultRGB[i].g = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eHEIGHTFIELD);
            if(m_ocean->isChoppy())
            {
                resultRGB[i].r = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPX);
                resultRGB[i].b = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPZ);
            }
            else
                resultRGB[i].r = resultRGB[i].b = 0.0f;
        }

        if(resultEigenvalue)
            resultEigenvalue[i] = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eFOAM);
    }

    return 0;
}

RIX_PATTERNCREATE
{
    return new aaOceanShader();
}

RIX_PATTERNDESTROY
{
    delete ((aaOceanShader*)pattern);
}
