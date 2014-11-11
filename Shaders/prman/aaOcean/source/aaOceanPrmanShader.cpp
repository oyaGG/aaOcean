#include "RixPattern.h"
#include "RixShadingUtils.h"
#include <string.h>

#include "aaOceanClass.cpp"
#ifdef WRITE_OPENEXR
#include "openEXROutput.h"
#endif /* WRITE_EXR */

#define getParam(param, ctx) (ctx->EvalParam(k_##param, -1, &param))

class aaOceanShader : public RixPattern
{
public:

    aaOceanShader();
    virtual ~aaOceanShader();

    virtual int Init(RixContext &, char const *pluginpath);
    virtual RixSCParamInfo const *GetParamTable();
    virtual void Finalize(RixContext &);

    virtual int CreateInstanceData(RixContext &,
                                   char const *handle,
                                   RixParameterList const *,
                                   InstanceData *id);

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
};


aaOceanShader::aaOceanShader() :
    m_writeFile(0),
    m_currentFrame(1),
    m_msg(NULL)
{
    strcpy (m_outputFolder, "");
    strcpy (m_postfix, "");
}

aaOceanShader::~aaOceanShader()
{
    
}

enum paramId
{
    k_outDisplacement=0,      // vector displacement output
    k_outEigenvalue, // eigenvalue (foam/spray) output
    k_uvCoords,
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
    k_currentTime,
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
        RixSCParamInfo("outDisplacement",k_RixSCColor, k_RixSCOutput),
        RixSCParamInfo("outEigenvalue",  k_RixSCFloat, k_RixSCOutput),

        // inputs - constant
        RixSCParamInfo("uvCoords",       k_RixSCVector),
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
        RixSCParamInfo("currentTime",    k_RixSCFloat),
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
    return 0;
}

// Finalize: Called when the plugin is unloaded from memory by the renderer.
void aaOceanShader::Finalize(RixContext &ctx)
{
 
}
// to clear shader instance-specific ocean data
static void cleanUp(void *oceandata)
{
     //delete static_cast<aaOcean*>(oceandata);
    aaOcean *pOcean = (aaOcean*)oceandata;
    delete pOcean;
    pOcean = 0;
}
int
aaOceanShader::CreateInstanceData(RixContext &ctx,
                                  char const *handle,
                                  RixParameterList const *plist,
                                  InstanceData *idata)
{
    // evaluate input parameters
    RtInt   resolution;
    RtFloat oceanScale;
    RtFloat oceanDepth;
    RtFloat surfaceTension;
    RtInt   seed;
    RtFloat waveHeight;
    RtFloat velocity;
    RtFloat waveSpeed;
    RtFloat chopAmount;
    RtFloat cutoff;
    RtFloat windDir;
    RtFloat damp;
    RtInt   windAlign;
    RtFloat currentTime;
    RtFloat repeatTime;
    RtInt   writeFile;
    RtInt   currentFrame;
    RtConstString outputFolder = NULL;
    RtConstString postfix = NULL;

    // retrieving critical ocean inputs first
    getParam(resolution, plist);
    getParam(oceanScale, plist);
    getParam(oceanDepth, plist);
    getParam(surfaceTension, plist);
    getParam(seed, plist);
    getParam(waveHeight, plist);
    getParam(velocity, plist);
    getParam(waveSpeed, plist);
    getParam(chopAmount, plist);
    getParam(cutoff, plist);
    getParam(windDir, plist);
    getParam(damp, plist);
    getParam(windAlign, plist);
    getParam(currentTime, plist);
    getParam(repeatTime, plist);

    //// openexr output params
    getParam(writeFile, plist);
    getParam(currentFrame, plist);
    getParam(outputFolder, plist);
    getParam(postfix, plist);
        
    // initialize_ocean
    aaOcean *pOcean = new aaOcean;
    pOcean->input(resolution, 
    seed,
    oceanScale, 
    oceanDepth,
    surfaceTension,
    velocity, 
    cutoff, 
    windDir, 
    windAlign, 
    damp, 
    waveSpeed, 
    waveHeight,
    chopAmount, 
    currentTime,
    repeatTime,
    TRUE,
    FALSE);

    int tileRes = (int)pow(2.0f, (4 + abs(resolution)));
    m_msg->WarningAlways("[aaOcean] Created new aaOcean %dx%d tile", tileRes, tileRes);

#ifdef WRITE_OPENEXR
    if(writeFile)
    {
        if(!dirExists(&outputFolder[0]))
            m_msg->WarningAlways("[aaOcean] Invalid folder path: %s", outputFolder);
        else
        {
            char outputFileName[255];
            sprintf(outputFileName, "none");
            oceanDataToEXR(pOcean, &outputFolder[0], &postfix[0], currentFrame, &outputFileName[0]);
            m_msg->WarningAlways("[aaOcean] Image written to %s", outputFileName);
        }
    }
#endif
    idata->datalen = sizeof(aaOcean*);
    idata->data = (void*) pOcean;
    idata->freefunc = cleanUp;
    return 0; // no error
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
    aaOcean *pOcean = (aaOcean*)instanceData;

    // working variables
    RtFloat2 const *uv2;
    RtVector3 const *uvCoords;
    RtInt const *useUVInput;
    RtFloat const *fade;
    RtFloat const *gamma;
    RtFloat const *brightness;
    RtFloat const *fMin;
    RtFloat const *fMax;
    RtInt const *invertFoam;

   
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

    RtColorRGB* outDisplacement = (RtColorRGB*) out[k_outDisplacement].value;
    RtFloat* outEigenvalue = (RtFloat*) out[k_outEigenvalue].value;

    if(!outDisplacement)
        outDisplacement = pool.AllocForPattern<RtColorRGB>(sctx->numPts);
    if(!outEigenvalue)
        outEigenvalue = pool.AllocForPattern<RtFloat>(sctx->numPts);

    // pulling UVs from primvars
    getParam(useUVInput, sctx);

    if(cInfo != k_RixSCNetworkValue && *useUVInput==0)
    {
        RtFloat const *uv2Width;
        RtFloat2 fill (0.f, 0.f);
        sctx->GetPrimVar("st", fill, &uv2, &uv2Width);
    }
    else
    {
    // TODO:
        getParam(uvCoords, sctx);
    }
    getParam(gamma, sctx);
    getParam(brightness, sctx);
    getParam(fMin, sctx);
    getParam(fMax, sctx);
    getParam(fade, sctx);
    getParam(invertFoam, sctx);

    // loop over our points and shade them
    for (int i = 0; i < sctx->numPts; i++)
    {
        if(outDisplacement)
        {
            outDisplacement[i].g = pOcean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eHEIGHTFIELD)* (1.0f - *fade);
            if(pOcean->isChoppy())
            {
                outDisplacement[i].r = pOcean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPX)* (1.0f - *fade);
                outDisplacement[i].b = pOcean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPZ)* (1.0f - *fade);
            }
            else
                outDisplacement[i].r = outDisplacement[i].b = 0.0f;
        }

        if(outEigenvalue)
        {
            outEigenvalue[i] = pOcean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eFOAM);

            if(*invertFoam)
                outEigenvalue[i]  = 1.0f - rescale(outEigenvalue[i], *fMin, *fMax, 0.0f, 1.0f);
            else
                outEigenvalue[i]  = rescale(outEigenvalue[i], *fMin, *fMax, 0.0f, 1.0f);

            outEigenvalue[i]  = maximum<float>(outEigenvalue[i], 0.0f);
            outEigenvalue[i]  = pow(outEigenvalue[i], *gamma);
            outEigenvalue[i] *= *brightness;
            outEigenvalue[i] *= (1.0f - *fade);
        }
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
