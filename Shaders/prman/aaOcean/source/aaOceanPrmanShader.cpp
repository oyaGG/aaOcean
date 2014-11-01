#include "RixPattern.h"
#include "RixShadingUtils.h"
#include "aaOceanClass.cpp"

class aaOceanShader : public RixPattern
{
public:

    aaOceanShader();
    virtual ~aaOceanShader();

    virtual int Init(RixContext &, char const *pluginpath);
    virtual RixSCParamInfo const *GetParamTable();
    virtual void Finalize(RixContext &);

    virtual int ComputeOutputParams(RixShadingContext const *,
                                    RtInt *n, RixPattern::OutputSpec **outputs,
                                    RtConstPointer instanceData,
                                    RixSCParamInfo const *);
private:
    RtVector3   m_uv_coords;     
    RtInt       m_use_uv_input;
    RtFloat     m_fade;  
    RtInt       m_resolution;
    RtFloat     m_oceanScale;
    RtFloat     m_oceanDepth;
    RtFloat     m_surfaceTension;
    RtInt       m_seed;
    RtFloat     m_waveHeight;
    RtFloat     m_velocity;
    RtFloat     m_waveSpeed;
    RtFloat     m_chopAmount;
    RtFloat     m_cutoff;
    RtFloat     m_windDir;
    RtFloat     m_damp;
    RtInt       m_windAlign;
    RtFloat     m_time;
    RtFloat     m_repeatTime;
    RtFloat     m_gamma;
    RtFloat     m_brightness;
    RtInt       m_raw;
    RtInt       m_invertFoam;
    RtFloat     m_fMin;
    RtFloat     m_fMax;
    RtInt       m_writeFile;
    RtString    m_outputFolder;
    RtString    m_postfix;
    RtInt       m_currentFrame;
    RtInt       m_rotateUV;
    RtMatrix    m_transform;

    RixMessages *m_msg;
    RixMutex *m_mutex;

    bool isOceanReady;

    aaOcean *m_ocean;

};


aaOceanShader::aaOceanShader() :
    m_uv_coords(1.0f, 1.0f, 1.0f),
    m_use_uv_input(0),
    m_fade(0.0f),
    m_resolution(4),
    m_oceanScale(100.0f),
    m_oceanDepth(10000.0f),
    m_surfaceTension(0.0f),
    m_seed(1),
    m_waveHeight(5.0f),
    m_velocity(4.5f),
    m_waveSpeed(1.0f),
    m_chopAmount(1.0f),
    m_cutoff(0.0f),
    m_windDir(45.0f),
    m_damp(0.985f),
    m_windAlign(0),
    m_time(0.1f),
    m_repeatTime(1000.f),
    m_gamma(1.0f),
    m_brightness(1.0f),
    m_raw(0),
    m_invertFoam(0),
    m_fMin(0.0f),
    m_fMax(0.0f),
    m_writeFile(0),
    m_outputFolder(""),
    m_postfix(""),
    m_currentFrame(1),
    m_rotateUV(0),
   // m_transform,
    m_msg(NULL),
    m_mutex(0L),
    isOceanReady(0)
{
}

aaOceanShader::~aaOceanShader()
{
    
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

enum paramId
{
    k_resultRGB=0,      // vector displacement output
    k_resultEigenvalue, // eigenvalue (foam/spray)
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

// Finalize: Called when the plugin is unloaded from memory by the renderer.
void aaOceanShader::Finalize(RixContext &ctx)
{
    delete m_mutex;
    delete m_ocean;
}

// ComputeOutputParams: This  reads the input \
// parameters and computes the output parameters. It is called once per graph execution. \
// All outputs must be computed during this one call. The renderer provides a list of the \
// outputs it expects the plugin to compute. Most often, this is exactly the same as the \
// outputs declared in the parameter table.
int aaOceanShader::ComputeOutputParams(RixShadingContext const *sctx,
                                RtInt *noutputs, OutputSpec **outputs,
                                RtConstPointer instanceData,
                                RixSCParamInfo const *ignored)
{
    RixSCType type;
    RixSCConnectionInfo cInfo;
    bool varying = 1;
    bool uniform = 0;

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
    RtBoolean const *raw;
    RtBoolean const *invertFoam;
    RtFloat   const *fMin;
    RtFloat   const *fMax;
    RtBoolean const *writeFile;
    RtString  const *outputFolder;
    RtString  const *postfix;
    RtInt     const *currentFrame;
    RtBoolean const *rotateUV;
    RtMatrix  const *transform;

    sctx->EvalParam(k_resolution, -1, &resolution, &m_resolution, uniform);
    sctx->EvalParam(k_oceanScale, -1, &oceanScale, &m_oceanScale, uniform);
    sctx->EvalParam(k_oceanDepth, -1, &oceanDepth, &m_oceanDepth, uniform);
    sctx->EvalParam(k_surfaceTension, -1, &surfaceTension, &m_surfaceTension, uniform);
    sctx->EvalParam(k_seed, -1, &seed, &m_seed, uniform);
    sctx->EvalParam(k_waveHeight, -1, &waveHeight, &m_waveHeight, uniform);
    sctx->EvalParam(k_velocity, -1, &velocity, &m_velocity, uniform);
    sctx->EvalParam(k_waveSpeed, -1, &waveSpeed, &m_waveSpeed, uniform);
    sctx->EvalParam(k_chopAmount, -1, &chopAmount, &m_chopAmount, uniform);
    sctx->EvalParam(k_cutoff, -1, &cutoff, &m_cutoff, uniform);
    sctx->EvalParam(k_windDir, -1, &windDir, &m_windDir, uniform);
    sctx->EvalParam(k_damp, -1, &damp, &m_damp, uniform);
    sctx->EvalParam(k_windAlign, -1, &windAlign, &m_windAlign, uniform);
    sctx->EvalParam(k_time, -1, &time, &m_time, uniform);
    sctx->EvalParam(k_repeatTime, -1, &repeatTime, &m_repeatTime, uniform);

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

        m_msg->WarningAlways("[aaOcean Shader] Initializing ocean. In Mutex zone" );
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
            resultRGB[i].r = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPX);
            resultRGB[i].g = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eHEIGHTFIELD);
            resultRGB[i].b = m_ocean->getOceanData(uv2[i][0], uv2[i][1], aaOcean::eCHOPZ);
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
