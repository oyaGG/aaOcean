#include "RixPattern.h"
#include "RixShadingUtils.h"

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
    //m_transform
    m_msg(NULL)
{
}

aaOceanShader::~aaOceanShader()
{
}

int
aaOceanShader::Init(RixContext &ctx, char const *pluginpath)
{
    m_msg = (RixMessages*)ctx.GetRixInterface(k_RixMessages);

    if (!m_msg)
        return 1;

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

RixSCParamInfo const *
aaOceanShader::GetParamTable()
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

void
aaOceanShader::Finalize(RixContext &ctx)
{
}


int
aaOceanShader::ComputeOutputParams(RixShadingContext const *sctx,
                                RtInt *noutputs, OutputSpec **outputs,
                                RtConstPointer instanceData,
                                RixSCParamInfo const *ignored)
{
    RixSCType type;
    RixSCConnectionInfo cInfo;
    bool varying = true;
    bool uniform = false;

    // TODO: evaluate input parameters

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

    // TODO: write ouput
    // loop over our points
    for (int i = 0; i < sctx->numPts; i++)
    {
        // compute
        resultRGB[i].r = 0.0f;
        resultRGB[i].g = 1.1f;
        resultRGB[i].b = 0.0f;

        if(resultEigenvalue)
        {
            resultEigenvalue[i] = resultRGB[i].g;
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
