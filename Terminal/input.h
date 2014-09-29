#include "log.h"
#include "optionparser.h"

class input
{
public:
	int resolution;
	int seed;
	float oceanScale;
	float oceanDepth;
	float surfaceTension;
	float velocity;
	float smooth;
	float windDir;
	int windAlign;
	float reflectedWaves;
	float waveSpeed;
	float waveHeight;
	float waveChop;
    float repeatTime;
    char outputFolder[512];
    char postfix[512];

    int startFrame;
    int endFrame;
    float fps;

    int logLevel;

	input();
};

input::input()
{
	resolution      = 5;
	seed            = 5;
	oceanScale      = 100.f;
	oceanDepth      = 10000.f;
	surfaceTension  = 0.0f;
	velocity        = 30.0f;
	smooth          = 0.0f;
	windDir         = 0.0f;
	windAlign       = 1;
	reflectedWaves  = 0.2f;
	waveSpeed       = 1.0f;
	waveHeight      = 1.0f;
	waveChop        = 1.0f;
	
    repeatTime      = 10000.f;
    sprintf(outputFolder, "/tmp");
    sprintf(postfix, "");
    startFrame      = 1;
    endFrame        = 1;
    fps = 24;

    logLevel = 0;
}

// define options and usage
enum  optionIndex { UNKNOWN, HELP, RESOLUTION, SEED, OCEANSCALE, \
	OCEANDEPTH, SURFACETENSION, VELOCITY, SMOOTH, \
	WINDDIRECTION, WINDALIGN, REFLECTEDWAVES, WAVESPEED, \
	WAVEHEIGHT, WAVECHOP, STARTFRAME, ENDFRAME, FPS, OUTPUTFOLDER, POSTFIX, LOGLEVEL};

const option::Descriptor usage[] =
{
	{UNKNOWN, 0, "", "",						option::Arg::None,	   "USAGE: aaOcean [options]\n\nOptions:" },
	{RESOLUTION, 0,"res","resolution",			option::Arg::Optional, "  --resolution=<arg>" },
	{SEED, 0,"seed","seed",						option::Arg::Optional, "  --seed=<arg>" },
	{OCEANSCALE, 0,"scale","oceanscale",		option::Arg::Optional, "  --oceanscale=<arg>" },
	{OCEANDEPTH, 0,"depth","oceandepth",		option::Arg::Optional, "  --oceandepth=<arg>" },
	{SURFACETENSION, 0,"tens","surfacetension",	option::Arg::Optional, "  --surfacetension=<arg>" },
	{VELOCITY, 0,"vel","velocity",				option::Arg::Optional, "  --velocity=<arg>" },
	{SMOOTH, 0,"sm","smooth",					option::Arg::Optional, "  --smooth=<arg>" },
	{WINDDIRECTION, 0,"dir","winddir",			option::Arg::Optional, "  --winddir=<arg>" },
	{WINDALIGN, 0,"align","windalign",			option::Arg::Optional, "  --windalign=<arg>" },
	{REFLECTEDWAVES, 0,"ref","reflectedwaves",	option::Arg::Optional, "  --reflectedwaves=<arg>" },
	{WAVESPEED, 0,"sp","speed",					option::Arg::Optional, "  --speed=<arg>" },
	{WAVEHEIGHT, 0,"ht","waveheight",			option::Arg::Optional, "  --waveheight=<arg>" },
	{WAVECHOP, 0,"chp","wavechop",				option::Arg::Optional, "  --wavechop=<arg>" },
	{STARTFRAME, 0,"start","startframe",        option::Arg::Optional, "  --startframe=<arg>" },
    {ENDFRAME, 0,"end","endframe",              option::Arg::Optional, "  --endframe=<arg>" },
    {FPS, 0,"fps","fps",                        option::Arg::Optional, "  --fps=<arg>" },
    {OUTPUTFOLDER, 0,"o","outputfolder",        option::Arg::Optional, "  --outputfolder=<arg>" },
    {POSTFIX, 0,"pfix","postfix",               option::Arg::Optional, "  --postfix=<arg>" },
    {LOGLEVEL, 0,"log","loglevel",              option::Arg::Optional, "  --loglevel=<arg>" },
	{HELP, 0,"h", "help",                       option::Arg::None, "-h or --help  \tUSAGE: aaOcean [options]\n\nOptions:" },
	{UNKNOWN, 0, "", "",                        option::Arg::None, "\nExamples:\n example --unknown -- --this_is_no_option\n"},
	{0,0,0,0,0,0}
};

bool processInput(int argc, char** argv, input &oceanInput)
{
	FILELog::ReportingLevel() = logINFO;
    int level;
	char msg[512];

	argc-=(argc>0); argv+=(argc>0); 

	option::Stats  stats(usage, argc, argv);
	option::Option* options = new option::Option[stats.options_max];
	option::Option* buffer  = new option::Option[stats.buffer_max];
	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 0;

	if (argc == 0)
	{
		option::printUsage(std::cout, usage);
		return 0;
	}

	for (int i = 0; i < parse.optionsCount(); ++i)
	{
		option::Option& opt = buffer[i];
		switch (opt.index())
		{
		case HELP:
			option::printUsage(std::cout, usage);
			break;
        case LOGLEVEL:
            level = atoi(opt.arg);
            if(level > 0)
            {
                sprintf(msg, "logging DEBUG and INFO messages %d\n", atoi(opt.arg));
                FILELog::ReportingLevel() = logDEBUG;
			    LOG(logDEBUG) << msg;
            }
			break;
		case RESOLUTION:
			oceanInput.resolution = atoi(opt.arg);
			sprintf(msg, "--resolution with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case SEED:
			oceanInput.seed = atoi(opt.arg);
			sprintf(msg, "--seed with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case OCEANSCALE:
			oceanInput.oceanScale = (float)atof(opt.arg);
			sprintf(msg, "--oceanscale with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case OCEANDEPTH:
			oceanInput.oceanDepth = (float)atof(opt.arg);
			sprintf(msg, "--oceandepth with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case SURFACETENSION:
			oceanInput.surfaceTension = (float)atof(opt.arg);
			sprintf(msg, "--surfacetension with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case VELOCITY:
			oceanInput.velocity = (float)atof(opt.arg);
			sprintf(msg, "--velocity with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case SMOOTH:
			oceanInput.smooth = (float)atof(opt.arg);
			sprintf(msg, "--smooth with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case WINDDIRECTION:
			oceanInput.windDir = (float)atof(opt.arg);
			sprintf(msg, "--winddir with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case WINDALIGN:
			oceanInput.windAlign = atoi(opt.arg);
			sprintf(msg, "--windalign with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case REFLECTEDWAVES:
			oceanInput.reflectedWaves = (float)atof(opt.arg);
			sprintf(msg, "--reflectedwaves with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case WAVESPEED:
			oceanInput.waveSpeed = (float)atof(opt.arg);
			sprintf(msg, "--speed with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case WAVEHEIGHT:
			oceanInput.waveHeight = (float)atof(opt.arg);
			sprintf(msg, "--waveheight with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case WAVECHOP:
			oceanInput.waveChop = (float)atof(opt.arg);
			sprintf(msg, "--wavechop with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		case STARTFRAME:
            oceanInput.startFrame = atoi(opt.arg);
			sprintf(msg, "--startframe with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
        case ENDFRAME:
            oceanInput.endFrame = atoi(opt.arg);
			sprintf(msg, "--endframe with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
        case FPS:
            oceanInput.fps = (float)atof(opt.arg);
			sprintf(msg, "--fps with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
        case OUTPUTFOLDER:
            sprintf(oceanInput.outputFolder, opt.arg);
			sprintf(msg, "--outputfolder with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
        case POSTFIX:
            sprintf(oceanInput.postfix, opt.arg);
			sprintf(msg, "--postfix with  argument '%s'\n", opt.arg);
			LOG(logDEBUG) << msg;
			break;
		}
	}

	for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
		LOG(logWARNING) << "Unknown option: " << std::string(opt->name,opt->namelen) << "\n";

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
		LOG(logWARNING) << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";

	delete[] options;
	delete[] buffer;

    return 1;
}