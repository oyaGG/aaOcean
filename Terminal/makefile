# Requires FFTW, OPENEXR
# To enable OpenEXR support uncomment the LDFLAGS and CXXFLAGS below, and comment out the existing ones

# This specific makefile links against single-precision, static version of FFTW

# fftw, openexr, ilmbase and libz header and libs (change this if required)
DEPLOC      := /development/sdks/deploy

# aaocean's repository location (change this if required)
AAOCEANLOC    := /development/Repos/aaocean

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

# the compiler: gcc for C program, define as g++ for C++
CXX     := $(GXXLOC)g++

# compiler flags:
CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden -pthread

# includes
INCLUDES := -I$(AAOCEANLOC)/externals/aaOcean/src -I$(AAOCEANLOC)/externals/helpers -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR

# libs where exr libs and zlib are installed:
LIBSLOC     := $(DEPLOC)/lib

# libs to link against 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a $(LIBSLOC)/libz.a

# this should not be changed
LDFLAGS  := -fPIC -fopenmp -pthread -fvisibility=hidden -w -O2  -Wl,-rpath -Wl,$(LIBSLOC)

# the build target executable:
TARGET = aaOceanTerminal

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cpp $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

clean:
    $(RM) $(TARGET)