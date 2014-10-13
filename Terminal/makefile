# Requires the following packages
# fftw, download source from http://fftw.org/download.html
# ./configure CFLAGS="-fPIC" --enable-float --enable-openmp --enable-sse --prefix=
#ILMBASE, OPENEXR, zlib
#
# for portability, it is advised to link against the static libraries of the above packages
# build instructions are included in <aaocean_repo_path>/externals

# fftw, openexr, ilmbase and libz header and libs (change this if required)
DEPLOC      := ../../../sdks/deploy

# aaocean's repository location (change this if required)
AAOCEANLOC    := ../../aaocean

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
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC)

# this should not be changed
LDFLAGS  := -fPIC -fopenmp -pthread -fvisibility=hidden -w -O2  

# the build target executable:
TARGET = aaOceanTerminal

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cpp $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

clean:
	$(RM) $(TARGET)
