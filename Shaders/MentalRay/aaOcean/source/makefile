# Makefile for mental ray aaOcean shader

# this makefile requires fftw, openexr, ilmbase and zlib
# please visit the following URL to see how to compile these
# https://bitbucket.org/amaanakram/aaocean/wiki/Home

# $(DEPLOC) is where you have installed 'dependencies' like fftw, openexr, ilmbase, and zlib
# the following folder structure is expected under $(DEPLOC)
# $(DEPLOC)/include
# $(DEPLOC)/lib
DEPLOC		:= ../../../../../../sdks/deploy

# this makefile expects the following folder structure under $(MENTALRAY) folder
# $(MENTALRAY)/include
# $(MENTALRAY)/lib
MENTALRAY 	:= $(DEPLOC)/softimage/softimage2013

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

LIBSLOC		:= $(DEPLOC)/lib
GXXLOC	 	:= 
TARGETLOC	:= 

INCLUDES := -I$(MENTALRAY)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR
LDFLAGS  := -fopenmp 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC)

CXXFLAGS	:= -O2 -fPIC -fopenmp -fvisibility=hidden -pthread -DWRITE_OPENEXR -ffast-math
CXX		:= $(GXXLOC)g++
SOURCES		:= shaders/aaOceanData.cpp
OBJECTS		:= $(SOURCES:.cpp=.o)
TARGET		:= $(TARGETLOC)aaOceanShaders.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
