# Requires FFTW, 
# To enable OpenEXR support uncomment the LDFLAGS and CXXFLAGS below, and comment out the existing ones

# This specific makefile links against single-precision, static version of FFTW. 
# Please make sure that threads are enabled for FFTW

DEPLOC		:= ../../../../../../sdks/deploy
LIBSLOC		:= $(DEPLOC)/lib
ARNOLD	 	:= $(DEPLOC)/arnold/Arnold-4.2.0.5-linux
GXXLOC	 	:= 
TARGETLOC	:= 


INCLUDES := -I$(ARNOLD)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR


LDFLAGS  := -fopenmp -L$(ARNOLD)/bin 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC) -lai

CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden -pthread -DWRITE_OPENEXR

CXX	 := $(GXXLOC)g++
SOURCES  := shaders/aaOceanArnoldShader.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanArnold.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
