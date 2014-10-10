# Requires FFTW, OPENEXR
# To enable OpenEXR support uncomment the LDFLAGS and CXXFLAGS below, and comment out the existing ones

# This specific makefile links against single-precision, static version of FFTW. 

FFTWLOC  	:= /development/sdks/deploy
OPENEXRLOC	:= /usr/local
AAOCEANLOC	:= /development/Repos/aaocean
TARGETLOC	:= 

INCLUDES := -I$(AAOCEANLOC)/externals/aaOcean/src -I$(AAOCEANLOC)/externals/helpers  -I$(FFTWLOC)/include -I$(OPENEXRLOC)/include/OpenEXR

LDFLAGS  := -fopenmp $(FFTWLOC)/lib/libfftw3f.a $(FFTWLOC)/lib/libfftw3f_omp.a -L/usr/local/lib $(OPENEXRLOC)/lib/libHalf.a $(OPENEXRLOC)/lib/libIex.a $(OPENEXRLOC)/lib/libIlmImf.a $(OPENEXRLOC)/lib/libIlmThread.a $(OPENEXRLOC)/lib/libImath.a /usr/local/lib/libIexMath.a -lz


#OPENEXR CXXFLAGS BELOW
CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden -DWRITE_OPENEXR

CXX	 := $(GXXLOC)g++
SOURCES  := aaOceanTerminal.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOcean

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
