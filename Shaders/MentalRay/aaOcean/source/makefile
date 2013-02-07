# Requires FFTW, 
# To enable OpenEXR support uncomment the LDFLAGS and CXXFLAGS below, and comment out the existing ones

# This specific makefile links against single-precision, static version of FFTW. 
# Please make sure that threads are enabled for FFTW

MENTALSDKLOC	:= 
GXXLOC	 		:= 
FFTWLOC  		:= ../../../../externals/fftw
TARGETLOC		:= 
OPENEXRLOC		:=

INCLUDES := -I$(MENTALSDKLOC) -I$(FFTWLOC) -I$(OPENEXRLOC) -I../../../../externals/aaOcean/src -I../../../../externals/helpers

#OPENEXR LDFLAGS BELOW
#LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a $(FFTWLOC)libfftw3f_omp.a -lHalf -lIex -lIlmImf -lIlmThread -lImath -lz
LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a $(FFTWLOC)libfftw3f_omp.a 

#OPENEXR CXXFLAGS BELOW
#CXXFLAGS := -w -O3 -fPIC -fopenmp -DWRITE_OPENEXR

CXXFLAGS := -w -O3 -fPIC -fopenmp
CXX	 	 := $(GXXLOC)g++
SOURCES  := shaders/aaOceanData.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanShaders.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
