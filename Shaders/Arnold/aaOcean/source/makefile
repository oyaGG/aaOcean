# Requires FFTW, 
# To enable OpenEXR support uncomment the LDFLAGS and CXXFLAGS below, and comment out the existing ones

# This specific makefile links against single-precision, static version of FFTW. 
# Please make sure that threads are enabled for FFTW

ARNOLD	 	:= 
FFTWLOC  	:= ../../../../externals/fftw
GXXLOC	 	:= 
TARGETLOC	:= 
OPENEXRLOC	:=

INCLUDES := -I$(ARNOLD)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(FFTWLOC) -I$(OPENEXRLOC)

#OPENEXR LDFLAGS BELOW
#LDFLAGS  := -fopenmp -L$(ARNOLD)/bin -lai $(FFTWLOC)/libfftw3f.a $(FFTWLOC)libfftw3f_omp.a -lHalf -lIex -lIlmImf -lIlmThread -lImath -lz
LDFLAGS  := -fopenmp -L$(ARNOLD)/bin -lai $(FFTWLOC)/libfftw3f.a $(FFTWLOC)libfftw3f_omp.a

#OPENEXR CXXFLAGS BELOW
#CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden -DWRITE_OPENEXR
CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden

CXX	 := $(GXXLOC)g++
SOURCES  := aaOceanArnold.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanArnold.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
