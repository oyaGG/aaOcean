# Requires FFTW, OpenEXR and Z-Lib libraries. Please make sure they are installed before running this makefile
# This specific makefile links against single-precision, static version of FFTW. 
# Please make sure that threads are enabled for FFTW
# If you are linking dynamically to FFTW, please use -lfftwf3_threads in LDFLAGS instead of $(FFTWLOC)libfftw3f.a

ARNOLD	 := 
FFTWLOC  := ../../../../externals/fftw
GXXLOC	 := 
TARGETLOC:= 
OPENEXRLOC:=

INCLUDES := -I$(ARNOLD)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(FFTWLOC) -I$(OPENEXRLOC)
LDFLAGS  := -fopenmp -L$(ARNOLD)/bin -lai $(FFTWLOC)/libfftw3f.a -lHalf -lIex -lIlmImf -lIlmThread -lImath -lz
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
