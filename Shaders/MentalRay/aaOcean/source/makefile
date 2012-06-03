# Requires FFTW, OpenEXR and Z-Lib libraries. Please make sure they are installed before running this makefile
# This specific makefile links against single-precision, static version of FFTW. 
# Please make sure that threads are enabled for FFTW
# If you are linking dynamically to FFTW, please use -lfftwf3_threads in LDFLAGS instead of $(FFTWLOC)libfftw3f.a

MENTALSDKLOC:= 
GXXLOC	 := 
FFTWLOC  := ../../../../external/fftw
TARGETLOC:= 
OPENEXRLOC:=

INCLUDES := -I$(MENTALSDKLOC) -I$(FFTWLOC) -I$(OPENEXRLOC) -I../../../../external/aaOcean/src -I../../../../external/helpers
LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a -lHalf -lIex -lIlmImf -lIlmThread -lImath -lz
CXXFLAGS := -w -O2 -fPIC -fopenmp
CXX	 	 := $(GXXLOC)g++
SOURCES  := aaOceanShaders.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanShaders.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
