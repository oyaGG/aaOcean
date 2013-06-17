## This makefile is for gcc-4.2.x, which is ABI-compatible with gcc-4.1.2
MAYAINC	 := 
MAYA_LIB :=
FFTWLOC  := ../../externals
GXXLOC	 := 
TARGETLOC:= 

INCLUDES := -I$(MAYAINC) -I$(FFTWLOC) -I../../externals/aaOcean/src -I../../externals/helpers 
LDFLAGS  := -fopenmp $(FFTWLOC)/libfftw3f.a $(FFTWLOC)/libfftw3f_omp.a -L$(MAYA_LIB) -lOpenMaya
CXXFLAGS := -w -O2 -fPIC -DLINUX -D_BOOL -DI386 -ffast-math -DREQUIRE_IOSTREAM
CXX	     := $(GXXLOC)g++
SOURCES  := aaOceanDeformer.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanDeformer.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
