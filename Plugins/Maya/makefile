MAYAINC	 := 
MAYA_LIB :=
FFTWLOC  := ../../externals
GXXLOC	 := 
TARGETLOC:= 

INCLUDES := -I$(MAYAINC) -I$(FFTWLOC) -I../../externals/aaOcean/src -I../../externals/helpers 
LDFLAGS  := -fopenmp $(FFTWLOC)/libfftw3f.a $(FFTWLOC)libfftw3f_omp.a -L$(MAYA_LIB)/lib -lOpenMaya
CXXFLAGS := -w -O2 -fPIC -DLINUX -D_BOOL -DI386 -ffast-math
CXX	     := $(GXXLOC)g++
SOURCES  := aaOceanMayaNode.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanDeformer.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
