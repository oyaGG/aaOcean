MENTALSDKLOC:= 
GXXLOC	 := 
FFTWLOC  := ../../../../external/fftw
TARGETLOC:= 

INCLUDES := -I$(MENTALSDKLOC) -I$(FFTWLOC)/libfftw3f.a -I../../../../external/aaOcean/src -I../../../../external/helpers
LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a
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
