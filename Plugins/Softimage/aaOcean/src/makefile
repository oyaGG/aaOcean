XSISDKLOC:= 
GXXLOC	 := 
FFTWLOC  := ../../../../externals/fftw
TARGETLOC:= 

INCLUDES := -I$(XSISDKLOC)/include -I$(FFTWLOC) -I../../../../external/aaOcean/src -I../../../../external/helpers
LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a
CXXFLAGS := -w -O2 -fPIC -fvisibility=hidden
CXX	     := $(GXXLOC)g++
SOURCES  := icenode.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanICE.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
