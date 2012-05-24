ARNOLD	 := 
FFTWLOC  := 
GXXLOC	 := 
TARGETLOC:= ../deploy/

INCLUDES := -I$(ARNOLD)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(FFTWLOC)
LDFLAGS  := -fopenmp -L$(ARNOLD)/bin -lai $(FFTWLOC)/libfftw3f.a
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
