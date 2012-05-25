MENTALSDKLOC:= 
GXXLOC	 := 
FFTWLOC  := 
TARGETLOC:= 

INCLUDES := -I$(MENTALSDKLOC) -I$(FFTWLOC)/libfftw3f.a -I../../../../aaOcean/src -I../../../../helpers
LDFLAGS  := -fopenmp $(FFTWLOC)libfftw3f.a
CXXFLAGS := -w -O2 -fPIC -fopenmp -fvisibility=hidden
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
