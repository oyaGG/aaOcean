XSISDK	 := 
FFTWLOC  := 
GXXLOC	 := 
TARGETLOC:= 

INCLUDES := -I$(XSISDK)
LDFLAGS  := 
CXXFLAGS := -w -O2 -fPIC
CXX      := $(GXXLOC)g++
SOURCES  := shaderDefs.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := aaOceanShaderDefs.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
