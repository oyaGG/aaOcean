# Makefile for mental ray aaOcean shader's Shader Definition Plugin for Softimage

DEPLOC		:= ../../../../../../../../sdks/deploy
XSISDK 		:= $(DEPLOC)/softimage/softimage2013
GXXLOC	 	:= 
TARGETLOC	:= 

# you should not need to change anything below this line unless aaOcean's repo
# has not been deployed as stated in the wiki page

INCLUDES := -I$(XSISDK)/include
LDFLAGS  := -Wl,-rpath="\\$$ORIGIN:\\$$ORIGIN/../lib:\\$$ORIGIN/../bin"
CXXFLAGS := -w -O2 -fPIC
CXX      := $(GXXLOC)g++
SOURCES  := shaderDefsPlugin.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := aaOceanShaderDefs.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
