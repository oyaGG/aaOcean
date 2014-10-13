# makefile for softimage aaocean deformer

DEPLOC		:= ../../../../../../sdks/deploy
SOFTIMAGE	:= $(DEPLOC)/softimage/softimage2013
LIBSLOC		:= $(DEPLOC)/lib
GXXLOC		:= 
TARGETLOC	:= 

INCLUDES 	:= -I$(SOFTIMAGE)/include -I$(DEPLOC)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers
LDFLAGS		:= -fopenmp
LIBS		:= $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a
CXXFLAGS	:= -w -O2 -fPIC
CXX		:= $(GXXLOC)g++
SOURCES  	:= icenode.cpp
OBJECTS  	:= $(SOURCES:.cpp=.o)
TARGET   	:= $(TARGETLOC)aaOceanICE.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
