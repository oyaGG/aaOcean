# This makefile is for gcc-4.2.x, which is ABI-compatible with gcc-4.1.2
# makefile will also work for subsequent gcc versions, but please check
# ABI compatibility with your maya version

DEPLOC		:= ../../../../sdks/deploy
MAYA		:= $(DEPLOC)/maya/maya2014
MAYALIB		:= $(MAYA)/lib
LIBSLOC		:= $(DEPLOC)/lib
GXXLOC		:= 
TARGETLOC	:= 

INCLUDES	:= -I$(MAYA)/include -I$(DEPLOC)/include -I../../externals/aaOcean/src -I../../externals/helpers
LDFLAGS		:= -fopenmp
LIBS		:= $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a  -L$(MAYALIB) -lOpenMaya -lOpenMayaAnim
CXXFLAGS	:= -w -O2 -fPIC -fvisibility=hidden -DLINUX -D_BOOL -DI386 -ffast-math -DREQUIRE_IOSTREAM
CXX		:= $(GXXLOC)g++
SOURCES		:= aaOceanDeformer.cpp
OBJECTS		:= $(SOURCES:.cpp=.o)
TARGET		:= $(TARGETLOC)aaOceanDeformer.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
