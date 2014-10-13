# Makefile for mental ray aaOcean shader

DEPLOC		:= ../../../../../../sdks/deploy
LIBSLOC		:= $(DEPLOC)/lib
MENTALRAY 	:= $(DEPLOC)/softimage/softimage2013
GXXLOC	 	:= 
TARGETLOC	:= 

# you should not need to change anything below this line unless aaOcean's repo
# has not been deployed as stated in the wiki page

INCLUDES := -I$(MENTALRAY)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR
LDFLAGS  := -fopenmp 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libfftw3f.a $(LIBSLOC)/libfftw3f_omp.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC)

CXXFLAGS	:= -O2 -fPIC -fopenmp -fvisibility=hidden -pthread -DWRITE_OPENEXR
CXX		:= $(GXXLOC)g++
SOURCES		:= shaders/aaOceanData.cpp
OBJECTS		:= $(SOURCES:.cpp=.o)
TARGET		:= $(TARGETLOC)aaOceanShaders.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
