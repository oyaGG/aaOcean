# hcustom build command on windows
# please make sure that you are using the compiler required for your version of houdini
# Launch Houdini Command Line Tools
# run the following

set MSVCDir="C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC"
cd %MSVCDir%
vcvarsall.bat

# now 'cd' to aaOceanSOP.c dir
set HCUSTOM_CFLAGS=-fp:fast -Ox -Oy -GL -D_OPENMP -openmp
hcustom -I ../../../externals/aaOcean/src -I ../../../externals/helpers aaOceanSOP.c

