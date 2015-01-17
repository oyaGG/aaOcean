# hcustom build command on linux
# cd to houdini's install directory and source the following environment
source houdini.setup

# now 'cd' to the folder containing aaOceanSOP.c and run
hcustom -d -I../../../externals/aaOcean/src -I../../../externals/helpers aaOceanSOP.c