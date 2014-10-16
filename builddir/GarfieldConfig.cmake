# - Config file for the Garfield package
# It defines the following variables
#  GARFIELD_INCLUDE_DIRS - include directories for Garfield
 
SET(GARFIELD_INCLUDE_DIRS
  "/pnfs/user/unclok/garfield/trunk/Include"
  "/pnfs/user/unclok/garfield/trunk/Heed/heed++/code")
  
SET( GARFIELD_SOURCE_DIR "/pnfs/user/unclok/garfield/trunk")
SET( GARFIELD_BINARY_DIR "/pnfs/user/unclok/garfield/builddir")
 
include("/pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake")
