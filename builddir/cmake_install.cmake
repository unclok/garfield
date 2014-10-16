# Install script for directory: /pnfs/user/unclok/garfield/trunk

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/pnfs/user/unclok/garfield/garfield")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0"
         RPATH "")
  ENDIF(EXISTS "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
  FILE(INSTALL DESTINATION "/pnfs/user/unclok/garfield/trunk/Library" TYPE SHARED_LIBRARY FILES
    "/pnfs/user/unclok/garfield/builddir/lib/libGarfield.so.0.1.0"
    "/pnfs/user/unclok/garfield/builddir/lib/libGarfield.so"
    )
  IF(EXISTS "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
    FILE(RPATH_REMOVE
         FILE "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF(EXISTS "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/libGarfield.so.0.1.0")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "cmakefiles")
  IF(EXISTS "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake")
    FILE(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake"
         "/pnfs/user/unclok/garfield/builddir/CMakeFiles/Export/_pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake")
    IF(EXPORT_FILE_CHANGED)
      FILE(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends-*.cmake")
      IF(OLD_CONFIG_FILES)
        MESSAGE(STATUS "Old export file \"$ENV{DESTDIR}/pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        FILE(REMOVE ${OLD_CONFIG_FILES})
      ENDIF(OLD_CONFIG_FILES)
    ENDIF(EXPORT_FILE_CHANGED)
  ENDIF()
  FILE(INSTALL DESTINATION "/pnfs/user/unclok/garfield/trunk/Library/cmake" TYPE FILE FILES "/pnfs/user/unclok/garfield/builddir/CMakeFiles/Export/_pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends.cmake")
  IF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    FILE(INSTALL DESTINATION "/pnfs/user/unclok/garfield/trunk/Library/cmake" TYPE FILE FILES "/pnfs/user/unclok/garfield/builddir/CMakeFiles/Export/_pnfs/user/unclok/garfield/trunk/Library/cmake/GarfieldLibraryDepends-debug.cmake")
  ENDIF("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "cmakefiles")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "cmakefiles")
  FILE(INSTALL DESTINATION "/pnfs/user/unclok/garfield/trunk/Library/cmake" TYPE FILE FILES
    "/pnfs/user/unclok/garfield/builddir/GarfieldConfig.cmake"
    "/pnfs/user/unclok/garfield/builddir/GarfieldConfigVersion.cmake"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "cmakefiles")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/pnfs/user/unclok/garfield/builddir/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/pnfs/user/unclok/garfield/builddir/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
