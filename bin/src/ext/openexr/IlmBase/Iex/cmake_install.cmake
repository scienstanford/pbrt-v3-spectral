# Install script for directory: /Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/Iex/libIex.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIex.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIex.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIex.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexBaseExc.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexMathExc.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexErrnoExc.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexMacros.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/Iex.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexNamespace.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexExport.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Iex/IexForward.h"
    )
endif()

