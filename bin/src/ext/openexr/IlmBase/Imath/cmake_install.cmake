# Install script for directory: /Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/Imath/libImath.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathBoxAlgo.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathBox.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathColorAlgo.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathColor.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathEuler.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathExc.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathExport.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathForward.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathFrame.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathFrustum.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathFrustumTest.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathFun.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathGL.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathGLU.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathHalfLimits.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathInt64.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathInterval.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathLimits.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathLineAlgo.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathLine.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathMath.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathMatrixAlgo.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathMatrix.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathNamespace.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathPlane.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathPlatform.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathQuat.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathRandom.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathRoots.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathShear.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathSphere.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathVecAlgo.h"
    "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase/Imath/ImathVec.h"
    )
endif()

