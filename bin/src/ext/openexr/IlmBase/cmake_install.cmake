# Install script for directory: /Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/openexr/IlmBase

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES "/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/config/IlmBaseConfig.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/Half/cmake_install.cmake")
  include("/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/Iex/cmake_install.cmake")
  include("/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/IexMath/cmake_install.cmake")
  include("/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/Imath/cmake_install.cmake")
  include("/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/openexr/IlmBase/IlmThread/cmake_install.cmake")

endif()

