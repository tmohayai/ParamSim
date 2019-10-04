# Install script for directory: /dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so"
         RPATH "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/lib:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/build/libCAFAnaModule.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so"
         OLD_RPATH "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/lib:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCAFAnaModule.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module"
         RPATH "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/lib:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/build/cafanatree_module")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module"
         OLD_RPATH "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/build:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib:"
         NEW_RPATH "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/lib:/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+2.6-2.12-e17-prof/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cafanatree_module")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/dune/app/users/mtanaz/GArSoft/srcs/garsoft/Ana/ParamSim_debug_9/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
