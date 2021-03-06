#################################################
cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
#################################################

# ----- project name -----
PROJECT( CAFAnaModule )

set( ${PROJECT_NAME}_VERSION_MAJOR 1 )
set( ${PROJECT_NAME}_VERSION_MINOR 0 )
set( ${PROJECT_NAME}_VERSION_PATCH 0 )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo
    MinSizeRel."
  FORCE )
endif()

if( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
  set( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
endif()

# write this variable to cache
set( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )

set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
mark_as_advanced( CMAKE_INSTALL_RPATH )

#append link pathes to rpath list
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
mark_as_advanced( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# Find packages
find_package (GSL REQUIRED )
find_package( GENIE REQUIRED )

IF( GENIE_FOUND )
  find_package( LIBXML2 REQUIRED )
  find_package( PYTHIA6 REQUIRED )
  find_package( LOG4CPP REQUIRED )
ENDIF()

find_package( ROOT REQUIRED COMPONENTS Core RIO Hist Tree Rint XMLIO Geom Physics EGPythia6 GenVector )
add_definitions( ${ROOT_DEFINITIONS} )

set(CMAKE_CXX_FLAGS "-Werror -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing -std=c++17 ${CMAKE_CXX_FLAGS}")

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
endif()

#make libs
include_directories( BEFORE include ${ROOT_INCLUDE_DIRS} ${GENIE_INCLUDE_DIRS} ${LOG4CPP_INCLUDE_DIRS} ${PYTHIA6_INCLUDE_DIRS} ${LIBXML2_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} )
aux_source_directory( src SRC_FILES )

add_library( ${PROJECT_NAME} SHARED ${SRC_FILES} )
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} ${GENIE_LIBRARIES} ${LOG4CPP_LIBRARIES} ${PYTHIA6_LIBRARIES} ${LIBXML2_LIBRARIES} ${GSL_LIBRARIES} )

install( TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib )

#skimmer
add_executable(mcp_skimmer main/mcp_skimmer.cxx)
target_link_libraries(mcp_skimmer ${PROJECT_NAME})

#caf module
add_executable(cafanatree_module main/cafanatree_module.cxx)
target_link_libraries(cafanatree_module ${PROJECT_NAME})

#install in the bin folder
install( TARGETS mcp_skimmer RUNTIME DESTINATION bin )
install( TARGETS cafanatree_module RUNTIME DESTINATION bin )
