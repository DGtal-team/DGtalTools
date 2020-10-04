
#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(DGtalTools)

cmake_minimum_required (VERSION 3.1)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")


# -----------------------------------------------------------------------------
# DGtal dependencies
# -----------------------------------------------------------------------------
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "Checking if DGtal is installed:")
FIND_PACKAGE(DGtal 1.0 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")



 # -----------------------------------------------------------------------------
 # CPP11 
 # -----------------------------------------------------------------------------
 set (CMAKE_CXX_STANDARD 11)
 set (CMAKE_CXX_STANDARD_REQUIRED TRUE)



# -------------------------------------------------------------------------
# This test is for instance used for ITK v3.x. ITK forces a limited
# template depth. We remove this option.
# --------------------------------------------------------------------------
if (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")
    message( "Warning: some package has enabled a limited template depth for the C++ compiler." )
    message( "         Disabling option -ftemplate-depth-xx in CMAKE_CXX_FLAGS." )
    set( CMAKE_CXX_FLAGS_TMP ${CMAKE_CXX_FLAGS} )
    STRING( REGEX REPLACE "-ftemplate-depth-[0-9]*" ""
    CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_TMP}" )
     message ("         CMAKE_CXX_FLAGS=" ${CMAKE_CXX_FLAGS} )
endif (CMAKE_CXX_FLAGS MATCHES "-ftemplate-depth-[0-9]*")



# CLI11
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )



# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(DGtalTools_VERSION_MAJOR 1)
SET(DGtalTools_VERSION_MINOR 1)
SET(DGtalTools_VERSION_PATCH beta)
SET(DGTALTOOLS_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(PROJECT_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(VERSION ${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH})





 # -----------------------------------------------------------------------------
 # Looking for boost
 # -----------------------------------------------------------------------------

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.50.0 REQUIRED)
if ( Boost_FOUND )
  ADD_DEFINITIONS(${BOOST_DEFINITIONS} -DBOOST_ALL_NO_LIB)
  # SYSTEM to avoid warnings from boost.
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies}
     ${Boost_LIBRAIRIES})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})  
endif( Boost_FOUND )




# -----------------------------------------------------------------------------
# Documentation
# -----------------------------------------------------------------------------
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "Checking if doxygen/dot is installed:")
set(INSTALL_DOC_PATH ${CMAKE_INSTALL_PREFIX}/share/DGtalTools CACHE PATH "Installation directory for DGtal documentation files.")
INCLUDE(doxygen)
INCLUDE(TargetDoxygenDoc OPTIONAL)
message(STATUS "-------------------------------------------------------------------------------")

#------Tools subdirectories ------
add_subdirectory(volumetric)
add_subdirectory(generators)
add_subdirectory(estimators)
add_subdirectory(converters)
add_subdirectory(visualisation)
add_subdirectory(geometricTransform)
add_subdirectory(imageProcessing)
