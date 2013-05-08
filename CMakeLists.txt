#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(DGtalTools)

cmake_minimum_required (VERSION 2.6) 


FIND_PACKAGE(DGtal 0.6 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(DGtalTools_VERSION_MAJOR 0)
SET(DGtalTools_VERSION_MINOR 6)
SET(DGtalTools_VERSION_PATCH 0)
SET(DGTALTOOLS_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(PROJECT_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(VERSION ${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH})

SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

#------Tools subdirectories ------

add_subdirectory(volumetric)
add_subdirectory(shapeGenerator)
add_subdirectory(estimators)
add_subdirectory(converters)
add_subdirectory(volumetric/voltools)
add_subdirectory(visualisation)
