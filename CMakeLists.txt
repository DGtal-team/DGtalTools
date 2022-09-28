
#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(DGtalTools)

cmake_minimum_required (VERSION 3.1)
cmake_policy(SET CMP0057 NEW) # the new interpretation of IN_LIST is by default not set (at least until  3.21.1).

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
SET(DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD "100" CACHE INTERNAL "Threshold for the random selection of tools to build.")
SET(DGTALTOOLS_RANDOMIZED_BUILD_WHITELIST "" CACHE INTERNAL "List of whitelisted tools to build.")

# -----------------------------------------------------------------------------
# DGtal dependencies
# -----------------------------------------------------------------------------
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "Checking if DGtal is installed:")
FIND_PACKAGE(DGtal 1.2 REQUIRED)
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



# -----------------------------------------------------------------------------
# CLI11
# -----------------------------------------------------------------------------
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )



# -----------------------------------------------------------------------------
# Define variables and cmake parameters, and display some information
# -----------------------------------------------------------------------------
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
SET(DGtalTools_VERSION_MAJOR 1)
SET(DGtalTools_VERSION_MINOR 3)
SET(DGtalTools_VERSION_PATCH beta)
SET(DGTALTOOLS_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(PROJECT_VERSION "${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH}")
SET(VERSION ${DGtalTools_VERSION_MAJOR}.${DGtalTools_VERSION_MINOR}.${DGtalTools_VERSION_PATCH})




 # -----------------------------------------------------------------------------
 # Randomized build.
 # -----------------------------------------------------------------------------

message(STATUS "Randomized build")
if (DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD EQUAL "100")
   message(STATUS "All tools will be compiled (cmake variable DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD)")
else()
   message(STATUS "Only ~${DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD}% of the (randomly selected) tools will be  compiled.")
   message(STATUS "  (you can adjust this with the cmake variable DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD)")
   message(STATUS " Whitelist: ${DGTALTOOLS_RANDOMIZED_BUILD_WHITELIST}")
   message(STATUS "  (you can adjust this with the cmake variable DGTALTOOLS_RANDOMIZED_BUILD_WHITELIST)")
 endif()
message(STATUS "-------------------------------------------------------------------------------")



#-----------------------------------------------------------------------------
# ITK setting
#-----------------------------------------------------------------------------
if (  WITH_ITK )
    include(${ITK_USE_FILE})
ENDIF( WITH_ITK )


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
include(cmake/AddToolsFunction.cmake)
add_subdirectory(volumetric)
add_subdirectory(generators)
add_subdirectory(estimators)
add_subdirectory(converters)
add_subdirectory(visualisation)
add_subdirectory(geometricTransform)
add_subdirectory(imageProcessing)
