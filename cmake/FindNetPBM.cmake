# Try to find the NetPBM librairies
#  NETPBM_FOUND - system has NetPBM lib
#  NETPBM_INCLUDE_DIR - the NetPBM include directory
#  NETPBM_LIBRARIES - Libraries needed to use NetPBM

# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (NETPBM_INCLUDE_DIR AND NETPBM_LIBRARIES)
  # Already in cache, be silent
  set(NETPBM_FIND_QUIETLY TRUE)
endif (NETPBM_INCLUDE_DIR AND NETPBM_LIBRARIES)

find_path(NETPBM_INCLUDE_DIR NAMES pm.h PATHS /usr/local/netpbm/include)
find_library(NETPBM_LIBRARIES NAMES netpbm PATHS /usr/local/netpbm/lib)
#find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NETPBM DEFAULT_MSG NETPBM_INCLUDE_DIR NETPBM_LIBRARIES)

mark_as_advanced(NETPBM_INCLUDE_DIR NETPBM_LIBRARIES)
