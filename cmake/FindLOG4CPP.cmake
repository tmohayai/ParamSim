#########################################################
# cmake module for finding LOG4CPP
#
#
# returns:
#   LOG4CPP_FOUND        : set to TRUE or FALSE
#   LOG4CPP_INCLUDE_DIRS : paths to gsl includes
#   LOG4CPP_LIBRARIES    : list of gsl libraries
#
# @author Eldwan Brianne, DESY
#########################################################

# find log4cpp-config (seems that there is a problem currently on cvmfs to execute it)
# SET( LOG4CPP_CONFIG_EXECUTABLE LOG4CPP_CONFIG_EXECUTABLE-NOTFOUND )
# MARK_AS_ADVANCED( LOG4CPP_CONFIG_EXECUTABLE )
# FIND_PROGRAM( LOG4CPP_CONFIG_EXECUTABLE log4cpp-config PATHS $ENV{LOG4CPP_FQ_DIR}/bin NO_DEFAULT_PATH )
# IF( NOT LOG4CPP_FQ_DIR )
#     FIND_PROGRAM( LOG4CPP_CONFIG_EXECUTABLE log4cpp-config )
# ENDIF()
#
# IF( LOG4CPP_CONFIG_EXECUTABLE )
#
#     # ==============================================
#     # ===          LOG4CPP_PREFIX                ===
#     # ==============================================
#
#     EXECUTE_PROCESS( COMMAND "${LOG4CPP_CONFIG_EXECUTABLE}" --prefix
#         OUTPUT_VARIABLE LOG4CPP_ROOT
#         OUTPUT_STRIP_TRAILING_WHITESPACE
#     )
#
#     # ==============================================
#     # ===          LIBXML2_LIBRARIES             ===
#     # ==============================================
#
#     EXECUTE_PROCESS( COMMAND "${LOG4CPP_CONFIG_EXECUTABLE}" --libs
#         OUTPUT_VARIABLE LOG4CPP_LIBRARIES
#         OUTPUT_STRIP_TRAILING_WHITESPACE
#     )
#
# ENDIF( LOG4CPP_CONFIG_EXECUTABLE )

SET( LOG4CPP_ROOT LOG4CPP_ROOT-NOTFOUND )
MARK_AS_ADVANCED( LOG4CPP_ROOT )
SET( LOG4CPP_ROOT $ENV{LOG4CPP_FQ_DIR} )

FIND_LIBRARY(LOG4CPP_LIBRARIES NAMES log4cpp PATHS ${LOG4CPP_ROOT} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)

# ---------- includes ---------------------------------------------------------
SET( LOG4CPP_INCLUDE_DIRS LOG4CPP_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( LOG4CPP_INCLUDE_DIRS )

SET( LOG4CPP_INCLUDE_DIRS NAMES log4cpp/config.h PATHS ${LOG4CPP_ROOT}/include NO_DEFAULT_PATH )
IF( NOT LOG4CPP_ROOT )
    FIND_PATH( LOG4CPP_INCLUDE_DIRS NAMES log4cpp/config.h )
ENDIF()

# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set GENIE_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( LOG4CPP DEFAULT_MSG LOG4CPP_ROOT LOG4CPP_LIBRARIES LOG4CPP_INCLUDE_DIRS )
