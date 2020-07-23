#########################################################
# cmake module for finding GENIE
#
#
# returns:
#   GENIE_FOUND        : set to TRUE or FALSE
#   GENIE_INCLUDE_DIRS : paths to gsl includes
#   GENIE_LIBRARIES    : list of gsl libraries
#
# @author Eldwan Brianne, DESY
#########################################################


# find genie-config
SET( GENIE_CONFIG_EXECUTABLE GENIE_CONFIG_EXECUTABLE-NOTFOUND )
MARK_AS_ADVANCED( GENIE_CONFIG_EXECUTABLE )
FIND_PROGRAM( GENIE_CONFIG_EXECUTABLE genie-config PATHS $ENV{GENIE_FQ_DIR}/bin NO_DEFAULT_PATH )
IF( NOT GENIE_FQ_DIR )
    FIND_PROGRAM( GENIE_CONFIG_EXECUTABLE genie-config )
ENDIF()

IF( GENIE_CONFIG_EXECUTABLE )

    # ==============================================
    # ===          GENIE_LIBRARIES               ===
    # ==============================================

    EXECUTE_PROCESS( COMMAND "${GENIE_CONFIG_EXECUTABLE}" --libs
        OUTPUT_VARIABLE GENIE_LIBRARIES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

ENDIF( GENIE_CONFIG_EXECUTABLE )

SET( GENIE_ROOT GENIE_ROOT-NOTFOUND )
MARK_AS_ADVANCED( GENIE_ROOT )
SET( GENIE_ROOT $ENV{GENIE_FQ_DIR} )

# ---------- includes ---------------------------------------------------------
SET( GENIE_INCLUDE_DIRS GENIE_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( GENIE_INCLUDE_DIRS )

SET( GENIE_INCLUDE_DIRS NAMES Framework/Utils/AppInit.h PATHS ${GENIE_ROOT}/include/GENIE NO_DEFAULT_PATH )
IF( NOT GENIE_ROOT )
    FIND_PATH( GENIE_INCLUDE_DIRS NAMES Framework/Utils/AppInit.h )
ENDIF()

# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set GENIE_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GENIE DEFAULT_MSG GENIE_ROOT GENIE_LIBRARIES GENIE_INCLUDE_DIRS )
