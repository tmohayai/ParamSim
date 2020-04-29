#########################################################
# cmake module for finding LIBXML2
#
#
# returns:
#   LIBXML2_FOUND        : set to TRUE or FALSE
#   LIBXML2_INCLUDE_DIRS : paths to gsl includes
#   LIBXML2_LIBRARIES    : list of gsl libraries
#
# @author Eldwan Brianne, DESY
#########################################################


# find genie-config
SET( LIBXML2_CONFIG_EXECUTABLE LIBXML2_CONFIG_EXECUTABLE-NOTFOUND )
MARK_AS_ADVANCED( LIBXML2_CONFIG_EXECUTABLE )
FIND_PROGRAM( LIBXML2_CONFIG_EXECUTABLE xml2-config PATHS $ENV{LIBXML2_FQ_DIR}/bin NO_DEFAULT_PATH )
IF( NOT LIBXML2_FQ_DIR )
    FIND_PROGRAM( LIBXML2_CONFIG_EXECUTABLE xml2-config )
ENDIF()

IF( LIBXML2_CONFIG_EXECUTABLE )

    # ==============================================
    # ===          LIBXML2_PREFIX                ===
    # ==============================================

    EXECUTE_PROCESS( COMMAND "${LIBXML2_CONFIG_EXECUTABLE}" --prefix
        OUTPUT_VARIABLE LIBXML2_ROOT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # ==============================================
    # ===          LIBXML2_LIBRARIES             ===
    # ==============================================

    EXECUTE_PROCESS( COMMAND "${LIBXML2_CONFIG_EXECUTABLE}" --libs
        OUTPUT_VARIABLE LIBXML2_LIBRARIES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

ENDIF( LIBXML2_CONFIG_EXECUTABLE )

# ---------- includes ---------------------------------------------------------
SET( LIBXML2_INCLUDE_DIRS LIBXML2_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( LIBXML2_INCLUDE_DIRS )

SET( LIBXML2_INCLUDE_DIRS NAMES libxml/xmlversion.h PATHS ${LIBXML2_ROOT}/include NO_DEFAULT_PATH )
IF( NOT LIBXML2_ROOT )
    FIND_PATH( LIBXML2_INCLUDE_DIRS NAMES libxml/xmlversion.h )
ENDIF()

# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set GENIE_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( LIBXML2 DEFAULT_MSG LIBXML2_ROOT LIBXML2_LIBRARIES LIBXML2_INCLUDE_DIRS )
