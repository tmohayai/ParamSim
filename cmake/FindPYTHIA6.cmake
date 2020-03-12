#########################################################
# cmake module for finding PYTHIA6
#
#
# returns:
#   PYTHIA6_FOUND        : set to TRUE or FALSE
#   PYTHIA6_INCLUDE_DIRS : paths to gsl includes
#   PYTHIA6_LIBRARIES    : list of gsl libraries
#
# @author Eldwan Brianne, DESY
#########################################################


SET( PYTHIA6_ROOT PYTHIA6_ROOT-NOTFOUND )
MARK_AS_ADVANCED( PYTHIA6_ROOT )
SET( PYTHIA6_ROOT $ENV{PYTHIA_FQ_DIR} )

FIND_LIBRARY(PYTHIA6_LIBRARIES NAMES Pythia6 pythia6 PATHS ${PYTHIA6_ROOT} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)

# ---------- includes ---------------------------------------------------------
SET( PYTHIA6_INCLUDE_DIRS PYTHIA6_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( PYTHIA6_INCLUDE_DIRS )

SET( PYTHIA6_INCLUDE_DIRS NAMES pyjets.inc PATHS ${PYTHIA6_ROOT}/include NO_DEFAULT_PATH )
IF( NOT PYTHIA6_ROOT )
    FIND_PATH( PYTHIA6_INCLUDE_DIRS NAMES pyjets.inc )
ENDIF()

# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set GENIE_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( PYTHIA6 DEFAULT_MSG PYTHIA6_ROOT PYTHIA6_LIBRARIES PYTHIA6_INCLUDE_DIRS )
