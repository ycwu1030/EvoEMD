# Finding EVOEMD

# Defines the following variables
# EVOEMD_INCLUDE_DIRS       - Location of EVOEMD's include directory.
# EVOEMD_LIBRARIES          - Location of EVOEMD's libraries.
# EVOEMD_FOUND              - True if EVOEMD has been found.

# One may provide a hint to where EVOEDM might be through EVOEMD_ROOT

FIND_PACKAGE(GSL 2.1 REQUIRED)
SET(_INCLUDES ${GSL_INCLUDE_DIRS})
SET(_LIBS ${GSL_LIBRARIES})

IF(EVOEMD_INCLUDES)
    SET(MULTINEST_FIND_QUIETLY TRUE)
ENDIF(EVOEMD_INCLUDES)

FIND_PATH(EVOEMD_ROOT_DIR
    NAMES include/EvoEMD/EvoEMD.h
    HINTS /usr/local ${EVOEMD_ROOT}
    DOC "EVOEMD root directory.")

FIND_PATH(_EVOEMD_INCLUDE_DIRS
    NAMES EvoEMD/EvoEMD.h
    HINTS ${EVOEMD_ROOT_DIR}/include
    DOC "EVOEMD include directory.")

SET(_EVOEMD_LIB_NAME "EVOEMD")
FIND_LIBRARY(_EVOEMD_LIBRARY
    NAMES ${_EVOEMD_LIB_NAME}
    HINTS ${EVOEMD_ROOT_DIR}/lib)

SET(EVOEMD_INCLUDE_DIRS ${_EVOEMD_INCLUDE_DIRS} ${_INCLUDES})
SET(EVOEMD_LIBRARIES ${_EVOEMD_LIBRARY} ${_LIBS})

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(EVOEMD DEFAULT_MSG EVOEMD_LIBRARIES EVOEMD_INCLUDE_DIRS)
MARK_AS_ADVANCED(EVOEMD_LIBRARIES EVOEMD_INCLUDE_DIRS)
