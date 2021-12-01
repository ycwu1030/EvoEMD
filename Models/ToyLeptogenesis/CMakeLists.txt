CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12 FATAL_ERROR)

# If this ToyLG example is used as standalone, please uncomment following command
# if (POLICY CMP0042)
#   cmake_policy(SET CMP0042 NEW)
# endif ()
# if (POLICY CMP0074)
#   cmake_policy(SET CMP0074 NEW)
# endif ()
# PROJECT(ToyLG)
# # SET(EVOEMD_ROOT "_PATH_TO_EVOEMD_")
# SET(CMAKE_MODULE_PATH "${EVOEMD_ROOT}/share/cmake/Modules;${CMAKE_MODULE_PATH}")
# FIND_PACKAGE(EVOEMD REQUIRED)
# INCLUDE_DIRECTORIES(${EVOEMD_INCLUDE_DIRS})
# SET(CMAKE_CXX_STANDARD 11)

add_executable(ToyLG Amplitudes.cpp Parameters.cpp ToyLG.cpp)
target_link_libraries(ToyLG ${GSL_LIBRARIES} ${EVOEMD_LIBRARIES})
