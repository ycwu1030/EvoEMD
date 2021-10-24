
add_executable(SeesawFI_LG SeesawFI_LG.cpp SeesawTypeI.cpp)
target_link_libraries(SeesawFI_LG ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})
