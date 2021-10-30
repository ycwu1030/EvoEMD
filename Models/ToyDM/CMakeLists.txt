
add_executable(ToyDM_FO Amplitudes.cpp ToyDM_FO.cpp)
target_link_libraries(ToyDM_FO ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})

add_executable(ToyDM_FI Amplitudes.cpp ToyDM_FI.cpp)
target_link_libraries(ToyDM_FI ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})
