
add_executable(ToyDM Amplitudes.cpp ToyDM.cpp)
target_link_libraries(ToyDM ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})

# add_executable(ToyDM_FI Amplitudes.cpp ToyDM_FI.cpp)
# target_link_libraries(ToyDM_FI ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})
