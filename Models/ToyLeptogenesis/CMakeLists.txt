
add_executable(ToyLG Amplitudes.cpp Parameters.cpp ToyLG.cpp)
target_link_libraries(ToyLG ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})
