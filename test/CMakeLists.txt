

add_executable(RK_Test RungeKutta_Test.cpp)
target_link_libraries(RK_Test ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})

add_executable(NU_Test Neutrino_Test.cpp)
target_link_libraries(NU_Test ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})


add_executable(Hubble_Test Hubble_in_EMD.cpp)
target_link_libraries(Hubble_Test ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})
