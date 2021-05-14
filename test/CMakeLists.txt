

add_executable(RK_Test RungeKutta_Test.cpp)
target_link_libraries(RK_Test ${GSL_LIBRARIES} ${EVOEMD_LIBRARY})