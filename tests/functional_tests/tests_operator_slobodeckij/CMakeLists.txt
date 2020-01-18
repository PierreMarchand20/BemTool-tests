#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#

add_executable(Tests_operator_slobodeckij_2D_P1 tests_operator_slobodeckij_2D_P1.cpp)
target_link_libraries(Tests_operator_slobodeckij_2D_P1 bemtool)
add_dependencies(build-tests-functional-operators-slobodeckij Tests_operator_slobodeckij_2D_P1)
add_test(NAME Tests_operator_slobodeckij_2D_P1 COMMAND  ${CMAKE_CURRENT_BINARY_DIR}/Tests_operator_slobodeckij_2D_P1)

add_executable(Tests_operator_slobodeckij_3D_P1 tests_operator_slobodeckij_3D_P1.cpp)
target_link_libraries(Tests_operator_slobodeckij_3D_P1 bemtool)
add_dependencies(build-tests-functional-operators-slobodeckij Tests_operator_slobodeckij_3D_P1)
add_test(NAME Tests_operator_slobodeckij_3D_P1 COMMAND  ${CMAKE_CURRENT_BINARY_DIR}/Tests_operator_slobodeckij_3D_P1)