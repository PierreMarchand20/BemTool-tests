#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
set(equations "LA")
list(APPEND equations "HE")
list(APPEND equations "YU")

set(operators "SL")
list(APPEND operators "DL")
list(APPEND operators "TDL")
list(APPEND operators "HS")

set(dimensions "2")
list(APPEND dimensions "3")

add_custom_target(build-tests-functional-operators-2D)
add_custom_target(build-tests-functional-operators-3D)

foreach(equation ${equations})
  foreach(operator ${operators})
    foreach(dimension ${dimensions})

    # Tests with hmats
    add_executable(Tests_${equation}_${operator}_${dimension}D_P1xP1_hmat tests_${equation}_${operator}_${dimension}D_P1xP1_hmat.cpp)
    target_link_libraries(Tests_${equation}_${operator}_${dimension}D_P1xP1_hmat bemtool)
    add_dependencies(build-tests-functional-operators-${dimension}D Tests_${equation}_${operator}_${dimension}D_P1xP1_hmat)
    add_test(NAME Tests_${equation}_${operator}_${dimension}D_P1xP1_hmat COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/Tests_${equation}_${operator}_${dimension}D_P1xP1_hmat)
    endforeach()
  endforeach()
endforeach()

add_executable(Tests_LA_SL_2D_P0xP0_hmat tests_LA_SL_2D_P0xP0_hmat.cpp)
target_link_libraries(Tests_LA_SL_2D_P0xP0_hmat bemtool)

add_executable(Tests_HE_SL_2D_P2xP2_hmat tests_HE_SL_2D_P2xP2_hmat.cpp)
target_link_libraries(Tests_HE_SL_2D_P2xP2_hmat bemtool)


add_executable(Tests_HE_SL_2D_P0xP0_hmat tests_HE_SL_2D_P0xP0_hmat.cpp)
target_link_libraries(Tests_HE_SL_2D_P0xP0_hmat bemtool)


add_dependencies(build-tests-functional-operators build-tests-functional-operators-2D)
add_dependencies(build-tests-functional-operators build-tests-functional-operators-3D)
