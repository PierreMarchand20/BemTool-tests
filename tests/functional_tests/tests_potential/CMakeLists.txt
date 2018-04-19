#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
set(equations "HE")
list(APPEND equations "LA")
list(APPEND equations "YU")


foreach(equation ${equations})
    # Tests with hmats
    add_executable(Tests_potential_${equation}_2D_P1_hmat tests_potential_${equation}_2D_P1_hmat.cpp)
    target_link_libraries(Tests_potential_${equation}_2D_P1_hmat bemtool)
    add_dependencies(build-tests-functional-potentials-${dimension}D Tests_potential_${equation}_2D_P1_hmat)
    add_test(NAME Tests_potential_${equation}_2D_P1_hmat COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/Tests_potential_${equation}_2D_P1_hmat)

endforeach()


add_dependencies(build-tests-functional-potentials build-tests-functional-potentials-2D)
# add_dependencies(build-tests-functional-potentials build-tests-functional-potentials-3D)
