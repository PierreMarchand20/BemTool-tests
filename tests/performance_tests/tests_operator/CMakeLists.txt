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

foreach(equation ${equations})
  foreach(operator ${operators})

    # Tests with dense matices
    add_executable(Tests_${equation}_${operator}_2D_P1xP1_dense tests_${equation}_${operator}_2D_P1xP1_dense.cpp)
    target_include_directories(Tests_${equation}_${operator}_2D_P1xP1_dense PUBLIC ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(Tests_${equation}_${operator}_2D_P1xP1_dense bemtool)

    # Tests with hmats
    add_executable(Tests_${equation}_${operator}_2D_P1xP1_hmat tests_${equation}_${operator}_2D_P1xP1_hmat.cpp)
    target_link_libraries(Tests_${equation}_${operator}_2D_P1xP1_hmat bemtool)
    add_dependencies(build-performance-tests Tests_${equation}_${operator}_2D_P1xP1_hmat)
    add_test(NAME Tests_${equation}_${operator}_2D_P1xP1_hmat COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/Tests_${equation}_${operator}_2D_P1xP1_hmat)
  endforeach()
endforeach()

add_executable(Tests_HE_SL_3D_P1xP1_hmat tests_HE_SL_3D_P1xP1_hmat.cpp)
target_link_libraries(Tests_HE_SL_3D_P1xP1_hmat bemtool)
add_dependencies(build-performance-tests Tests_HE_SL_3D_P1xP1_hmat)
add_test(NAME Tests_HE_SL_3D_P1xP1_hmat COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/Tests_HE_SL_3D_P1xP1_hmat)
