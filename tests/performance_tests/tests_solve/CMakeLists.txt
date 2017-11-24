#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
set(equations "HE")
# list(APPEND equations "LA")
# list(APPEND equations "YU")

set(dimensions "2D")
list(APPEND dimensions "3D")

set(formulations "direct")
list(APPEND formulations "indirect")

foreach(equation ${equations})
  foreach(dimension ${dimensions})
    foreach(formulation ${formulations})
  # # Tests with dense matices
  # add_executable(Tests_potential_${equation}_2D_P1_dense tests_potential_${equation}_2D_P1_dense.cpp)
  # target_include_directories(Tests_potential_${equation}_2D_P1_dense PUBLIC ${EIGEN3_INCLUDE_DIR})
  # target_link_libraries(Tests_potential_${equation}_2D_P1_dense bemtool)

#   # Tests with hmats
#   add_executable(Tests_potential_${equation}_2D_P1_hmat tests_potential_${equation}_2D_P1_hmat.cpp)
#   target_link_libraries(Tests_potential_${equation}_2D_P1_hmat bemtool)
#   add_dependencies(build-performance-tests Tests_potential_${equation}_2D_P1_hmat)
#   add_test(NAME Tests_potential_${equation}_2D_P1_hmat COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/Tests_potential_${equation}_2D_P1_hmat)
#
# add_executable(Tests_solve_${equation}_2D_P1_dir_hmat tests_solve_${equation}_2D_P1_dir_hmat.cpp)
# target_link_libraries(Tests_solve_${equation}_2D_P1_dir_hmat bemtool)
#
# add_executable(Tests_solve_${equation}_3D_P1_dir_hmat tests_solve_${equation}_3D_P1_dir_hmat.cpp)
# target_link_libraries(Tests_solve_${equation}_3D_P1_dir_hmat bemtool)

      add_executable(Tests_ddm_${equation}_${dimension}_P1_${formulation}_dir_hmat tests_ddm_${equation}_${dimension}_P1_${formulation}_dir_hmat.cpp)
      target_link_libraries(Tests_ddm_${equation}_${dimension}_P1_${formulation}_dir_hmat bemtool)

    endforeach()
  endforeach()
endforeach()
