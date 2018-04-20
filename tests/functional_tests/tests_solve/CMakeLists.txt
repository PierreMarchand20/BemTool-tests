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

      add_executable(Tests_ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat tests_ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat.cpp)
      target_link_libraries(Tests_ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat bemtool)

    endforeach()
  endforeach()
endforeach()


add_executable(Tests_proto_YU_2D_P1_indirect_hmat tests_proto_YU_2D_P1_indirect_hmat.cpp)
target_link_libraries(Tests_proto_YU_2D_P1_indirect_hmat bemtool)

# add_executable(Tests_ddm_ff tests_ddm_ff.cpp)
# target_link_libraries(Tests_ddm_ff bemtool)
#
# add_executable(Tests_ddm_dir_YU_2D_P1_indirect_hmat tests_ddm_dir_YU_2D_P1_indirect_hmat.cpp)
# target_link_libraries(Tests_ddm_dir_YU_2D_P1_indirect_hmat bemtool)
#
# add_executable(Tests_ddm_dir_YU_3D_P1_indirect_hmat tests_ddm_dir_YU_3D_P1_indirect_hmat.cpp)
# target_link_libraries(Tests_ddm_dir_YU_3D_P1_indirect_hmat bemtool)
#
# add_executable(Tests_ddm_neu_YU_3D_P1_indirect_hmat tests_ddm_neu_YU_3D_P1_indirect_hmat.cpp)
# target_link_libraries(Tests_ddm_neu_YU_3D_P1_indirect_hmat bemtool)
#
# add_executable(Tests_ddm_dir_LA_2D_P1_indirect_hmat tests_ddm_dir_LA_2D_P1_indirect_hmat.cpp)
# target_link_libraries(Tests_ddm_dir_LA_2D_P1_indirect_hmat bemtool)
#
# add_executable(Tests_ddm_neu_YU_2D_P1_indirect_hmat tests_ddm_neu_YU_2D_P1_indirect_hmat.cpp)
# target_link_libraries(Tests_ddm_neu_YU_2D_P1_indirect_hmat bemtool)
