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

      add_executable(Ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat.cpp)
      target_link_libraries(Ddm_dir_${equation}_${dimension}_P1_${formulation}_hmat bemtool)

    endforeach()
  endforeach()
endforeach()

add_executable(Ddm_neu_HE_2D_P1_indirect_hmat ddm_neu_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_neu_HE_2D_P1_indirect_hmat bemtool)

add_executable(Ddm_neu_HE_3D_P1_indirect_hmat ddm_neu_HE_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_neu_HE_3D_P1_indirect_hmat bemtool)