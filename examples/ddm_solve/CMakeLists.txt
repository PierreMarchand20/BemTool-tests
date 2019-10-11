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

      add_executable(Ddm_first_kind_dir_${equation}_${dimension}_P1_${formulation}_hmat ddm_first_kind_dir_${equation}_${dimension}_P1_${formulation}_hmat.cpp)
      target_link_libraries(Ddm_first_kind_dir_${equation}_${dimension}_P1_${formulation}_hmat bemtool)
      target_compile_definitions(Ddm_first_kind_dir_${equation}_${dimension}_P1_${formulation}_hmat PRIVATE "-DWITH_HPDDM")
    endforeach()
  endforeach()
endforeach()

add_executable(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat ddm_first_kind_neu_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM")

add_executable(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat ddm_first_kind_neu_HE_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM")

add_executable(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat ddm_second_kind_dir_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM")