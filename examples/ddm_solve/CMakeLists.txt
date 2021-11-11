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
      target_compile_definitions(Ddm_first_kind_dir_${equation}_${dimension}_P1_${formulation}_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")
    endforeach()
  endforeach()
endforeach()

add_executable(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat ddm_first_kind_neu_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_neu_HE_2D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat ddm_first_kind_neu_HE_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_neu_HE_3D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_first_kind_neu_LA_3D_P1_indirect_hmat ddm_first_kind_neu_LA_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_neu_LA_3D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_neu_LA_3D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_first_kind_dir_LA_3D_P1_indirect_hmat ddm_first_kind_dir_LA_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_first_kind_dir_LA_3D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_first_kind_dir_LA_3D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat ddm_second_kind_dir_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_second_kind_dir_HE_2D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_second_kind_neu_HE_2D_P1_direct_hmat ddm_second_kind_neu_HE_2D_P1_direct_hmat.cpp)
target_link_libraries(Ddm_second_kind_neu_HE_2D_P1_direct_hmat bemtool)
target_compile_definitions(Ddm_second_kind_neu_HE_2D_P1_direct_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_combined_formulation_dir_HE_2D_P1_indirect_hmat ddm_combined_formulation_dir_HE_2D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_combined_formulation_dir_HE_2D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_combined_formulation_dir_HE_2D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")

add_executable(Ddm_combined_formulation_neu_LA_3D_P1_indirect_hmat ddm_combined_formulation_neu_LA_3D_P1_indirect_hmat.cpp)
target_link_libraries(Ddm_combined_formulation_neu_LA_3D_P1_indirect_hmat bemtool)
target_compile_definitions(Ddm_combined_formulation_neu_LA_3D_P1_indirect_hmat PRIVATE "-DWITH_HPDDM" "-DHPDDM_MKL" "-DHTOOL_MKL")