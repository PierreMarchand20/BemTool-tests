add_subdirectory(ddm_solve)

add_executable(Compute_err compute_err.cpp)
target_link_libraries(Compute_err bemtool)

add_executable(Generate generate.cpp)
target_link_libraries(Generate bemtool)
target_compile_definitions(Generate PRIVATE "-DWITH_HPDDM -DHPDDM_MKL")

add_executable(Generate_sym generate_sym.cpp)
target_link_libraries(Generate_sym bemtool)
target_compile_definitions(Generate_sym PRIVATE "-DWITH_HPDDM -DHPDDM_MKL")
