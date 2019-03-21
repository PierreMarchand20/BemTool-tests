add_subdirectory(ddm_solve)

add_executable(Compute_err compute_err.cpp)
target_link_libraries(Compute_err bemtool)
