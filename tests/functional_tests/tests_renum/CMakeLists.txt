#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_boundary_edges tests_boundary_edges.cpp)
target_link_libraries(Tests_boundary_edges bemtool)

add_executable(Tests_renum_wo_boundary tests_renum_wo_boundary.cpp)
target_link_libraries(Tests_renum_wo_boundary bemtool)

add_executable(Tests_is_boundary_nodes tests_is_boundary_nodes.cpp)
target_link_libraries(Tests_is_boundary_nodes bemtool)
