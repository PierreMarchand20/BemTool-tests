#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#

# #=== lrmat_SVD
# if (EIGEN3_FOUND)
#   add_executable(Test_lrmat_SVD test_lrmat_SVD.cpp)
#   target_link_libraries(Test_lrmat_SVD bemtool
#                             ${EIGEN3_LIBRARIES})
#   target_include_directories(Test_lrmat_SVD PUBLIC
#                             ${EIGEN3_INCLUDE_DIR})
#   add_dependencies(build-tests Test_lrmat_SVD)
#   add_test(Test_lrmat_SVD Test_lrmat_SVD)
# endif()


set(compressions "partialACA")
# list(APPEND compressions "partialACA")
foreach(compression ${compressions})
    add_executable(Test_lrmat_${compression} test_lrmat_${compression}.cpp)
    target_link_libraries(Test_lrmat_${compression} bemtool)
    add_dependencies(build-tests Test_lrmat_${compression})
    add_test(Test_lrmat_${compression} Test_lrmat_${compression})
endforeach()
