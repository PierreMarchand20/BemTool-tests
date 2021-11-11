#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(BemTool-tests)

# To force c++11
if (${CMAKE_VERSION} VERSION_LESS 3.1)
	add_compile_options(-std=c++11)
elseif (${CMAKE_VERSION} VERSION_LESS 3.6.3 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	add_compile_options(-std=c++11)
else()
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Files to do find_package for some module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

# Options
option(BEMTOOL_WITH_GUI      "Build bemtool visualization tools ?" OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#=============================================================================#
#========================== External Libraries ===============================#
#=============================================================================#

# Boost package - special functions
find_package(Boost REQUIRED)

# HPDDM
find_package(HPDDM)

# MPI
find_package (MPI REQUIRED)
message("-- MPI libraries found in " "${MPI_LIBRARIES}")
message("-- MPI include files found in " "${MPI_INCLUDE_PATH}")

# OPENMP
find_package(OpenMP)

# BLAS
find_package(BLAS REQUIRED)
message("-- Found Blas implementation:" "${BLAS_LIBRARIES}")

# LAPACK
find_package(LAPACK)
message("-- Found Lapack:" "${LAPACK_LIBRARIES}")


#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#



add_library(bemtool INTERFACE)
target_include_directories(bemtool INTERFACE
					"${CMAKE_CURRENT_SOURCE_DIR}/include"
					"external/bemtool/"
					${Boost_INCLUDE_DIRS}
					"external/htool/include/"
					${MPI_INCLUDE_PATH}
					${HPDDM_INCLUDE_DIR}
                    )

target_link_libraries(bemtool INTERFACE
					${MPI_LIBRARIES}
					${BLAS_LIBRARIES}
                    ${LAPACK_LIBRARIES}
					${OpenMP_CXX_LIBRARIES}
                    )


#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#
add_subdirectory(examples EXCLUDE_FROM_ALL)

enable_testing()
add_custom_target(build-tests-functional)
add_custom_target(build-tests-functional-operators)
add_custom_target(build-tests-functional-operators-slobodeckij)
add_custom_target(build-tests-functional-potentials)
add_subdirectory(tests EXCLUDE_FROM_ALL)
