#
# Try to find BEMTOOL library and include path.
# Once done this will define
#
# BEMTOOL_FOUND
# BEMTOOL_INCLUDE_DIR
#

FIND_PATH(
  BEMTOOL_INCLUDE_DIR
  NAMES bemtool/tools.hpp
  PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/../BemTool
    )

# Handle the QUIETLY and REQUIRED arguments and set the BEMTOOL_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BemTool DEFAULT_MSG
BEMTOOL_INCLUDE_DIR)

mark_as_advanced(BEMTOOL_INCLUDE_DIR)
