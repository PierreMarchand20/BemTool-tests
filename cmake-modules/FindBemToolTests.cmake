#
# Try to find BEMTOOL library and include path.
# Once done this will define
#
# BEMTOOLTESTS_FOUND
# BEMTOOLTESTS_INCLUDE_DIRS
#

FIND_PATH(
  BEMTOOLTESTS_INCLUDE_DIR
  NAMES bemtool-tests/tools.hpp
  PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/../BemTool-tests/include
    )

# Handle the QUIETLY and REQUIRED arguments and set the BEMTOOLTESTS_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BemToolTests DEFAULT_MSG
BEMTOOLTESTS_INCLUDE_DIR)

mark_as_advanced(BEMTOOLTESTS_INCLUDE_DIR)
if (BEMTOOLTESTS_FOUND)
    set(BEMTOOLTESTS_INCLUDE_DIRS ${BEMTOOLTESTS_INCLUDE_DIR} )
endif()
# message("${BEMTOOLTESTS_INCLUDE_DIRS}")
