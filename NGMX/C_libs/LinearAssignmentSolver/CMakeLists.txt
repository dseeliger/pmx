set(CORELIBS ${GROMACS_LIBRARY_DIRS} m)

#set(CMAKE_C_FLAGS "-g -Wall")

include_directories(${GROMACS_INCLUDE_DIRS} ${G_PERMUTE_SOURCE_DIR}/src)
LINK_directories(${G_PERMUTE_BINARY_DIR})

add_library( lap lap.c )
target_link_libraries( lap ${CORELIBS} )
add_library( memory memory.c )
target_link_libraries( memory ${CORELIBS} )

add_executable(g_permute g_permute.c )

target_link_libraries( g_permute lap memory ${CORELIBS} )

