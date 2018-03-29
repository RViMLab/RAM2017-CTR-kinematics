cmake_minimum_required(VERSION 2.8.9)
project(ConcentricTube)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
set(CMAKE_LINK_FLAG "${CMAKE_LINK_LIBRARY_FLAG} -std=gnu++14")


option (CTR_DOUBLE "Compile Robot As Double" ON) 
option (CTR_FLOAT  "Compile Robot As Float" ON) 


if(CTR_DOUBLE)
  add_definitions(-DCT_ROBOT_DOUBLE)
endif(CTR_DOUBLE)
if(CTR_FLOAT)
  add_definitions(-DCT_ROBOT_FLOAT)
endif(CTR_FLOAT)

set(CTR_COMMON_ALIGNMENT "16" CACHE STRINGS "alignment parameter for classes")
set_property(CACHE CTR_COMMON_ALIGNMENT PROPERTY STRINGS 16 32 64 128)
add_definitions(-DERL_COMMON_ALIGNMENT=${CTR_COMMON_ALIGNMENT})
add_definitions(-DCTR_ALIGNMENT=${CTR_COMMON_ALIGNMENT})

set(CTR_MAXSECTION "3" CACHE STRINGS "maximal number of ctr sections")
add_definitions(-DCTR_MAXSECTION=${CTR_MAXSECTION})

# FLAGS
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(COMPILER_IS_CLANG 1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
    set(CMAKE_LINK_LIBRARY_FLAG "${CMAKE_LINK_LIBRARY_FLAG} -lc++ -lc++abi -lgomp")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(COMPILER_IS_GXX 1)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    option (COMPILE_XEON_PHI "Compile for the Xeon Phi (k1om)" OFF) 
    set(COMPILER_IS_IXX 1)
endif()

set(CTR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ctr_source" CACHE PATH "${CMAKE_CURRENT_SOURCE_DIR}")
set(CTR_EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ctr_external" CACHE PATH "${CMAKE_CURRENT_SOURCE_DIR}")


message("Start Adding Source Files")
set(SOURCES 
   ${CTR_DIR}/ctr_robot.cpp
)
set(HEADERS 
    ${CTR_DIR}/ctr_robot.h
)

include_directories("${CTR_DIR}")
include_directories("${CTR_EXTERNAL_DIR}/eigen")
include_directories("${CTR_EXTERNAL_DIR}/erl")


# EXECUTABLE
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(COMPILER_IS_CLANG 1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
    set(CMAKE_LINK_LIBRARY_FLAG "${CMAKE_LINK_LIBRARY_FLAG} -lc++ -lc++abi -lgomp")
    
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ")
    set(LIBRARY_NAME "ctr_kinematics_c")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(COMPILER_IS_GXX 1)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
    set(LIBRARY_NAME "ctr_kinematics_g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(COMPILER_IS_IXX 1)
    set(LIBRARY_NAME "ctr_kinematics_i")
endif()

message("Library Name: ${LIBRARY_NAME}")
#message("SOURCES: ${SOURCES}")
add_library(${LIBRARY_NAME} SHARED ${SOURCES})




set(ADD_LIBS ${ADD_LIBS} "pthread")      



############### EXECUTABLE ############################
set(EXECUTABLE_SOURCES ${CTR_DIR}/ctr_test.cpp )
set(EXECUTABLE_NAME ${LIBRARY_NAME}_test.bin )
add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCES})
target_link_libraries(${EXECUTABLE_NAME} ${LIBRARY_NAME})
##########################################


install (TARGETS ${LIBRARY_NAME} ${EXECUTABLE_NAME}
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)

install(FILES ${HEADERS}
        DESTINATION include
        )
