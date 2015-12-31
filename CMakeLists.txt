cmake_minimum_required( VERSION 3.0.0 )

set( CMAKE_CXX_COMPILER /usr/bin/g++ )
set( CMAKE_C_COMPILER /usr/bin/gcc )

project( EHApplication )

#set( CMAKE_BUILD_TYPE Release )
set( CMAKE_BUILD_TYPE Debug )

#set( BUILD_SHARED_LIBS ON )


set( DEFAULT_FLAGS "-std=c++14 -Wall -lglfw -lGLEW -lGL" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEFAULT_FLAGS}" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEFAULT_FLAGS} -Werror -O0" )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DEFAULT_FLAGS} -O2" )

set( SOURCE_ROOT_DIR "~/Documents/Workspace/EHApplication" )

list( APPEND SOURCES
    ${SOURCE_ROOT_DIR}/GL/EHBuffer.cpp
    ${SOURCE_ROOT_DIR}/GL/EHDrawing.cpp
    ${SOURCE_ROOT_DIR}/GL/EHProgram.cpp
    ${SOURCE_ROOT_DIR}/GL/EHGL_Context.cpp
    )
list( APPEND SOURCES
    ${SOURCE_ROOT_DIR}/EHApplication.cpp
    ${SOURCE_ROOT_DIR}/EHApplication_window_glfw.cpp
    )

add_executable( EHApplication ${SOURCES} main.cpp engine.cpp )
