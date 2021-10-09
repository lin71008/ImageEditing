cmake_minimum_required(VERSION 3.1)

project(ImageEditing)
set(CMAKE_FIND_ROOT_PATH /c/MinGW/)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAEK_EXE_LINKER_FLAGS_INIT "-static-libgcc -static-libstdc++")
set(CMAKE_CREATE_WIN32_EXE  "/subsystem:windowsce -mwindows")

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/)

add_Definitions("-D_XKEYCHECK_H")
add_definitions(-DPROJECT_DIR="${PROJECT_SOURCE_DIR}")

add_library(targa ${SRC_DIR}libtarga.h ${SRC_DIR}libtarga.c)

add_executable(ImageEditing 
    ${SRC_DIR}Main.cpp
    ${SRC_DIR}Globals.h
    ${SRC_DIR}Globals.inl
    ${SRC_DIR}ImageWidget.h
    ${SRC_DIR}ImageWidget.cpp
    ${SRC_DIR}ScriptHandler.h
    ${SRC_DIR}ScriptHandler.cpp
    ${SRC_DIR}TargaImage.h
    ${SRC_DIR}TargaImage.cpp)

target_link_libraries(ImageEditing targa fltk fltk_images comctl32)
