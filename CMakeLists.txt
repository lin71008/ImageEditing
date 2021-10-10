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

enable_testing()

add_test(NAME test COMMAND mkdir -p test_result)
add_test(NAME gray COMMAND ImageEditing -headless ../test/gray.txt)
add_test(NAME quant-unif COMMAND ImageEditing -headless ../test/quant-unif.txt)
add_test(NAME quant-pop COMMAND ImageEditing -headless ../test/quant-pop.txt)
add_test(NAME dither-thresh COMMAND ImageEditing -headless ../test/dither-thresh.txt)
add_test(NAME dither-rand COMMAND ImageEditing -headless ../test/dither-rand.txt)
add_test(NAME dither-fs COMMAND ImageEditing -headless ../test/dither-fs.txt)
add_test(NAME dither-bright COMMAND ImageEditing -headless ../test/dither-bright.txt)
add_test(NAME dither-cluster COMMAND ImageEditing -headless ../test/dither-cluster.txt)
add_test(NAME dither-color COMMAND ImageEditing -headless ../test/dither-color.txt)
add_test(NAME filter-box COMMAND ImageEditing -headless ../test/filter-box.txt)
add_test(NAME filter-bartlett COMMAND ImageEditing -headless ../test/filter-bartlett.txt)
add_test(NAME filter-gauss COMMAND ImageEditing -headless ../test/filter-gauss.txt)
add_test(NAME filter-gauss-n COMMAND ImageEditing -headless ../test/filter-gauss-n.txt)
add_test(NAME filter-edge COMMAND ImageEditing -headless ../test/filter-edge.txt)
add_test(NAME filter-enhance COMMAND ImageEditing -headless ../test/filter-enhance.txt)
add_test(NAME npr-paint COMMAND ImageEditing -headless ../test/npr-paint.txt)
add_test(NAME half COMMAND ImageEditing -headless ../test/half.txt)
add_test(NAME double COMMAND ImageEditing -headless ../test/double.txt)
add_test(NAME scale COMMAND ImageEditing -headless ../test/scale.txt)
add_test(NAME rotate COMMAND ImageEditing -headless ../test/rotate.txt)
