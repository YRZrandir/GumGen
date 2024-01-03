cmake_minimum_required(VERSION 3.0)
project(Elasticity)
set(CMAKE_CXX_STANDARD 20)
set(EMROOT D:/dependencies/emsdk/upstream/emscripten)
set(CMAKE_TOOLCHAIN_FILE ${EMROOT}/cmake/Modules/Platform/Emscripten.cmake)
set(CMAKE_CXX_COMPILER ${EMROOT}/em++.bat)
set(CMAKE_C_COMPILER ${EMROOT}/emcc.bat)
set(CMAKE_BUILD_TYPE Release)
set(EMLIBS "D:/dev/emlibs")
#GLM
include_directories( ../dep/glm )
#Eigen
include_directories(../dep/eigen-3.4.0/include/eigen3)
set(Eigen3_DIR ../dep/eigen-3.4.0/share/eigen3/cmake)
#cgal
include_directories(../dep/CGAL)
#ceres
#set(MINIGLOG ON)
#include_directories(${EMLIBS}/em_ceres/include)
#boost
include_directories(../dep/boost)
#gmp
include_directories(${EMLIBS}/em_gmp/include)
#mpfr
include_directories(${EMLIBS}/em_mpfr/include)
#json
include_directories(../dep/json/include)
#emscripten
include_directories(${EMROOT}/cache/sysroot/include)

include_directories(${EMLIBS}/em_bezier/include)
include_directories(../dep/libigl/include)

aux_source_directory(../src/ SRCS)
aux_source_directory(../src/util UTIL_SRCS)
aux_source_directory(../src/gumgen GUM_SRCS)
aux_source_directory(../src/Bezier BEZIER_SRC)

add_executable(GumGen
        ${SRCS} ${UTIL_SRCS} ${GUM_SRCS} ${BEZIER_SRC})

set_target_properties(GumGen
        PROPERTIES SUFFIX ".mjs"
        COMPILE_FLAGS "-O3 --no-entry "
        LINK_FLAGS " -O3\
        -s \"EXPORTED_RUNTIME_METHODS=['cwrap','ccall', 'getValue', 'setValue']\" \
        -s \"EXPORTED_FUNCTIONS=['_CreateGumInit', '_CreateStep', '_malloc', '_free']\"  \
        -sINITIAL_MEMORY=2684354560 -sSTACK_SIZE=83886080 \
        -sSINGLE_FILE=1 -sASSERTIONS=1\
        --no-entry -s ENVIRONMENT='web' -sMODULARIZE=1 -s EXPORT_NAME='GumGenModule'")

#target_compile_definitions(GumGen PRIVATE CGAL_PMP_USE_CERES_SOLVER)
target_link_libraries(
        GumGen 
        ${EMLIBS}/em_gmp/lib/libgmp.a 
        ${EMLIBS}/em_mpfr/lib/libmpfr.a 
        #${EMLIBS}/em_ceres/lib/libceres.a
        ${EMLIBS}/em_bezier/lib/libbezier.a
        )
