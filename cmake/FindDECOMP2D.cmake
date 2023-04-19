# - Find the 2decomp-fft library

if( NOT D2D_ROOT AND DEFINED ENV{D2D_DIR} )
    set( D2D_ROOT $ENV{D2D_DIR} )
endif()
message(STATUS "D2D_ROOT ${D2D_ROOT}")

# find libs
find_library(
    DECOMP2D_FFT_LIB
    NAMES "decomp2d" libdecomp2d
    PATHS ${D2D_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
    REQUIRED
    )
set(2decomp_INCLUDE_DIR "${D2D_ROOT}/include")
message(STATUS "D2D found at: ${D2D_ROOT}")


