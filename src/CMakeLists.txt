include_directories(${2decomp_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

set(X3DFILES
    bc_tgv2d.f90
    bc_tgv.f90
    case.f90
    CMakeLists.txt
    module_param.f90
    navier.f90
    parameters.f90
    poisson.f90
    thomas.f90
    time_integrators.f90
    tools.f90
    transeq.f90
    var.f90
    x3d_derive.f90
    x3d_filters.f90
    x3d_operator_1d.f90
    x3d_operator_x_data.f90
    x3d_operator_y_data.f90
    x3d_operator_z_data.f90
    x3d_precision.f90
    x3d_staggered.f90
    x3d_tools.f90
    x3d_transpose.f90)
add_library(x3d STATIC ${X3DFILES})
#target_include_directories(x3d PRIVATE ${2decomp_INCLUDE_DIR})
#target_link_libraries(x3d PRIVATE decomp2d)

add_executable(xcompact3d xcompact3d.f90)
target_include_directories(xcompact3d PRIVATE ${2decomp_INCLUDE_DIR})
target_link_libraries(x3d PRIVATE ${DECOMP2D_FFT_LIB})
target_link_libraries(xcompact3d PRIVATE x3d)
target_link_libraries(xcompact3d PRIVATE ${DECOMP2D_FFT_LIB})

if (MPI_FOUND)
  target_link_libraries(x3d PRIVATE MPI::MPI_Fortran)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if (FFTW_FOUND)
  target_link_libraries(fft_c2c_x PRIVATE ${FFTW_DOUBLE_LIB})
  target_link_libraries(fft_c2c_x PRIVATE ${FFTW_FLOAT_LIB})
  target_link_libraries(fft_r2c_x PRIVATE ${FFTW_DOUBLE_LIB})
  target_link_libraries(fft_r2c_x PRIVATE ${FFTW_FLOAT_LIB})
  target_link_libraries(fft_grid_x PRIVATE ${FFTW_DOUBLE_LIB})
  target_link_libraries(fft_grid_x PRIVATE ${FFTW_FLOAT_LIB})
elseif (MKL_FOUND)
  target_link_libraries(fft_c2c_x PUBLIC $<LINK_ONLY:MKL::MKL>)
  target_link_libraries(fft_r2c_x PUBLIC $<LINK_ONLY:MKL::MKL>)
  target_link_libraries(fft_grid_x PUBLIC $<LINK_ONLY:MKL::MKL>)
endif (FFTW_FOUND)

install(TARGETS xcompact3d x3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

