
set(MLF_UNSUPERVISED unsupervised/mlf_kmeans.f90 unsupervised/mlf_emgmm.f90 unsupervised/mlf_gaussian.f90)
set(MLF_CSRC cfuns/qsort.c)
set(MLF_CPPINTF cpp_intf/mlfcpp.cpp cpp_intf/mlf_dlopen.cpp)
set(MLF_REDUCE reduce/mlf_funbasis.f90 reduce/mlf_gridmodel.f90)
set(MLF_MULTI_OBJ multi_objective/mlf_pareto.f90)
set(MLF_MONO_OBJ mono_objective/mlf_cmaes.f90)
set(MLF_MAIN_SRC mlf_cfuns.f90 mlf_errors.f90 mlf_intf.f90 mlf_rsc_array.f90 mlf_fun_intf.f90 mlf_step_algo.f90 mlf_models.f90 mlf_hdf5.f90 mlf_optim.f90 mlf_utils.f90 mlf_rand.f90 mlf_matrix.f90 mlf_optim_intf.f90)
set(MLF_SRC ${MLF_CSRC} ${MLF_MAIN_SRC} ${MLF_UNSUPERVISED} ${MLF_REDUCE} ${MLF_MONO_OBJ} ${MLF_MULTI_OBJ} ${MLF_CPPINTF})

set(CMAKE_Fortran_MODULE_DIRECTORY ${MLF_INCLUDES})
add_library(${MLF_LIB} SHARED ${MLF_SRC})
#add_library(${MLF_LIB}_s STATIC ${MLF_SRC})
target_link_libraries(${MLF_LIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

install(TARGETS ${MLF_LIB} DESTINATION lib)
install(DIRECTORY ${MLF_INCLUDES} DESTINATION include)


