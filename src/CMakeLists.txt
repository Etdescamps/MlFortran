
set(MLF_UNSUPERVISED unsupervised/mlf_kmeans.f90 unsupervised/mlf_kmeans_naive.f90 unsupervised/mlf_emgmm.f90 unsupervised/mlf_gaussian.f90)
set(MLF_CSRC cfuns/qsort.c)
set(MLF_STOCHASTIC stochastic/mlf_kmc.f90 stochastic/mlf_hybrid_kmc.f90)
set(MLF_CPPINTF cpp_intf/mlfcpp.cpp cpp_intf/mlf_dlopen.cpp)
set(MLF_REDUCE reduce/mlf_funbasis.f90 reduce/mlf_gridmodel.f90)
set(MLF_MULTI_OBJ multi_objective/mlf_pareto.f90)
set(MLF_MONO_OBJ mono_objective/mlf_cmaes.f90)
set(MLF_ODE_OBJ ode/mlf_ode_class.f90 ode/mlf_ode45.f90)
set(MLF_MODEL_OBJ mlf_models.f90 mlf_supervised_model.f90)
set(MLF_SAMPLERS_OBJ samplers/mlf_inverse_cdf.f90 statistics/mlf_histogram.f90 samplers/mlf_radius_sampler.f90 samplers/mlf_sobol1111.f90 samplers/mlf_sobol1111_table.f90 statistics/mlf_gamma_dist.f90 statistics/mlf_beta_dist.f90)
set(MLF_MAIN_SRC mlf_cfuns.f90 mlf_errors.f90 mlf_intf.f90 mlf_rsc_array.f90 mlf_fun_intf.f90 mlf_step_algo.f90 mlf_hdf5.f90 mlf_optim.f90 mlf_utils.f90 mlf_rand.f90 mlf_matrix.f90 mlf_optim_intf.f90 mlf_poly.f90 mlf_linalg.f90)
set(MLF_SRC ${MLF_CSRC} ${MLF_MAIN_SRC} ${MLF_UNSUPERVISED} ${MLF_REDUCE} ${MLF_MONO_OBJ} ${MLF_MULTI_OBJ} ${MLF_CPPINTF} ${MLF_ODE_OBJ} ${MLF_STOCHASTIC} ${MLF_MODEL_OBJ} ${MLF_SAMPLERS_OBJ})

set(CMAKE_Fortran_MODULE_DIRECTORY ${MLF_INCLUDES})
add_library(${MLF_LIB} SHARED ${MLF_SRC})
#add_library(${MLF_LIB}_s STATIC ${MLF_SRC})
target_link_libraries(${MLF_LIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

install(TARGETS ${MLF_LIB} DESTINATION lib)
install(DIRECTORY ${MLF_INCLUDES} DESTINATION include)


