
add_executable(test_rsc test_rsc.f90 test_rsc_main.c)
add_executable(test_kmeans test_kmeans.f90)
add_executable(test_polynomial test_polynomial.f90)
add_executable(test_emgmm test_common.f90 test_emgmm.f90)
add_executable(test_kmeans_c test_kmeans_c.c)
add_executable(test_pareto test_common.f90 test_pareto.f90)
add_executable(test_funbasis test_functions.f90 test_funbasis.f90)
add_executable(test_funbasis_c test_funbasis_c.c)
add_executable(test_cmaes test_functions.f90 test_cmaes.f90)
add_executable(test_ode45 test_functions.f90 test_ode45.f90)
add_executable(test_ode45_dense test_functions.f90 test_ode45_dense.f90)
add_executable(test_inverse_cdf test_inverse_cdf.f90)
add_executable(test_kmc test_kmc.f90 test_kmc_model.f90)
add_executable(test_hybrid_kmc test_hybrid_kmc.f90 test_hybrid_kmc_model.f90)
add_executable(test_cstr_kmc test_cstr_kmc.f90 test_cstr_kmc_model.f90)
add_executable(test_sobol test_sobol.f90)
add_executable(test_linalg test_linalg.f90)
add_executable(test_beta test_beta.f90)
add_executable(test_gamma test_gamma.f90)
add_executable(test_utils test_utils.f90)
add_executable(test_2dplane test_2dplane.f90 test_functions.f90)
add_executable(test_sumtree test_sumtree.f90)

target_link_libraries(test_rsc ${MLF_LIB})
target_link_libraries(test_kmeans ${MLF_LIB})
target_link_libraries(test_polynomial ${MLF_LIB})
target_link_libraries(test_emgmm ${MLF_LIB})
target_link_libraries(test_kmeans_c ${MLF_LIB})
target_link_libraries(test_pareto ${MLF_LIB})
target_link_libraries(test_funbasis ${MLF_LIB})
#target_link_libraries(test_funbasis ${MLF_LIB}_s ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(test_funbasis_c ${MLF_LIB})
target_link_libraries(test_cmaes ${MLF_LIB})
target_link_libraries(test_ode45 ${MLF_LIB})
target_link_libraries(test_ode45_dense ${MLF_LIB})
target_link_libraries(test_inverse_cdf ${MLF_LIB})
target_link_libraries(test_kmc ${MLF_LIB})
target_link_libraries(test_hybrid_kmc ${MLF_LIB})
target_link_libraries(test_cstr_kmc ${MLF_LIB})
target_link_libraries(test_sobol ${MLF_LIB})
target_link_libraries(test_linalg ${MLF_LIB})
target_link_libraries(test_beta ${MLF_LIB})
target_link_libraries(test_gamma ${MLF_LIB})
target_link_libraries(test_2dplane ${MLF_LIB})
target_link_libraries(test_utils ${MLF_LIB})
target_link_libraries(test_sumtree ${MLF_LIB})

