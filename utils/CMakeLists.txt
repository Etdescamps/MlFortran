
add_executable(mlf_optim mlf_optim.cpp tools_dlopen.cpp)
target_link_libraries(mlf_optim ${MLF_LIB} ${CMAKE_DL_LIBS})
install(TARGETS mlf_optim DESTINATION bin)

