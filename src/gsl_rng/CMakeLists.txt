find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})
add_executable(gsl_rng gsl_rng.cpp)
target_link_libraries(gsl_rng ${LIBS})
