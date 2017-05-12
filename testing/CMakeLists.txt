cmake_minimum_required(VERSION 2.6)

project(tests)

add_definitions("-fopenmp")
add_definitions("-lomp")
add_definitions("-std=c++11")
add_definitions("-Wall")
add_definitions("-Wno-reorder")
add_definitions("-msse3")

set(COXITER_PATH ../../)

add_executable(tests tests.cpp ${COXITER_PATH}/lib/numbers/number_template.cpp ${COXITER_PATH}/lib/numbers/mpz_rational.cpp ${COXITER_PATH}/lib/string.cpp ${COXITER_PATH}/lib/regexp.cpp  ${COXITER_PATH}/lib/polynomials.cpp ${COXITER_PATH}/graphs.product.set.cpp ${COXITER_PATH}/graphs.product.cpp ${COXITER_PATH}/graphs.list.n.cpp ${COXITER_PATH}/graphs.list.iterator.cpp ${COXITER_PATH}/graphs.list.cpp ${COXITER_PATH}/graph.cpp ${COXITER_PATH}/coxiter.cpp ${COXITER_PATH}/arithmeticity.cpp main.cpp)

set( CMAKE_BUILD_TYPE Release )

find_package( OpenMP )
if( OPENMP_FOUND )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	
	add_definitions("-fopenmp")
	add_definitions("-lomp")
	add_definitions(-D_COMPILE_WITH_OPENMP_)
endif()

find_library(PCRE_LIBRARY pcre)
target_link_libraries(tests ${PCRE_LIBRARY})

find_library(GMP_LIBRARY gmp)
target_link_libraries(tests ${GMP_LIBRARY})

find_library(GMP_LIBRARYXX gmpxx)
target_link_libraries(tests ${GMP_LIBRARYXX})

find_library( PARI_LIBRARY pari )
add_library( coxiter_growthrate STATIC ${COXITER_PATH}/lib/paripolynomials.cpp ${COXITER_PATH}/growthrate.cpp ${COXITER_PATH}/signature.cpp )
target_link_libraries( tests coxiter_growthrate )
target_link_libraries( tests ${PARI_LIBRARY} )

install(TARGETS tests RUNTIME DESTINATION bin)
