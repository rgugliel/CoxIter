project(multitools)

add_executable(multitools arithmeticity.number.cpp multitools.cpp main.cpp bruteforce.cpp ../../../tools/numbers/number_template.cpp ../../../tools/numbers/number_rational.cpp  ../../../tools/string.cpp ../../../tools/regexp.cpp ../../../tools/polynomials.cpp ../../gbd.cpp ../../graphs.product.set.cpp ../../graphs.product.cpp ../../graphs.list.n.cpp ../../graphs.list.iterator.cpp ../../graphs.list.cpp ../../graph.cpp ../../coxiter.cpp ../../arithmeticity.cpp)

add_definitions("-fopenmp")
add_definitions("-lomp")
add_definitions("-std=c++11")
add_definitions("-Wall")
add_definitions("-Wno-reorder")

set( CMAKE_BUILD_TYPE Release )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_library(PCRE_LIBRARY pcre)
target_link_libraries(multitools ${PCRE_LIBRARY})

find_library(GMP_LIBRARY gmp)
target_link_libraries(multitools ${GMP_LIBRARY})

find_library(GMP_LIBRARYXX gmpxx)
target_link_libraries(multitools ${GMP_LIBRARYXX})

install(TARGETS multitools RUNTIME DESTINATION bin)
