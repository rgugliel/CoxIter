cmake_minimum_required(VERSION 2.4)

project(temp)

add_executable(temp temp.cpp ../../../../tools/string.cpp ../../../../Vino/rationalinteger_invariantsqf.cpp ../../../../tools/regexp.cpp ../../gbd.cpp  ../../coxiter.cpp ../../arithmeticity.cpp main.cpp)

# Maths: polynomials, fractions, varia functions
add_library( coxiter_maths STATIC ../../../tools/maths_tools.cpp ../../../tools/polynomials.cpp ../../../tools/numbers/number_template.cpp ../../../tools/numbers/number_rational.cpp )
target_link_libraries( temp coxiter_maths )

# Graphs: graphs, products of graphs and iterators
add_library( coxiter_graphs STATIC ../../graphs.product.set.cpp ../../graphs.product.cpp ../../graphs.list.n.cpp ../../graphs.list.iterator.cpp ../../graphs.list.cpp ../../graph.cpp  )
target_link_libraries( temp coxiter_graphs )

add_library( vino_numbers STATIC ../../../Vino/quadraticinteger.cpp ../../../Vino/rationalinteger.cpp ../../../Vino/algebraicinteger.cpp )
target_link_libraries( temp vino_numbers )

add_definitions("-fopenmp")
add_definitions("-lomp")
add_definitions("-std=c++11")
add_definitions("-Wall")
add_definitions("-Wno-reorder")

set( CMAKE_BUILD_TYPE Release )

# ------------------------------------------------------------------
# External libraries
find_package( OpenMP )
if( OPENMP_FOUND )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_library( PCRE_LIBRARY pcre )
target_link_libraries( temp ${PCRE_LIBRARY} )

find_library( GMP_LIBRARY gmp )
target_link_libraries( temp ${GMP_LIBRARY} )

find_library( GMP_LIBRARYXX gmpxx )
target_link_libraries( temp ${GMP_LIBRARYXX} )

if( COMPILE_WITH_PARI ) # Used for the growth rate and the Perron/Salem/Pisot tests
	find_library( PARI_LIBRARY pari )
	if( PARI_LIBRARY )
		add_library( coxiter_growthrate STATIC ../../../tools/paripolynomials.cpp ../../growthrate.cpp )
		target_link_libraries( temp coxiter_growthrate )
		
		target_link_libraries( temp ${PARI_LIBRARY} )
		add_definitions(-D_COMPILE_WITH_PARI_)
	else()
		message( WARNING "Warning: libpari not found. Some functions won't be avaible." )
	endif()
endif( COMPILE_WITH_PARI )

install(TARGETS temp RUNTIME DESTINATION bin)
