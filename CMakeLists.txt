cmake_minimum_required(VERSION 2.4)

project(coxiter)

add_definitions("-std=c++11")
if( CYGWIN )
	add_definitions("-std=gnu++11")
endif()

add_definitions("-Wall")
add_definitions("-Wno-reorder")
add_definitions("-Wno-unknown-pragmas")

set( CMAKE_BUILD_TYPE Release )

# Main files
add_executable(coxiter index2.cpp lib/string.cpp lib/regexp.cpp coxiter.cpp arithmeticity.cpp app.cpp main.cpp)

# Maths: polynomials, fractions, varia functions
add_library( coxiter_maths STATIC lib/math_tools.cpp lib/polynomials.cpp lib/numbers/number_template.cpp lib/numbers/mpz_rational.cpp )
target_link_libraries( coxiter coxiter_maths )

# Graphs: graphs, products of graphs and iterators
add_library( coxiter_graphs STATIC graphs.product.set.cpp graphs.product.cpp graphs.list.n.cpp graphs.list.iterator.cpp graphs.list.cpp graph.cpp  )
target_link_libraries( coxiter coxiter_graphs )

# ------------------------------------------------------------------
# External libraries
find_package( OpenMP )
if( OPENMP_FOUND )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	
	add_definitions("-fopenmp")
	add_definitions("-lomp")
endif()

find_library( PCRE_LIBRARY pcre )
target_link_libraries( coxiter ${PCRE_LIBRARY} )

find_library( GMP_LIBRARY gmp )
target_link_libraries( coxiter ${GMP_LIBRARY} )
# find <gmpxx.h> in APPLE
if( APPLE )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/include" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include" )
endif()

find_library( GMP_LIBRARYXX gmpxx )
target_link_libraries( coxiter ${GMP_LIBRARYXX} )

find_library( PARI_LIBRARY pari )
if( PARI_LIBRARY )
	add_library( coxiter_pari STATIC lib/paripolynomials.cpp growthrate.cpp signature.cpp )
	target_link_libraries( coxiter_pari ${GMP_LIBRARY} )
	target_link_libraries( coxiter_pari ${GMP_LIBRARYXX} )
	target_link_libraries( coxiter coxiter_pari )
	
	target_link_libraries( coxiter ${PARI_LIBRARY} )
	add_definitions(-D_COMPILE_WITH_PARI_)
else()
	message( WARNING "Warning: PARI library was not found. CoxIter won't be able to compute growth rate and signature." )
endif()

# ------------------------------------------------------------------
# If the dot program is present, we will be able to create images
find_program( PROGRAM_DOT dot )
if( PROGRAM_DOT )
	add_definitions(-D_DOT_PROGRAM_FOUND_)
endif()

install(TARGETS coxiter RUNTIME DESTINATION bin)
