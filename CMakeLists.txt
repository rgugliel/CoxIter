cmake_minimum_required(VERSION 3.12)

# Taken from https://blog.kitware.com/cmake-and-the-default-build-type/
# Set a default build type if none was specified.
# This must precede the project() line, which would set the CMAKE_BUILD_TYPE
# to 'Debug' with single-config generators on Windows.
set(default_build_type "Release")
get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT isMultiConfig AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

project(coxiter)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT MSVC)
	add_compile_options(-Wall -Wno-reorder -Wno-unknown-pragmas)
endif()

# ------------------------------------------------------------------
# External libraries
find_package( OpenMP )

find_library( PCRE_LIBRARY pcre REQUIRED )
find_path( PCRE_INCLUDE_DIR pcre.h REQUIRED )

find_library( GMP_LIBRARY gmp REQUIRED )

find_library( GMPXX_LIBRARY gmpxx REQUIRED )
find_path( GMPXX_INCLUDE_DIR gmpxx.h REQUIRED )

find_library( PARI_LIBRARY pari )
find_path( PARI_INCLUDE_DIR pari/pari.h )

# ------------------------------------------------------------------
# Set up targets

# Maths: polynomials, fractions, varia functions
add_library( coxiter_maths STATIC lib/math_tools.cpp lib/polynomials.cpp lib/numbers/number_template.cpp lib/numbers/mpz_rational.cpp )
target_link_libraries( coxiter_maths PUBLIC ${GMP_LIBRARY} ${GMPXX_LIBRARY} )
target_include_directories( coxiter_maths PUBLIC ${GMPXX_INCLUDE_DIR} )

# Graphs: graphs, products of graphs and iterators
add_library( coxiter_graphs STATIC graphs.product.set.cpp graphs.product.cpp graphs.list.n.cpp graphs.list.iterator.cpp graphs.list.cpp graph.cpp  )

# Main files
add_executable(coxiter index2.cpp lib/string.cpp lib/regexp.cpp coxiter.cpp arithmeticity.cpp app.cpp main.cpp)
target_link_libraries( coxiter PUBLIC ${PCRE_LIBRARY} )
target_include_directories( coxiter PUBLIC ${PCRE_INCLUDE_DIR} )
target_link_libraries( coxiter PRIVATE coxiter_maths coxiter_graphs )

if( OpenMP_FOUND )
	target_link_libraries( coxiter PRIVATE OpenMP::OpenMP_CXX )
endif()

if( PARI_LIBRARY )
	add_library( coxiter_pari STATIC lib/paripolynomials.cpp growthrate.cpp signature.cpp )
	target_link_libraries( coxiter_pari PUBLIC ${GMP_LIBRARY} ${GMPXX_LIBRARY} ${PARI_LIBRARY} )
	target_include_directories( coxiter_pari PUBLIC ${GMPXX_INCLUDE_DIR} ${PARI_INCLUDE_DIR} )
	
	target_link_libraries( coxiter PRIVATE coxiter_pari )
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
