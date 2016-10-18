cmake_minimum_required(VERSION 2.8.5)

find_package(OpenMP)
if (OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(target_name "${PROJECT_NAME}app")

file(GLOB src_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.c")

include_directories(
	${CMAKE_SOURCE_DIR}/lib/include
)

add_executable(
	${target_name}
	${src_files}
)

set(link_libraries
	em
	gomp
	gmp
)

target_link_libraries(
	${target_name} -Wl,-E
	${link_libraries}
)

