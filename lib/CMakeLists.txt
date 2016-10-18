cmake_minimum_required(VERSION 2.8.5)

set(target_name ${PROJECT_NAME})

file(GLOB src_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.c")

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(
	${target_name} SHARED
	${src_files}
)

set(link_libraries
)

target_link_libraries(
	${target_name} -Wl,-E
	${link_libraries}
)

