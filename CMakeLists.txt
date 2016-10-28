cmake_minimum_required(VERSION 2.8.5)

project(em)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include_directories(
	${CMAKE_SOURCE_DIR}/include
)

file(GLOB libs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*")

foreach(dir ${libs})
	message(STATUS "Adding subdirectory: ${dir}")
	add_subdirectory("${CMAKE_SOURCE_DIR}/${dir}")
endforeach()
