cmake_minimum_required(VERSION 2.8.5)

set(target_name ${PROJECT_NAME})

set(link_libraries
	m
	matrix
)

find_package(CUnit)
if(CUNIT_FOUND)
	include_directories(
		${CUNIT_INCLUDE_DIRS}
	)
	set(link_libraries
		${link_libraries}
		${CUNIT_LIBRARIES}
	)
endif()

find_package(GMP)
if(GMP_FOUND)
	include_directories(
		${GMP_INCLUDE_DIR}
	)
	set(link_libraries
		${link_libraries}
		${GMP_LIBRARIES}
	)
endif()

file(GLOB src_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.c")

add_library(
	${target_name} SHARED
	${src_files}
)

target_link_libraries(
	${target_name} -Wl,-E
	${link_libraries}
)

if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
	file(GLOB test_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "test/*.c")

	set(test_name ${target_name}_test)

	add_executable(
		${test_name}
		${test_files}
	)

	target_link_libraries(
		${test_name} -Wl,-E
		${target_name}
	)
endif()
