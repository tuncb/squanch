project(squanch-tester CXX)
cmake_minimum_required(VERSION 2.8.12)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

file(GLOB_RECURSE SRC_TESTING
    "./src/*.cpp"
    "./src/*.h"
)
file(GLOB_RECURSE SRC_SQUANCH
    "../include/squanch/*.cpp"
    "../include/squanch/*.h"
)

add_executable(squanch-tester ${SRC_TESTING} ${SRC_SQUANCH})

target_link_libraries(squanch-tester ${CONAN_LIBS})
target_include_directories(squanch-tester PRIVATE "../include/squanch/")

source_group("squanch" FILES ${SRC_SQUANCH})

