cmake_minimum_required(VERSION 3.22)
project(complex_predict)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include_directories(./include)

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB SOURCES ${SOURCE_DIR}/*.cpp) # This line is not needed if you list your sources explicitly

add_executable(randomwalk
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/randomwalk.cpp
)

add_executable(randomwalk2
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/randomwalk2.cpp
)

add_executable(random_core
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/random_core.cpp
)

add_executable(new_delete
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        src/new_delete.cpp
)
