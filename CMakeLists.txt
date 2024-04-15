cmake_minimum_required(VERSION 3.22)
project(complex_predict)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include_directories(./include)

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB SOURCES ${SOURCE_DIR}/*.cpp) # This line is not needed if you list your sources explicitly

# List your source files explicitly here
add_executable(bops
        ${SOURCE_DIR}/lib/graph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/bops.cpp
        src/lib/main.cpp
)

add_executable(cns
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/cns.cpp
)

add_executable(new_cns
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/new_cns.cpp
)

add_executable(max_matched
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/max_matched.cpp
)

add_executable(randomwalk
        ${SOURCE_DIR}/lib/ungraph.cpp
        ${SOURCE_DIR}/lib/bio_information.cpp
        ${SOURCE_DIR}/lib/dag.cpp
        ${SOURCE_DIR}/randomwalk.cpp
)