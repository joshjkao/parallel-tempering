cmake_minimum_required(VERSION 3.23) 
project(parallel_tempering)
set(CMAKE_CXX_STANDARD 17)

#include_directories(partem)

# export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
# export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

# Add main.cpp file of the project root directory as a source file
set(SOURCE_FILES 
    main.cpp
    src/PT.cpp
    src/Replica.cpp
    src/RNG.cpp
    src/util.cpp
)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(run ${SOURCE_FILES})

target_include_directories(run
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)