cmake_minimum_required(VERSION 3.28)
project(projectC C)

set(CMAKE_C_STANDARD 11)


add_executable(projectC
        main.c
        includes/vec.h
        lib/vec.c
        lib/arr.c
        includes/arr.h
        includes/phisics.h
)


target_link_libraries(projectC m)