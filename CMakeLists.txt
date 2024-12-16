cmake_minimum_required(VERSION 3.22)
project(projectC C)

set(CMAKE_C_STANDARD 11)


add_executable(projectC
        main.c
        includes/vec.h
        lib/vec.c
        lib/arr.c
        includes/arr.h
        includes/phisics.h
        lib/until.c
        includes/until.h
)


target_link_libraries(projectC m)