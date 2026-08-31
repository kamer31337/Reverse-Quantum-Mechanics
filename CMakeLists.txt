cmake_minimum_required(VERSION 3.15)
project(ReverseQuantumMechanics C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(include)

set(SOURCES
    src/main.c
    src/rqm_math.c
    src/rqm_classical.c
    src/rqm_hilbert_pathology.c
    src/rqm_ensembles.c
    src/rqm_born_entropy.c
    src/rqm_unitary.c
    src/rqm_equilibration.c
    src/rqm_classical_limit.c
    src/rqm_triple_equilibria.c
    src/rqm_ui.c
    src/rqm_tests.c
)

add_executable(rqm ${SOURCES})

if(NOT MSVC)
    target_link_libraries(rqm m)
endif()

if(MSVC)
    target_compile_options(rqm PRIVATE /W4 /O2)
else()
    target_compile_options(rqm PRIVATE -Wall -Wextra -O3)
endif()

enable_testing()
add_test(NAME rqm_all_theorems COMMAND rqm --test)
