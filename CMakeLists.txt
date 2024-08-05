cmake_minimum_required(VERSION 3.20)
project(fch)


# CONFIG ================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
        message("RELEASE MODE: ENABLING ALL OPTIMIZATION FLAGS.")
        set(OPTIM_FLAGS
                "-Ofast"
                "-march=native"
                "-mtune=native"
                #"-mfpmath=sse"
        )
else ()
        set(OPTIM_FLAGS "")
endif ()

if(TEST_CASE)
        message("TEST CASE SELECTED: ${TEST_CASE}.")
else ()
        message("TEST CASE NOT SELECTED. USING DEFAULT.")
endif ()


# DEPENDENCIES ==========================================
add_subdirectory(eigen)
add_subdirectory(argparse)

set(DEPENDENCIES
        Eigen3::Eigen
        argparse::argparse
)


# COMMON ================================================
set(COMMON_SRC
        common/bmi2.h
        common/bmi2.cpp
        common/utility.h
        common/utility.cpp
        common/definitions.h
        common/definitions.cpp
        common/test_cases/test_case_template.h
        common/test_cases/all_cases.h
        common/test_cases/klein_5to3.h
        common/test_cases/hypersphere_5to1.h
        common/test_cases/circle_9to8.h
)


# PTA ===================================================
set(PTA_SRC
        pta/permutahedron.h
        pta/permutahedron.cpp
)
#---
add_executable(pta
        ${COMMON_SRC}
        ${PTA_SRC}
        pta/main_pta.cpp
)
target_link_libraries(pta
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(pta
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(pta 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()
#---
add_executable(skeleton_pta
        ${COMMON_SRC}
        ${PTA_SRC}
        pta/main_skeleton_pta.cpp
)
target_link_libraries(skeleton_pta
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(skeleton_pta
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(skeleton_pta 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()


# GCCH ==================================================
set(GCCH_SRC
        gcch/gch.h
        gcch/gch.cpp
)
#---
add_executable(gcch
        ${COMMON_SRC}
        ${GCCH_SRC}
        gcch/main_gcch.cpp
)
target_link_libraries(gcch
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(gcch
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(gcch 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()
#---
add_executable(reorder_hypercubes
        ${COMMON_SRC}
        ${GCCH_SRC}
        gcch/main_reorder_hypercubes.cpp
)
target_link_libraries(reorder_hypercubes
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(reorder_hypercubes
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(reorder_hypercubes 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()
#---
add_executable(skeleton_hypercubes
        ${COMMON_SRC}
        ${GCCH_SRC}
        gcch/main_skeleton_hypercubes.cpp
)
target_link_libraries(skeleton_hypercubes
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(skeleton_hypercubes
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(skeleton_hypercubes 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()


# FCH ===================================================
set(FCH_SRC
        fch/fch.h
        fch/fch.cpp
)
#---
add_executable(fch
        ${COMMON_SRC}
        ${GCCH_SRC}
        ${FCH_SRC}
        fch/main_fch.cpp
)
target_link_libraries(fch
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(fch
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(fch 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()
#---
add_executable(fch_to_hypercubes
        ${COMMON_SRC}
        ${GCCH_SRC}
        ${FCH_SRC}
        fch/main_fch_to_hypercubes.cpp
)
target_link_libraries(fch_to_hypercubes
        PRIVATE
        ${DEPENDENCIES}
)
target_compile_options(fch_to_hypercubes
        PRIVATE
        ${OPTIM_FLAGS}
)
if(TEST_CASE)
        target_compile_definitions(fch_to_hypercubes 
                PRIVATE 
                -DSELECTED_CASE=${TEST_CASE}
        )
endif ()
