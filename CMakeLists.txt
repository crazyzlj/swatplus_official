# SWAT+ cmake, req min v3.22 (ifx, etc. support)
###############################################################
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" OR
   "${CMAKE_SOURCE_DIR}/src" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "'in-source' builds not permitted, e.g. use 'cmake -B build' instead.")
endif()

cmake_minimum_required(VERSION 3.22)

# Query the host system information
cmake_host_system_information(RESULT os_id QUERY DISTRIB_ID)
cmake_host_system_information(RESULT os_version QUERY DISTRIB_RELEASE)
cmake_host_system_information(RESULT os_like QUERY DISTRIB_ID_LIKE)
cmake_host_system_information(RESULT pretty_name QUERY DISTRIB_PRETTY_NAME)

project(swatplus VERSION 0.2 LANGUAGES Fortran)

include(CTest)

# Resolve SWAT+ version used in executable names and generated source files.
# Priority:
#   1. Git tag information, when available.
#   2. A manually maintained VERSION file in the project root.
#   3. The short hash of the latest Git commit.
# A user-provided -DTAG=... still overrides this auto-detection block.
if(NOT TAG)
    set(SWAT_VERSION_SOURCE "unknown")
    find_package(Git QUIET)

    # 1) Prefer Git tag information when the repository has tags.
    if(GIT_FOUND)
        execute_process(
          COMMAND ${GIT_EXECUTABLE} describe --tags
          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          OUTPUT_VARIABLE TAG
          RESULT_VARIABLE TAG_RESULT
          OUTPUT_STRIP_TRAILING_WHITESPACE
          ERROR_QUIET
        )

        if(TAG_RESULT EQUAL 0 AND NOT TAG STREQUAL "")
            set(SWAT_VERSION_SOURCE "git tag")
        else()
            unset(TAG)
        endif()
    endif()

    # 2) If no tag is available, read the first non-empty line from VERSION.
    if(NOT TAG)
        set(VERSION_FILE "${PROJECT_SOURCE_DIR}/VERSION")
        if(EXISTS "${VERSION_FILE}")
            file(STRINGS "${VERSION_FILE}" VERSION_LINES LIMIT_COUNT 1)
            if(VERSION_LINES)
                list(GET VERSION_LINES 0 TAG)
                string(STRIP "${TAG}" TAG)
            endif()

            if(TAG)
                set(SWAT_VERSION_SOURCE "VERSION file")
            else()
                unset(TAG)
            endif()
        endif()
    endif()

    # 3) If neither tags nor VERSION exist, fall back to the latest commit short hash.
    if(NOT TAG AND GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE TAG
            RESULT_VARIABLE HASH_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )

        if(HASH_RESULT EQUAL 0 AND TAG)
            set(SWAT_VERSION_SOURCE "git commit")
        else()
            unset(TAG)
        endif()
    endif()

    if(NOT TAG)
        set(TAG "unknown")
        set(SWAT_VERSION_SOURCE "unknown")
    endif()
else()
    set(SWAT_VERSION_SOURCE "manual TAG")
endif()

# user-configurable optimization level for Release.
# default to O1, can be overridden using -DOPT_LEVEL=O2
set(OPT_LEVEL "O1" CACHE STRING "Optimization level for non-Debug configs: O0, O1, O2, O3")
set_property(CACHE OPT_LEVEL PROPERTY STRINGS O0 O1 O2 O3)
string(TOLOWER "${OPT_LEVEL}" OPT_LEVEL_TAG)

# Normalize Linux distribution strings for robust matching.
string(TOLOWER "${os_id}" os_id_lc)
string(TOLOWER "${os_like}" os_like_lc)

# define compiler/link option lists
set(F_COMMON "")
set(F_DEBUG "")
set(F_NONDEBUG "")
set(F_LINK "")
set(FFC "unknown")

# cross-platform / cross-compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set(FFC "ifx")
    else()
        set(FFC "ifo")
    endif()

    if(WIN32)
        list(APPEND F_COMMON /fpp /free /fpe0 /traceback)
        if(FFC STREQUAL "ifo")
            # suppress noisy diagnostics/remarks.
            list(APPEND F_COMMON /Qdiag-disable:all /Qdiag-disable:remarks)
        endif()
        list(APPEND F_DEBUG  /warn:all /check:bounds /Od)
        list(APPEND F_NONDEBUG /${OPT_LEVEL})
        list(APPEND F_LINK -static)
    else()
        list(APPEND F_COMMON -fpp -free -fpe0 -traceback)
        list(APPEND F_DEBUG  "SHELL:-warn all" "SHELL:-check bounds" -O0)
        list(APPEND F_NONDEBUG -${OPT_LEVEL})
        list(APPEND F_LINK -static)

        if(FFC STREQUAL "ifo")
            list(APPEND F_COMMON -diag-disable=10448)
        endif()
    endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(FFC "gnu")
    # keep only dialect and semantic/numeric-stability flags in all configurations;
    # put runtime diagnostic traps in Debug only. This avoids Release crashes caused by
    # benign FP flags raised inside optimized system math routines.
    list(APPEND F_COMMON
        -cpp
        -fimplicit-none
        -ffree-line-length-none
        -fbacktrace
        -finit-local-zero
        -fno-unsafe-math-optimizations
        -frounding-math
    )

    list(APPEND F_DEBUG
        -Wall
        -fcheck=all
        -ffpe-trap=invalid,zero,overflow,underflow
        -fsignaling-nans
        -O0
    )

    list(APPEND F_NONDEBUG
        -${OPT_LEVEL}
    )

    if(APPLE)
        # macOS cannot fully statically link libSystem. Static-link the GNU Fortran
        # runtime libraries to avoid depending on Homebrew gcc dylibs where possible.
        list(APPEND F_LINK -static-libgfortran -static-libquadmath -static-libgcc)
    elseif(WIN32 OR (UNIX AND NOT ("${os_id_lc}" STREQUAL "arch" OR
                                   "${os_like_lc}" MATCHES "arch" OR
                                   "${os_id_lc}" STREQUAL "fedora" OR
                                   "${os_like_lc}" MATCHES "fedora")))
        # Avoid fully static GNU builds on Arch/Fedora-like Linux systems, where
        # static runtime/library availability is commonly incomplete or discouraged.
        list(APPEND F_LINK -static)
    endif()
endif()

# To avoid MSVC automatically add the default /O2 to override the user-specific OPT_LEVEL
if(CMAKE_GENERATOR MATCHES "Visual Studio" AND CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE "" CACHE STRING "Release flags" FORCE)
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "" CACHE STRING "RelWithDebInfo flags" FORCE)
    set(CMAKE_Fortran_FLAGS_MINSIZEREL "" CACHE STRING "MinSizeRel flags" FORCE)
endif()

# Get linker properties for later output.
get_property(CURRENT_LINK_LIBS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_LIBRARIES)


string(TOLOWER "${CMAKE_HOST_SYSTEM_PROCESSOR}" ARCH)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(AR "lin_${ARCH}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(AR "win_${ARCH}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(AR "mac_${ARCH}")
else()
    set(AR "unknown")
endif()

# SWAT Version number
set(SWAT_VERSION ${TAG})
set(SWATPLUS_EXE_BASE "swatplus-${SWAT_VERSION}-${FFC}-${AR}")


# Set this to 'TRUE' to see the fortran command on compile
set(CMAKE_VERBOSE_MAKEFILE FALSE)

#############################################################################
# Build
if(EXISTS "${PROJECT_SOURCE_DIR}/src/main.f90.in")
    string(TIMESTAMP ISO   "%Y-%m-%d %H:%M:%S")
    string(TIMESTAMP TODAY "%Y-%m-%d")  # e.g. produces 2023-12-07
    string(TIMESTAMP YEAR  "%Y")

    configure_file(
        "${PROJECT_SOURCE_DIR}/src/main.f90.in"
        "${PROJECT_SOURCE_DIR}/src/main.f90"
    )
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/ford.md.in")
    configure_file(
        "${PROJECT_SOURCE_DIR}/ford.md.in"
        "${PROJECT_SOURCE_DIR}/ford.md"
    )   
endif()

### Set default or specified installation directory.
IF (INSTALL_PREFIX)
    SET(INSTALL_DIR ${INSTALL_PREFIX} CACHE PATH "Set specified install path" FORCE)
ELSE()
    IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        SET(INSTALL_DIR ${PROJECT_BINARY_DIR} CACHE PATH "Set specified install path" FORCE)
    ELSE ()
        SET(INSTALL_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Set default install path" FORCE)
    ENDIF ()
ENDIF ()

file(GLOB sources src/*.f90)
add_executable(swatplus_target ${sources})

target_compile_options(swatplus_target PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:${F_COMMON}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${F_DEBUG}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<NOT:$<CONFIG:Debug>>>:${F_NONDEBUG}>"
)

if(F_LINK)
    target_link_options(swatplus_target PRIVATE ${F_LINK})
endif()

set_target_properties(swatplus_target PROPERTIES
    OUTPUT_NAME_DEBUG          "${SWATPLUS_EXE_BASE}-Dbg"
    OUTPUT_NAME_RELEASE        "${SWATPLUS_EXE_BASE}-Rel-${OPT_LEVEL_TAG}"
    OUTPUT_NAME_RELWITHDEBINFO "${SWATPLUS_EXE_BASE}-RelWithDebInfo-${OPT_LEVEL_TAG}"
    OUTPUT_NAME_MINSIZEREL     "${SWATPLUS_EXE_BASE}-MinSizeRel-${OPT_LEVEL_TAG}"
)

install(TARGETS swatplus_target
        RUNTIME DESTINATION ${INSTALL_DIR}/bin)

find_package(Python3 QUIET)
if(Python3_FOUND)
    set(spcheck "${PROJECT_SOURCE_DIR}/test/spcheck.py")
    set(exe_path "$<TARGET_FILE:swatplus_target>")
    set(test_dir "${PROJECT_BINARY_DIR}/data")
    set(ref_dir  "${PROJECT_SOURCE_DIR}/data")

    set(rerr "0.01")
    set(aerr "1e-8")

    add_test(NAME Ames_sub1
        COMMAND ${Python3_EXECUTABLE} ${spcheck} ctest ${exe_path}
                ${ref_dir}/Ames_sub1 ${test_dir}
                --abserr ${aerr} --relerr ${rerr})
endif()


#############################################################################
# Summary
message( "-- Verbose Makefile      ${CMAKE_VERBOSE_MAKEFILE}")
message( "-- Generator             ${CMAKE_GENERATOR_TOOLSET}")
message( "-- VS Platform           ${CMAKE_VS_PLATFORM_TOOLSET}")
message( "-- Detected OS ID:       ${os_id}")
message( "-- Detected OS Version:  ${os_version}")
message( "-- Detected OS Like:     ${os_like}")
message( "-- Detected Pretty Name: ${pretty_name}")
message( "-- Build Type            ${CMAKE_BUILD_TYPE}")
message( "-- Compiler              ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
message( "-- Optimization Level    ${OPT_LEVEL}")
message( "-- Target COMMON Flags   ${F_COMMON}")
message( "-- Target DEBUG Flags    ${F_DEBUG}")
message( "-- Target NONDEBUG Flags ${F_NONDEBUG}")
message( "-- Target LINK Flags     ${F_LINK}")
message( "-- PROJECT_DIR           ${PROJECT_SOURCE_DIR}")
message( "-- BINARY_DIR            ${PROJECT_BINARY_DIR}")
message( "-- Install Prefix        ${CMAKE_INSTALL_PREFIX}")
message( "-- Executable Base       ${SWATPLUS_EXE_BASE}, rev/tag: ${SWAT_VERSION}")
message( "-- Version Source        ${SWAT_VERSION_SOURCE}")
message( "--  Debug Name           ${SWATPLUS_EXE_BASE}-Dbg")
message( "--  Release Name         ${SWATPLUS_EXE_BASE}-Rel-${OPT_LEVEL_TAG}")
message( "--  RelWithDebInfo Name  ${SWATPLUS_EXE_BASE}-RelWithDebInfo-${OPT_LEVEL_TAG}")
message( "--  MinSizeRel Name      ${SWATPLUS_EXE_BASE}-MinSizeRel-${OPT_LEVEL_TAG}")
message( "-- FORD project          ${PROJECT_SOURCE_DIR}/ford.md")
if(Python3_FOUND)
    message( "-- Test Data Dir         ${ref_dir}")
    message( "--  Abs Error            ${aerr}")
    message( "--  Rel Error            ${rerr}")
else()
    message( "-- Python3               not found; tests not configured")
endif()