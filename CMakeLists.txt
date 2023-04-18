cmake_minimum_required(VERSION 3.0)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
endif()

# Set the project name and basic settings
project(BCNUTAU LANGUAGES CXX VERSION 1.0.0)
include(CMake/StandardProjectSettings.cmake)

file(RELATIVE_PATH relativeRpath
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${relativeRpath})

# Link this 'library' to set the c++ standard / compile-time options requested
# Additionally, link to get include and external dependencies
add_library(BCNUTAU_project_options INTERFACE)
target_compile_features(BCNUTAU_project_options INTERFACE cxx_std_17)
target_include_directories(BCNUTAU_project_options INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(BCNUTAU_project_warnings INTERFACE)

# Standard compiler warnings
include(CMake/CompilerWarnings.cmake)
set_project_warnings(BCNUTAU_project_warnings)

# Sanitizer options if supported by compiler
include(CMake/Sanitizers.cmake)
enable_sanitizers(BCNUTAU_project_options)

# Allow for static analysis options
include(CMake/StaticAnalyzers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)
option(BCNUTAU_ENABLE_TESTING "Enable Test Builds" OFF)

# Add dependencies
include(CMake/CPM.cmake)
add_subdirectory(external)

# Testing
if(BCNUTAU_ENABLE_TESTING)
    if(BCNUTAU_COVERAGE)
        target_compile_options(BCNUTAU_project_options INTERFACE --coverage)
        target_link_libraries(BCNUTAU_project_options INTERFACE gcov)
    endif()
    target_compile_definitions(BCNUTAU_project_options INTERFACE TESTING)
    enable_testing()
    message(STATUS "Building Tests.")
    add_subdirectory(test)
endif()

# Main code
add_subdirectory(src)
