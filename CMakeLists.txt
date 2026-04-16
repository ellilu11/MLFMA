# CMakeList.txt : CMake project for MLFMA, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo,Release>:ProgramDatabase>")
endif()

project ("MLFMA" CXX)

set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/DEBUG")

# Add eigen library
add_library(eigen INTERFACE)
target_include_directories(eigen
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0
)

# Add OpenMP support
find_package(OpenMP REQUIRED)

# Add source executable
add_executable (MLFMA "src/main.cpp" "src/main.h")

target_link_libraries(MLFMA PRIVATE eigen)
target_link_libraries(MLFMA PRIVATE OpenMP::OpenMP_CXX)

# target_compile_options(MLFMA PRIVATE $<$<CONFIG:Release>:/Wall>)

# Add test executable
add_executable (test "test/tRwg.cpp" "src/main.h")

target_link_libraries(test PRIVATE eigen)

# Specify C++20 standard
if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET MLFMA PROPERTY CXX_STANDARD 20)
  set_property(TARGET test PROPERTY CXX_STANDARD 20)
endif()