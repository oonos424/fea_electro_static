cmake_minimum_required(VERSION 3.16)
project(FEM_Electrostatics LANGUAGES CXX)

# ---------------- Language & compile_commands ----------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------- Build type (enforce for single-config) -------------
# VS/Xcode have multi-config, NMake/Ninja/Makefiles are single-config.
if(NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "" CACHE STRING "Choose build type: Debug, Release, RelWithDebInfo, MinSizeRel")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
  if(NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR
      "CMAKE_BUILD_TYPE is not set.\n"
      "Choose one of: Debug, Release, RelWithDebInfo, MinSizeRel.\n"
      "Example (NMake): cmake -G \"NMake Makefiles\" -S . -B build-nmake -DCMAKE_BUILD_TYPE=Debug")
  endif()
endif()

# ------------------- Platform tags (explicit) -------------------
if(WIN32)
  add_compile_definitions(PLATFORM_WINDOWS)
elseif(APPLE)
  add_compile_definitions(PLATFORM_MACOS)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    add_compile_definitions(PLATFORM_APPLE_ARM64)
  endif()
elseif(UNIX)
  add_compile_definitions(PLATFORM_LINUX)
else()
  message(FATAL_ERROR "Unsupported platform. Expected WIN32, APPLE, or UNIX.")
endif()

# -------------------------- Target --------------------------
add_executable(fem_solver
  main.cpp
  electrostatic.cpp
)

# Nice output dirs (works for VS multi-config & single-config generators)
set_target_properties(fem_solver PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/$<CONFIG>"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/$<CONFIG>"
)

# Visual Studio QoL: startup project + debugger working dir
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT fem_solver)
set_target_properties(fem_solver PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

# MSVC: enable parallel compilation (/MP)
if(MSVC)
  target_compile_options(fem_solver PRIVATE /MP)
endif()

# Option: select MSVC runtime (default: dynamic MD/MDd). Change in cache if desired.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL" CACHE STRING
    "MSVC runtime: MultiThreaded[Debug] or MultiThreaded[Debug]DLL")

# --------------- Compiler families (explicit branches) ---------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # Warnings & conformance
  target_compile_options(fem_solver PRIVATE /W4 /permissive-)
  # Config flags
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} /Zi /Od /RTC1")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(fem_solver PRIVATE -Wall -Wextra -Wpedantic)
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  target_compile_options(fem_solver PRIVATE -Wall -Wextra -Wpedantic)
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  if(WIN32)
    # Intel on Windows: MSVC-like flags
    target_compile_options(fem_solver PRIVATE /W4 /permissive-)
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} /Zi /Od")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
  else()
    # Intel on Linux/macOS: GCC-like flags
    target_compile_options(fem_solver PRIVATE -Wall -Wextra -Wpedantic)
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
  endif()

else()
  message(FATAL_ERROR
    "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}\n"
    "Supported: MSVC, GNU, Clang, AppleClang, Intel, IntelLLVM.\n"
    "Hint: set compilers explicitly, e.g.\n"
    "  -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl (MSVC)\n"
    "  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ (GCC)\n"
    "  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ (Clang)\n"
    "  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icx|icpx (Intel oneAPI)")
endif()

# ------------------------ Linking rules ------------------------
if(WIN32)
  # No libm on Windows/MSVC
else()
  target_link_libraries(fem_solver PRIVATE m)
endif()
