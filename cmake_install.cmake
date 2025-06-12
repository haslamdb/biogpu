# Install script for directory: /home/david/Documents/Code/biogpu

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/david/Documents/Code/biogpu/runtime/kernels/resistance/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/david/Documents/Code/biogpu/runtime/kernels/profiler/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/home/david/Documents/Code/biogpu/src/python/build_integrated_resistance_db.py"
    "/home/david/Documents/Code/biogpu/src/python/build_clean_resistance_db.py"
    "/home/david/Documents/Code/biogpu/src/python/validate_integrated_db.py"
    "/home/david/Documents/Code/biogpu/src/python/resistance_profile_analyzer.py"
    "/home/david/Documents/Code/biogpu/src/python/clinical_interpreter.py"
    "/home/david/Documents/Code/biogpu/src/python/benchmark_integrated.py"
    "/home/david/Documents/Code/biogpu/src/python/read_hdf5_alignments.py"
    "/home/david/Documents/Code/biogpu/src/python/enhanced_kmer_builder.py"
    "/home/david/Documents/Code/biogpu/src/python/index_validator.py"
    "/home/david/Documents/Code/biogpu/src/python/build_protein_resistance_db.py"
    "/home/david/Documents/Code/biogpu/src/python/diagnostic_analyzer.py"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/biogpu/scripts" TYPE DIRECTORY FILES "/home/david/Documents/Code/biogpu/scripts/" FILES_MATCHING REGEX "/[^/]*\\.py$" REGEX "/[^/]*\\.sh$")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/david/Documents/Code/biogpu/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
