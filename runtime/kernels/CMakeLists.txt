# runtime/kernels/CMakeLists.txt
# Kernels directory CMakeLists.txt

# Add shared components subdirectory first (needed by other subdirectories)
add_subdirectory(shared)

# Add profiler kernels subdirectory
add_subdirectory(profiler)

# Add resistance detection kernels subdirectory  
add_subdirectory(resistance)

# Add AMR gene detection kernels subdirectory
add_subdirectory(genes)