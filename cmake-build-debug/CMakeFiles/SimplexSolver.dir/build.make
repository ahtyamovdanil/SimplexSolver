# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/mirror/soft/clion-2020.1/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/mirror/soft/clion-2020.1/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mirror/university/diploma/SimplexSolver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mirror/university/diploma/SimplexSolver/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SimplexSolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SimplexSolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SimplexSolver.dir/flags.make

CMakeFiles/SimplexSolver.dir/main.cpp.o: CMakeFiles/SimplexSolver.dir/flags.make
CMakeFiles/SimplexSolver.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mirror/university/diploma/SimplexSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SimplexSolver.dir/main.cpp.o"
	/usr/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SimplexSolver.dir/main.cpp.o -c /home/mirror/university/diploma/SimplexSolver/main.cpp

CMakeFiles/SimplexSolver.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SimplexSolver.dir/main.cpp.i"
	/usr/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mirror/university/diploma/SimplexSolver/main.cpp > CMakeFiles/SimplexSolver.dir/main.cpp.i

CMakeFiles/SimplexSolver.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SimplexSolver.dir/main.cpp.s"
	/usr/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mirror/university/diploma/SimplexSolver/main.cpp -o CMakeFiles/SimplexSolver.dir/main.cpp.s

# Object files for target SimplexSolver
SimplexSolver_OBJECTS = \
"CMakeFiles/SimplexSolver.dir/main.cpp.o"

# External object files for target SimplexSolver
SimplexSolver_EXTERNAL_OBJECTS =

SimplexSolver: CMakeFiles/SimplexSolver.dir/main.cpp.o
SimplexSolver: CMakeFiles/SimplexSolver.dir/build.make
SimplexSolver: CMakeFiles/SimplexSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mirror/university/diploma/SimplexSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SimplexSolver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SimplexSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SimplexSolver.dir/build: SimplexSolver

.PHONY : CMakeFiles/SimplexSolver.dir/build

CMakeFiles/SimplexSolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SimplexSolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SimplexSolver.dir/clean

CMakeFiles/SimplexSolver.dir/depend:
	cd /home/mirror/university/diploma/SimplexSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mirror/university/diploma/SimplexSolver /home/mirror/university/diploma/SimplexSolver /home/mirror/university/diploma/SimplexSolver/cmake-build-debug /home/mirror/university/diploma/SimplexSolver/cmake-build-debug /home/mirror/university/diploma/SimplexSolver/cmake-build-debug/CMakeFiles/SimplexSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SimplexSolver.dir/depend
