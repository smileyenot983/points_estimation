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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/user/points_estimation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/points_estimation/build

# Include any dependencies generated for this target.
include CMakeFiles/estimate_points.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/estimate_points.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/estimate_points.dir/flags.make

CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o: CMakeFiles/estimate_points.dir/flags.make
CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o: ../src/estimate_points.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/points_estimation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o -c /home/user/points_estimation/src/estimate_points.cpp

CMakeFiles/estimate_points.dir/src/estimate_points.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/estimate_points.dir/src/estimate_points.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/points_estimation/src/estimate_points.cpp > CMakeFiles/estimate_points.dir/src/estimate_points.cpp.i

CMakeFiles/estimate_points.dir/src/estimate_points.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/estimate_points.dir/src/estimate_points.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/points_estimation/src/estimate_points.cpp -o CMakeFiles/estimate_points.dir/src/estimate_points.cpp.s

# Object files for target estimate_points
estimate_points_OBJECTS = \
"CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o"

# External object files for target estimate_points
estimate_points_EXTERNAL_OBJECTS =

estimate_points: CMakeFiles/estimate_points.dir/src/estimate_points.cpp.o
estimate_points: CMakeFiles/estimate_points.dir/build.make
estimate_points: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
estimate_points: /usr/lib/x86_64-linux-gnu/libpthread.so
estimate_points: CMakeFiles/estimate_points.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/user/points_estimation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable estimate_points"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/estimate_points.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/estimate_points.dir/build: estimate_points

.PHONY : CMakeFiles/estimate_points.dir/build

CMakeFiles/estimate_points.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/estimate_points.dir/cmake_clean.cmake
.PHONY : CMakeFiles/estimate_points.dir/clean

CMakeFiles/estimate_points.dir/depend:
	cd /home/user/points_estimation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/user/points_estimation /home/user/points_estimation /home/user/points_estimation/build /home/user/points_estimation/build /home/user/points_estimation/build/CMakeFiles/estimate_points.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/estimate_points.dir/depend

