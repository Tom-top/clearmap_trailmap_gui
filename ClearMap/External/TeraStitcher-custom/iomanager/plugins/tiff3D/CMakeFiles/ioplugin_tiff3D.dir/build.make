# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nicolas.renier/Programs/TeraStitcher-custom

# Include any dependencies generated for this target.
include iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/depend.make

# Include the progress variables for this target.
include iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/progress.make

# Include the compile flags for this target's objects.
include iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/flags.make

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/flags.make
iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o: /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/iomanager/plugins/tiff3D/tiff3D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas.renier/Programs/TeraStitcher-custom/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o -c /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/iomanager/plugins/tiff3D/tiff3D.cpp

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.i"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/iomanager/plugins/tiff3D/tiff3D.cpp > CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.i

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.s"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/iomanager/plugins/tiff3D/tiff3D.cpp -o CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.s

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.requires:

.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.requires

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.provides: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.requires
	$(MAKE) -f iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/build.make iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.provides.build
.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.provides

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.provides.build: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o


# Object files for target ioplugin_tiff3D
ioplugin_tiff3D_OBJECTS = \
"CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o"

# External object files for target ioplugin_tiff3D
ioplugin_tiff3D_EXTERNAL_OBJECTS =

iomanager/plugins/tiff3D/libioplugin_tiff3D.a: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o
iomanager/plugins/tiff3D/libioplugin_tiff3D.a: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/build.make
iomanager/plugins/tiff3D/libioplugin_tiff3D.a: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nicolas.renier/Programs/TeraStitcher-custom/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libioplugin_tiff3D.a"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && $(CMAKE_COMMAND) -P CMakeFiles/ioplugin_tiff3D.dir/cmake_clean_target.cmake
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ioplugin_tiff3D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/build: iomanager/plugins/tiff3D/libioplugin_tiff3D.a

.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/build

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/requires: iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/tiff3D.cpp.o.requires

.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/requires

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/clean:
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D && $(CMAKE_COMMAND) -P CMakeFiles/ioplugin_tiff3D.dir/cmake_clean.cmake
.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/clean

iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/depend:
	cd /home/nicolas.renier/Programs/TeraStitcher-custom && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/iomanager/plugins/tiff3D /home/nicolas.renier/Programs/TeraStitcher-custom /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D /home/nicolas.renier/Programs/TeraStitcher-custom/iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : iomanager/plugins/tiff3D/CMakeFiles/ioplugin_tiff3D.dir/depend

