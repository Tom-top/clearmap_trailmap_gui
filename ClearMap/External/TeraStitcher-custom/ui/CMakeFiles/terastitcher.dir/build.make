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
include ui/CMakeFiles/terastitcher.dir/depend.make

# Include the progress variables for this target.
include ui/CMakeFiles/terastitcher.dir/progress.make

# Include the compile flags for this target's objects.
include ui/CMakeFiles/terastitcher.dir/flags.make

ui/CMakeFiles/terastitcher.dir/main.cpp.o: ui/CMakeFiles/terastitcher.dir/flags.make
ui/CMakeFiles/terastitcher.dir/main.cpp.o: /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas.renier/Programs/TeraStitcher-custom/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ui/CMakeFiles/terastitcher.dir/main.cpp.o"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/terastitcher.dir/main.cpp.o -c /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/main.cpp

ui/CMakeFiles/terastitcher.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/terastitcher.dir/main.cpp.i"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/main.cpp > CMakeFiles/terastitcher.dir/main.cpp.i

ui/CMakeFiles/terastitcher.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/terastitcher.dir/main.cpp.s"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/main.cpp -o CMakeFiles/terastitcher.dir/main.cpp.s

ui/CMakeFiles/terastitcher.dir/main.cpp.o.requires:

.PHONY : ui/CMakeFiles/terastitcher.dir/main.cpp.o.requires

ui/CMakeFiles/terastitcher.dir/main.cpp.o.provides: ui/CMakeFiles/terastitcher.dir/main.cpp.o.requires
	$(MAKE) -f ui/CMakeFiles/terastitcher.dir/build.make ui/CMakeFiles/terastitcher.dir/main.cpp.o.provides.build
.PHONY : ui/CMakeFiles/terastitcher.dir/main.cpp.o.provides

ui/CMakeFiles/terastitcher.dir/main.cpp.o.provides.build: ui/CMakeFiles/terastitcher.dir/main.cpp.o


ui/CMakeFiles/terastitcher.dir/CLI.cpp.o: ui/CMakeFiles/terastitcher.dir/flags.make
ui/CMakeFiles/terastitcher.dir/CLI.cpp.o: /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/CLI.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas.renier/Programs/TeraStitcher-custom/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object ui/CMakeFiles/terastitcher.dir/CLI.cpp.o"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/terastitcher.dir/CLI.cpp.o -c /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/CLI.cpp

ui/CMakeFiles/terastitcher.dir/CLI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/terastitcher.dir/CLI.cpp.i"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/CLI.cpp > CMakeFiles/terastitcher.dir/CLI.cpp.i

ui/CMakeFiles/terastitcher.dir/CLI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/terastitcher.dir/CLI.cpp.s"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui/CLI.cpp -o CMakeFiles/terastitcher.dir/CLI.cpp.s

ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.requires:

.PHONY : ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.requires

ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.provides: ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.requires
	$(MAKE) -f ui/CMakeFiles/terastitcher.dir/build.make ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.provides.build
.PHONY : ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.provides

ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.provides.build: ui/CMakeFiles/terastitcher.dir/CLI.cpp.o


# Object files for target terastitcher
terastitcher_OBJECTS = \
"CMakeFiles/terastitcher.dir/main.cpp.o" \
"CMakeFiles/terastitcher.dir/CLI.cpp.o"

# External object files for target terastitcher
terastitcher_EXTERNAL_OBJECTS =

bin/terastitcher: ui/CMakeFiles/terastitcher.dir/main.cpp.o
bin/terastitcher: ui/CMakeFiles/terastitcher.dir/CLI.cpp.o
bin/terastitcher: ui/CMakeFiles/terastitcher.dir/build.make
bin/terastitcher: stitcher/libstitcher.a
bin/terastitcher: volumemanager/libvolumemanager.a
bin/terastitcher: crossmips/libcrossmips.a
bin/terastitcher: iomanager/libiomanager.a
bin/terastitcher: 3rdparty/tinyxml/libtinyxml.a
bin/terastitcher: common/libcommon.a
bin/terastitcher: iomanager/plugins/tiff3D/libioplugin_tiff3D.a
bin/terastitcher: iomanager/plugins/tiff2D/libioplugin_tiff2D.a
bin/terastitcher: 3rdparty/libtiff/libtiff.a
bin/terastitcher: 3rdparty/zlib/libzlib.a
bin/terastitcher: imagemanager/libimagemanager.a
bin/terastitcher: ui/CMakeFiles/terastitcher.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nicolas.renier/Programs/TeraStitcher-custom/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../bin/terastitcher"
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/terastitcher.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ui/CMakeFiles/terastitcher.dir/build: bin/terastitcher

.PHONY : ui/CMakeFiles/terastitcher.dir/build

ui/CMakeFiles/terastitcher.dir/requires: ui/CMakeFiles/terastitcher.dir/main.cpp.o.requires
ui/CMakeFiles/terastitcher.dir/requires: ui/CMakeFiles/terastitcher.dir/CLI.cpp.o.requires

.PHONY : ui/CMakeFiles/terastitcher.dir/requires

ui/CMakeFiles/terastitcher.dir/clean:
	cd /home/nicolas.renier/Programs/TeraStitcher-custom/ui && $(CMAKE_COMMAND) -P CMakeFiles/terastitcher.dir/cmake_clean.cmake
.PHONY : ui/CMakeFiles/terastitcher.dir/clean

ui/CMakeFiles/terastitcher.dir/depend:
	cd /home/nicolas.renier/Programs/TeraStitcher-custom && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src /home/nicolas.renier/Programs/TeraStitcher-19190f8f56698079b03d6313a69650488b42df77/src/ui /home/nicolas.renier/Programs/TeraStitcher-custom /home/nicolas.renier/Programs/TeraStitcher-custom/ui /home/nicolas.renier/Programs/TeraStitcher-custom/ui/CMakeFiles/terastitcher.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ui/CMakeFiles/terastitcher.dir/depend

