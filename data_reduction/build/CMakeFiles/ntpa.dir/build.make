# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /cvmfs/ams.cern.ch/Offline/dbar/public/installed/bin/cmake

# The command to remove a file.
RM = /cvmfs/ams.cern.ch/Offline/dbar/public/installed/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build

# Include any dependencies generated for this target.
include CMakeFiles/ntpa.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ntpa.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ntpa.dir/flags.make

G__ntp.cxx: ../include/LinkDef.h
G__ntp.cxx: ../include/Ntp.h
G__ntp.cxx: ../include/VPSArchive.h
G__ntp.cxx: ../include/RichOccupancy.h
G__ntp.cxx: ../include/Classifier.h
G__ntp.cxx: ../include/RichBDT.h
G__ntp.cxx: ../include/MSpline.h
G__ntp.cxx: ../include/IsoPoly.h
G__ntp.cxx: ../include/PDF1.h
G__ntp.cxx: ../include/PDF2.h
G__ntp.cxx: ../include/PDF2DB.h
G__ntp.cxx: ../include/LikelihoodVar.h
G__ntp.cxx: ../include/Likelihood.h
G__ntp.cxx: ../include/Ntp.h
G__ntp.cxx: ../include/VPSArchive.h
G__ntp.cxx: ../include/RichOccupancy.h
G__ntp.cxx: ../include/Classifier.h
G__ntp.cxx: ../include/RichBDT.h
G__ntp.cxx: ../include/MSpline.h
G__ntp.cxx: ../include/IsoPoly.h
G__ntp.cxx: ../include/PDF1.h
G__ntp.cxx: ../include/PDF2.h
G__ntp.cxx: ../include/PDF2DB.h
G__ntp.cxx: ../include/LikelihoodVar.h
G__ntp.cxx: ../include/Likelihood.h
G__ntp.cxx: ../include/LinkDef.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating G__ntp.cxx, libntp_rdict.pcm, libntp.rootmap"
	/cvmfs/ams.cern.ch/Offline/dbar/public/installed/bin/cmake -E env LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/lib:/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/lib64:/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/lib64:/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/lib64:/afs/cern.ch/work/f/fdimicco/private/Deutons/data_reduction/lib:/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/lib64:/afs/cern.ch/work/f/fdimicco/private/Deutons/data_reduction/lib:/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/lib64:/afs/cern.ch/ams/local2/opt/xrootd-gcc64-41/lib64:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/lib:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/lib:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/lib:/cvmfs/ams.cern.ch/Offline/dbar/public/installed/root-6.14.04/lib:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/lib /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/bin/rootcling -v2 -f G__ntp.cxx -s /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/libntp.so -rml libntp.so -rmf /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/libntp.rootmap -D__ROOTSHAREDLIBRARY__ -D_PGTRACK_ -D_NTPLIB_ -I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/include -I/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include -I/data1/home/pzuccon/work/dbar/AMS/include /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Ntp.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/VPSArchive.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichOccupancy.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Classifier.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/RichBDT.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/MSpline.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/IsoPoly.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF1.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/PDF2DB.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/LikelihoodVar.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/Likelihood.h /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/include/LinkDef.h

libntp_rdict.pcm: G__ntp.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libntp_rdict.pcm

libntp.rootmap: G__ntp.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libntp.rootmap

CMakeFiles/ntpa.dir/src/Ntp.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/Ntp.cxx.o: ../src/Ntp.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ntpa.dir/src/Ntp.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/Ntp.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Ntp.cxx

CMakeFiles/ntpa.dir/src/Ntp.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/Ntp.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Ntp.cxx > CMakeFiles/ntpa.dir/src/Ntp.cxx.i

CMakeFiles/ntpa.dir/src/Ntp.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/Ntp.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Ntp.cxx -o CMakeFiles/ntpa.dir/src/Ntp.cxx.s

CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o: ../src/VPSArchive.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/VPSArchive.cxx

CMakeFiles/ntpa.dir/src/VPSArchive.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/VPSArchive.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/VPSArchive.cxx > CMakeFiles/ntpa.dir/src/VPSArchive.cxx.i

CMakeFiles/ntpa.dir/src/VPSArchive.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/VPSArchive.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/VPSArchive.cxx -o CMakeFiles/ntpa.dir/src/VPSArchive.cxx.s

CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o: ../src/RichOccupancy.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichOccupancy.cxx

CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichOccupancy.cxx > CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.i

CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichOccupancy.cxx -o CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.s

CMakeFiles/ntpa.dir/src/RichBDT.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/RichBDT.cxx.o: ../src/RichBDT.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/ntpa.dir/src/RichBDT.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/RichBDT.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichBDT.cxx

CMakeFiles/ntpa.dir/src/RichBDT.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/RichBDT.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichBDT.cxx > CMakeFiles/ntpa.dir/src/RichBDT.cxx.i

CMakeFiles/ntpa.dir/src/RichBDT.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/RichBDT.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/RichBDT.cxx -o CMakeFiles/ntpa.dir/src/RichBDT.cxx.s

CMakeFiles/ntpa.dir/G__ntp.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/G__ntp.cxx.o: G__ntp.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/ntpa.dir/G__ntp.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/G__ntp.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/G__ntp.cxx

CMakeFiles/ntpa.dir/G__ntp.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/G__ntp.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/G__ntp.cxx > CMakeFiles/ntpa.dir/G__ntp.cxx.i

CMakeFiles/ntpa.dir/G__ntp.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/G__ntp.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/G__ntp.cxx -o CMakeFiles/ntpa.dir/G__ntp.cxx.s

CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o: ../src/TrackerBDT.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/TrackerBDT.cxx

CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/TrackerBDT.cxx > CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.i

CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/TrackerBDT.cxx -o CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.s

CMakeFiles/ntpa.dir/src/Classifier.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/Classifier.cxx.o: ../src/Classifier.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/ntpa.dir/src/Classifier.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/Classifier.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Classifier.cxx

CMakeFiles/ntpa.dir/src/Classifier.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/Classifier.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Classifier.cxx > CMakeFiles/ntpa.dir/src/Classifier.cxx.i

CMakeFiles/ntpa.dir/src/Classifier.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/Classifier.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Classifier.cxx -o CMakeFiles/ntpa.dir/src/Classifier.cxx.s

CMakeFiles/ntpa.dir/src/MSpline.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/MSpline.cxx.o: ../src/MSpline.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/ntpa.dir/src/MSpline.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/MSpline.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/MSpline.cxx

CMakeFiles/ntpa.dir/src/MSpline.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/MSpline.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/MSpline.cxx > CMakeFiles/ntpa.dir/src/MSpline.cxx.i

CMakeFiles/ntpa.dir/src/MSpline.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/MSpline.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/MSpline.cxx -o CMakeFiles/ntpa.dir/src/MSpline.cxx.s

CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o: ../src/IsoPoly.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/IsoPoly.cxx

CMakeFiles/ntpa.dir/src/IsoPoly.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/IsoPoly.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/IsoPoly.cxx > CMakeFiles/ntpa.dir/src/IsoPoly.cxx.i

CMakeFiles/ntpa.dir/src/IsoPoly.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/IsoPoly.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/IsoPoly.cxx -o CMakeFiles/ntpa.dir/src/IsoPoly.cxx.s

CMakeFiles/ntpa.dir/src/PDF1.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/PDF1.cxx.o: ../src/PDF1.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/ntpa.dir/src/PDF1.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/PDF1.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF1.cxx

CMakeFiles/ntpa.dir/src/PDF1.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/PDF1.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF1.cxx > CMakeFiles/ntpa.dir/src/PDF1.cxx.i

CMakeFiles/ntpa.dir/src/PDF1.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/PDF1.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF1.cxx -o CMakeFiles/ntpa.dir/src/PDF1.cxx.s

CMakeFiles/ntpa.dir/src/PDF2.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/PDF2.cxx.o: ../src/PDF2.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/ntpa.dir/src/PDF2.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/PDF2.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2.cxx

CMakeFiles/ntpa.dir/src/PDF2.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/PDF2.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2.cxx > CMakeFiles/ntpa.dir/src/PDF2.cxx.i

CMakeFiles/ntpa.dir/src/PDF2.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/PDF2.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2.cxx -o CMakeFiles/ntpa.dir/src/PDF2.cxx.s

CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o: ../src/PDF2DB.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2DB.cxx

CMakeFiles/ntpa.dir/src/PDF2DB.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/PDF2DB.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2DB.cxx > CMakeFiles/ntpa.dir/src/PDF2DB.cxx.i

CMakeFiles/ntpa.dir/src/PDF2DB.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/PDF2DB.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/PDF2DB.cxx -o CMakeFiles/ntpa.dir/src/PDF2DB.cxx.s

CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o: ../src/LikelihoodVar.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/LikelihoodVar.cxx

CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/LikelihoodVar.cxx > CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.i

CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/LikelihoodVar.cxx -o CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.s

CMakeFiles/ntpa.dir/src/Likelihood.cxx.o: CMakeFiles/ntpa.dir/flags.make
CMakeFiles/ntpa.dir/src/Likelihood.cxx.o: ../src/Likelihood.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/ntpa.dir/src/Likelihood.cxx.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ntpa.dir/src/Likelihood.cxx.o -c /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Likelihood.cxx

CMakeFiles/ntpa.dir/src/Likelihood.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ntpa.dir/src/Likelihood.cxx.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Likelihood.cxx > CMakeFiles/ntpa.dir/src/Likelihood.cxx.i

CMakeFiles/ntpa.dir/src/Likelihood.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ntpa.dir/src/Likelihood.cxx.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/src/Likelihood.cxx -o CMakeFiles/ntpa.dir/src/Likelihood.cxx.s

# Object files for target ntpa
ntpa_OBJECTS = \
"CMakeFiles/ntpa.dir/src/Ntp.cxx.o" \
"CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o" \
"CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o" \
"CMakeFiles/ntpa.dir/src/RichBDT.cxx.o" \
"CMakeFiles/ntpa.dir/G__ntp.cxx.o" \
"CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o" \
"CMakeFiles/ntpa.dir/src/Classifier.cxx.o" \
"CMakeFiles/ntpa.dir/src/MSpline.cxx.o" \
"CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o" \
"CMakeFiles/ntpa.dir/src/PDF1.cxx.o" \
"CMakeFiles/ntpa.dir/src/PDF2.cxx.o" \
"CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o" \
"CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o" \
"CMakeFiles/ntpa.dir/src/Likelihood.cxx.o"

# External object files for target ntpa
ntpa_EXTERNAL_OBJECTS =

libntpa.a: CMakeFiles/ntpa.dir/src/Ntp.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/VPSArchive.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/RichOccupancy.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/RichBDT.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/G__ntp.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/TrackerBDT.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/Classifier.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/MSpline.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/IsoPoly.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/PDF1.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/PDF2.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/PDF2DB.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/LikelihoodVar.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/src/Likelihood.cxx.o
libntpa.a: CMakeFiles/ntpa.dir/build.make
libntpa.a: CMakeFiles/ntpa.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking CXX static library libntpa.a"
	$(CMAKE_COMMAND) -P CMakeFiles/ntpa.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ntpa.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ntpa.dir/build: libntpa.a

.PHONY : CMakeFiles/ntpa.dir/build

CMakeFiles/ntpa.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ntpa.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ntpa.dir/clean

CMakeFiles/ntpa.dir/depend: G__ntp.cxx
CMakeFiles/ntpa.dir/depend: libntp_rdict.pcm
CMakeFiles/ntpa.dir/depend: libntp.rootmap
	cd /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build /afs/cern.ch/user/f/fdimicco/Work/Deutons/data_reduction/build/CMakeFiles/ntpa.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ntpa.dir/depend
