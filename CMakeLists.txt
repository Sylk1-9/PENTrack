cmake_minimum_required(VERSION 3.18)
# cmake_minimum_required (VERSION 3.1)
project (PENTrack CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("alglib-3.15.0/cpp/src")
add_library(alglib OBJECT alglib-3.15.0/cpp/src/alglibinternal.cpp
  alglib-3.15.0/cpp/src/integration.cpp
  alglib-3.15.0/cpp/src/optimization.cpp
  alglib-3.15.0/cpp/src/alglibmisc.cpp
  alglib-3.15.0/cpp/src/interpolation.cpp
  alglib-3.15.0/cpp/src/solvers.cpp
  alglib-3.15.0/cpp/src/ap.cpp
  alglib-3.15.0/cpp/src/linalg.cpp
  alglib-3.15.0/cpp/src/specialfunctions.cpp)

include_directories("tricubic-1.0/src/libtricubic")
add_library(libtricubic OBJECT tricubic-1.0/src/libtricubic/libtricubic.cpp tricubic-1.0/src/libtricubic/tricubic_utils.cpp)

include_directories("exprtk")
include_directories("include")

set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Tell CGAL not to override compiler flags")
find_package(CGAL HINTS "cgal")
if (CGAL_FOUND)
  message(STATUS "CGAL version: " ${CGAL_VERSION})
  if (CGAL_VERSION VERSION_LESS 4.14)
    message(SEND_ERROR "CGAL versions older than 4.14 are not supported!")
  endif()
  include_directories(${CGAL_INCLUDE_DIRS})
  include(${CGAL_USE_FILE})
else()
  message(SEND_ERROR "CGAL library not found!")
endif()

set(Boost_REALPATH ON)
if (BUILD_TESTS)
  message(STATUS "Tests will be built, requires Boost 1.59.0 or newer with unit test framework")
  # find_package(Boost 1.59.0 REQUIRED system filesystem iostreams unit_test_framework)
  find_package(Boost REQUIRED python system filesystem iostreams unit_test_framework)
else()
  # find_package(Boost 1.53.0 REQUIRED system filesystem iostreams)
    find_package(Boost REQUIRED python system filesystem iostreams)
endif()
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  link_directories(${Boost_LIBRARY_DIRS})
endif()


set(pentrack_SOURCES src/globals.cpp src/trianglemesh.cpp src/geometry.cpp src/mc.cpp src/field.cpp src/edmfields.cpp src/tracking.cpp src/logger.cpp src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/harmonicfields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp src/config.cpp src/analyticFields.cpp src/threadpool.cpp)
 

# ind_package(Python3 REQUIRED Development)
find_package (Python3 REQUIRED COMPONENTS Interpreter Development)
if (Python3_FOUND)
  message (STATUS "Found Python3 " ${Python3_VERSION})
  message (STATUS "Found Python3 in " ${Python3_EXECUTABLE})
  if (Python3_INCLUDE_DIRS)
    message (STATUS "Found Python3 include directories: ")
    include_directories(${Python3_INCLUDE_DIRS})
    # add_library(pythonFields_src OBJECT src/pythonFields.cpp)
    set(pentrack_SOURCES ${pentrack_SOURCES} src/pythonFields.cpp)
  else()
    message (STATUS "Couldn’t find Python3 include directories!")
  endif()
else()
  message (STATUS "Couldn’t find Python3!")
endif()

if (Python3_Interpreter_FOUND)
  message (STATUS "Found Python3 Interpreter")
else()
  message (STATUS "Couldn’t find Python3 Interpreter!")
endif()

if (Python3_Development_FOUND)
  message (STATUS "Found Python3 Development")
else()
  message (STATUS "Couldn’t find Python3 Development!")
endif()


# add_library(PENTrack_src OBJECT src/globals.cpp src/trianglemesh.cpp src/geometry.cpp src/mc.cpp src/field.cpp src/edmfields.cpp src/tracking.cpp src/logger.cpp
#   src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/harmonicfields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp
#   src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp src/config.cpp src/analyticFields.cpp src/threadpool.cpp)

add_library(PENTrack_src OBJECT ${pentrack_SOURCES})

add_executable(PENTrack src/main.cpp $<TARGET_OBJECTS:PENTrack_src> $<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:libtricubic>)


# target_link_libraries (PENTrack ${Boost_LIBRARIES} ${CGAL_LIBRARIES} Python3::Python)
target_link_libraries (PENTrack ${Boost_LIBRARIES} ${CGAL_LIBRARIES})


# Copy the Python module to the build directory
# configure_file(src/MagpyField.py MagpyField.py COPYONLY)
# configure_file(in/magpytest.py magpytest.py COPYONLY)

if (Python3_FOUND)
  message (STATUS "Found Python3 " ${Python3_VERSION})
  message (STATUS "Found Python3 in " ${Python3_EXECUTABLE})
  if (Python3_INCLUDE_DIRS)
    message (STATUS "Found Python3 include directories: ")
    # include_directories(${Python3_INCLUDE_DIRS})
    target_compile_definitions(PENTrack_src PUBLIC USEPYTHON=1)
    target_compile_definitions(PENTrack PUBLIC USEPYTHON=1)
    target_link_libraries (PENTrack_src Python3::Python)
    target_link_libraries (PENTrack Python3::Python)
  endif()
else()
  message (STATUS "Couldn’t find Python3!")
endif()



find_package(ROOT)
if (ROOT_FOUND)
  include(${ROOT_USE_FILE})
  target_compile_definitions(PENTrack_src PUBLIC USEROOT=1)
  target_link_libraries (PENTrack ${ROOT_LIBRARIES})
else()
  message(STATUS "Could not find ROOT, you won't be able to log to ROOT files")
endif()
find_package(HDF5 COMPONENTS HL)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_compile_definitions(PENTrack_src PUBLIC USEHDF5=1)
  target_link_libraries (PENTrack ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  # set(HDF5_USE_FILE_LOCKING FALSE)
else()
  message(STATUS "Could not find HDF5 library, you won't be able to log to HDF5 files")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  target_compile_options(PENTrack_src PUBLIC -Wall)
endif()



if (BUILD_TESTS)
  enable_testing()
  add_executable(runTests test/test.cpp test/fieldTests.cpp test/microroughnessTests.cpp $<TARGET_OBJECTS:PENTrack_src> $<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:libtricubic>)
  target_link_libraries(runTests ${Boost_LIBRARIES} ${CGAL_LIBRARIES})
  if (ROOT_FOUND)
    target_link_libraries(runTests ${ROOT_LIBRARIES})
  endif()
  if (HDF5_FOUND)
    target_link_libraries(runTests ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    set(HDF5_USE_FILE_LOCKING FALSE)
  endif()
  target_compile_definitions(runTests PRIVATE "BOOST_TEST_DYN_LINK=1")
  add_test(COMMAND runTests)
endif()




