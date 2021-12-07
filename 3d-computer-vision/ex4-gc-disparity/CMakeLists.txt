cmake_minimum_required(VERSION 3.0)
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?

project(GraphCuts)
find_package(Imagine REQUIRED)

set(CMAKE_CXX_STANDARD 11)

add_executable(exampleGC exampleGC.cpp maxflow/graph.cpp)

add_executable(stereoGC stereoGC.cpp maxflow/graph.cpp)
ImagineUseModules(stereoGC Images)
