cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_STANDARD 23)

project(GravitationalSearch LANGUAGES CXX)

# find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-shared) 
find_package (Eigen3 REQUIRED NO_MODULE)

set(EXEC_NAME "grav-search")

file(GLOB SRC_FILES "src/*.cpp")
add_executable(${EXEC_NAME} ${SRC_FILES})

target_link_libraries (${EXEC_NAME} Eigen3::Eigen)

# target_link_libraries(${EXEC_NAME} PRIVATE SDL3::SDL3)
# set_target_properties(${EXEC_NAME} PROPERTIES
#   CXX_STANDARD 23
#   CXX_STANDARD_REQUIRED FALSE
# )
