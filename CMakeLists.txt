project(BTA_service_tools)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#aux_source_directory(. SRC_LIST)

#add_executable(${PROJECT_NAME} ${SRC_LIST})

find_package(Boost REQUIRED COMPONENTS system program_options filesystem)


set(ROTCEN_APP rotation_center)
add_executable(${ROTCEN_APP} rotation_center.cpp ascii_file.cpp)
target_link_libraries(${ROTCEN_APP} ${Boost_LIBRARIES})

message(STATUS ${Boost_LIBRARIES})
