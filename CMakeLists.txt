project(BTA_service_tools)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

#add_executable(${PROJECT_NAME} ${SRC_LIST})

find_package(Boost REQUIRED COMPONENTS program_options)


set(ROTCEN_APP rotation_center)
add_executable(${ROTCEN_APP} ./rotation_center.cpp)
target_link_libraries(${ROTCEN_APP} ${Boost_LIBRARIES})

message(STATUS ${Boost_LIBRARIES})
