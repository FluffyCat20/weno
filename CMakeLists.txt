cmake_minimum_required(VERSION 3.5)

project(weno LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (nlohmann_json_DIR json_parser/)
find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(weno weno.hpp WenoSimpleTransfer.hpp printer.hpp boundary_conditions.hpp init.hpp data.hpp main.cpp)

target_link_libraries(weno nlohmann_json::nlohmann_json)
