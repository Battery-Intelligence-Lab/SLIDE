cmake_minimum_required(VERSION 3.16)
add_library(cooling "")

target_sources(cooling
  PRIVATE
    CoolSystem.cpp
    CoolSystem_open.cpp
    CoolSystem_HVAC.cpp
    data_storage/CoolSystemData.cpp
    data_storage/CoolSystem_HVACData.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/CoolSystem.hpp
    ${CMAKE_CURRENT_LIST_DIR}/CoolSystem_open.hpp
    ${CMAKE_CURRENT_LIST_DIR}/CoolSystem_HVAC.hpp
    ${CMAKE_CURRENT_LIST_DIR}/data_storage/cool_data.hpp
    ${CMAKE_CURRENT_LIST_DIR}/data_storage/CoolSystemData.hpp
    ${CMAKE_CURRENT_LIST_DIR}/data_storage/CoolSystem_HVACData.hpp
    ${CMAKE_CURRENT_LIST_DIR}/data_storage/cool_data.hpp
  )

target_include_directories(cooling 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )




