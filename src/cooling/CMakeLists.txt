cmake_minimum_required(VERSION 3.17)
add_library(cooling STATIC)

target_sources(cooling
  PRIVATE
    CoolSystem.cpp
    CoolSystem_open.cpp
    CoolSystem_HVAC.cpp
    data_storage/CoolSystemData.cpp
    data_storage/CoolSystem_HVACData.cpp
  PUBLIC
    CoolSystem.hpp
    CoolSystem_open.hpp
    CoolSystem_HVAC.hpp
    data_storage/cool_data.hpp
    data_storage/CoolSystemData.hpp
    data_storage/CoolSystem_HVACData.hpp
    data_storage/cool_data.hpp
  )

target_include_directories(cooling PUBLIC .)




