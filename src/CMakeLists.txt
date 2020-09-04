set(SOURCE_FILES main.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)

#---------------------------------------------------[Main Executable]---------------------------------------------------
add_executable(ProteinElectricField ${SOURCE_FILES})
if(CMAKE_BUILD_TYPE MATCHES Debug)
    enable_sanitizers(ProteinElectricField)
    enable_warnings(ProteinElectricField)
endif()
target_link_libraries(ProteinElectricField spdlog)
target_link_libraries(ProteinElectricField cxxopts)