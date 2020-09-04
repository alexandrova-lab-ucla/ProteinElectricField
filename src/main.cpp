//
// Created by Matthew Hennefarth on 6/30/20.
//

#include <filesystem>

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "spdlog/spdlog.h"
#include "cxxopts.hpp"
#include "System.h"

int main(int argc, char** argv){

    cxxopts::Options options("ProteinElectricField",
        "Calculates the Electric Field and Gradient of the Electric Field");

    options.add_options()
        ("p,protein", "PDB", cxxopts::value<std::string>())
        ("o,options", "Option file", cxxopts::value<std::string>())
        ("h,help", "Print usage")
        ;

    SPDLOG_DEBUG("Parsing options");
    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        SPDLOG_WARN(options.help());
        exit(0);
    }

    std::string proteinFile;
    if (result.count("protein")){
        proteinFile = result["protein"].as<std::string>();
        if(!std::filesystem::exists(proteinFile)){
            SPDLOG_ERROR("{} file does not exist", proteinFile);
            exit(1);
        }
    }
    else{
        SPDLOG_ERROR("No protein file specified");
        SPDLOG_WARN(options.help());
        exit(1);
    }

    std::string optionsFile;
    if (result.count("options")){
        optionsFile = result["options"].as<std::string>();
        if(!std::filesystem::exists(optionsFile)){
            SPDLOG_ERROR("{} file does not exist", optionsFile);
            exit(1);
        }
    }
    else{
        SPDLOG_ERROR("No options file specified");
        SPDLOG_WARN(options.help());
        exit(1);
    }

    System s(proteinFile, optionsFile);
    s.calculate(true);


    return 0;
}