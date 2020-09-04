//
// Created by Matthew Hennefarth on 6/30/20.
//

#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <array>
#include <fstream>
#include <cmath>

//For V/Ang
#define PERM_SPACE 0.0055263495
#define PI 3.1415926535
#define QUAD_TO_ANG (1.602176634 * std::pow(10, -19) * 0.529177 * 0.529177)
#define JOULE_TO_EV (1.602176634 *std::pow(10, -19))

#include "Vector.h"
#include "spdlog/spdlog.h"

struct PointCharge{
    Vector coordinate;
    double charge;
};

std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> out;

    std::string::size_type start = 0;
    for (std::string::size_type i = 0; i < str.size(); i++) {
        if (str[i] == delim && start != i) {
            out.emplace_back(str.substr(start, i - start));
            start = i + 1;
        }
        else if(str[i] == delim && start == i){
            start++;
        }
    }
    out.push_back(str.substr(start, str.size()));

    //Now we filter out
    for(std::vector<std::string>::size_type i = 0; i < out.size(); i++) {
        if(out[i].empty()){
            out.erase(out.begin() + static_cast<long>(i));
        }
    }

    return out;
}

class System{
    public:
        System(const std::string& proteinFile, const std::string& optionsFile);

        [[nodiscard]] constexpr Vector electricField(const Vector& position) const;
        [[nodiscard]] std::array<std::array<double, 3>, 3>  electricFieldGradient(const Vector& position) const;
        constexpr void addPointCharge(const PointCharge& pc);
        constexpr void addPosition(const Vector& pos);
        void calculate(const bool efg) const;

    private:
        std::vector<PointCharge> _pointCharges;
        double _dielectric;
        std::vector<Vector> _points;
        std::array<std::array<double, 3>, 3> _quadrapole;
};


System::System(const std::string& proteinFile, const std::string& optionsFile) : _dielectric(1), _quadrapole(){
    std::fstream inFile(proteinFile);
    std::string line;
    if(inFile.is_open()){
        while(std::getline(inFile, line)){
            if(line.substr(0, 4) == "ATOM" || line.substr(0,6) == "HETATM"){
                PointCharge pc = {{std::stod(line.substr(31, 8)),
                                            std::stod(line.substr(39,8)),
                                            std::stod(line.substr(47, 8))},
                                  std::stod(line.substr(55,8))};
                this->addPointCharge(pc);
            }
        }
        inFile.close();
    }
    else{
        SPDLOG_ERROR("Could not open protein file to read!");
    }

    inFile.open(optionsFile);
    if(inFile.is_open()){
        while(std::getline(inFile, line)){
            if (line.substr(0,2) == "xx"){
                _quadrapole[0][0] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if (line.substr(0,2) == "yy"){
                _quadrapole[1][1] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if(line.substr(0,2) == "zz"){
                _quadrapole[2][2] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if(line.substr(0,2) == "xy"){
                _quadrapole[0][1] = std::stod(line.substr(2)) * QUAD_TO_ANG;
                _quadrapole[1][0] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if(line.substr(0,2) == "xz"){
                _quadrapole[0][2] = std::stod(line.substr(2)) * QUAD_TO_ANG;
                _quadrapole[2][0] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if(line.substr(0,2) == "yz"){
                _quadrapole[1][2] = std::stod(line.substr(2)) * QUAD_TO_ANG;
                _quadrapole[2][1] = std::stod(line.substr(2)) * QUAD_TO_ANG;
            }
            else if(line.substr(0,10) == "dielectric"){
                _dielectric = std::stod(line.substr(10));
            }
            else if (line.size() > 2){
                auto pos = split(line, ' ');
                if (pos.size() == 3){
                    addPosition({std::stod(pos[0]), std::stod(pos[1]), std::stod(pos[2])});
                }
                else{
                    SPDLOG_WARN("Skipping position: {}", line);
                }
            }
        }
        inFile.close();
    }
    else{
        SPDLOG_ERROR("Could not open options file to read!");
    }

}

constexpr void System::addPointCharge(const PointCharge& pc){
    _pointCharges.push_back(pc);
}

constexpr void System::addPosition(const Vector& pos){
    _points.push_back(pos);
}

constexpr Vector System::electricField(const Vector &position) const {
    Vector result;

    for (const auto& pc : _pointCharges){
        Vector d = position - pc.coordinate;
        double dNorm = d.norm();
        result += ((pc.charge * d) / (dNorm*dNorm*dNorm));
    }

    return result/(4.0 * PI * PERM_SPACE * _dielectric);
}

std::array<std::array<double, 3>, 3> System::electricFieldGradient (const Vector &position) const {
    std::array<std::array<double, 3>, 3> result{};

    for(std::array<std::array<double, 3>, 3>::size_type i = 0; i < 3; i++){
        for(std::array<std::array<double, 3>, 3>::size_type j = 0; j < 3; j++){
            for(const auto& pc : _pointCharges){
                Vector d = position - pc.coordinate;
                double dNorm = d.norm();
                result[i][j] += pc.charge * -3.0 * (d[i]) * (d[j])/(dNorm*dNorm*dNorm*dNorm*dNorm);
                if (i == j){
                    result[i][j] += (pc.charge/(dNorm*dNorm*dNorm));
                }
            }
            // TM in atomic units electron charge * bohr^2, so should convert to C*Ang^2...maybe easier than this
            result[i][j] /= (4.0 * PI * PERM_SPACE * _dielectric);
        }
    }
    // This is in V/Ang^2, or J/(C*Ang^2), so make sure the quadrapole is in correct units
    return result;
}

void System::calculate (const bool efg) const {

    SPDLOG_INFO("Calculating electric field at points (V/Ang)");
    SPDLOG_INFO("[x, y, z] [Ex, Ey, Ez] Mag");
    for (const auto& point: _points){
        auto result = this->electricField(point);
        SPDLOG_INFO("{} {} {}", point, result, result.norm());
    }
    if (efg){
        SPDLOG_INFO("Calculating electric field gradient! (V/Ang^2)");
        for(const auto& point: _points){
            auto result = this->electricFieldGradient(point);

            SPDLOG_INFO("Position: {}", point);
            std::string output;
            for(const auto& r : result){
                output = "";
                for(const auto& v: r){
                    output += std::to_string(v) + ", ";
                }
                SPDLOG_INFO(output);
            }

            double energy = 0;
            for(size_t i = 0; i < 3; i++){
                for(size_t j = 0; j < 3; j++){
                    energy += _quadrapole[i][j]*result[j][i];
                }
            }
            energy *= (-0.5);

            energy /= JOULE_TO_EV;

            SPDLOG_INFO("Energy: {} eV", energy);

        }
    }
}

#endif //SYSTEM_H
