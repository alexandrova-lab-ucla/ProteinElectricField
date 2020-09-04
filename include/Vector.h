//
// Created by Matthew Hennefarth on 6/30/20.
//

#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <cmath>
#include <iostream>

#include "spdlog/fmt/ostr.h"

#define VECTOR_SIZE 3

class Vector{
    public:
        typedef std::array<double,3>::size_type size_t;
        typedef size_t size_type;

        constexpr Vector();
        constexpr Vector(const double& x, const double& y, const double& z);
        constexpr Vector(const std::array<double,VECTOR_SIZE>& position);
        constexpr Vector(const Vector& vec);

        [[nodiscard]] double norm() const;
        [[nodiscard]] constexpr double dot(const Vector& vec) const;

        constexpr Vector& operator=(const Vector &vec);
        constexpr bool operator==(const Vector& vec) const;
        constexpr Vector operator+(const Vector &vec) const;
        constexpr Vector& operator+=(const Vector &vec);
        constexpr Vector operator-(const Vector &vec) const;
        constexpr Vector& operator-=(const Vector &vec);
        constexpr Vector operator*(const Vector &vec) const;
        constexpr Vector& operator*=(const Vector &vec);
        constexpr Vector operator/(const Vector &vec) const;
        constexpr Vector& operator/=(const Vector &vec);

        constexpr Vector operator+(const double& val) const;
        constexpr Vector& operator+=(const double& val);
        constexpr Vector operator*(const double& val) const;
        constexpr Vector& operator*=(const double& val);
        constexpr Vector operator-(const double& val) const;
        constexpr Vector& operator-=(const double& val);
        constexpr Vector operator/(const double& val) const;
        constexpr Vector& operator/=(const double& val);

        constexpr const double& operator[](const size_t index) const;
        constexpr double& operator[](const size_t index);

        constexpr friend Vector operator+(const double& val, const Vector &vec);
        constexpr friend Vector operator*(const double& val, const Vector &vec);
        constexpr friend Vector operator-(const double& val, const Vector &vec);
        constexpr friend Vector operator/(const double& val, const Vector &vec);

        //For spdlog library
        template<typename OStream>
        inline friend OStream& operator<<(OStream &os, const Vector& vec);

        //for normal std::cout
        inline friend std::ostream& operator<<(std::ostream &os, const Vector& vec){return operator<<<std::ostream>(os, vec);}

    private:
        std::array<double,VECTOR_SIZE> _values;

};

constexpr Vector::Vector(): _values({0.0, 0.0, 0.0}){}

constexpr Vector::Vector(const double& x, const double& y, const double& z) : _values({x, y, z}){}

constexpr Vector::Vector(const std::array<double,VECTOR_SIZE>& position) : _values(position){}

constexpr Vector::Vector(const Vector& vec) : _values(vec._values){}

double Vector::norm() const{
    return sqrt(this->dot(*this));
}

constexpr double Vector::dot(const Vector& vec) const{
    double sum = 0.0;

    for(size_type i = 0; i < VECTOR_SIZE; i++){
        sum += (_values[i] * vec._values[i]);
    }

    return sum;
}

constexpr Vector& Vector::operator=(const Vector &vec){
    if (this != &vec){
        _values = vec._values;
    }

    return *this;
}

constexpr bool Vector::operator==(const Vector& vec) const{
    return _values == vec._values;
}

constexpr Vector Vector::operator+(const Vector &vec) const{
    Vector result;

    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i] + vec._values[i];
    }

    return result;
}

constexpr Vector& Vector::operator+=(const Vector &vec){
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        _values[i] += vec._values[i];
    }

    return *this;
}

constexpr Vector Vector::operator-(const Vector&vec) const{
    return (*this + (vec * -1.0));
}

constexpr Vector& Vector::operator-=(const Vector& vec) {
    return (*this += (vec * -1.0));
}

constexpr Vector Vector::operator*(const Vector &vec) const{
    Vector result;
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i] * vec._values[i];
    }

    return result;
}

constexpr Vector& Vector::operator*=(const Vector &vec){
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        _values[i] *= vec._values[i];
    }
    return *this;
}

constexpr Vector Vector::operator/(const Vector &vec) const {
    return (*this * (1.0 / vec));
}

constexpr Vector& Vector::operator/=(const Vector &vec){
    return (*this *= (1.0 / vec));
}

constexpr Vector Vector::operator+(const double& val) const{
    Vector result;
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i] + val;
    }

    return result;
}

constexpr Vector& Vector::operator+=(const double& val) {
    for(auto& value: _values){
        value += val;
    }
    return *this;
}

constexpr Vector Vector::operator*(const double& val) const{
    Vector result;
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i]*val;
    }

    return result;
}

constexpr Vector& Vector::operator*=(const double& val){
    for(auto& value: _values){
        value *= val;
    }
    return *this;
}

constexpr Vector Vector::operator-(const double& val) const{
    Vector result;
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i] - val;
    }
    return result;
}

constexpr Vector& Vector::operator-=(const double& val){
    for(auto& value : _values){
        value -= val;
    }
    return *this;
}

constexpr Vector Vector::operator/(const double& val) const{
    Vector result;
    for(size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = _values[i] / val;
    }
    return result;
}

constexpr Vector& Vector::operator/=(const double& val){
    for(auto& value: _values){
        value /= val;
    }
    return *this;
}

constexpr Vector operator+(const double& val, const Vector &vec){
    return (vec + val);
}

constexpr Vector operator*(const double& val, const Vector &vec){
    return (vec * val);
}

constexpr Vector operator-(const double& val, const Vector &vec){
    return (val + (vec * -1.0));
}

constexpr Vector operator/(const double& val, const Vector &vec){
    Vector result;
    for (size_t i = 0; i < VECTOR_SIZE; i++){
        result._values[i] = (val/vec._values[i]);
    }

    return result;
}

constexpr const double& Vector::operator[](const size_t index) const{
    return _values[index];
}

constexpr double& Vector::operator[](const size_t index){
    return const_cast<double &>(static_cast<const Vector&>(*this)[index]);
}

template<typename OStream>
OStream &operator<< (OStream &os, const Vector &vec) {
    os << "[ ";
    for(size_t i = 0; i < VECTOR_SIZE-1; i++){
        os << vec._values[i] << ", ";
    }
    os << vec._values[VECTOR_SIZE-1] << " ]";

    return os;
}



#endif //VECTOR_H
