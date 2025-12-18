#pragma once

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>

enum class Mode {
    READ,
    WRITE
};

enum class Dist {
    UNIFORM,
    GAUSSIAN,
    SPHERE,
    CYLINDER
};

enum class Precision {
    LOW,
    MEDIUM,
    HIGH
};

void getDigit(std::istringstream& iss, char ch) {
    while (iss.get(ch)) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            iss.unget();
            break;
        }
    }
}

template <typename T>
std::ifstream& operator>>(std::ifstream& is, T& val) {
    std::string line;
    if (std::getline(is, line)) {
        std::istringstream iss(line);

        char ch = '\0';
        getDigit(iss, ch);

        if constexpr (std::is_enum_v<T>) {
            typename std::underlying_type<T>::type eval;

            while (iss >> eval) {
                val = static_cast<T>(eval);
                getDigit(iss, ch);
            }

        } else
            while (iss >> val)
                getDigit(iss, ch);
    }

    return is;
}

struct Config {
    Config() = default;
    Config(const std::string& fileName) {
        std::ifstream is(fileName);
        is >> mode >> pdist >> quadPrec // enums first
            >> digits >> interpOrder
            >> rootLeng >> maxNodeSrcs  >> evalDirect
            >> nsrcs;
    }

    Precision quadPrec;
    int digits;
    int interpOrder;
    double rootLeng;
    int maxNodeSrcs;
    bool evalDirect;

    // For point sources only
    Mode mode;
    Dist pdist;
    int nsrcs;
};