#pragma once

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>

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
        is >> prec
           >> order >> nsrcs >> rootLeng >> maxNodeSrcs  >> evalDirect;
    }

    Precision prec;
    int order;
    int nsrcs;
    double rootLeng;
    int maxNodeSrcs;
    bool evalDirect;
};

/*
std::filesystem::path makePath(const Config& config) {
    std::string distStr =
        [&]() -> std::string {
        switch (config.dist) {
            case Dist::UNIFORM:    return "uniform";
            case Dist::GAUSSIAN:   return "gauss";
            case Dist::SPHERE:     return "sphere";
            case Dist::CYLINDER:   return "cyl";
        }
        }();
    std::string cdistStr =
        [&]() -> std::string {
        switch (config.qdist) {
            case ChargeDist::PLUS:  return "plus";
            case ChargeDist::MINUS: return "minus";
            case ChargeDist::DIP:   return "dip";
            case ChargeDist::QUAD:  return "quad";
            case ChargeDist::OCT:   return "oct";
            case ChargeDist::RAND:  return "rand";
        }
        }();

    return
        std::filesystem::path("config") /
        (distStr + "_" + cdistStr + ".txt");
}*/