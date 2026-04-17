#pragma once

#include "fileio.h"
#include "phys.h"

enum class Mode { FMM, DIR, FMMDIR };

enum class IE { EFIE, MFIE, CFIE };

enum class Precision { VERYLOW, LOW, MEDLOW, MEDIUM, MEDHIGH, HIGH, VERYHIGH};

struct Config {
    Config() = default;
    
    Config(const std::filesystem::path&);
    
    static std::array<double, 16> importConfig(const std::filesystem::path&);

    // General
    Mode mode;          // FMM, direct, or FMM+direct
    Precision quadPrec; // Triangle quadrature precision
    double alpha;       // CFIE parameter
    double k;           // Wavenumber
    double leps;        // Length error tolerance for nearfield

    // FMM
    int maxNodeSrcs;    // Max # sources in leaf
    int digits;         // Digit precision
    int interpOrder;    // Interpolation order
    double overInterp;  // Translation operator angular oversampling factor

    // GMRES
    double epsIter;     // Error tolerance for GMRES
    int maxIter;        // Max GMRES iterations (0 for one MVM)
    double iluTol;      // Drop tolerance for ILU preconditioner
    int iluFactor;      // Fill factor for ILU preconditioner

    // OpenMP
    int numThreads;     // Number of OpenMP threads
    int pivotLvl;       // Level at which to switch from parallel nodes to parallel angles

    // CFIE (from alpha)
    IE ie;              // Integral equation type
    cmplx C_efie;       // -eta * alpha * ik
    cmplx C_mfie;       // eta * (1-alpha)

    // File I/O
    int nsrcs;          
    double eleng;
    std::string lengStr;
};