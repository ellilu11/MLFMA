#include "config.h"

std::string getModeStr(Mode mode) {
    return [&]() {
        switch (mode) {
            case Mode::FMM:    return "FMM";
            case Mode::DIR:    return "DIRECT";
            case Mode::FMMDIR: return "FMM+DIRECT";
        };
        } ();
}

std::string getIEStr(IE ie) {
    return [&]() {
        switch (ie) {
            case IE::EFIE: return "EFIE";
            case IE::MFIE: return "MFIE";
            case IE::CFIE: return "CFIE";
        };
        } ();
}

size_t getNumQuads(Precision prec) {
    return [&]() {
        switch (prec) {
            case Precision::VERYLOW:  return 1;
            case Precision::LOW:      return 3;
            case Precision::MEDLOW:   return 4;
            case Precision::MEDIUM:   return 6;
            case Precision::MEDHIGH:  return 7;
            case Precision::HIGH:     return 12;
            case Precision::VERYHIGH: return 13;
        };
        } ();
}

Config::Config(const std::filesystem::path& fileName) {
    auto configVals = importConfig(fileName);

    mode = static_cast<Mode>(configVals[0]);
    quadPrec = static_cast<Precision>(configVals[1]);
    alpha = configVals[2];
    k = configVals[3];
    leps = configVals[4];
    nsrcs = static_cast<int>(configVals[5]);

    maxNodeSrcs = static_cast<int>(configVals[6]);
    digits = static_cast<int>(configVals[7]);
    interpOrder = static_cast<int>(configVals[8]);
    overInterp = configVals[9];

    epsIter = configVals[10];
    maxIter = static_cast<int>(configVals[11]);
    iluTol = configVals[12];
    iluFactor = static_cast<int>(configVals[13]);

    numThreads = configVals[14];
    pivotLvl = configVals[15];

    if (alpha < 0.0 || alpha > 1.0)
        throw std::runtime_error("Alpha must be in [0,1]");
    ie = (alpha == 0.0) ? IE::MFIE : (alpha == 1.0) ? IE::EFIE : IE::CFIE;
    C_efie = -Phys::eta * alpha * iu * k;
    C_mfie = Phys::eta * (1.0 - alpha);

    eleng = configVals[16];
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << eleng;
    lengStr = ss.str();

    std::cout << " *************************** \n";
    std::cout << " ***** Helmholtz-MLFMM ***** \n";
    std::cout << " *************************** \n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "   Mode:            " << getModeStr(mode) << '\n';
    std::cout << "   IE (alpha):      " << getIEStr(ie) << " (" << alpha << ")\n";
    std::cout << "   Wavenumber:      " << k << " /m\n";
    std::cout << "   Length tol:      " << leps << " m\n";
    std::cout << "   Max in node:     " << maxNodeSrcs << '\n';
    std::cout << "   Digit precision: " << digits << '\n';
    std::cout << "   Interp order:    " << interpOrder << '\n';
    std::cout << "   Overinterp:      " << overInterp << '\n';
    std::cout << "   Tri quad rule:   " << getNumQuads(quadPrec) << "-point\n";
    std::cout << "   ILU drop tol:    " << iluTol << '\n';
    std::cout << "   ILU fill factor: " << iluFactor << "\n";
    std::cout << "   Num threads:     " << numThreads << "\n\n";
}

std::array<double, 17> Config::importConfig(const std::filesystem::path& path) {
    std::array<double, 17> arr;
    std::ifstream file(path);
    std::string line;
    size_t i = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        // Skip empty and non-numeric lines
        if (line.empty() || !std::any_of(line.begin(), line.end(),
            [](unsigned char c) { return std::isdigit(c); }))
            continue;

        char ch = '\0';
        getDigit(iss, ch);

        while (iss >> arr[i]) {
            getDigit(iss, ch);
            ++i;
        }
    }

    return arr;
}