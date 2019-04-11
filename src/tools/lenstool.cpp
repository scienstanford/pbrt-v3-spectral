//
// imgtool.cpp
//
// Various useful operations on images.
//

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include "floatfile.h"
#include "fileutil.h"
#include "imageio.h"
#include "pbrt.h"
#include "spectrum.h"
#include "parallel.h"
#include "ext/json.hpp"
#include <glog/logging.h>
using namespace nlohmann;
using namespace pbrt;

static void usage(const char *msg = nullptr, ...) {
    if (msg) {
        va_list args;
        va_start(args, msg);
        fprintf(stderr, "lenstool: ");
        vfprintf(stderr, msg, args);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, R"(usage: lenstool <command> [options] <filenames...>

commands: convert

convert options:
    --inputscale <n>   Input units per mm (which are used in the output). Default: 1.0
)");
    exit(1);
}

int convert(int argc, char *argv[]) {
    float scale = 1.f;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, double> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;

        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, atof(value)};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-') break;
        std::pair<std::string, double> arg = parseArg();
        if (std::get<0>(arg) == "inputscale") {
            scale = std::get<1>(arg);
            if (scale == 0) usage("--inputscale value must be non-zero");
        } else
            usage();
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"convert\"");
    else if (i >= argc)
        usage("missing filenames for \"convert\"");

    const char *inFilename = argv[i], *outFilename = argv[i + 1];



    auto endsWith = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
    };
    std::string lensFile = inFilename;
    if (endsWith(lensFile, ".dat")) {
        // Load element data from lens description file
        std::vector<Float> lensData;
        std::string name;
        std::string description;
        if (!ReadFloatFile(lensFile.c_str(), &lensData)) {
            Error("Error reading lens specification file \"%s\".",
                lensFile.c_str());
            return -1;
        } else {
            // TODO: More robust extraction of names, descriptions
            std::ifstream infile(lensFile);
            if (infile.good()) {
                std::string nameLine, descriptionLine;
                getline(infile, nameLine);
                if (nameLine[0] == '#') {
                    name = nameLine.substr(1);
                }
                getline(infile, descriptionLine);
                while (descriptionLine[0] == '#') {
                    description += descriptionLine.substr(1);
                    getline(infile, descriptionLine);
                } 
            }
        }
        if (lensData.size() % 4 != 0) {
            // Trisha: If the size has an extra value, it's possible this lens type was meant for pbrt-v2-spectral and has an extra focal length value at the top. In this case, let's automatically convert it by removing this extra value.
            if (lensData.size() % 4 == 1) {
                Warning("Extra value in lens specification file, this lens file may be for pbrt-v2-spectral. Removing extra value to make it compatible with pbrt-v3-spectral...");
                lensData.erase(lensData.begin());
            }
            else {
                Error(
                    "Excess values in lens specification file \"%s\"; "
                    "must be multiple-of-four values, read %d.",
                    lensFile.c_str(), (int)lensData.size());
                return -1;
            }
        }
        
        std::vector<json> surfaces;
        for (int i = 0; i < lensData.size(); i += 4) {
            json jsurf;
            jsurf["radius"] = lensData[i+0];
            jsurf["thickness"] = lensData[i+1];
            jsurf["ior"] = lensData[i+2];
            jsurf["semi_aperture"] = lensData[i+3] / 2.0;
            surfaces.push_back(jsurf);
        }

        json j;
        j["name"] = name;
        j["description"] = description;
        j["surfaces"] = surfaces;
        std::ofstream outfile(outFilename);
        outfile << std::setw(4) << j << std::endl;

        // Covert to outfile here
        printf("Input file: %s, Ouput file: %s; %zd surfaces\n", inFilename, outFilename, surfaces.size());
        
    } else {
        Error(
            "Input to lenstool conver must be a .dat file, but given \"%s\"; ",
            lensFile.c_str());
    }

    return 0;
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 1; // Warning and above.

    if (argc < 2) usage();

    if (!strcmp(argv[1], "convert"))
        return convert(argc - 2, argv + 2);
    else
        usage("unknown command \"%s\"", argv[1]);

    return 0;
}
