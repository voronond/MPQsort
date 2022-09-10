#include <cxxopts.hpp>
#include <iostream>
#include <random>
#include <string>

auto main(int argc, char** argv) -> int {
    cxxopts::Options options(*argv, "A program to generate input for benchmarking.");

    int32_t num = 10;
    int32_t from = 0, to = 10;
    std::string distribution;

    // clang-format off
    options.add_options()
    ("h,help", "Show help")

    ("n,num", "How many numbers to generate", cxxopts::value<int32_t>(num))
    ("f,from", "Smallest number in a sequence", cxxopts::value<int32_t>(from))
    ("t,to", "The biggest number in a sequence", cxxopts::value<int32_t>(to))
    ("d,distribution", "Selected distribution [unif]", cxxopts::value<std::string>(distribution)->default_value("unif"))
  ;
    // clang-format on

    auto result = options.parse(argc, argv);

    if (result["help"].as<bool>()) {
        std::cout << options.help() << std::endl;
    }

    // Seed with a real radom value
    std::random_device random_dev;
    std::mt19937 en(random_dev());

    if (distribution == "unif") {
        std::uniform_int_distribution<int32_t> uniform_dist(from, to);

        std::cout << num << std::endl;
        for (int32_t i = 0; i < num; ++i) std::cout << uniform_dist(en) << " ";
        std::cout << std::endl;
    }

    return 0;
}
