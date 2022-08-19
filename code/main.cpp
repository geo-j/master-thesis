#include "arrangement.hpp"


int main() {
    Arrangement arrangement;
    double learning_rate, pull_attraction;

    std::cin >> learning_rate >> pull_attraction >> arrangement;
    arrangement.read_guards(std::cin, learning_rate, pull_attraction);
    // arrangement.print_reflex_intersections();
    std::cout << arrangement;
    // arrangement.print_guards(std::cout);
    arrangement.optimise();

    // auto visible_vertices = arrangement.full_visibility();
    // for (auto v : visible_vertices)
        // std::cout << v << std::endl;
    // std::cout << arrangement.is_completely_visible(visible_vertices) << std::endl;

    return 0;
}