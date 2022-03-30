#include "arrangement.hpp"


int main() {
    Arrangement arrangement;
    float learning_rate;

    std::cin >> learning_rate >> arrangement;
    arrangement.read_guards(std::cin);
    // arrangement.print_reflex_intersections();
    std::cout << arrangement;
    // arrangement.print_guards(std::cout);
    arrangement.optimise(learning_rate);

    // auto visible_vertices = arrangement.compute_full_visibility();
    // for (auto v : visible_vertices)
        // std::cout << v << std::endl;
    // std::cout << arrangement.is_completely_visible(visible_vertices) << std::endl;

    return 0;
}