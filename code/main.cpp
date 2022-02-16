#include "arrangement.hpp"


int main() {
    Arrangement arrangement;
    // TODO: piping multiple files in doesn't work. probably a fish problem?
    std::cin >> arrangement;
    arrangement.read_guards(std::cin);

    // std::cout << arrangement;
    // arrangement.print_guards(std::cout);

    auto visible_vertices = arrangement.compute_full_visibility();
    // for (auto v : visible_vertices)
        // std::cout << v << std::endl;
    std::cout << arrangement.is_completely_visible(visible_vertices) << std::endl;

    return 0;
}