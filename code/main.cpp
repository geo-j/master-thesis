#include "arrangement.hpp"


int main() {
    std::ifstream f("pentagram.in");
    std::ofstream g("pentagram.out");
    Arrangement arrangement;
    // TODO: maybe leave cin/cout, and pipe files in?
    std::cin >> arrangement;

    // TODO: for now guards are manually added/hard-coded. See how it evolves in the future. Maybe only keep guards for output, since it's the algorithm's job to place them.
    Point_2 p_crown(4, 5.5), p_pentagram1(12.5, 4.5), p_pentagram2(14, 6), p_triangle1(5.5, 10.5), p_triangle2(2.5, 10.5);
    arrangement.add_guard(p_pentagram1);
    arrangement.add_guard(p_pentagram2);

    std::cout << arrangement;
    arrangement.print_guards(std::cout);

    auto visible_vertices = arrangement.compute_visibility();
    // for (auto v : visible_vertices)
        // std::cout << v << std::endl;
    // std::cout << arrangement.is_completely_visible(visible_vertices) << std::endl;

    return 0;
}