#include "arrangement_io.hpp"


int main() {
    std::ifstream f("regularised_visibility_in");
    std::ofstream g("regularised_visibility_out");
    Arrangement arrangement;
    f >> arrangement;

    // TODO: for now guards are manually added/hard-coded. See how it evolves in the future. Maybe only keep guards for output, since it's the algorithm's job to place them.
    Point_2 p1(0.5, 2), p2(3.5, 2.5);
    arrangement.add_guard(p1);
    arrangement.add_guard(p2);

    g << arrangement;
    arrangement.print_guards(g);

    auto visible_vertices = arrangement.compute_visibility();
    for (auto v : visible_vertices)
        std::cout << v << std::endl;
    std::cout << arrangement.is_completely_visible(visible_vertices) << std::endl;

    return 0;
}