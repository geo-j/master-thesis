#include "arrangement_io.hpp"


int main() {
    std::ifstream f("regularised_visibility_in");
    std::ofstream g("regularised_visibility_out");
    Arrangement arrangement;
    f >> arrangement;

    // TODO: for now guards are manually added/hard-coded. See how it evolves in the future. Maybe only keep guards for output, since it's the algorithm's job to place them.
    Point_2 p(0.5, 2);
    arrangement.add_guard(p);

    g << arrangement;
    arrangement.print_guards(g);

    return 0;
}