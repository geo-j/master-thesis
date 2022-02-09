#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>

#include "arrangement_io.hpp"


int main() {
    std::ifstream f("arrangement_main_in");
    std::ofstream g("arrangement_main_out");
    Arrangement arrangement;
    f >> arrangement;

    // TODO: for now guards are manually added/hard-coded. See how it evolves in the future
    Point_2 p(0.5, 2);
    arrangement.add_guard(p);

    g << arrangement;
    arrangement.print_guards(g);

    return 0;
}