#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>

#include "arrangement_io.hpp"


int main() {
    std::ofstream g("arrangement_main_out");
    Arrangement arrangement;

    arrangement.print(g);

    return 0;
}