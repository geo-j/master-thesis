#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <stack>


typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
// typedef typename CGAL::Geometry_traits_2::Kernel                     Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef typename CGAL::Arrangement_2<Traits_2>                          Arrangement_2;
typedef typename Arrangement_2::Geometry_traits_2                       Geometry_traits_2;

typedef Kernel::Point_2                                                 Point_2;



Arrangement_2 compute_visibility_polygon(Point_2 p) {

}

int main() {


    return 0;
}