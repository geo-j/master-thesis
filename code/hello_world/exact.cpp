// If you must ensure that your numbers get interpreted at their full precision you can use a CGAL kernel that performs exact predicates and extract constructions.
// you will have floating point numbers that are "exact", in the sense that they were computed by some application or obtained from a sensor. They are not the string "0.1" or computed on the fly as "1.0/10.0", but a full precision floating point number. If they are input to an algorithm that makes no constructions, you can use a kernel that provides exact predicates but inexact constructions

#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <sstream>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
int main()
{
  //  the coordinates you see as text get turned into floating point numbers. When they are then turned into arbitrary precision rationals, they exactly represent the floating point number, but not the text
  Point_2 p(0, 0.3), q, r(2, 0.9);
  {
    q  = Point_2(1, 0.6);
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");
  }
  // reading numbers from a file. The arbitrary precision rationals are then directly constructed from a string so that they represent exactly the text
  {
    std::istringstream input("0 0.3   1 0.6   2 0.9");
    input >> p >> q >> r;
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");
  }
  // constructions as midpoint constructions are exact, just as the name of the kernel type suggests
  {
    q = CGAL::midpoint(p,r);
    std::cout << (CGAL::collinear(p,q,r) ? "collinear\n" : "not collinear\n");
  }
  return 0;
}