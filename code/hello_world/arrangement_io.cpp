// Using the arrangement I/O operators.
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <fstream>
#include "point_location_utils.h"
typedef CGAL::Cartesian<CGAL::Exact_rational>         Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;
  construct_segments_arr (arr);
  std::cout << "Writing an arrangement of size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;
  // Write the arrangement to a file.
  std::ofstream    out_file ("arr_ex_io.dat"), f("arrangement");
  out_file << arr;
  out_file.close();
  // Read the arrangement from the file.
  Arrangement_2    arr2;
  std::ifstream    in_file ("arr_ex_io.dat");
  in_file >> arr2;
  in_file.close();
  // std::cout << "Read an arrangement of size:" << std::endl
  //           << "   V = " << arr2.number_of_vertices()
  //           << ",  E = " << arr2.number_of_edges()
  //           << ",  F = " << arr2.number_of_faces() << std::endl;
  Arrangement_2::Vertex_iterator            vit;
  Arrangement_2::Edge_iterator              eit;
  std::cout << "The arrangement vertices: " << std::endl;
  f << arr2.number_of_vertices() << std::endl;
  for (vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++ vit) {
    // std::cout << '(' << vit->point() << ") - " << std::endl;
    f << vit->point() << std::endl;
  }
  
  f << arr2.number_of_edges() << std::endl;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // std::cout << "edge from " << eit->source()->point() << " to "  << eit->target()->point() << std::endl;
    f << eit->source()->point() << ' ' << eit->target()->point() << std::endl;
  }
  return (0);
}