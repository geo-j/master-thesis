#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>


// TODO: think about how the kernel would need to be changed
typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef CGAL::Polygon_2<Kernel>                                         Polygon_2;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Face_handle                                      Face_handle;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                          Ccb_halfedge_circulator;


class Arrangement {
    public:

        /*
        * Arrangement() class constructor
        * Initialises a hard-coded arrangement
        */
        Arrangement() {
            // first create the polygon and its boundaring points
            Polygon_2 P;
            Point_2 p1(0, 4), p2(0, 0), p3(3, 2), p4(4, 0), p5(4, 4), p6(1, 2);
            P.push_back(p1);
            P.push_back(p2);
            P.push_back(p3);
            P.push_back(p4);
            P.push_back(p5);
            P.push_back(p6);

            // then draw the segments between the points
            std::vector<Segment_2> segments;
            segments.push_back(Segment_2(p1, p2));
            segments.push_back(Segment_2(p2, p3));
            segments.push_back(Segment_2(p3, p4));
            segments.push_back(Segment_2(p4, p5));
            segments.push_back(Segment_2(p5, p6));
            segments.push_back(Segment_2(p6, p1));
            CGAL::insert_non_intersecting_curves(arrangement, segments.begin(), segments.end());
        }

        /* print method
        * :param stream f: data stream where the arrangement should be output
        *
        * The format of the output file is:
        * E                     * number of edges
        * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        * p3.x p3.y p4.x p4.y
        * ...
        */
        template<typename stream>
        void print(stream &f) {
            f << this->arrangement.number_of_edges() << std::endl;

            for (auto eit = this->arrangement.edges_begin(); eit != arrangement.edges_end(); ++ eit) {
                f << eit->source()->point() << ' ' << eit->target()->point() << std::endl;
            }
        }
    
    private:
        Arrangement_2 arrangement;
};