#include <fstream>
#include <algorithm>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>

#include "utils.hpp"


// TODO: think about how the kernel would need to be changed
typedef CGAL::Exact_predicates_exact_constructions_kernel                   Kernel;
typedef CGAL::Polygon_2<Kernel>                                             Polygon_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Segment_2                                                   Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                                  Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                       Arrangement_2;
typedef Arrangement_2::Face_handle                                          Face_handle;
typedef Arrangement_2::Edge_const_iterator                                  Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                              Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>    NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_true>  RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>              TEV;


// TODO: make function that checks whether the whole polygon is visible
class Arrangement {
    public:
        /*
        * Arrangement() class constructor
        */
        // TODO: should it also have hardcoded stuff?
        Arrangement() {
            // Point_2 p1(0, 4), p2(0, 0), p3(3, 2), p4(4, 0), p5(4, 4), p6(1, 2), p(0.5, 2);
            // // draw the segments between the points
            // std::vector<Segment_2> segments;
            // segments.push_back(Segment_2(p1, p2));
            // segments.push_back(Segment_2(p2, p3));
            // segments.push_back(Segment_2(p3, p4));
            // segments.push_back(Segment_2(p4, p5));
            // segments.push_back(Segment_2(p5, p6));
            // segments.push_back(Segment_2(p6, p1));
            // CGAL::insert_non_intersecting_curves(arrangement, seif (it == a.boundary_vertices.end())
        }

       /* overloaded input operator
        * The format of the input file is:
        * E                     * number of edges
        * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        * p3.x p3.y p4.x p4.y
        * ...
        */
        friend std::istream &operator>>(std::istream &f, Arrangement &a) {
            std::size_t E, x1, y1, x2, y2;
            std::vector<Segment_2> segments;

            f >> E;

            for (auto i = 0; i < E; i ++) {
                f >> x1 >> y1 >> x2 >> y2;
                Point_2 p1(x1, y1), p2(x2, y2);

                segments.push_back(Segment_2(p1, p2));
                push_back_unique(a.boundary_vertices, p1);
                push_back_unique(a.boundary_vertices, p2);
            }

            CGAL::insert_non_intersecting_curves(a.input_arrangement, segments.begin(), segments.end());

            return f;
        }

        /* overloaded output operator
        * The format of the output file is:
        * E                     * number of edges
        * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        * p3.x p3.y p4.x p4.y
        * ...
        */
        friend std::ostream &operator<<(std::ostream &f, const Arrangement &a) {
            f << a.input_arrangement.number_of_edges() << std::endl;

            for (auto eit = a.input_arrangement.edges_begin(); eit != a.input_arrangement.edges_end(); ++ eit) {
                f << eit->source()->point() << ' ' << eit->target()->point() << std::endl;
            }

            return f;
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
            f << this->input_arrangement.number_of_edges() << std::endl;

            for (auto eit = this->input_arrangement.edges_begin(); eit != this->input_arrangement.edges_end(); ++ eit) {
                f << eit->source()->point() << ' ' << eit->target()->point() << std::endl;
            }
        }

        // TODO: think about whether I should also print guards in the arrangement

        /* print_guards method
        * :param stream f: data stream where the arrangement should be output
        *
        * The format of the output file is:
        * IV                    * number of isolated vertices (intended to be guards)
        * q1.x q1.y             * isolated vertex with coordinates q1(x, y) separated by spaces
        * q2.x q2.y
        * ...
        */
        template<typename stream>
        void print_guards(stream &f) {
            f << guards.size() << std::endl;

            for (auto guard : guards) {
                f << guard.x() << ' ' << guard.y() << std::endl;
            }
        }

        /* add_guard method
        *  :param Point_2 p: point that would guard the arrangement
        */
        void add_guard(const Point_2 p) {
            this->guards.push_back(p);
        }

        // TODO: either take the visible vertices vector as input, or use the guards
        bool is_completely_visible(std::vector<Point_2> visible_points) {
            std::size_t count = 0;

            for (auto point : visible_points) {
                auto it = std::find(this->boundary_vertices.begin(), this->boundary_vertices.end(), point);

                if (it != this->boundary_vertices.end())
                    count ++;
            }

            return count == boundary_vertices.size();
        }

        /* compute_visibility method
        *  :param Point_2 p                :point whose visibility region needs to be computed
        *  :return Arrangement_2           :arrangement visible from point p
        * 
        *  This method computes the visibility region arrangement given a point p in the arrangement 
        */
        Arrangement_2 compute_visibility() {
            std::vector<Point_2> visible_points;
            Arrangement_2 visibility_region;

            for (auto guard : this->guards) {
                // find the face of the guard
                Arrangement_2::Face_const_handle *face;
                CGAL::Arr_naive_point_location<Arrangement_2> pl(this->input_arrangement);
                CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(guard);

                // The query point is located in the interior of a face
                face = boost::get<Arrangement_2::Face_const_handle> (&obj);

                // define type of visibility algorithm used
                // TODO: should make it changeable at invocation time
                TEV visibility(this->input_arrangement);
                visibility.compute_visibility(guard, *face, visibility_region);

                // add the visible points only once
                // TODO: probably optimise this
                // for (auto eit = visibility_region.edges_begin(); eit != visibility_region.edges_end(); ++ eit) {
                //     push_back_unique(visible_points, eit->source()->point());
                //     push_back_unique(visible_points, eit->target()->point());
                // }
            }

            return visibility_region;
        }
    
    private:
        Arrangement_2 input_arrangement;
        std::vector<Point_2> guards;
        std::vector<Point_2> boundary_vertices;
};