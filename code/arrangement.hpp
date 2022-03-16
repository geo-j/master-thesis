#include <fstream>
#include <algorithm>

#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/draw_polygon_2.h>

#include "utils.hpp"


typedef Arrangement_2::Face_handle                                          Face_handle;
typedef Arrangement_2::Edge_const_iterator                                  Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                              Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>    NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_true>  RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>              TEV;



class Arrangement {
    public:
       /* overloaded input operator
        * The format of the input file is:
        * E                     * number of edges
        * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        * p3.x p3.y p4.x p4.y
        * ...
        */
        friend std::istream &operator>>(std::istream &f, Arrangement &a) {
            std::size_t E;
            double x1, y1, x2, y2;
            std::vector<Segment_2> segments;

            f >> E;

            for (auto i = 0; i < E; i ++) {
                f >> x1 >> y1 >> x2 >> y2;
                Point_2 p1(x1, y1), p2(x2, y2);

                // create segments for arrangement
                segments.push_back(Segment_2(p1, p2));

                // create points for polygon
                a.input_polygon.push_back(p1);
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

        /* read_guards method
        * :param stream f: data stream from where the guards are input
        * The format of the input guards file is:
        * G             * number of guards
        * q1.x q1.y     * coordinates of the guard in the format of q(x, y) separated by spaces
        * q2.x q2.y     
        */
        template<typename stream>
        void read_guards(stream &f) {
            std::size_t n_guards;
            f >> n_guards;

            for (auto i = 0; i < n_guards; i ++) {
                double x, y;
                f >> x >> y;
                Point_2 q(x, y);
                this->add_guard(q);
            }
        }

        /* print_polygon method
        * :param stream f: data stream where the arrangement should be output
        *
        * The format of the output file is:
        * V                     * number of vertices of the polygon
        * p1.x p1.y             * vertex with coordinates p1(x, y), separated by spaces, in clockwise order
        * p2.x p2.y
        * ...
        * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
        * p3.x p3.y p4.x p4.y
        * ...
        */
        template<typename stream>
        void print_polygon(stream &f) {
            f << this->input_polygon.size() << std::endl;

            for (auto vit = this->input_polygon.vertices_begin(); vit != this->input_polygon.vertices_end(); ++ vit) {
                f << *vit << std::endl;
            }
            f << std::cout;
            for (auto eit = this->input_polygon.edges_begin(); eit != this->input_polygon.edges_end(); ++ eit) {
                f << *eit << std::endl;
            }
        }

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
        void add_guard(const Point_2 q) {
            this->guards.push_back(q);
        }


        // TODO: probably some edge cases are missing based on the arrangement to polygon conversion
        /* is_completely_visible method
        *  :param Arrangement_2 visibility_arrangement: the visibility arrangement resulted from computing the visibility of one or multiple guards
        *  :return bool                               : whether the whole input arrangement is seen by the visibility arrangement
        *       1                                     : whole arrangement is seen
        *       0                                     : otherwise
        * 
        *  This method computes whether the whole input arrangement is completely seen by the visibility arrangement of one or multiple guards
        *  This is achieved by converting the visibility arrangement to a polygon and using the built-in == operator.
        */
        bool is_completely_visible(Arrangement_2 visibility_arrangement) {
            Polygon_2 visibility_polygon = arrangement_to_polygon(visibility_arrangement);
            // CGAL::draw(visibility_polygon);
            // CGAL::draw(this->input_polygon);
            // for (auto vit = visibility_polygon.vertices_begin(); vit != visibility_polygon.vertices_end(); ++ vit) {
            //     std::cout << *vit << std::endl;
            // }
            // std::cout << std::endl;

            return visibility_polygon == this->input_polygon;
        }

        /* compute_full_visibility method
        *  :return Arrangement_2           :arrangement visible from all guards
        * 
        *  This method computes the visibility region arrangement of all the guards
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 compute_full_visibility() {
            std::vector<Point_2> visible_points;
            Arrangement_2 prev_visibility_polygon, cur_visibility_polygon, joined_visibility_polygon;

            for (auto i = 0; i < this->guards.size(); i ++) {
                auto guard = this->guards.at(i);

                // find the face of the guard
                // Arrangement_2::Face_const_handle *face;
                CGAL::Arr_naive_point_location<Arrangement_2> pl(this->input_arrangement);
                auto obj = pl.locate(guard);

                // The query point is located in the interior of a face
                auto *face = boost::get<Arrangement_2::Face_const_handle> (&obj);

                // define type of visibility algorithm used
                // TODO: should make it changeable at invocation time
                TEV visibility(this->input_arrangement);
                Arrangement_2 visibility_arrangement;
                visibility.compute_visibility(guard, *face, visibility_arrangement);

                if (i == 0)
                    // initialise first visibility polygon
                    prev_visibility_polygon = visibility_arrangement;
                else {
                    // if there are multiple guards, join the previously computed visibility polygon with the current one
                    cur_visibility_polygon = visibility_arrangement;
                    CGAL::overlay(prev_visibility_polygon, cur_visibility_polygon, joined_visibility_polygon);
                    prev_visibility_polygon = joined_visibility_polygon;
                }
                // add the visible points only once
                // for (auto eit = visibility_region.edges_begin(); eit != visibility_region.edges_end(); ++ eit) {
                //     push_back_unique(visible_points, eit->source()->point());
                //     push_back_unique(visible_points, eit->target()->point());
                // }
            }

            return prev_visibility_polygon;
        }
    
    private:
        Arrangement_2 input_arrangement;
        // TODO: should probably adapt it to polygons with holes
        Polygon_2 input_polygon;
        // TODO: maybe guards class in the future?
        std::vector<Point_2> guards;
};