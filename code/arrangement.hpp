#include <fstream>
#include <algorithm>
#include <utility>

#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Vector_2.h>

#include "utils.hpp"


typedef Arrangement_2::Face_handle                                          Face_handle;
typedef Arrangement_2::Edge_const_iterator                                  Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                              Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>    NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_true>  RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>              TEV;

typedef Kernel::Ray_2                                                       Ray_2;
typedef Kernel::Intersect_2                                                 Intersect_2;
typedef Kernel::Object_2                                                    Object_2;
typedef Kernel::Line_2                                                      Line_2;
typedef Kernel::FT                                                          FT;
typedef Kernel::Vector_2                                                    Vector_2;




class Arrangement {
    public:
        /* *****************
           *     I/O       *
           *****************
        */

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
        * :out param stream f: data stream from where the guards are input
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
        * :out param stream f: data stream where the arrangement should be output
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
        * :out param stream f: data stream where the arrangement should be output
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
        *  :in param Point_2 p: point that would guard the arrangement
        */
        void add_guard(const Point_2 q) {
            this->guards.push_back(q);
        }



        /* *****************
           *  Visibility   *
           *****************
        */

        // TODO: probably some edge cases are missing based on the arrangement to polygon conversion
        /* is_completely_visible method
        *  :in param Arrangement_2 visibility_arrangement: the visibility arrangement resulted from computing the visibility of one or multiple guards
        *  :return bool                               : whether the whole input arrangement is seen by the visibility arrangement
        *       1                                     : whole arrangement is seen
        *       0                                     : otherwise
        * 
        *  This method computes whether the whole input arrangement is completely seen by the visibility arrangement of one or multiple guards
        *  This is achieved by converting the visibility arrangement to a polygon and using the built-in == operator.
        */
        bool is_completely_visible(Arrangement_2 visibility_arrangement) {
            Polygon_2 visibility_polygon = arrangement_to_polygon(visibility_arrangement);

            return visibility_polygon == this->input_polygon;
        }

        /* compute_guard_visibility method
        *  :in param Point_2 guard         :guard whose visibility region needs to be computed
        *  :return Arrangement_2           :arrangement visible from the guard
        * 
        *  This method computes the visibility region arrangement of a guards
        */
        Arrangement_2 compute_guard_visibility(Point_2 guard) {
            // find the face of the guard
            CGAL::Arr_naive_point_location<Arrangement_2> pl(this->input_arrangement);
            auto obj = pl.locate(guard);

            // The query point is located in the interior of a face
            auto *face = boost::get<Arrangement_2::Face_const_handle> (&obj);

            // define type of visibility algorithm used
            // TODO: should make it changeable at invocation time
            TEV visibility(this->input_arrangement);
            Arrangement_2 visibility_arrangement;
            visibility.compute_visibility(guard, *face, visibility_arrangement);

            return visibility_arrangement;
        }

        /* compute_full_visibility method
        *  :return Arrangement_2           :arrangement visible from all guards
        * 
        *  This method computes the visibility region arrangement of all the guards
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 compute_full_visibility() {
            std::vector<Point_2> visible_points;
            Arrangement_2 prev_visibility_arrangement, cur_visibility_arrangement, joined_visibility_arrangement;

            for (auto i = 0; i < this->guards.size(); i ++) {
                auto guard = this->guards.at(i);

                Arrangement_2 visibility_arrangement = this->compute_guard_visibility(guard);

                if (i == 0)
                    // initialise first visibility polygon
                    prev_visibility_arrangement = visibility_arrangement;
                else {
                    // if there are multiple guards, join the previously computed visibility polygon with the current one
                    cur_visibility_arrangement = visibility_arrangement;
                    CGAL::overlay(prev_visibility_arrangement, cur_visibility_arrangement, joined_visibility_arrangement);
                    prev_visibility_arrangement = joined_visibility_arrangement;
                }
            }

            return prev_visibility_arrangement;
        }

        /* is_visible_from method
        * :in param Point_2 p:  the viewpoint
        * :in param Point_2 r:  the point that needs to be checked for visibility from viewpoint p
        * :return bool:         true if r is visible from p,
        *                           false otherwise
        * 
        * This method checks whether a point r is visible from p by checking whether r is in the visibility region of p.
        */
        bool is_visible_from(Point_2 p, Point_2 r) {
            // first compute the visibility region of the guard
            auto visibility_region = this->compute_guard_visibility(p);

            // find where in the visibility region the guard is placed
            CGAL::Arr_naive_point_location<Arrangement_2> pl(visibility_region);
            auto obj = pl.locate(r);

            // The query point may be a vertex on the visibility region boundary
            // TODO: handle degenerate collinear case
            auto *vertex = boost::get<Arrangement_2::Vertex_const_handle> (&obj);
            
            
            // if d is in the visibility region of the guard, then it can be counted as an intersection point
            if (vertex) {
                // std::cout << '\t' << p << " can see " << r << std::endl;
                return true;
            }

            return false;
        }

        /* is_inside_arrangement method
        * :in param Point_2 p:  the point that needs to be checked whether it is inside the arrangement
        * :return bool:         true if r is inside the arrangement
        *                           false otherwise
        * 
        * This method checks whether a point p is inside the input arrangement
        */
        bool is_inside_arrangement(Point_2 p) {
            // find where in the visibility region the guard is placed
            CGAL::Arr_naive_point_location<Arrangement_2> pl(this->input_arrangement);
            auto obj = pl.locate(p);

            // The query point may be a face on the visibility region boundary
            auto *face = boost::get<Arrangement_2::Face_const_handle> (&obj);
            
            if (!(*face)->is_unbounded()) {
                // std::cout << '\t' << p << " can see " << r << std::endl;
                return true;
            }

            return false;
        }

        /* *****************
           *   Gradient    *
           *****************
        */

    
        /* get_guard_reflex_arrangement_intersection method
        * 
        * :in param Point_2    guard:             the guard whose reflex-arrangement intersection we want to calculate
        * :in param Point_2    reflex_vertex:     the reflex vertex past which the guard's visibility intersects the arrangement
        * :out param Point_2   d:                 the intersection point between the ray [guard, reflex_vertex] and the arrangement as output
        * :return bool:                           true if an intersection point between the ray [guard, reflex_vertex] and the arrangement was found;
        *                                            false otherwise
        * 
        * This method computes the intersection point between a given guard and a reflex vertex in the input arrangement.
        */
        bool get_guard_reflex_arrangement_intersection(Point_2 guard, Point_2 reflex_vertex, Point_2 &d) {
            // std::cout << "Ray [" << guard << ';' << reflex_vertex << "] " << std::endl;
            Ray_2 pr(guard, reflex_vertex);

            // only look at reflex vertices that are actually visible from the guard
            if (this->is_visible_from(guard, reflex_vertex)) {
                auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();

                do {
                    Segment_2 edge = Segment_2(eit->source()->point(), eit->target()->point());
                    // std::cout << "- tries to intersect boundary segment [" << eit->source()->point() << ", " << eit->target()->point() << "]" << std::endl;

                    Object_2 intersection = CGAL::intersection(pr, edge);
                    const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&intersection);

                    // check if the ray intersects any of the edges of the visibility region
                    if (intersection_point) {
                        d = *intersection_point;
                        // d = *boost::get<Point_2 >(&*intersection);

                        if (d != eit->target()->point() && d != eit->source()->point() // intersection point should be different than the segment endpoints it is on (collinear guard - reflex vertex - intersection point case)
                            && this->is_visible_from(guard, d)) {  // intersection point should be visible from the guard
                                // std::cout << "\tmanages to at " << d << std::endl;
                                return true;
                        }
                    }
                } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());
            }

            return false;
        }

        std::vector<std::pair<Point_2, Point_2>> get_reflex_intersection_pairs(Point_2 guard) {
            Arrangement_2 visibility_arrangement = this->compute_guard_visibility(guard);
            std::vector<std::pair<Point_2, Point_2>> boundary_intersections; //= get_reflex_vertices(visibility_arrangement);

            auto eit = *visibility_arrangement.unbounded_face()->inner_ccbs_begin();
            // Identify reflex vertices
            // use do/while for circular loop
            // The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
            do {
                //Left turn, because the boundary is clockwise...
                if (CGAL::orientation(eit->prev()->source()->point(), eit->prev()->target()->point(), eit->target()->point()) == CGAL::LEFT_TURN) {
                    // identify reflex vertex
                    Point_2 reflex_vertex = eit->source()->point();

                    // check the clockwise orientation of the reflex vertex and the guard
                    // whether we're in the order intersection point - reflex vertex - guard
                    if (CGAL::collinear(guard, reflex_vertex, eit->prev()->source()->point()))
                        boundary_intersections.push_back(std::make_pair(reflex_vertex, eit->prev()->source()->point()));

                    // or whether we're in the order guard - reflex vertex - intersection point
                    else if (CGAL::collinear(guard, reflex_vertex, eit->target()->point()))
                        boundary_intersections.push_back(std::make_pair(reflex_vertex, eit->target()->point()));
                }

            } while (++ eit != *visibility_arrangement.unbounded_face()->inner_ccbs_begin());

            return boundary_intersections;
        }
        /* gradient method
        * :in param Point_2 guard:  guard point whose gradient needs to be computed
        * :return Vector_2:         gradient of the guard as a vector
        * 
        * This method computes the gradient of a guard around all the reflex vertices it sees
        */
        Vector_2 gradient(Point_2 guard) {
            Vector_2 Df;
            
            for (auto i = 0; i < this->reflex_vertices.size(); i ++) {
                Point_2 arrangement_intersection_point;

                if (this->get_guard_reflex_arrangement_intersection(guard, this->reflex_vertices.at(i), arrangement_intersection_point)) {
                    auto alpha = CGAL::squared_distance(guard, this->reflex_vertices.at(i));
                    auto beta = CGAL::squared_distance(arrangement_intersection_point, this->reflex_vertices.at(i));

                    // get the line of the guard p and the arrangement intersection point d
                    Line_2 pd(guard, arrangement_intersection_point);
                    // get the orthogonal line on pd through p
                    Line_2 orthogonal_pd = pd.perpendicular(guard);

                    // calculate the y-coord of the gradient
                    auto y = (beta * beta) / (2 * alpha);
                    // calculate the x-coord of the gradient as rotated on the orthogonal line
                    auto x = orthogonal_pd.x_at_y(y);

                    // create a vector with the partial gradient to be able to add them for each reflex vertex
                    Vector_2 df(guard, Point_2(x, y));

                    // std::cout << "Guard " << guard << " for reflex vertex " << this->reflex_vertices.at(i) << " has gradient " << x << ", " << y << std::endl;

                    // initialise gradient if first reflex vertex,
                    //      otherwise just add all gradient vectors
                    if (i == 0) 
                        Df = df;
                    else
                        Df += df;
                    // std::cout << ">>> Ray [" << guard << "; " << reflex_vertex << "] intersects boundary point " << intersection_point << " and has gradient " << gradient << std::endl;
                }
            }
            return Df;
        }

        void optimise() {
            float learning_rate = 0.01;
            for (auto i = 0; i < this->guards.size(); i ++) {
                Vector_2 gradient;
                Point_2 cur_guard_position = this->guards.at(i), prev_guard_position;
                std::cout << cur_guard_position << std::endl;

                do {
                    prev_guard_position = cur_guard_position;

                    gradient = this->gradient(prev_guard_position);
                    
                    cur_guard_position = Point_2(prev_guard_position.x() + learning_rate * gradient.x(), prev_guard_position.y() + learning_rate * gradient.y());

                    std::cout << prev_guard_position << ';' << cur_guard_position << std::endl;
                    std::cout << this->is_inside_arrangement(cur_guard_position) << std::endl;
                    
                    // std::cout << ">>>> Guard " << this->guards.at(i) << " has gradient (" << this->guards.at(i).x() + learning_rate * gradient.x() << ", " << this->guards.at(i).y() + learning_rate * gradient.y() << ")" << std::endl;
                } while(prev_guard_position != cur_guard_position && this->is_inside_arrangement(cur_guard_position));
            }
        }

        void print_reflex_intersections() {
            for (auto guard : this->guards) {
                auto pairs = this->get_reflex_intersection_pairs(guard);

                for (auto pair : pairs) {
                    std::cout << "Guard " << guard << " reflex vertex " << pair.first << " intersection point " << pair.second << std::endl;
                }

            }
        }


    private:
        Arrangement_2 input_arrangement;
        // TODO: should probably adapt it to polygons with holes
        Polygon_2 input_polygon;
        // TODO: maybe guards class in the future?
        std::vector<Point_2> guards;

        std::vector<Point_2> reflex_vertices;
};