#include <fstream>
#include <algorithm>
#include <tuple>

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
        *  This method computes the visibility region arrangement of a guard
        */
        Arrangement_2 compute_guard_visibility(Point_2 guard) {
            // define type of visibility algorithm used
            // TODO: should make it changeable at invocation time
            NSPV visibility(this->input_arrangement);
            Arrangement_2 visibility_arrangement;

            // find the face of the guard
            CGAL::Arr_naive_point_location<Arrangement_2> pl(this->input_arrangement);
            auto obj = pl.locate(guard);

            // check if the query point is located in the interior of a face
            auto *face = boost::get<Arrangement_2::Face_const_handle>(&obj);

            if (face)
                visibility.compute_visibility(guard, *face, visibility_arrangement);
            // if the guard is on the arrangement boundary, then get the inner face of that edge to compute its visibility
            else {
                auto *edge = boost::get<Arrangement_2::Halfedge_const_handle>(&obj);

                if (edge) {
                    try {
                        visibility.compute_visibility(guard, (*edge)->twin()->ccb(), visibility_arrangement);
                    } catch(CGAL::Assertion_exception) {
                        visibility.compute_visibility(guard, (*edge)->ccb(), visibility_arrangement);
                        std::cout << "here\n";

                    }
                }
            }

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
            if (vertex)
                return true;

            return false;
        }


        /* *****************
           *   Gradient    *
           *****************
        */

        /* get_reflex_intersection_pairs method
        * :in param Arrangement_2 visibility_arrangement:                           the visibility region of the guard
        * :in param Point_2 guard:                                                  the guard whose reflex intersection point we want to compute
        * :return std::vector<std::tuple<Point_2, Point_2, CGAL:Oriented_side>>:    vector of tuples containing 
        *                                                                               (all reflex vertices,
        *                                                                                their arrangement boundary intersection points,
        *                                                                                the orientation of the guard w.r.t. the boundary line of the reflex vertex)
        * 
        *  This method computes the tuples between all the reflex vertices a guard sees, their intersection points with the input arrangement boundaries and the orientation of the guard in relation to the boundary of the polygon and a specific reflex vertex
        */
        std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> get_reflex_intersection_pairs(Arrangement_2 visibility_arrangement, Point_2 guard) {
            std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> boundary_intersections;

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
                    if (CGAL::collinear(guard, reflex_vertex, eit->prev()->source()->point())) {
                        Line_2 boundary(reflex_vertex, eit->target()->point());

                        CGAL::Oriented_side orientation = boundary.oriented_side(guard);
                        boundary_intersections.push_back(std::make_tuple(reflex_vertex, eit->prev()->source()->point(), orientation));
                    }

                    // or whether we're in the order guard - reflex vertex - intersection point
                    else if (CGAL::collinear(guard, reflex_vertex, eit->target()->point())) {
                        Line_2 boundary(reflex_vertex, eit->prev()->source()->point());

                        CGAL::Oriented_side orientation = boundary.oriented_side(guard);
                        boundary_intersections.push_back(std::make_tuple(reflex_vertex, eit->target()->point(), orientation));
                    }
                }

            } while (++ eit != *visibility_arrangement.unbounded_face()->inner_ccbs_begin());

            return boundary_intersections;
        }

        /* gradient method
        * :in param Arrangement_2 visibility_arrangement:   the visibility region of the guard
        * :in param Point_2 guard:                          guard point whose gradient needs to be computed
        * :return Vector_2:                                 gradient of the guard as a vector
        * 
        * This method computes the gradient of a guard around all the reflex vertices it sees
        */
        Vector_2 gradient(Arrangement_2 visibility_arrangement, Point_2 guard) {
            Vector_2 Df;
            // get all (reflex vertex, boundary intersection point, orientation) tuples for the guard
            auto reflex_intersections = this->get_reflex_intersection_pairs(visibility_arrangement, guard);

            for (auto i = 0; i < reflex_intersections.size(); i ++) {
                // unpack the tuple
                Point_2 reflex_vertex = std::get<0>(reflex_intersections.at(i));
                Point_2 intersection = std::get<1>(reflex_intersections.at(i));
                CGAL::Oriented_side orientation = std::get<2>(reflex_intersections.at(i));

                // compute distances between guard - reflex vertex - intersection point
                auto alpha = CGAL::squared_distance(guard, reflex_vertex);
                auto beta = CGAL::squared_distance(reflex_vertex, intersection);

                // compute guard-reflex vertex vector
                Vector_2 v = Vector_2(guard, reflex_vertex);

                // compute orthogonal vector to the guard-reflex vector
                // if the guard is on the positive side of the one of the edges of the arrangement the reflex vertex is on, the vector needs to be clockwise perpendicular,
                //      otherwise, counterclockwise
                Vector_2 vp;
                if (orientation == CGAL::ON_POSITIVE_SIDE)
                    vp = v.perpendicular(CGAL::CLOCKWISE);
                else
                    vp = v.perpendicular(CGAL::COUNTERCLOCKWISE);

                // compute Df for reflex vertex r
                Vector_2 Dfr = vp * ((beta * beta) / (2 * alpha)) * (1 / alpha);

                // initialise total gradient Df if first reflex vertex,
                //      otherwise just add all gradient vectors
                if (i == 0) 
                    Df = Dfr;
                else
                    Df += Dfr;
            }
        
            return Df;
        }

        /* optimise method
        *
        * This method optimises the position of all guards using gradient descent.
        * The optimisation process stops when the guard position cannot be changed, or the guard is moved outside of the polygon
        */
       // TODO: see where to move the learning rate; probably in main?
        void optimise() {
            float learning_rate = 0.2;

            for (auto i = 0; i < this->guards.size(); i ++) {
                Vector_2 gradient;
                Point_2 cur_guard_position = this->guards.at(i), prev_guard_position;
                Arrangement_2 visibility_arrangement;
                std::cout << cur_guard_position << std::endl;

                int j = 10;
                // try to update the guard position until there are no more changes, or it goes outside the arrangement
                do {
                    // compute visibility arrangement of each guard position
                    visibility_arrangement = this->compute_guard_visibility(cur_guard_position);
                    prev_guard_position = cur_guard_position;

                    // compute gradient of current guard position
                    gradient = this->gradient(visibility_arrangement, prev_guard_position);

                    // update current guard position
                    cur_guard_position = Point_2(prev_guard_position.x() + learning_rate * gradient.x(), prev_guard_position.y() + learning_rate * gradient.y());

                    // std::cout << prev_guard_position << ';' << cur_guard_position << std::endl;
                    // if the current guard position is not inside the arrangement, then it means the gradient requires it to be outside; so place it on the boundary
                    if (!is_inside_arrangement(this->input_arrangement, cur_guard_position)) {
                        // std::cout << "not inside\n";
                        cur_guard_position = this->place_guard_on_boundary(prev_guard_position, cur_guard_position);
                    }
                    std::cout << cur_guard_position << std::endl;

                    j ++;
                    // std::cout << is_inside_arrangement(this->input_arrangement, cur_guard_position) << std::endl;
                    
                } while (prev_guard_position != cur_guard_position && is_inside_arrangement(this->input_arrangement, cur_guard_position) && j < 17);

            }
        }

    /* place_guard_on_boundary method
    * :in param Point_2 prev_guard:     the previous position of the guard
    * :in param Point_2 guard:          the guard position it should have according to the gradient
    * :return Point_2:                  the boundary position of the guard
    * 
    * This method computes the guard's position on the arrangement's boundary in the case when the gradient requires it to be outside of the polygon
    */
    Point_2 place_guard_on_boundary(Point_2 prev_guard, Point_2 guard) {
        // std::cout << "prev guard " << prev_guard << " wished guard " << guard << std::endl;
        auto guard_movement = Segment_2(prev_guard, guard);
        auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
        Point_2 new_guard;
        Segment_2 edge;

        do {
            // std::cout<<"here";
            edge = Segment_2(eit->source()->point(), eit->target()->point());
            // compute the intersection between the guard's gradient direction and the arrangement boundary
            auto intersection = CGAL::intersection(edge, guard_movement);

            if (intersection) {
                new_guard = *boost::get<Point_2>(&*intersection);
                break;
            }

        } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());

        if (prev_guard == new_guard) {
            auto edge_line = edge.supporting_line();
            new_guard = edge_line.projection(guard);
            // std::cout << "boundary guard " << new_guard << std::endl;
        }

        return new_guard;
    }


    private:
        Arrangement_2 input_arrangement;
        // TODO: should probably adapt it to polygons with holes
        Polygon_2 input_polygon;
        // TODO: maybe guards class in the future?
        std::vector<Point_2> guards;

        std::vector<Point_2> reflex_vertices;
};