#include <fstream>
#include <algorithm>
#include <tuple>
#include <chrono>

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
#include <CGAL/Polygon_with_holes_2.h>

#include "utils.hpp"


typedef Arrangement_2::Face_handle                                                  Face_handle;
typedef Arrangement_2::Edge_const_iterator                                          Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                                      Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>            NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_true>          RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_true>      TEV;

typedef Kernel::Ray_2                                                               Ray_2;
typedef Kernel::Intersect_2                                                         Intersect_2;
typedef Kernel::Object_2                                                            Object_2;
typedef Kernel::Line_2                                                              Line_2;
typedef Kernel::Vector_2                                                            Vector_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                                          Polygon_with_holes_2;



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
            // skip newline character that isn't read by f, so that newline can properly read the next lines
            f.ignore();

            for (auto i = 0; i < E; i ++) {
                // keep count of which one of the 4 coords we are reading
                int coord = 0;

                // read each line with points
                std::string line, number; 
                std::getline(f, line);
                // stringstream tokenises the string by spaces
                std::stringstream s(line);

                // for each string between spaces, check which one of the 4 coords it is, and get its double value, based on whether it's a fraction or not
                while (s >> number) {
                    switch (coord) {
                        case 0:
                            x1 = get_number(number);
                            break;
                        case 1:
                            y1 = get_number(number);
                            break;
                        case 2: 
                            x2 = get_number(number);
                            break;
                        case 3:
                            y2 = get_number(number);
                            break;
                    }

                    coord ++;
                }

                Point_2 p1(x1, y1), p2(x2, y2);

                // create segments for arrangement
                segments.push_back(Segment_2(p1, p2));

                // create points for polygon
                a.input_polygon.push_back(p1);
            }

            CGAL::insert_non_intersecting_curves(a.input_arrangement, segments.begin(), segments.end());

            a.visibility.attach(a.input_arrangement);
            a.pl = CGAL::Arr_naive_point_location<Arrangement_2>(a.input_arrangement);

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
                Arrangement_2 visibility_region = this->compute_guard_visibility(q);

                this->add_guard(q, visibility_region);

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
        void add_guard(const Point_2 q, const Arrangement_2 visibility_region) {
            this->guards.push_back(q);
            this->visibility_regions.push_back(visibility_region);
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
        bool is_completely_visible(Arrangement_2 &visibility_arrangement) {
            Polygon_2 visibility_polygon = arrangement_to_polygon(visibility_arrangement);

            return visibility_polygon == this->input_polygon;
        }

        /* compute_guard_visibility method
        *  :in param Point_2 guard         :guard whose visibility region needs to be computed
        *  :return Arrangement_2           :arrangement visible from the guard
        * 
        *  This method computes the visibility region arrangement of a guard
        */
        Arrangement_2 compute_guard_visibility(const Point_2 guard) {
            Arrangement_2 visibility_arrangement;

            auto obj = this->pl.locate(guard);

            // check if the query point is located in the interior of a face
            auto *face = boost::get<Arrangement_2::Face_const_handle>(&obj);

            if (face) {
                // std::cout << "before face\n";
                this->visibility.compute_visibility(guard, *face, visibility_arrangement);
                // std::cout << (*face)->number_of_isolated_vertices() << std::endl;

                // std::cout << "after face\n";

            }
            // if the guard is on the arrangement boundary, then get the inner face of that edge to compute its visibility
            else {                
                auto *edge = boost::get<Arrangement_2::Halfedge_const_handle>(&obj);

                if (edge) {
                    // std::cout << "before edge\n";

                    if ((*edge)->is_on_inner_ccb())
                        this->visibility.compute_visibility(guard, (*edge)->twin()->ccb(), visibility_arrangement);
                    else
                        this->visibility.compute_visibility(guard, (*edge)->ccb(), visibility_arrangement);
                    
                    // std::cout << "after edge\n";
                } 
                // TODO: find one of the halfedges where the vertex is located on
                // else {
                //     auto *vertex = boost::get<Arrangement_2::Vertex_const_handle>(&obj);

                //     if (vertex) {
                //         std::cout << "vertex " << (*vertex)->incident_halfedges()->source()->point() << std::endl;
                //         this->visibility.compute_visibility(guard, (*vertex)->incident_halfedges()->twin()->ccb(), visibility_arrangement);
                //     }
                // }

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

            for (auto i = 0; i < this->visibility_regions.size(); i ++) {
                auto visibility_arrangement = this->visibility_regions.at(i);

                // Arrangement_2 visibility_arrangement = this->compute_guard_visibility(guard);

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
        bool is_visible_from(const Point_2 p, const Point_2 r) {
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
        std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> get_reflex_intersection_pairs(Arrangement_2 &visibility_arrangement, const Point_2 guard) {
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
        Vector_2 gradient(Arrangement_2 &visibility_arrangement, const Point_2 guard) {
            // std::cout << "start gradient\n";
            Vector_2 Df;
            // get all (reflex vertex, boundary intersection point, orientation) tuples for the guard
            auto reflex_intersections = this->get_reflex_intersection_pairs(visibility_arrangement, guard);
            // std::cout << reflex_intersections.size() << " intersections\n";

            for (auto i = 0; i < reflex_intersections.size(); i ++) {
                // std::cout << i << ' ';
                // unpack the tuple
                Point_2 reflex_vertex = std::get<0>(reflex_intersections.at(i));
                Point_2 intersection = std::get<1>(reflex_intersections.at(i));
                CGAL::Oriented_side orientation = std::get<2>(reflex_intersections.at(i));

                // compute distances between guard - reflex vertex - intersection point
                // auto alpha = CGAL::squared_distance(guard, reflex_vertex);
                auto alpha = distance(guard, reflex_vertex);
                auto beta = distance(reflex_vertex, intersection);
                // std::cout << "beta = " << beta << std::endl;
                // bool first = true;

                auto visible_segment = Segment_2(reflex_vertex, intersection);

                // compute intersection points with other guards' visibility regions
                for (auto i = 0; i < this->guards.size(); i ++) {
                    // auto left = (guard < (this->guards.at(i))) ? 1 : -1;
                    std::vector<Point_2> intersection_points;

                    // look at all the other guards but the current one
                    if (this->guards.at(i) != guard) {
                        auto eit = *this->visibility_regions.at(i).unbounded_face()->inner_ccbs_begin();
                        do {
                            auto edge = Segment_2(eit->source()->point(), eit->target()->point());
                            const auto visibility_intersection = CGAL::intersection(edge, visible_segment);

                            // check if the reflex vertex - boundary intersection segment of the guard intersects the visibility region of any of the existing guards
                            if (visibility_intersection) {
                                // std::cout << "intersect with " << edge << std::endl;
                                auto *intersection_point = boost::get<Point_2>(&*visibility_intersection);

                                // // if it completely overlaps in a segment, add the 2 segment edges
                                if (!intersection_point) {
                                    auto intersection_segment = *boost::get<Segment_2>(&*visibility_intersection);
                                    // std::cout << "segment\n";
                                }
                                // otherwise just add the point it intersects
                                else {
                                    // std::cout << *intersection_point << std::endl;
                                    // intersection_points.push_back(*intersection_point);
                                    push_back_unique(intersection_points, *intersection_point);
                                }
                            }
                        } while (++ eit != *this->visibility_regions.at(i).unbounded_face()->inner_ccbs_begin());
                    }

                    if (intersection_points.size() == 1) {
                        // beta = 0;
                        // std::cout << "here, point: " << intersection_points.at(0) << ", reflex_vertex = " << reflex_vertex << ", intersection boundary = " << intersection << std::endl;
                    }
                    else if (intersection_points.size() > 1) {
                        std::sort(intersection_points.begin(), intersection_points.end());
                        
                        // TODO: I think I need to also rotate?
                        // case where the reflex vertex is seen
                        if (intersection_points.at(0) == reflex_vertex && intersection_points.at(1) != intersection) {
                            // if the reflex vertex and another part of the reflex vertex - boundary intersection segment is seen, but not the whole segment
                                // std::cout << "\treflex vertex seen\n";
                                // std::cout << "size = " << intersection_points.size() << std::endl;
                                // std::cout << reflex_vertex << ',' << intersection << ',' << intersection_points.at(1) << std::endl;
                                beta -= distance(intersection_points.at(1), reflex_vertex);
                        }
                        // case where the boundary intersection point is seen, but not the reflex vertex
                        else if (intersection_points.at(0) != reflex_vertex && intersection_points.at(1) == intersection) {
                            // whole segment seen case already tackled (do nothing)
                            // std::cout << "\tboundary intersection seen\n";
                            beta -= distance(intersection_points.at(0), intersection);
                        }
                        //rotated case where the boundary is on the left side
                        else if (intersection_points.at(0) != intersection && intersection_points.at(1) == reflex_vertex) {
                                // std::cout << "\trotated seen reflex vertex\n";
                                beta -= distance(intersection_points.at(0), reflex_vertex);
                        }
                        // rotated case where the intersection point is seen first, but not the rest of the segment up to the reflex vertex
                        else if (intersection_points.at(1) != reflex_vertex && intersection_points.at(0) == intersection) {
                            // std::cout << "rotated seen intersection\n";
                            beta -= distance(intersection_points.at(1), intersection);
                        }
                        // if the whole segment is seen (overlapping visibility regions), don't move the guard
                        else if (intersection_points.at(0) == reflex_vertex && intersection_points.at(1) == intersection) {
                            // std::cout << "\there\n";
                            // beta = -beta;
                            beta = 0;
                        }
                        // if only a small part of the segment is seen
                        else {
                            // std::cout << "\tpart of segment seen\n";
                            // beta -= distance(reflex_vertex, intersection_points.at(0)) + distance(reflex_vertex, intersection_points.at(1));
                            beta -= distance(intersection_points.at(0), intersection_points.at(1));
                        }
                        // else
                        //     beta = 0;
                    }

                    // std::cout << "beta after = " << beta << std::endl;
                    // auto minus = -1 * left;
                    // for (auto j = 0; j < intersection_points.size(); j ++) {
                    //     auto b = distance(reflex_vertex, intersection_points.at(j));
                    //     // TODO: check how to find out whether the first region is seen or not
                    //     if (j == 0)
                    //         beta = b;
                    //     else {
                    //         beta += minus * b;
                    //         minus = -minus;
                    //     }
                    // }
                }

                // std::cout << "beta = " << beta << std::endl;
                // std::sort(intersection_points.begin(), intersection_points.end());

                // compute guard-reflex vertex vector
                Vector_2 v = Vector_2(guard, reflex_vertex);
                // Vector_2 w = Vector_2(reflex_vertex, intersection);

                // compute orthogonal vector to the guard-reflex vector
                // if the guard is on the positive side of the one of the edges of the arrangement the reflex vertex is on, the vector needs to be clockwise perpendicular,
                //      otherwise, counterclockwise
                Vector_2 vp;
                if (orientation == CGAL::ON_POSITIVE_SIDE)
                    vp = v.perpendicular(CGAL::CLOCKWISE);
                else
                    vp = v.perpendicular(CGAL::COUNTERCLOCKWISE);

                // compute Df for reflex vertex r
                Vector_2 Dfr = vp * (beta / (2 * alpha));

                // initialise total gradient Df if first reflex vertex,
                //      otherwise just add all gradient vectors
                if (i == 0) 
                    Df = Dfr;
                else
                    Df += Dfr;
            }
        
            // std::cout << "end gradient\n";
            return Df;
        }

        /* optimise method
        *
        * This method optimises the position of all guards using gradient descent.
        * The optimisation process stops when the guard position cannot be changed, or the guard is moved outside of the polygon
        */
        void optimise(double learning_rate) {
            auto full_arrangement = this->compute_full_visibility();
            do {
                for (auto i = 0; i < this->guards.size(); i ++) {
                    Vector_2 gradient, prev_gradient;
                    Point_2 cur_guard_position = this->guards.at(i), prev_guard_position;
                    std::vector<Vector_2> gradients;
                    Arrangement_2 visibility_arrangement;
                    std::cout << 'g' << i << '=' << cur_guard_position << std::endl;

                    // int j = 0;
                    // try to update the guard position until there are no more changes, or it goes outside the arrangement
                    // do {
                        this->visibility_regions[i].clear();
                        // std::cout << "here";

                        // compute visibility arrangement of each guard position
                        this->visibility_regions[i] = this->compute_guard_visibility(cur_guard_position);

                        // prev_guard_position.reset();
                        prev_guard_position = cur_guard_position;
                        // prev_gradient = gradient;

                        // compute gradient of current guard position
                        gradient = this->gradient(this->visibility_regions.at(i), prev_guard_position);
                        // std::cout << "Df = " << gradient << std::endl;

                        // gradient smoothening
                        // if (gradients.size() < 3)
                        // gradient = prev_gradient * 0.3 + 0.7 * gradient;
                        // if (gradients.size() < 2)
                        //     gradients.push_back(gradient);
                        // else {
                        //     for (auto g : gradients)
                        //         gradient += g;
                        //     gradient /= gradients.size();
                        //     gradients.erase(gradients.begin());
                        //     gradients.push_back(gradient);
                        // }

                        // cur_guard_position.reset();
                        // update current guard position
                        cur_guard_position = Point_2(prev_guard_position.x() + learning_rate * gradient.x(), prev_guard_position.y() + learning_rate * gradient.y());

                        // std::cout << prev_guard_position << ';' << cur_guard_position << std::endl;
                        // if the current guard position is not inside the arrangement, then it means the gradient requires it to be outside; so place it on the boundary
                        if (this->input_polygon.bounded_side(cur_guard_position) == CGAL::ON_UNBOUNDED_SIDE) {
                            // std::cout << "not inside\n";
                            Point_2 new_guard_position;
                            if (this->place_guard_on_boundary(prev_guard_position, cur_guard_position, new_guard_position))
                                cur_guard_position = new_guard_position;
                            // std::cout << "now inside\n";
                        }

                        if (this->input_polygon.bounded_side(cur_guard_position) != CGAL::ON_UNBOUNDED_SIDE) {
                            std::cout << 'g' << i << '=' << cur_guard_position << std::endl;
                            this->guards[i] = cur_guard_position;
                            this->visibility_regions[i].clear();
                            this->visibility_regions[i] = this->compute_guard_visibility(cur_guard_position);
                        }
                        

                        // j ++;
                        
                    // } while (!this->is_completely_visible(visibility_arrangement) && (prev_guard_position != cur_guard_position && this->input_polygon.bounded_side(cur_guard_position) != CGAL::ON_UNBOUNDED_SIDE));
                }
                full_arrangement = this->compute_full_visibility();
            } while(!this->is_completely_visible(full_arrangement));
        }

    /* place_guard_on_boundary method
    * :in param Point_2 prev_guard:     the previous position of the guard
    * :in param Point_2 guard:          the guard position it should have according to the gradient
    * :out param Point_2 new_guard:     the new guard boundary position
    * :return bool:                     true if the guard can be placed on the boundary,
    *                                       false otherwise
    * 
    * This method computes the guard's position on the arrangement's boundary in the case when the gradient requires it to be outside of the polygon
    */
   // TODO: no square roots; compute closest segment to point
    bool place_guard_on_boundary(Point_2 prev_guard, Point_2 guard, Point_2 &new_guard) {
        // std::cout << "prev guard " << prev_guard << " wished guard " << guard << std::endl;
        auto guard_movement = Segment_2(prev_guard, guard);
        auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
        Segment_2 edge;
        // std::cout << guard_movement << std::endl;
        bool placed = false;

        do {
            // std::cout<<"here";
            edge = Segment_2(eit->source()->point(), eit->target()->point());
            // compute the intersection between the guard's gradient direction and the arrangement boundary
            auto intersection = CGAL::intersection(edge, guard_movement);

            if (intersection) {
                new_guard = *boost::get<Point_2>(&*intersection);
                placed = true;
                // std::cout << "hello???\n";
                break;
            }

            // edge.reset();
        } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());

        // if the gradient is perpendicular with the arrangement edge, hence no move, then project the guard on the boundary and start from there
        if (prev_guard == new_guard && placed) {
            auto edge_line = edge.supporting_line();
            new_guard = edge_line.projection(guard);
            // std::cout << "boundary guard " << new_guard << std::endl;
        }

        // std::cout << "placed? " << placed << std::endl;
        return placed;
    }


    private:
        Arrangement_2 input_arrangement;
        // TODO: should probably adapt it to polygons with holes
        Polygon_2 input_polygon;
        // define type of visibility algorithm used
        // TODO: should make it changeable at invocation time
        TEV visibility;
        // find the face of the guard
        CGAL::Arr_naive_point_location<Arrangement_2> pl;
        // TODO: maybe guards class in the future?
        std::vector<Point_2> guards;
        std::vector<Arrangement_2> visibility_regions;
        std::vector<Point_2> reflex_vertices;
};