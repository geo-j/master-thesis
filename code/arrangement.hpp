#include <fstream>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <set>
#include <iterator>

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
#include <CGAL/Polygon_2_algorithms.h>

#include "guard.hpp"


typedef Arrangement_2::Face_handle                                                  Face_handle;
typedef Arrangement_2::Edge_const_iterator                                          Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                                      Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_false>           NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_false>         RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_false>     TEV;

typedef Kernel::Ray_2                                                               Ray_2;
typedef Kernel::Intersect_2                                                         Intersect_2;
typedef Kernel::Object_2                                                            Object_2;
typedef Kernel::Line_2                                                              Line_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                                          Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                                             Pwh_list_2;



class Arrangement {
    public:
        // TODO: move to .cpp file
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
            a.visibility_algo.attach(a.input_arrangement);
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
        void read_guards(stream &f, double learning_rate, double pull_attraction) {
            std::size_t n_guards;
            f >> n_guards;

            for (auto i = 0; i < n_guards; i ++) {
                double x, y;
                f >> x >> y;
                Point_2 q(x, y);

                this->add_guard(q, learning_rate, pull_attraction);

            }
        }

        template<typename stream>
        void init_guards(stream &f, double learning_rate, double pull_attraction) {
            std::size_t n_guards;
            f >> n_guards;

            for (auto i = 0; i < n_guards; i ++) {
                if (i == 0) {
                    auto eit = *(this->input_arrangement).unbounded_face()->inner_ccbs_begin();
                    auto edge = Segment_2(eit->source()->point(), eit->target()->point());

                    Point_2 q((eit->source()->point().x() + eit->target()->point().x()) / 2, (eit->source()->point().y() + eit->target()->point().y()) / 2);
                    this->add_guard(q, learning_rate, pull_attraction);
                } else {
                    auto full_visibility_arrangement = this->full_visibility();
                    auto full_visibility_polygon = arrangement_to_polygon(full_visibility_arrangement);

                    if (this->input_polygon.is_clockwise_oriented())
                        this->input_polygon.reverse_orientation();
                    if (full_visibility_polygon.is_clockwise_oriented())
                        full_visibility_polygon.reverse_orientation();

                    Pwh_list_2 symmR;
                    CGAL::difference(this->input_polygon, full_visibility_polygon, std::back_inserter(symmR));

                    auto v = symmR.begin()->outer_boundary().edge(0);

                    Point_2 q((v.source().x() + v.target().x()) / 2, (v.source().y() + v.target().y()) / 2);
                    this->add_guard(q, learning_rate, pull_attraction);

                }
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
                f << guard << std::endl;
            }
        }

        /* add_guard method
        *  :in param Point_2 q:                         point that would guard the arrangement
        *  :in param Arrangement_2 visibility_region:   visibility region of the guard
        * 
        * This method adds the guard with its corresponding visibility region to the guard vector
        */
        void add_guard(const Point_2 q, const double learning_rate, double pull_attraction) {
            Arrangement_2 visibility_region = this->visibility(q);
            this->guards.push_back(Guard(q, visibility_region, learning_rate, pull_attraction));
        }

        void add_reflex_vertices() {
            // std::cout << "here\n";
            auto eit = *(this->input_arrangement).unbounded_face()->inner_ccbs_begin();

            // Identify reflex vertices
            // use do/while for circular loop
            // The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
            do {
                //Left turn, because the boundary is clockwise...
                // std::cout << "++++++ checking left turn for " << eit->prev()->source()->point() << ',' << eit->prev()->target()->point() << ',' << eit->target()->point() << std::endl;
                if (CGAL::orientation(eit->prev()->source()->point(), eit->prev()->target()->point(), eit->target()->point()) == CGAL::LEFT_TURN) {
                    // identify reflex vertex
                    Point_2 reflex_vertex = eit->source()->point();
                    // std::cout << "found reflex vertex " << reflex_vertex << std::endl;

                    this->reflex_vertices.push_back(reflex_vertex);
                }
            } while (++ eit != *(this->input_arrangement).unbounded_face()->inner_ccbs_begin());
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

            // TODO: deal with inner holes of whole visibility arrangement
            // visibility_arrangement.number_of_faces() == 2 && 
            return visibility_polygon == this->input_polygon;
        }

        /* visibility method
        *  :in param Point_2 guard         :guard whose visibility region needs to be computed
        *  :return Arrangement_2           :arrangement visible from the guard
        * 
        *  This method computes the visibility region arrangement of a guard
        */
        Arrangement_2 visibility(const Point_2 guard) {
            Arrangement_2 visibility_arrangement;

            auto obj = this->pl.locate(guard);

            // check if the query point is located in the interior of a face
            auto *face = boost::get<Arrangement_2::Face_const_handle>(&obj);

            if (face) {
                // std::cout << "before face\n";
                this->visibility_algo.compute_visibility(guard, *face, visibility_arrangement);
                // std::cout << (*face)->number_of_isolated_vertices() << std::endl;

                // std::cout << "after face\n";

            }
            // if the guard is on the arrangement boundary, then get the inner face of that edge to compute its visibility
            else {                
                auto *edge = boost::get<Arrangement_2::Halfedge_const_handle>(&obj);

                if (edge) {
                    if ((*edge)->is_on_inner_ccb())
                        this->visibility_algo.compute_visibility(guard, (*edge)->twin()->ccb(), visibility_arrangement);
                    else
                        this->visibility_algo.compute_visibility(guard, (*edge)->ccb(), visibility_arrangement);
                } else {
                    auto *vertex = boost::get<Arrangement_2::Vertex_const_handle>(&obj);
                    // std::cout << "vertex\n";

                    if (vertex) {
                        auto he = this->input_arrangement.halfedges_begin();
                        while (he->target()->point() != (*vertex)->point()
                               || he->face()->is_unbounded()) {
                            he ++;
                            // std::cout << Segment_2(he->source()->point(), he->target()->point()) << std::endl;
                        }

                        // std::cout << (he->target()->point() == (*vertex)->point()) << std::endl;
                        
                        this->visibility_algo.compute_visibility(guard, he, visibility_arrangement);
                    }
                }

            }

            return visibility_arrangement;
        }

        /* full_visibility method
        *  :return Arrangement_2           :arrangement visible from all guards
        * 
        *  This method computes the visibility region arrangement of all the guards
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 full_visibility() {
            std::vector<Point_2> visible_points;
            Arrangement_2 prev_visibility_arrangement, cur_visibility_arrangement, joined_visibility_arrangement;

            for (auto i = 0; i < this->guards.size(); i ++) {
                auto visibility_arrangement = this->guards.at(i).get_visibility_region();

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

        /* exclusive_visibility method
        *  :in param Guard guard         : guard whose visibility region should be excluded from the computation
        *  :return Arrangement_2           :arrangement visible from all guards except current guard
        * 
        *  This method computes the visibility region arrangement of all the guards except the current guard
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 exclusive_visibility(Guard guard, std::vector<Guard> zero_df_guards) {
            std::vector<Point_2> visible_points;
            Arrangement_2 prev_visibility_arrangement, cur_visibility_arrangement, joined_visibility_arrangement;
            bool computed_visibility = false;

            for (auto i = 0; i < zero_df_guards.size(); i ++) {
                if (zero_df_guards.at(i) != guard) {
                    // std::cout << "different guards " << this->guards.at(i) << " and " << guard << std::endl;
                    auto visibility_arrangement = zero_df_guards.at(i).get_visibility_region();
                    computed_visibility = true;

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
            }

            if (!computed_visibility)
                return guard.get_visibility_region();
            
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
            auto visibility_region = this->visibility(p);

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

        /* reflex_vertex_pairs method
        * :in param Arrangement_2 visibility_arrangement:                           the visibility region of the guard
        * :in param Guard guard:                                                  the guard whose reflex intersection point we want to compute
        * :return std::vector<std::tuple<Point_2, Point_2, CGAL:Oriented_side>>:    vector of tuples containing 
        *                                                                               (all reflex vertices,
        *                                                                                their arrangement boundary intersection points,
        *                                                                                the orientation of the guard w.r.t. the boundary line of the reflex vertex)
        * 
        *  This method computes the tuples between all the reflex vertices a guard sees, their intersection points with the input arrangement boundaries and the orientation of the guard in relation to the boundary of the polygon and a specific reflex vertex
        */
        std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> reflex_vertex_pairs(const Guard g) {
            std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> boundary_intersections;
            Point_2 guard = g.get_coords();
            auto visibility_arrangement = g.get_visibility_region();

            auto eit = *visibility_arrangement.unbounded_face()->inner_ccbs_begin();

            // Identify reflex vertices
            // use do/while for circular loop
            // The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
            do {

                for (auto reflex_vertex : this->reflex_vertices) {
                    // std::cout << "collinear " << eit->prev()->source()->point() << '-' << reflex_vertex << '-' << eit->prev()->target()->point() << '?' << CGAL::collinear(eit->prev()->source()->point(), reflex_vertex, eit->prev()->target()->point()) << std::endl;
                }

                //Left turn, because the boundary is clockwise...
                Point_2 reflex_vertex = eit->source()->point();
                // std::cout << "++++++ checking left turn in " << eit->prev()->source()->point() << ',' << eit->prev()->target()->point() << ',' << eit->target()->point() << " for reflex vertex " << reflex_vertex << std::endl;
                auto it = std::find(this->reflex_vertices.begin(), this->reflex_vertices.end(), reflex_vertex);

                if (CGAL::orientation(eit->prev()->source()->point(), eit->prev()->target()->point(), eit->target()->point()) == CGAL::LEFT_TURN
                    ||
                    it != this->reflex_vertices.end()
                ) {
                    // identify reflex vertex
                    // std::cout << "found reflex vertex " << reflex_vertex << std::endl;

                    // check the clockwise orientation of the reflex vertex and the guard
                    // whether we're in the order intersection point - reflex vertex - guard
                    if (CGAL::collinear_are_ordered_along_line(eit->prev()->source()->point(), reflex_vertex, guard)) {
                        // std::cout << eit->prev()->source()->point() << ' ' << reflex_vertex << ' ' << guard << " collinear\n";
                        Line_2 boundary(reflex_vertex, eit->target()->point());

                        CGAL::Oriented_side orientation = boundary.oriented_side(guard);
                        boundary_intersections.push_back(std::make_tuple(reflex_vertex, eit->prev()->source()->point(), orientation));
                    }
                    // or whether we're in the order guard - reflex vertex - intersection point
                    else if (CGAL::collinear_are_ordered_along_line(guard, reflex_vertex, eit->target()->point())) {
                        // std::cout << eit->target()->point() << ' ' << reflex_vertex << ' ' << guard << " collinear\n";

                        Line_2 boundary(reflex_vertex, eit->prev()->source()->point());

                        CGAL::Oriented_side orientation = boundary.oriented_side(guard);
                        boundary_intersections.push_back(std::make_tuple(reflex_vertex, eit->target()->point(), orientation));
                    }
                }
            } while (++ eit != *visibility_arrangement.unbounded_face()->inner_ccbs_begin());

            return boundary_intersections;
        }

        /* exclusive_beta method
        * :in param Point_2 guard:                      the guard whose overlapping beta needs to be computed
        * :in param Point_2 intersection:               the intersection point on the polygon's boundary that the guard sees behind the reflex vertex
        * :in param Point_2 reflex_vertex:              the reflex vertex around which the guard needs to move
        * :in param double beta:                        the initial beta without any overlapping seen regions
        * :return double:                               the value of beta given that the guard's visibility region is also seen by other guards
        * 
        * This method computes the value of beta for the area seen exclusively by the guard.
        */
        double exclusive_beta(Guard guard, Point_2 intersection, Point_2 reflex_vertex, std::vector<Guard> zero_df_guards) {
            auto reflex_boundary_segment = Segment_2(reflex_vertex, intersection);
            std::vector<Point_2> shared_visibility_intersection_points;
            double beta = 0;


            // if (zero_df_guards.size() > 1) {
            auto guard_exclusive_visibility_region = this->exclusive_visibility(guard, zero_df_guards);

            auto eit = *guard_exclusive_visibility_region.unbounded_face()->inner_ccbs_begin();

            // compute the intersection points of all the guards' visibility regions (excluding the current one) with the reflex vertex - polygon boundary segment
            // aka compute the intersection points of the reflex vertex - polygon boundary segment that are seen by other guards
            do {
                auto edge = Segment_2(eit->source()->point(), eit->target()->point());
                auto visibility_intersection = CGAL::intersection(edge, reflex_boundary_segment);

                if (visibility_intersection) {
                    auto *intersection_point = boost::get<Point_2>(&*visibility_intersection);

                    if (intersection_point) {
                        push_back_unique(shared_visibility_intersection_points, *intersection_point);
                        // std::cout << "intersection with visibility edge " << *intersection_point << std::endl;
                    }
                    else {
                        auto intersection_segment = *boost::get<Segment_2>(&*visibility_intersection);
                        // std::cout << "intersection with visibility edge " << intersection_segment << std::endl;
                        push_back_unique(shared_visibility_intersection_points, intersection_segment.source());
                        push_back_unique(shared_visibility_intersection_points, intersection_segment.target());
                    }
                }
            } while (++ eit != *guard_exclusive_visibility_region.unbounded_face()->inner_ccbs_begin());

            auto minus = 1;

            // check if the polygon boundary intersection point is seen. If it is seen, then it is removed. Otherwise, it is added to the intersection points vertex. 
            // This ensures that the intersection points vector is of the "shape" s.t. every pair of vertices starting from the rightmost signifies the end and the beginning of a triangle, respectively
            // e.g.: beta2 - beta1 <-> correspond to triangle 1
            auto it = std::find(shared_visibility_intersection_points.begin(), shared_visibility_intersection_points.end(), intersection);
            if (it == shared_visibility_intersection_points.end()) {
                shared_visibility_intersection_points.push_back(intersection);
                // std::cout << "intersection point " << intersection << " not found and pushed\n";
            }
            else {
                shared_visibility_intersection_points.erase(it);
                // std::cout << "intersection point " << intersection << " found and erased\n";
            }

            std::sort(shared_visibility_intersection_points.begin(), shared_visibility_intersection_points.end());

            if (shared_visibility_intersection_points.size() >= 1) {

                // if guard on the left-hand side of the visible segment
                if (guard.get_coords() < shared_visibility_intersection_points.at(0)) {
                    // if reflex vertex seen, switch the signs
                    // std::cout << "on left-side of visible segment " << shared_visibility_intersection_points.at(0) << ' ' << shared_visibility_intersection_points.at(shared_visibility_intersection_points.size() - 1)<< std::endl;

                    // starting from the rightmost point, sequentially add and subtract every beta s.t. the beta pairs resembling the triangles are created
                    for (int i = shared_visibility_intersection_points.size() - 1; i >= 0; i --) {
                        // std::cout << i << std::endl;
                        beta += distance(shared_visibility_intersection_points.at(i), reflex_vertex) * minus;
                        // std::cout << "adding distance up to " << shared_visibility_intersection_points.at(i) << " = " << CGAL::squared_distance(shared_visibility_intersection_points.at(i), reflex_vertex) * minus << std::endl;
                        minus = -minus;
                    }
                } else {
                    // std::cout << "on right-side of visible segment " << shared_visibility_intersection_points.at(0) << ' ' << shared_visibility_intersection_points.at(shared_visibility_intersection_points.size() - 1) << std::endl;


                    for (auto i = 0; i <= shared_visibility_intersection_points.size() - 1; i ++) {
                        // std::cout << i << std::endl;
                        // std::cout << "adding distance up to " << shared_visibility_intersection_points.at(i) << " = " << CGAL::squared_distance(shared_visibility_intersection_points.at(i), reflex_vertex) * minus << std::endl;
                        beta += distance(shared_visibility_intersection_points.at(i), reflex_vertex) * minus;
                        minus = -minus;
                    }
                }
            }
            // }

            // std::cout << "beta " << beta << std::endl;
            return beta;
        }

        /* point_behind_reflex_vertex method
        :in param Guard guard:             the guard who sees the reflex vertex
        :in param Point_2 reflex_vertex:   the reflex vertex behind which the angle is of interest
        :return double:                    the end of the segment behind the reflex vertex (unseen by the guard)

        This method computes the end of the segment (unseen by the guard) on the other side of the reflex vertex seen by a guard. 
        */
        Point_2 point_behind_reflex_vertex(const Guard guard, Point_2 reflex_vertex) {
            auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();

            // Identify reflex vertices
            // use do/while for circular loop
            // The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
            do {
                // std::cout << "------ checking segment " << eit->source()->point() << ',' << eit->target()->point() << std::endl;                
                
                // if we're on the segment with the reflex vertex
                if (eit->source()->point() == reflex_vertex) {
                    // get boundary segment whose end is the reflex vertex
                    // Line_2 boundary(eit->target()->point(), reflex_vertex);
                    // auto orientation = boundary.oriented_side(guard.get_coords());

                    // if the guard sees the other end of the segment, we need to take the beginning of the previous segment
                    if (this->is_visible_from(guard.get_coords(), eit->target()->point())) {
                        // std::cout << "guard sees end " << eit->target()->point() << " so point behind vertex is " << eit->prev()->source()->point() << std::endl;
                        return eit->prev()->source()->point();
                    }
                    // otherwise take current segment end
                    else {
                        // std::cout << "guard doesn't see end " << eit->target()->point() << " so point behind vertex is " << eit->target()->point() << std::endl;
                        return eit->target()->point();

                    }
                }
            } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());
        }

        /* gradient method
        * :in param Arrangement_2 visibility_arrangement:   the visibility region of the guard
        * :in param Guard guard:                            guard point whose gradient needs to be computed
        * :return Vector_2:                                 gradient of the guard as a vector
        * 
        * This method computes the gradient of a guard around all the reflex vertices it sees
        */
        std::tuple<std::vector<Vector_2>, std::vector<Vector_2>, std::vector<Point_2>, Line_2> gradient(const Guard g, std::vector<Guard> zero_df_guards) {
            std::vector<Vector_2> Dfs, hs;
            std::vector<Point_2> reflex_vertices;
            Line_2 segment_behind_reflex_vertex;
            Vector_2 Df(0, 0), h(0, 0);
            auto guard = g.get_coords();

            // get all (reflex vertex, boundary intersection point, orientation) tuples for the guard
            auto reflex_vertex_pairs = this->reflex_vertex_pairs(g);

            for (auto j = 0; j < reflex_vertex_pairs.size(); j ++) {
                // unpack the tuple
                auto reflex_intersection_orientation_pair = reflex_vertex_pairs.at(j);
                Point_2 reflex_vertex = std::get<0>(reflex_intersection_orientation_pair);
                Point_2 intersection = std::get<1>(reflex_intersection_orientation_pair);
                CGAL::Oriented_side orientation = std::get<2>(reflex_intersection_orientation_pair);
                auto point_behind_reflex_vertex = this->point_behind_reflex_vertex(g, reflex_vertex);

                if (guard == reflex_vertex)
                    segment_behind_reflex_vertex = Line_2(reflex_vertex, point_behind_reflex_vertex);

                if (guard != reflex_vertex) {
                    reflex_vertices.push_back(reflex_vertex);
                    // compute distances between guard - reflex vertex - intersection point
                    auto alpha = distance(guard, reflex_vertex);
                    auto beta = distance(reflex_vertex, intersection);
                    auto new_beta = beta;

                    // std::cout << "\tlooking at reflex vertex " << reflex_vertex << std::endl;

                    // if there's other guards to be considered, compute the exclusive beta
                    if (zero_df_guards.size() > 1) {
                        // std::cout << "multiple guards to be considered\n";
                        new_beta = this->exclusive_beta(g, intersection, reflex_vertex, zero_df_guards);
                    }
                    // std::cout << beta << ' ' << new_beta << std::endl;

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


                    Vector_2 v1(reflex_vertex, intersection), v2(reflex_vertex, point_behind_reflex_vertex);
                    double angle = 0.001 + acos(CGAL::to_double(v1 * v2 / (sqrt(CGAL::to_double(v1.squared_length().exact())) * sqrt(CGAL::to_double(v2.squared_length().exact())))));
                    // std::cout << "angle between " << point_behind_reflex_vertex << ' ' << reflex_vertex << ' ' << intersection << " is " << angle << std::endl;
                    
                    // only take the angle if no overlapping visibility regions
                    // if (new_beta != beta)
                        angle = 1;
                    // compute partial Df and h for reflex vertex r
                    Vector_2 Dfr = angle * vp * (new_beta / (2 * alpha));
                    Vector_2 hr = angle * v * (new_beta / (2 * alpha * sqrt(alpha)));
                    // std::cout << Dfr << std::endl;
                    Dfs.push_back(Dfr);
                    hs.push_back(hr);

                    // compute pull for reflex vertex r
                    Df += Dfr;
                    h += hr;
                }
            }
            // the sum of the partials is on the last index of the array
            Dfs.push_back(Df);
            hs.push_back(h);
        

            return std::make_tuple(Dfs, hs, reflex_vertices, segment_behind_reflex_vertex);
        }

        /* optimise method
        *
        * This method optimises the position of all guards using gradient descent.
        * 
        * The algorithm steps are:
        *   - as long as the polygon is not fully seen:
        *       - for each guard:
        *           - compute its visibility region for its current coordinates
        *           - update its visibility region for its current coordinates
        *           - compute the gradient for its current coordinates
        *           - update its coordinates based on its gradient
        *           - if the new position of the guard is outside of the polygon, place it on the boundary
        *           - if it's overshooting, increase the learning rate of the guard; otherwise, decrease it
        *           - update guard in the guards vector
        */
        void optimise() {
            auto full_arrangement = this->full_visibility();
            auto l = 0;
            this->add_reflex_vertices();
            std::cout << "total area=" << compute_area(this->input_arrangement) << std::endl;

            do {                
                std::vector<Guard> zero_df_guards(this->guards);
                std::vector<Guard> new_guards(this->guards);
                std::cout << "i=" << l << std::endl;
                std::cout << "area=" << compute_area(full_arrangement) << std::endl;

                while (!zero_df_guards.empty()) {
                    std::cout << "+++++" << zero_df_guards.size() << " guards with gradient 0\n";

                    for (auto j = 0; j < zero_df_guards.size(); j ++) {
                        auto itr = std::find(this->guards.begin(), this->guards.end(), zero_df_guards.at(j));
                        auto i = std::distance(this->guards.begin(), itr);

                        // }
                    // }
                    // for (auto i = 0; i < this->guards.size(); i ++) {
                        Vector_2 gradient, pull;
                        Guard cur_guard = Guard(this->guards.at(i)), prev_guard = Guard(cur_guard);
                        std::cout << "\t======== guard " << cur_guard << std::endl;
                        std::cout << "\twas guard reflex vertex? " << cur_guard.is_reflex_vertex() << std::endl;



                        // compute gradient of current guard position
                        auto cur_visibility = cur_guard.get_visibility_region();

                        auto gradient_pair = this->gradient(cur_guard, zero_df_guards);
                        auto gradients = std::get<0>(gradient_pair);
                        auto pulls = std::get<1>(gradient_pair);
                        auto reflex_vertices = std::get<2>(gradient_pair);
                        auto segment_behind_reflex_vertex = std::get<3>(gradient_pair);
                        gradient = gradients.at(gradients.size() - 1);
                        pull = pulls.at(pulls.size() - 1);
                        std::cout << "\tcomputed gradient " << gradient << " and pull " << pull << std::endl;

                        if (gradient != Vector_2(0, 0) || pull != Vector_2(0, 0)) {
                            std::cout << "\t------- removing guard " << cur_guard << std::endl;
                            zero_df_guards.erase(zero_df_guards.begin() + j);
                            j --;
                        // }

                            std::cout << 'g' << i << '=' << cur_guard << std::endl;
                            std::cout << "area" << i << '=' << cur_guard.get_area() << std::endl;
                            std::cout << "alpha" << i << '=' << cur_guard.get_learning_rate() << std::endl;

                            // update current guard position
                            cur_guard.update_coords(gradients, pulls, reflex_vertices);
                            auto ray = Segment_2(Point_2(prev_guard.get_coords()), Point_2(cur_guard.get_coords()));


                            std::cout << "\t>>>> prev guard " << prev_guard << " cur guard " << cur_guard << std::endl;
                            std::cout << "\twas guard reflex vertex? " << prev_guard.is_reflex_vertex() << std::endl;
                            std::cout << "\twas guard in reflex area? " << prev_guard.is_in_reflex_area() << std::endl;
                            std::cout << "\tdoes guard intersect boundary ray " << ray << '?' << this->intersects_boundary(ray) << std::endl;
                            // if the current guard position is not inside the arrangement, then it means the gradient requires it to be outside; so place it on the boundary
                            // exception for when the guard is on a reflex vertex
                            // TODO: if guard was reflex vertex, don't let it move away from the reflex vertex line

                            if (prev_guard.is_in_reflex_area()) {
                                Point_2 new_guard_position;
                                if (this->place_guard_on_reflex_line(prev_guard, cur_guard, new_guard_position) && this->input_polygon.has_on_bounded_side(new_guard_position) > 0) {
                                    std::cout << "\tnew reflex coords who dis\n";
                                    cur_guard.set_coords(new_guard_position);
                                    // cur_guard.set_reflex_vertex();
                                } else {
                                    std::cout << "\told reflex coords\n";
                                    cur_guard.set_coords(prev_guard.get_coords());
                                }
                            }

                            if (
                                // this->input_polygon.has_on_unbounded_side(cur_guard.get_coords())

                                !cur_guard.is_reflex_vertex() &&
                                this->intersects_boundary(ray)) {
                                    Point_2 new_guard_position;
                                    if (this->place_guard_on_boundary(prev_guard.get_coords(), cur_guard.get_coords(), new_guard_position)) {
                                        cur_guard.set_coords(new_guard_position);
                                        // std::cout << cur_guard << std::endl;
                                    } else {
                                        // std::cout << "YO\n";
                                        cur_guard.set_coords(prev_guard.get_coords());
                                    }
                            }

                            // if the guard is now inside the arrangement, update the guard position in the vector
                            if (!this->input_polygon.has_on_unbounded_side(cur_guard.get_coords())) {
                                // std::cout << "i=" << l << std::endl;
                                // std::cout << "area=" << compute_area(full_arrangement) << std::endl;
                                auto visibility_region = this->visibility(cur_guard.get_coords());
                                // auto new_visibility_arrangement = this->full_visibility();

                                cur_guard.update_visibility(visibility_region);
                                new_guards[i] = Guard(cur_guard);
                            }

                            // std::cout << 'g' << i << '=' << cur_guard << std::endl;

                            break;
                        }
                    }

                    full_arrangement = this->full_visibility();
                    std::cout << "is completely visible? " << this->is_completely_visible(full_arrangement) << std::endl;
                    std::cout << "are areas equal? " <<  compute_area(this->input_arrangement) << '=' << compute_area(full_arrangement) << '?' << (compute_area(this->input_arrangement) == compute_area(full_arrangement)) << std::endl;
                    
                    if (zero_df_guards.size() == this->guards.size() || this->is_completely_visible(full_arrangement) || (compute_area(this->input_arrangement) - compute_area(full_arrangement)) <= 0.005)
                        break;
                }
                l ++;

                this->guards = std::vector<Guard>(new_guards);
            } while(!this->is_completely_visible(full_arrangement) && compute_area(this->input_arrangement) - compute_area(full_arrangement) > 0.005);

            // print fully visible details, as they don't get printed if the while loop finishes (when everything is seen)
            std::cout << "i=" << l + 1 << std::endl;
            std::cout << "area=" << compute_area(full_arrangement) << std::endl;
            for (auto i = 0; i < this->guards.size(); i ++) {
                std::cout << 'g' << i << '=' << this->guards[i] << std::endl;
                std::cout << "alpha" << i << '=' << this->guards[i].get_learning_rate() << std::endl;
                std::cout << "area" << i << '=' << this->guards[i].get_area() << std::endl;
            }
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
    bool place_guard_on_boundary(Point_2 prev_guard, Point_2 guard, Point_2 &new_guard) {
        // std::cout << "prev guard " << prev_guard << std::endl;
        auto guard_movement = Segment_2(prev_guard, guard);
        auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
        Segment_2 min_edge;
        bool placed = false;

        do {
            auto edge = Segment_2(eit->source()->point(), eit->target()->point());
            // compute the intersection between the guard's gradient direction and the arrangement boundary
            auto intersection = CGAL::intersection(edge, guard_movement);

            if (intersection) {
                auto tmp_guard = *boost::get<Point_2>(&*intersection);
                // std::cout<<"here\n";
                // std::cout << "guard " << prev_guard << " wants to move to " << guard << " intersection with edge " << edge << " in " << tmp_guard << std::endl;

                if (!placed || distance(prev_guard, tmp_guard) < distance(prev_guard, new_guard)) {
                    new_guard = Point_2(tmp_guard);
                    min_edge = Segment_2(edge);
                    placed = true;
                }
            }

        } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());


        // if the gradient is perpendicular with the arrangement edge, hence no move, then project the guard on the boundary and start from there
        if (prev_guard == new_guard && placed) {
            auto edge_line = min_edge.supporting_line();
            new_guard = edge_line.projection(guard);

            // if the projection is outside the polygon, place the guard on the closest edge vertex
            if (this->input_polygon.has_on_unbounded_side(new_guard)) {
                if (distance(new_guard, min_edge.source()) <= distance(new_guard, min_edge.target()))
                    new_guard = min_edge.source();
                else
                    new_guard = min_edge.target();
            }
        }

        return placed;
    }

    int place_guard_on_reflex_line(Guard prev_guard, Guard guard, Point_2 &new_guard) {
        std::cout << "reflex vertex guard " << prev_guard << " tries to move to " << guard << std::endl;
        int n_placed = 0;

        auto guard_movement = Segment_2(prev_guard.get_coords(), guard.get_coords());

        auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
        do {
            auto edge = Segment_2(eit->source()->point(), eit->target()->point());
            std::cout << "intersection with reflex line " << edge << " ? " << CGAL::do_intersect(edge, guard_movement) << std::endl;

            if (CGAL::do_intersect(edge, guard_movement)) {
                auto intersection = CGAL::intersection(edge, guard_movement);
                if (intersection) {
                    auto tmp_guard = *boost::get<Point_2>(&*intersection);

                    if (n_placed == 0 || distance(prev_guard.get_coords(), tmp_guard) < distance(prev_guard.get_coords(), new_guard)) {
                        if (tmp_guard != prev_guard.get_coords()) {
                            std::cout << "placed at intersection on reflex line at " << tmp_guard << std::endl; 
                            new_guard = Point_2(tmp_guard);
                        } else {
                            new_guard = Point_2(edge.supporting_line().projection(guard.get_coords()));
                            std::cout << "projected on reflex line at " << new_guard << std::endl;
                        }
                        n_placed ++;
                        // std::cout << "wot\n";
                    }
                }
            }

        } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin() && n_placed < 2);

        return n_placed;
    }

    bool intersects_boundary(Segment_2 ray) {
        auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
        
        do {
            // std::cout << "bhere\n";
            auto edge = Segment_2(eit->source()->point(), eit->target()->point());
            // compute the intersection between the guard's gradient direction and the arrangement boundary
            if (ray.source() != eit->source()->point() && ray.source() != eit->target()->point() && ray.target() != eit->source()->point() && ray.target() != eit->target()->point() && CGAL::do_intersect(edge, ray))
                return true;

        } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());

        return false;
    }

    private:
        Arrangement_2 input_arrangement;
        // TODO: should probably adapt it to polygons with holes
        Polygon_2 input_polygon;
        // define type of visibility algorithm used
        // TODO: should make it changeable at invocation time
        TEV visibility_algo;
        // find the face of the guard
        CGAL::Arr_naive_point_location<Arrangement_2> pl;
        std::vector<Guard> guards;
        std::vector<Point_2> reflex_vertices;
};
