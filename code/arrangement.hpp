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
         friend std::istream &operator>>(std::istream &f, Arrangement &a);

         
         /* overloaded output operator
         * The format of the output file is:
         * E                     * number of edges
         * p1.x p1.y p2.x p2.y   * edge with endpoints coordinates separated by spaces p1(x, y)p2(x, y)
         * p3.x p3.y p4.x p4.y
         * ...
         */
         friend std::ostream &operator<<(std::ostream &f, const Arrangement &a);

        /* read_guards method
        * :out param stream f: data stream from where the guards are input
        *
        * This method reads G guards and their coordinates from a file of the format
        * 
        * G             * number of guards
        * q1.x q1.y     * coordinates of the guard in the format of q(x, y) separated by spaces
        * q2.x q2.y
        * 
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

        /* init_guards method
        * :out param stream f: data stream from where the guards are input
        *
        * This method reads G guards from a file and initialises their positions greedily. This means that every guard is sequentially placed in an unseen area.
        * 
        */
        template<typename stream>
        void init_guards(stream &f, double learning_rate, double pull_attraction) {
            std::size_t n_guards;
            f >> n_guards;

            for (auto i = 0; i < n_guards; i ++) {
                if (i == 0) {
                    auto eit = *(this->input_arrangement).unbounded_face()->inner_ccbs_begin();
                    auto edge = Segment_2(eit->source()->point(), eit->target()->point());

                    Point_2 q((eit->source()->point().x() + eit->target()->point().x()) / 2, (eit->source()->point().y() + eit->target()->point().y()) / 2);
                    // std::cout << "inserting " << q << std::endl;
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

                    auto v = symmR.begin()->outer_boundary().edge(2);

                    Point_2 q((v.source().x() + v.target().x()) / 2, (v.source().y() + v.target().y()) / 2);

                    q = Point_2(symmR.begin()->outer_boundary().bottom_vertex()->x(), symmR.begin()->outer_boundary().bottom_vertex()->y());
                    // std::cout << "inserting " << q << std::endl;
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
        void print_polygon(stream &f);

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
        void print_guards(stream &f);

        /* add_guard method
        *  :in param Point_2 q:                         point that would guard the arrangement
        *  :in param Arrangement_2 visibility_region:   visibility region of the guard
        * 
        * This method adds the guard with its corresponding visibility region to the guard vector
        */
        void add_guard(const Point_2 q, const double learning_rate, double pull_attraction);

         /* add_reflex_vertices method
         *
         * This method populates the reflex vertices vector with the reflex vertices of the arrangement
         */
        void add_reflex_vertices();

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
        bool is_completely_visible(Arrangement_2 &visibility_arrangement);

        /* visibility method
        *  :in param Point_2 guard         :guard whose visibility region needs to be computed
        *  :return Arrangement_2           :arrangement visible from the guard
        * 
        *  This method computes the visibility region arrangement of a guard
        */
        Arrangement_2 visibility(const Point_2 guard);

        /* full_visibility method
        *  :return Arrangement_2           :arrangement visible from all guards
        * 
        *  This method computes the visibility region arrangement of all the guards
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 full_visibility();

        /* exclusive_visibility method
        *  :in param Guard guard         : guard whose visibility region should be excluded from the computation
        *  :return Arrangement_2           :arrangement visible from all guards except current guard
        * 
        *  This method computes the visibility region arrangement of all the guards except the current guard
        *  This is achieved by computing the visibility region arrangement for each of the guards placed in the arrangement, and overlaying them
        */
        Arrangement_2 exclusive_visibility(Guard guard, std::vector<Guard> zero_df_guards);

        /* is_visible_from method
        * :in param Point_2 p:  the viewpoint
        * :in param Point_2 r:  the point that needs to be checked for visibility from viewpoint p
        * :return bool:         true if r is visible from p,
        *                           false otherwise
        * 
        * This method checks whether a point r is visible from p by checking whether r is in the visibility region of p.
        */
        bool is_visible_from(const Point_2 p, const Point_2 r);


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
        std::vector<std::tuple<Point_2, Point_2, CGAL::Oriented_side>> reflex_vertex_pairs(const Guard g);

        /* exclusive_beta method
        * :in param Point_2 guard:                      the guard whose overlapping beta needs to be computed
        * :in param Point_2 intersection:               the intersection point on the polygon's boundary that the guard sees behind the reflex vertex
        * :in param Point_2 reflex_vertex:              the reflex vertex around which the guard needs to move
        * :in param double beta:                        the initial beta without any overlapping seen regions
        * :return double:                               the value of beta given that the guard's visibility region is also seen by other guards
        * 
        * This method computes the value of beta for the area seen exclusively by the guard.
        */
        double exclusive_beta(Guard guard, Point_2 intersection, Point_2 reflex_vertex, std::vector<Guard> zero_df_guards);

        /* point_behind_reflex_vertex method
        :in param Guard guard:             the guard who sees the reflex vertex
        :in param Point_2 reflex_vertex:   the reflex vertex behind which the angle is of interest
        :return double:                    the end of the segment behind the reflex vertex (unseen by the guard)

        This method computes the end of the segment (unseen by the guard) on the other side of the reflex vertex seen by a guard. 
        */
        Point_2 point_behind_reflex_vertex(const Guard guard, Point_2 reflex_vertex);

        /* gradient method
        * :in param Guard guard:                            guard point whose gradient needs to be computed
        * :in param vector<Guard> zero_df_guards:           the vector of guards who have a gradient of 0                  
        * :return Gradient:                                 object containing all gradient information for the guard
        * 
        * This method computes the gradient of a guard around all the reflex vertices it sees
        */
        Gradient gradient(const Guard g, std::vector<Guard> zero_df_guards);

        /* compute_new_coords method
        * :in param Guard prev_guard:               the previous coordinates of the guard
        * :in param Guard cur_guard:                the current coordinates of the guard to be computed
        * :in Gradient gradient:                    object containing all gradient information for computing the new coordinates of the current guard
        * :in bool placed:                          check whether the vertex can be placed on a reflex vertex
        * :return Guard:                            returns the guard with its updated coordinates
        *
        * This method computes the new coordinates of the current guard based on its gradient. After the theoretically computed coordinates, the method also checks whether the reflex area heuristic applies (and if so, does it). Another check is done for whether the new coordinates are still within the polygon boundaries. If they are not and it's possible, they are placed on the polygon boundary.
        */
        Guard compute_new_coords(Guard prev_guard, Guard cur_guard, Gradient gradient,  bool placed);
        /* optimise method
        *
        * This method optimises the position of all guards using gradient descent.
        * 
        * The algorithm steps are:
        *   - as long as the polygon is not fully seen:
        *       - as long as not all guards have a gradient, or the polygon is not completely seen
        *           - for each guard:
        *               - compute its visibility region for its current coordinates
        *               - update its visibility region for its current coordinates
        *               - compute the gradient and the pull for its current coordinates
        *               - if both the gradient and the pull are 0, skip the guard, otherwise
        *                   - update its coordinates based on its gradient
        *                   - if the new position of the guard is outside of the reflex area (when it's supposed to be inside it), place it on the boundary of the reflex area; if not possible, don't move the guard
        *                   - if the new position of the guard is outside of the polygon, place it on the boundary; if not possible, don't move the guard
        *                   - update guard in the guards vector
        *                   - if the guard has the same coords as another guard, don't move the guard and recompute its gradient without the pull (edge-case for when multiple guards are placed on top of the same reflex vertex and they cannot escape the reflex region)
        */
        void optimise();

        /* place_guard_on_boundary method
        * :in param Point_2 prev_guard:     the previous position of the guard
        * :in param Point_2 guard:          the guard position it should have according to the gradient
        * :out param Point_2 new_guard:     the new guard boundary position
        * :return bool:                     true if the guard can be placed on the boundary,
        *                                       false otherwise
        * 
        * This method computes the guard's position on the arrangement's boundary in the case when the gradient requires it to be outside of the polygon.
        * In most cases, the guard is placed at the intersection between its movement and the polygon's boundary. 
        * Sometimes it can happen that the gradient is exactly perpendicular to the boundary, which results in no movement. To tackle that case, the new position of the guard is projected onto the polygon segment. 
        * If the project falls outside of the polygon, then the guard is placed onto the closest segment end.
        */
        bool place_guard_on_boundary(Point_2 prev_guard, Point_2 guard, Point_2 &new_guard);

        /* place_guard_inside_reflex_area method
        * :in param Guard prev_guard:       the previous position of the guard
        * :in param Guard guard:            the guard position it should have according to the gradient
        * :out param Point_2 new_guard:     the new guard boundary position
        * :return bool:                     true if the guard can be placed inside the reflex line
        *                                       false otherwise
        * 
        * This method computes the guard's position inside the reflex area in the case when the gradient requires it to be outside of it
        * If the new guard's position is inside the reflex area, then its position is unchanged. Otherwise, it is projected onto the closest reflex line.
        */
        bool place_guard_inside_reflex_area(Guard prev_guard, Guard guard, Point_2 &new_guard);

         /* intersects_boundary method
         * :in param Segment_2 ray:       the segment corresponding to the guard's movement
         * :return bool:                  true if the guard is moved outside of the polygon (intersects its boundary)
         *                                   false otherwise
         * 
         * This method checks whether the guard is to be moved outside of the polygon and thus intersects the polygon boundary.
         */
        bool intersects_boundary(Segment_2 ray);

    private:
        Arrangement_2 input_arrangement;
        Polygon_2 input_polygon;
        // define type of visibility algorithm used
        TEV visibility_algo;
        // find the face of the guard
        CGAL::Arr_naive_point_location<Arrangement_2> pl;
        std::vector<Guard> guards;
        std::vector<Point_2> reflex_vertices;
        bool reflex_area = true, hidden_gradient = true, angle = true;
};

