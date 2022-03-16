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
#include <CGAL/Ray_2.h>
#include <CGAL/intersections.h>

#include "utils.hpp"


typedef Arrangement_2::Face_handle                                          Face_handle;
typedef Arrangement_2::Edge_const_iterator                                  Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                              Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>    NSPV;
typedef CGAL::Rotational_sweep_visibility_2<Arrangement_2, CGAL::Tag_true>  RSV;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>              TEV;

// typedef CGAL::Ray_2<Traits_2>                                               Ray_2
typedef Kernel::Ray_2                                                       Ray_2;;
typedef Kernel::Intersect_2                                                 Intersect_2;



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



        /* *****************
           *   Gradient    *
           *****************
        */

        /* get_reflex_vertices method, as adapted from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArtGallery.cpp
        * 
        *  This method adds all the reflex vertices of the input arrangement in the reflex_vertices vector.
        */
        void get_reflex_vertices() {
            //Identify reflex vertices
            //use do/while for circular loop
            //The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
            auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();
            do {
                //Left turn, because the boundary is clockwise...
                if (CGAL::orientation(eit->prev()->source()->point(), eit->prev()->target()->point(), eit->target()->point()) == CGAL::LEFT_TURN) {
                    //eit->source() is a reflex vertex...
                    this->reflex_vertices.push_back(eit->source()->point());
                    // std::cout << "Reflex vertex (" << eit->source()->point().x() << ',' << eit->source()->point().y() << ')' << std::endl;
                    // eit->source()->data().isReflex = true;
                    // reflexHandles.push_back(eit->source());
                    // reflexNeighbours.push_back(
                    //     make_pair(eit->prev()->source(), eit->target()
                    //     ));

                }
            } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());
        }
    
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
        // TODO: account for intersection past the boundary of the polygon
        //       should I then first compute the visibility of a point?
        bool get_guard_reflex_arrangement_intersection(Point_2 guard, Point_2 reflex_vertex, Point_2 &d) {
            std::cout << "Ray [" << guard << ';' << reflex_vertex << "] " << std::endl;
            Ray_2 pr(guard, reflex_vertex);

            // only look at reflex vertices that are actually visible from the guard
            if (this->is_visible_from(guard, reflex_vertex)) {
                auto eit = *this->input_arrangement.unbounded_face()->inner_ccbs_begin();

                do {
                    Segment_2 edge = Segment_2(eit->source()->point(), eit->target()->point());
                    std::cout << "- tries to intersect boundary segment [" << eit->source()->point() << ", " << eit->target()->point() << "]" << std::endl;

                    auto intersection = CGAL::intersection(pr, edge);
                    if (intersection) {
                        d = *boost::get<Point_2 >(&*intersection);

                        // intersection point d should be visible from the guard, 
                        //      but also different than the reflex vertices (obviously the ray through the reflex vertex containts the reflex vertex itself)
                        if (d != reflex_vertex && this->is_visible_from(guard, d)) {
                                std::cout << "\tmanages to at " << d << std::endl;
                                return true;
                        }
                    }
                } while (++ eit != *this->input_arrangement.unbounded_face()->inner_ccbs_begin());
            }

            return false;
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
                std::cout << '\t' << p << " can see " << r << std::endl;
                return true;
            }

            return false;
        }


        void get_all_guard_reflex_boundary_intersections() {
            for (auto guard : this->guards) {
                for (auto reflex_vertex : this->reflex_vertices) {
                    Point_2 d;
                    if (this->get_guard_reflex_arrangement_intersection(guard, reflex_vertex, d))
                        std::cout << ">>> Ray [" << guard << "; " << reflex_vertex << "] intersects boundary point " << d << std::endl;
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