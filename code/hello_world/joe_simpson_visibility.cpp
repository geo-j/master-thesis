#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include <stack>


typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
// typedef typename CGAL::Geometry_traits_2                     Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef typename CGAL::Arrangement_2<Traits_2>                          Arrangement_2;
typedef typename Arrangement_2::Geometry_traits_2                       Geometry_traits_2;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Polygon_2<Kernel>                                         Polygon_2;
typedef typename Geometry_traits_2::Ray_2                               Ray_2;
typedef typename Geometry_traits_2::Object_2                            Object_2;
typedef typename Kernel::Intersect_2                                    Intersect_2;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Ccb_halfedge_circulator                          Ccb_halfedge_circulator;
typedef typename Arrangement_2::Ccb_halfedge_const_circulator           Ccb_halfedge_const_circulator;
typedef typename Arrangement_2::Face_const_handle                       Face_const_handle;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2>         Arr_point_location;
typedef typename Arr_point_location::result_type                        Location_result;
typedef typename Arrangement_2::Halfedge_const_handle                   Halfedge_const_handle;
typedef typename Geometry_traits_2::Line_2                              Line_2;
typedef typename Arrangement_2::Vertex_const_handle                     Vertex_const_handle;
typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;

typedef std::vector<Point_2>                                            Vertex_container;
typedef typename Vertex_container::size_type                            Size_type;


Vertex_container vertices;
std::stack<Point_2> S;
Geometry_traits_2 *traits;
Arr_point_location point_location;


enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} oper;

/*! Scan edges v_i,v_{i+1},...,v_n, until find an edge intersecting given ray
    or given segment. is_ray = true -> ray, false -> segment.
    The intersection point is returned by u */
Size_type scan_edges( Size_type i, Point_2 ray_begin, Point_2 ray_end, Point_2& u, bool is_ray) {
    CGAL::Orientation old_orientation = CGAL::RIGHT_TURN;
    Ray_2 ray(ray_begin, ray_end);
    Segment_2 ray_segment(ray_begin, ray_end);
    Size_type k;
    Object_2 result;
    for (k = i; k + 1 < vertices.size(); k ++) {
        CGAL::Orientation curr_orientation = traits->orientation_2_object()(ray_begin, ray_end, vertices[k + 1]);
        if (curr_orientation != old_orientation) {
            // Orientation switch, an intersection may occur
            Segment_2 seg(vertices[k], vertices[k + 1]);
            if (is_ray) {
                result = Intersect_2()(seg, ray);
                if(result)
                    break;
            } else {
                result = Intersect_2()(seg, ray_segment);
                if(result)
                    break;
            }
        }
        old_orientation = curr_orientation;
    }
    const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&result);
    if (intersection_point) {
            u = *intersection_point;
    } else {
            u = vertices[k + 1];
    }
    return k;
}

void left(Size_type &i, Point_2 &w, Point_2 p) {
    if (i == vertices.size() - 1) {
        oper = FINISH;
    } else {
        Point_2 tos = S.top();
        S.pop();
        Point_2 prev_tos = S.top();
        S.push(tos);
        CGAL::Orientation orientation_p = traits->orientation_2_object()(p, vertices[i], vertices[i + 1]);

        if (orientation_p != CGAL::RIGHT_TURN) {
            oper = LEFT;
            w = vertices[i + 1];
            S.push(w);
            i ++;
        } else {
            CGAL::Orientation orientation_prev_tos = traits->orientation_2_object()(prev_tos, vertices[i], vertices[i + 1]);

            if (orientation_prev_tos == CGAL::RIGHT_TURN) {
                oper = SCANA;
                w = vertices[i + 1];
                i ++;
            } else {
                oper = RIGHT;
                w = vertices[i];
                i ++;
            }
        }
    }
}

void right(Size_type &i, Point_2 w, Point_2 p) {
    Point_2 u, tos, prev_tos = S.top();
    CGAL::Orientation orientation_tos, orientation_prev_tos = traits->orientation_2_object()(p, prev_tos, vertices[i]);
    int mode = 0;

    while (S.size() > 1) {
        tos = prev_tos;
        orientation_tos = orientation_prev_tos;
        S.pop();
        prev_tos = S.top();
        orientation_prev_tos = traits->orientation_2_object()(p, prev_tos, vertices[i]);

        if (orientation_tos != CGAL::LEFT_TURN && orientation_prev_tos != CGAL::RIGHT_TURN) {
            mode = 1;
            break;
        }

        Segment_2 seg_v(vertices[i - 1], vertices[i]), seg_tos(prev_tos, tos);

        if (vertices[i - 1] != tos) {
            Object_2 intersection = Intersect_2()(seg_v, seg_tos);

            if (intersection) {
                const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&intersection);
                u = *intersection_point;
                mode = 2;
                break;
            }
        }
    }

    if (mode == 1) {
        orientation_tos = traits->orientation_2_object()(p, vertices[i], vertices[i + 1]);
        orientation_prev_tos = traits->orientation_2_object()(vertices[i - 1], vertices[i], vertices[i + 1]);

        if (orientation_tos == CGAL::RIGHT_TURN) {
            oper = RIGHT;
            S.push(tos);
            w = vertices[i];
            i ++;
        } else if (orientation_prev_tos == CGAL::RIGHT_TURN) {
            Ray_2 ray(p, vertices[i]);
            Segment_2 seg(prev_tos, tos);
            Object_2 intersection = Intersect_2()(seg, ray);
            const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&intersection);
            u = *intersection_point;

            if (S.top() != u) {
                S.push(u);
            }

            oper = LEFT;
            w = vertices[i + 1];
            S.push(vertices[i]);
            S.push(w);
            i ++;
        } else {
            Ray_2 ray(p, vertices[i]);
            Segment_2 seg(prev_tos, tos);
            Object_2 intersection = Intersect_2()(seg, ray);
            const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&intersection);
            u = *intersection_point;

            if (S.top() != u) {
                S.push(u);
            }

            oper = SCANC;
            w = vertices[i];
            i ++;
        }
    } else if (mode == 2) {
            // Case R4
            oper = SCAND;
            w = u;
    }
}

void scana(Size_type& i, Point_2& w, const Point_2& q) {
// Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    Point_2 u;
    Size_type k = scan_edges( i, q, S.top(), u, true);

    CGAL::Orientation orientation = traits->orientation_2_object()(q, vertices[k], vertices[k + 1]);

    if (orientation == CGAL::RIGHT_TURN) {
        bool forward = traits->collinear_are_ordered_along_line_2_object()(q, S.top(), u);

        if (!forward) {
        // Case A1
            oper = RIGHT;
            i = k + 1;
            w = u;
        } else {
        // Case A2
            oper = SCAND;
            i = k+1;
            w = u;
        }
    } else {
        // Case A3
        oper = LEFT;
        i = k + 1;
        S.push(u);
        w = vertices[k + 1];
        if (u != w) {
            S.push(w);
        }
    }
}

/*! Find the first edge interecting the segment (v_0, s_t) */
void scanb(Size_type& i, Point_2& w) {
    if (i == vertices.size() - 1) {
        oper = FINISH;
        return;
    }
    Point_2 u;
    Size_type k = scan_edges(i, S.top(), vertices[0], u, false);
    if ((k + 1 == vertices.size() - 1) && (vertices[0] == u)) {
        // Case B1
        oper = FINISH;
        S.push(vertices[0]);
    } else {
        // Case B2
        oper = RIGHT;
        i = k + 1;
        w = u;
    }
}

/*! Finds the exit from a general front hidden window by finding the first
vertex to the right of the ray defined by the query_point and w*/
void scanc(Size_type& i, Point_2& w) {
    Point_2 u;
    Size_type k = scan_edges(i, S.top(), w, u, false);
    oper = RIGHT;
    i = k + 1;
    w = u;
}

  /*! find the first edge intersecting the given window (s_t, w) */
void scand(Size_type& i, Point_2& w) {
    Point_2 u;
    Size_type k = scan_edges(i, S.top(), w, u, false);
    oper = LEFT;
    i = k + 1;
    S.push(u);
    if (u != vertices[k + 1]) {
        S.push(vertices[k + 1]);
    }
    w = vertices[k + 1];
}
  
void output(const Point_2& p, Arrangement_2& out_arr) {

    // if(inserted_artificial_starting_vertex)
    //   stack.pop();

    std::vector<Point_2> points;
    while(!S.empty()) {
        Point_2& tos = S.top();
        if (tos != p /*|| query_pt_is_vertex*/) {
            points.push_back(tos);
        }
        S.pop();
    }

    // if(inserted_artificial_starting_vertex) {
    //     points.back() = points[0];
    //     inserted_artificial_starting_vertex = false;
    // }
    points.pop_back();
      //std::cout << " ordanary " << std::endl; 	
    typename Arrangement_2::Vertex_handle v_last, v_first;
    v_last = v_first = out_arr.insert_in_face_interior(points[0], out_arr.unbounded_face());
	
    for(unsigned int i = 0; i < points.size() - 1; i ++){
        if(points[i] < points[i + 1]){
            v_last = out_arr.insert_from_left_vertex (Segment_2(points[i], points[i + 1]), v_last)->target();
        } else {
            v_last = out_arr.insert_from_right_vertex(Segment_2(points[i], points[i + 1]), v_last)->target();
        }        
    }
    out_arr.insert_at_vertices(Segment_2(points.front(), points.back()), v_last, v_first);
}

  /*! Finds a visible vertex from the query point 'q' in 'face' 
    to start the algorithm from*/
Ccb_halfedge_const_circulator find_visible_start(Face_const_handle face, Point_2 &q) {
    Location_result result = point_location.ray_shoot_up(q);

    if (const Halfedge_const_handle* e = boost::get<Halfedge_const_handle>(&(result))) {
	    Point_2 p(q.x(), traits->compute_y_at_x_2_object()(Line_2((*e)->source()->point(), (*e)->target()->point()), q.x()));

        vertices.push_back(p);
        // inserted_artificial_starting_vertex = true;

        return (*e)->next()->ccb();
    } else if (const Vertex_const_handle* v = boost::get<Vertex_const_handle>(&(result))) {
	    Halfedge_around_vertex_const_circulator cir = (*v)->incident_halfedges();

	while(face != cir->face()) {
	  ++cir;
	}
	return cir->next()->ccb();
      }
    else
      {
	CGAL_assertion_msg(false, "Should not be reachable.");
	return Ccb_halfedge_const_circulator();
      }
  }

Arrangement_2 compute_visibility_polygon(Point_2 p) {
    /*
    Orientation_2 operator()(Point_2 p, Point_2 q, Point_2 r) that returns CGAL::LEFT_TURN, if r lies to the left of the oriented line l defined by p and q,
                                                                           CGAL::RIGHT_TURN if r lies to the right of l, and 
                                                                           CGAL::COLLINEAR if r lies on l. 
   */
    Ccb_halfedge_const_circulator circ = find_visible_start(face, q);
    Ccb_halfedge_const_circulator curr = circ;

    do {
      vertices.push_back(curr->source()->point());
    } while(++curr != circ);

    vertices.push_back(vertices[0]);

    CGAL::Orientation orientation = traits->orientation_2_object()(p, vertices[0], vertices[1]);  // check where vertices[1] (v_1) lies w.r.t. the line pvertices[0] (pv_0)
    Point_2 w;
    Size_type i = 0;
    Arrangement_2 out_arr;

    // if (orientation != CGAL::RIGHT_TURN) {
    //     oper = LEFT;
    //     i = 1;
    //     S.push(vertices[0]);
    //     S.push(vertices[1]);
    // } else {
    //     oper = SCANA;
    //     i = 1;
    //     w = vertices[1];
    //     S.push(vertices[0]);
    // }

    // do {
    //     switch(oper) {
    //         case LEFT:
    //             left(i, w, p);
    //             break;
    //         case RIGHT:
    //             right(i, w, p);
    //             break;
    //         case SCANA:
    //             scana(i, w, p);
    //             break;
    //         case SCANB:
    //             scanb(i, w);
    //             break;
    //         case SCANC:
    //             scanc(i, w);
    //             break;
    //         case SCAND:
    //             scand(i, w);
    //             break;
    //     }

    //     Point_2 prev_tos = S.top();
    //     S.pop();

    //     if (oper == LEFT && 
    //         // check if s_{t - 1}s_t intersect pv_0
    //         traits->orientation_2_object()(p, vertices[0], S.top()) == CGAL::RIGHT_TURN && // check if s_t is on the right of v0
    //         traits->orientation_2_object()(p, vertices[0], prev_tos) == CGAL::LEFT_TURN) { // check if s_{t - 1} is on the left of v0
    //             Point_2 tos = S.top();
    //             S.pop();

    //             Segment_2 seg(S.top(), tos);
    //             Ray_2 ray_origin(p, vertices[0]);

    //             if (Object_2 result = Intersect_2()(seg, ray_origin)) {
    //                 const Point_2 *intersection_point = CGAL::object_cast<Point_2>(&result);
    //                 tos = *intersection_point;

    //                 oper = SCANB;
    //             }

    //         S.push(tos);
    //     }

    // } while(oper != FINISH);

    // output(p, out_arr);

    return out_arr;
}

int main() {
    Polygon_2 P;
    Point_2 p1(0, 4), p2(0, 0), p3(3, 2), p4(4, 0), p5(4, 4), p6(1, 2);
    P.push_back(p1);
    P.push_back(p2);
    P.push_back(p3);
    P.push_back(p4);
    P.push_back(p5);
    P.push_back(p6);

    std::vector<Segment_2> segments;
    segments.push_back(Segment_2(p1, p2));
    segments.push_back(Segment_2(p2, p3));
    segments.push_back(Segment_2(p3, p4));
    segments.push_back(Segment_2(p4, p5));
    segments.push_back(Segment_2(p5, p6));
    segments.push_back(Segment_2(p6, p1));
    Arrangement_2 env;
    CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());

    // find the face of the query point
    // (usually you may know that by other means)
    Point_2 p(0.5, 2);
    Arrangement_2::Face_const_handle * face;
    CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
    CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(p);

    // The query point locates in the interior of a face
    face = boost::get<Arrangement_2::Face_const_handle> (&obj);

    std::cout << "here";
    Arrangement_2 output = compute_visibility_polygon(p);

    // for (Edge_const_iterator eit = output.edges_begin(); eit != output.edges_end(); ++ eit)
    //     std::cout << "[" << eit->source()->point() << " -> " << eit->target()->point() << "]" << std::endl;

    return 0;
}