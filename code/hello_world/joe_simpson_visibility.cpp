#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>

#include <stack>


typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
// typedef typename CGAL::Geometry_traits_2                     Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef typename CGAL::Arrangement_2<Traits_2>                          Arrangement_2;
typedef typename Arrangement_2::Geometry_traits_2                       Geometry_traits_2;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef typename Geometry_traits_2::Ray_2                               Ray_2;
typedef typename Geometry_traits_2::Object_2                            Object_2;
typedef typename Kernel::Intersect_2                                    Intersect_2;



typedef std::vector<Point_2>                                            Vertex_container;
typedef typename Vertex_container::size_type                            Size_type;


Vertex_container vertices;
std::stack<Point_2> S;
Geometry_traits_2 *traits;


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
    Point_2 *intersection_point = object_cast<Point_2>(&result);
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
                Point_2 *intersection_point = object_cast<Point_2>(&intersection);
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
            Point_2 *intersection_point = object_cast<Point_2>(&intersection);
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
            Point_2 *intersection_point = object_cast<Point_2>(&intersection);
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

Arrangement_2 compute_visibility_polygon(Point_2 p) {
    /*
    Orientation_2 operator()(Point_2 p, Point_2 q, Point_2 r) that returns CGAL::LEFT_TURN, if r lies to the left of the oriented line l defined by p and q,
                                                                           CGAL::RIGHT_TURN if r lies to the right of l, and 
                                                                           CGAL::COLLINEAR if r lies on l. 
   */
    CGAL::Orientation orientation = traits->orientation_2_object()(p, vertices[0], vertices[1]);  // check where vertices[1] (v_1) lies w.r.t. the line pvertices[0] (pv_0)
    Point_2 w;
    Size_type i = 0;

    if (orientation != CGAL::RIGHT_TURN) {
        oper = LEFT;
        i = 1;
        S.push(vertices[0]);
        S.push(vertices[1]);
    } else {
        oper = SCANA;
        i = 1;
        w = vertices[1];
        S.push(vertices[0]);
    }

    do {
        switch(oper) {
            case LEFT:
                left(i, w, p);
                break;
            case RIGHT:
                right(i, w, p);
                break;
            case SCANA:
                scana(i, w, p);
                break;
            case SCANB:
                scanb(i, w, p);
                break;
            case SCANC:
                scanc(i, w, p);
                break;
            case SCAND:
                scand(i, w, p);
                break;
        }

        Point_2 prev_tos = S.top();
        S.pop();

        if (oper == LEFT && 
            // check if s_{t - 1}s_t intersect pv_0
            traits->orientation_2_object()(p, vertices[0], S.top()) == CGAL::RIGHT_TURN && // check if s_t is on the right of v0
            traits->orientation_2_object()(p, vertices[0], prev_tos) == CGAL::LEFT_TURN) { // check if s_{t - 1} is on the left of v0
                Point_2 tos = S.top();
                S.pop();

                Segment_2 seg(S.top(), tos);
                Ray_2 ray_origin(p, vertices[0]);

                if (Object_2 result = Intersect_2()(seg, ray_origin)) {
                    Point_2 *intersection_point = object_cast<Point_2>(&result);
                    tos = *intersection_point;

                    oper = SCANB;
                }

            S.push(tos);
        }

    } while(oper != FINISH);
}

int main() {


    return 0;
}