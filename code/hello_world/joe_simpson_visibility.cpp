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

enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} oper;


Arrangement_2 compute_visibility_polygon(Point_2 p) {
    Geometry_traits_2 *traits;
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