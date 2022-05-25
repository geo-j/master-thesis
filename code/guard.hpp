#include <CGAL/Vector_2.h>

#include "utils.hpp"


typedef Kernel::Vector_2                                                            Vector_2;


class Guard {
    public:
        /* *****************
           *     I/O       *
           *****************
        */
        Guard(Point_2 g, Arrangement_2 vis) : prev_coords(g), cur_coords(g), visibility_region(vis) {
            this->area = compute_area(vis);
        }

        Guard(Point_2 g, Arrangement_2 vis, double alpha) : prev_coords(g), cur_coords(g), visibility_region(vis), learning_rate(alpha) {
            this->area = compute_area(vis);
        }

        // copy constructor
        Guard(const Guard &g) {
            this->prev_coords = g.prev_coords;
            this->cur_coords = g.cur_coords;
            this->visibility_region = g.visibility_region;
            this->area = g.area;
        }

        /* *****************
           *    update     *
           *****************
        */

       /* update_visibility method
       * :in Arrangement_2 visibility_region:  newly computed visibility region of current guard
       *
       * This method updates the visibility region and its area with the newly computed visibility region arrangement
       */
        void update_visibility(const Arrangement_2 visibility_region) {
            this->visibility_region.clear();
            this->visibility_region = visibility_region;
            this->area = compute_area(this->visibility_region);
        }

        /* update_coords method
        * :in Vector_2 gradient:    gradient based on which the new position of the guard has to be updated
        *
        * This method updates the guard position based on the gradient, and saves the previous position
        */
        void update_coords(Vector_2 gradient) {
            this->prev_coords = this->cur_coords;
            this->cur_coords = Point_2(this->prev_coords.x() + this->learning_rate * gradient.x(), this->prev_coords.y() + this->learning_rate * gradient.y());
        }

    private:
        Point_2 prev_coords, cur_coords;
        Arrangement_2 visibility_region;
        double area, learning_rate{0.5};
};