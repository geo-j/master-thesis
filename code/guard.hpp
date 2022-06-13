#include <CGAL/Vector_2.h>

#include "utils.hpp"


typedef Kernel::Vector_2                                                            Vector_2;


class Guard {
    public:
        /* *****************
           *     I/O       *
           *****************
        */
        Guard(Point_2 g, Arrangement_2 vis) : cur_coords(g), visibility_region(vis) {
            this->area = compute_area(vis);
        }

        Guard(Point_2 g, Arrangement_2 vis, double alpha) : cur_coords(g), visibility_region(vis), learning_rate(alpha) {
            this->area = compute_area(vis);
        }

        // copy constructor
        Guard(const Guard &g) {
            // this->prev_coords = g.prev_coords;
            this->cur_coords = g.cur_coords;
            this->visibility_region = g.visibility_region;
            this->area = g.area;
            this->learning_rate = g.learning_rate;
            this->momentum = g.momentum;
        }

        // copy constructor
        Guard(const Guard &g, double alpha) {
            // this->prev_coords = g.prev_coords;
            this->cur_coords = g.cur_coords;
            this->visibility_region = g.visibility_region;
            this->area = g.area;
            this->learning_rate = alpha;
            this->momentum = g.momentum;
        }

        // visibility region getter
        Arrangement_2 get_visibility_region() const {
            return this->visibility_region;
        }

        // current coordinates getter
        Point_2 get_cur_coords() const {
            return this->cur_coords;
        }

        // area getter
        double get_area() const {
            return this->area;
        }

        double get_learning_rate() const {
            return this->learning_rate;
        }

        Vector_2 get_momentum() const {
            return this->momentum;
        }

        // current coordinates setter
        void set_cur_coords(Point_2 p) {
            this->cur_coords = Point_2(p);
        }

        // learning rate setter
        void set_learning_rate(double alpha) {
            this->learning_rate = alpha;
        }

        /* overloaded equality operator
        *
        * The equality is established if the previous, and current coordinates, as well as the area seen are the same. 
        * The visibility regions equality is not checked due to a lack of arrangement equality check from CGAL.
        */
        bool operator ==(Guard g) {
            if (this->cur_coords == g.cur_coords && this->area == g.area)
                return true;
            
            return false;
        }

        // overloaded inequality operator
        bool operator !=(Guard g) {
            return !operator==(g);
        }

        // TODO: should probably output all guard info
        /* overloaded output operator
        *
        * The output operator only outputs the current position of the guard
        */
        friend std::ostream &operator<<(std::ostream &f, const Guard &g) {
            f << g.cur_coords;

            return f;
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
        void update_coords(Vector_2 gradient, std::vector<Point_2> reflex_vertices) {
            if (reflex_vertices.size() > 1) {
                this->cur_coords = Point_2(reflex_vertices.at(1).x() - 0.001, reflex_vertices.at(1).y());
            } else {
                this->momentum = this->gamma * this->momentum + (1 - this->gamma) * gradient;

                this->cur_coords = Point_2(this->cur_coords + this->learning_rate * this->momentum);
                // std::cout << this->momentum << std::endl;
            }
        }

    private:
        Point_2 cur_coords;
        Arrangement_2 visibility_region;
        double area, learning_rate{0.5}, gamma{0.9};
        Vector_2 momentum{0, 0};
};