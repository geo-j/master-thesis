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
        void update_coords(std::vector<Vector_2> gradients, std::vector<Vector_2> pulls, std::vector<Point_2> reflex_vertices) {
            // print gradients
            for (auto i = 0; i < gradients.size() - 1; i ++) {
                std::cout << "Df=" << (this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(i) + this->pull_attraction * pulls.at(i))) * this->learning_rate << std::endl;
            }

            bool placed = false;


            // compute the min distance between all reflex vertices seen by the guard
            double D = min_dist_reflex_vertices(reflex_vertices);

            // if the guard is close enough (less than the min distance between 2 vertices) to a reflex vertex, then save the reflex vertex to place the guard later on it
            for (auto i = 0; i < reflex_vertices.size(); i ++) {
                if (distance(this->cur_coords, reflex_vertices.at(i)) < D / 2
                    && pulls.at(i).squared_length() > (2 / 3) * D
                ) {
                    this->momentum = Vector_2(0, 0);

                    this->cur_coords = Point_2(reflex_vertices.at(i));
                    placed = true;
                    std::cout << "event=placed on reflex vertex " << this->cur_coords << std::endl;

                    break;

                }
            }

            // place the guard on the reflex vertex, if the pull makes it go past it.
            // if (reflex_vertices.size() > 0) {
            //     auto j = -1;
            //     for (auto reflex_vertex : reflex_vertices) {
            //         while (j < pulls.size()) {
            //             j ++;

            //             if (distance(this->cur_coords, Point_2(this->cur_coords + pulls.at(j))) >= distance(this->cur_coords, reflex_vertex) 
            //             // &&  (pulls.at(j) * pulls.at(pulls.size() - 1)) / (pulls.at(j).squared_length() * pulls.at(pulls.size() - 1).squared_length()) > 0.8
            //             ) {
            //                 // update the momentum based on the move to the reflex vertex
            //                 this->momentum = Vector_2(0, 0);

            //                 // std::cout << "Df=" << this->momentum * this->learning_rate << std::endl;

            //                 // create a vector between the guard and the reflex vertex, s.t. we can get a coordinate close enough to the reflex vertex that is not the reflex vertex; otherwise visibility doesn't work.
            //                 auto new_reflex = this->momentum * this->learning_rate;
            //                 // std::cout << new_reflex << std::endl;

            //                 this->cur_coords = Point_2(this->cur_coords.x() + new_reflex.x(), this->cur_coords.y() + new_reflex.y());
            //                 placed = true;
            //                 std::cout << "event=placed on reflex vertex " << reflex_vertex << " with actual coords " << this->cur_coords << std::endl;
            //                 break;
            //             }
            //         }
            //     }
            // }

            // if the guard wasn't placed on a reflex vertex, place it normally based on its momentum
            if (!placed) {
                this->momentum = this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(gradients.size() - 1) + this->pull_attraction * pulls.at(pulls.size() - 1));
                // std::cout << "Df=" << this->momentum * this->learning_rate << std::endl;

                this->cur_coords = Point_2(this->cur_coords + this->learning_rate * this->momentum);
            }

            // print last gradient
            std::cout << "Df=" << this->momentum * this->learning_rate << std::endl;

        }

    private:
        Point_2 cur_coords;
        Arrangement_2 visibility_region;
        double area, learning_rate{0.5}, gamma{0.9}, pull_attraction{4};
        Vector_2 momentum{0, 0};
};