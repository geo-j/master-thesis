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

        Guard(Point_2 g, Arrangement_2 vis, double alpha, double beta) : cur_coords(g), visibility_region(vis), learning_rate(alpha), pull_attraction(beta) {
            this->area = compute_area(vis);
        }

        // copy constructor
        Guard(const Guard &g) {
            this->cur_coords = g.cur_coords;
            this->visibility_region = g.visibility_region;
            this->area = g.area;
            this->learning_rate = g.learning_rate;
            this->pull_attraction = g.pull_attraction;
            this->momentum = g.momentum;
            this->reflex_vertex = g.reflex_vertex;
            this->reflex_area = g.reflex_area;
            this->reflex_area_vertex = g.reflex_area_vertex;
        }

        // copy constructor
        Guard(const Guard &g, double alpha) {
            this->cur_coords = g.cur_coords;
            this->visibility_region = g.visibility_region;
            this->area = g.area;
            this->learning_rate = alpha;
            this->pull_attraction = g.pull_attraction;
            this->momentum = g.momentum;
            this->reflex_vertex = g.reflex_vertex;
            this->reflex_area = g.reflex_area;
            this->reflex_area_vertex = g.reflex_area_vertex;
        }

        // visibility region getter
        Arrangement_2 get_visibility_region() const {
            return this->visibility_region;
        }

        // current coordinates getter
        Point_2 get_coords() const {
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

        Point_2 get_reflex_area_vertex() const {
            return this->reflex_area_vertex;
        }

        bool is_reflex_vertex() const {
            return this->reflex_vertex;
        }

        bool is_in_reflex_area() const {
            return this->reflex_area;
        }

        // current coordinates setter
        void set_coords(Point_2 p) {
            this->cur_coords = Point_2(p);
        }

        void set_reflex_vertex(bool b) {
            this->reflex_vertex = b;
        }

        void set_in_reflex_area(bool b) {
            this->reflex_area = b;
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
            if (this->cur_coords == g.cur_coords) // && this->area == g.area)
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
        void update_coords(std::vector<Vector_2> gradients, std::vector<Vector_2> pulls, std::vector<Point_2> reflex_vertices, bool place_on_reflex_vertex) {
            // print gradients
            std::cout << gradients.size() << " " << (int) pulls.size() - 1 << " " << reflex_vertices.size() << std::endl;
            for (int i = 0; i < (int) pulls.size() - 1; i ++) {
                // std::cout << "h=" << this->pull_attraction * pulls.at(i) << std::endl;
                // std::cout << "Df=" << (this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(i) + this->pull_attraction * pulls.at(i))) * this->learning_rate << std::endl;
                std::cout << "Df=" << gradients.at(i) * this->learning_rate << std::endl;
                std::cout << "h=" << pulls.at(i) * this->learning_rate << std::endl;
            }

            bool placed = false;
            // compute the min distance between all reflex vertices seen by the guard
            double D = min_dist_reflex_vertices(reflex_vertices);
            // std::cout << D << std::endl;

            // if the guard is close enough (less than the min distance between 2 vertices) to a reflex vertex and far enough from the other one then save the reflex vertex to place the guard later on it
            for (int i = 0; i < reflex_vertices.size() && place_on_reflex_vertex; i ++) {
                auto d = distance(this->cur_coords, reflex_vertices.at(i));
                // std::cout << d << std::endl;
                // std::cout << pulls.at(i).squared_length() << " " << (2 / 3.0) * d << std::endl;
                if ((d < D / 2 || D == -1)
                    && pulls.at(i).squared_length() * this->learning_rate > (2 / 3.0) * d
                ) {
                    this->momentum = Vector_2(0, 0);

                    this->cur_coords = Point_2(reflex_vertices.at(i));
                    placed = true;
                    this->reflex_vertex = true;
                    this->reflex_area = true;
                    this->reflex_area_vertex = Point_2(this->cur_coords);
                    std::cout << "event=placed on reflex vertex " << this->cur_coords << std::endl;
                    // std::cout << "reflex vertex? " << this->is_reflex_vertex() << std::endl;

                    break;

                }
            }

            // if the guard wasn't placed on a reflex vertex, place it normally based on its momentum
            if (!placed) {
                // if (pulls.size() > 0)
                // if (CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) > CGAL::to_double(this->momentum.squared_length()) + 0.1) {
                //     std::cout << "pull reduced from " << CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) << " to ";
                //     pulls[pulls.size() - 1] = pulls.at(pulls.size() - 1) * CGAL::to_double(this->momentum.squared_length()) / CGAL::to_double(pulls.at(pulls.size() - 1).squared_length());
                //     std::cout << CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) << std::endl;

                // }
                    this->momentum = this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(gradients.size() - 1) + this->pull_attraction * pulls.at(pulls.size() - 1));
                // else
                //      this->momentum = this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(gradients.size() - 1));
                // std::cout << "Df=" << this->momentum * this->learning_rate << std::endl;

                this->cur_coords = Point_2(this->cur_coords + this->learning_rate * this->momentum);
                this->reflex_vertex = false;
            }

            // print last gradient
            if (pulls.size() > 0) {
                std::cout << "Df=" << this->momentum * this->learning_rate << std::endl;
                std::cout << "h=" << pulls.at(pulls.size() - 1) * this->learning_rate << std::endl;
            }
                
        }

    private:
        Point_2 cur_coords, reflex_area_vertex;
        Arrangement_2 visibility_region;
        double area, learning_rate{0.5}, gamma{0.5}, pull_attraction{1};
        Vector_2 momentum{0, 0};
        bool reflex_vertex = false, reflex_area = false;
};