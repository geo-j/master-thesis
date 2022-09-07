#include "gradient.hpp"



class Guard {
    public:
        /* *****************
           *     I/O       *
           *****************
        */

        Guard(){};
        
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
            this->prev_coords = g.prev_coords;
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
            this->prev_coords = g.prev_coords;

        }

        // visibility region getter
        Arrangement_2 get_visibility_region() const {
            return this->visibility_region;
        }

        Point_2 get_prev_coords() const {
            return this->prev_coords;
        }
        // current coordinates getter
        Point_2 get_coords() const {
            return this->cur_coords;
        }

        // area getter
        double get_area() const {
            return this->area;
        }

        // learning rate getter
        double get_learning_rate() const {
            return this->learning_rate;
        }

        // momentum getter
        Vector_2 get_momentum() const {
            return this->momentum;
        }

        // reflex area vertex getter
        Point_2 get_reflex_area_vertex() const {
            return this->reflex_area_vertex;
        }

        // whether guard is on reflex vertex checker
        bool is_reflex_vertex() const {
            return this->reflex_vertex;
        }

        // whether guard is in the reflex area checker
        bool is_in_reflex_area() const {
            return this->reflex_area;
        }

        // current coordinates setter
        void set_coords(Point_2 p) {
            this->cur_coords = Point_2(p);
        }

        // reflex vertex boolean setter
        void set_reflex_vertex(bool b) {
            this->reflex_vertex = b;
        }

        // reflex area boolean setter
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
        * :in vector<Vector_2> gradients:       gradients based on which the new position of the guard has to be updated; the last index is the full gradient, whereas the previous indices are the partial gradients for each reflex vertex the guard sees
        * :in vector<Vector_2> pulls:           pulls based on which the new position of the guard has to be updated; the last index is the full pull, whereas the previous indices are the partial pulls for each reflex vertex the guard sees
        * :in vector<Point_2> reflex_vertices:  the reflex vertices the guard sees; their index corresponds to the gradient and pull index in the other vertices
        * :in bool place_on_reflex_vertex:      check whether a guard is allowed to be placed on a reflex vertex in case of a pull (tackles the edge-case for multiple guards being pulled towards the same reflex vertex)
        *
        * This method updates the guard position based on the gradient, and saves the previous position
        */
        void update_coords(Gradient gradient, bool place_on_reflex_vertex) {
            if (!place_on_reflex_vertex)
                std::cout << "\tguard not allowed to use pull\n";
            // unpack the gradient
            auto gradients = gradient.get_gradients();
            auto pulls = gradient.get_pulls();
            auto reflex_vertices = gradient.get_reflex_vertices();
            // print gradients
            // std::cout << gradients.size() << " " << (int) pulls.size() - 1 << " " << reflex_vertices.size() << std::endl;


            bool placed = false;
            // compute the min distance between all reflex vertices seen by the guard

            double D;
            if (place_on_reflex_vertex && this->pull_onto)
                D = min_dist_reflex_vertices(reflex_vertices);
            // std::cout << D << std::endl;

            // if the guard is close enough (less than the min distance between 2 vertices) to a reflex vertex and far enough from the other one then save the reflex vertex to place the guard later on it
            for (int i = 0; i < reflex_vertices.size() && place_on_reflex_vertex && this->pull_onto; i ++) {
                auto d = distance(this->cur_coords, reflex_vertices.at(i));
                // std::cout << d << std::endl;
                // std::cout << pulls.at(i).squared_length() << " " << (2 / 3.0) * d << std::endl;
                if ((d < D / 2 || D == -1)
                    && pulls.at(i).squared_length() * this->learning_rate > (2 / 3.0) * d
                ) {
                    this->momentum = Vector_2(0, 0);
                    this->prev_coords = Point_2(cur_coords);
                    this->cur_coords = Point_2(reflex_vertices.at(i));
                    placed = true;
                    this->reflex_vertex = true;
                    this->reflex_area = true;
                    this->reflex_area_vertex = Point_2(this->cur_coords);
                    std::cout << "event=placed on reflex vertex " << this->cur_coords << std::endl;

                    break;

                }
            }

            // if the guard wasn't placed on a reflex vertex, place it normally based on its momentum
            if (!placed) {
                // cap the pull if too large compared to the momentum
                if (this->pull_capping &&
                    CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) > CGAL::to_double(this->momentum.squared_length()) + 0.1) {
                    std::cout << "pull reduced from " << CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) << " to ";
                    pulls[pulls.size() - 1] = pulls.at(pulls.size() - 1) * CGAL::to_double(this->momentum.squared_length()) / CGAL::to_double(pulls.at(pulls.size() - 1).squared_length());
                    // std::cout << CGAL::to_double(pulls.at(pulls.size() - 1).squared_length()) << std::endl;

                }
                    this->momentum = this->gamma * this->momentum + (1 - this->gamma) * (gradients.at(gradients.size() - 1) + this->pull_attraction * pulls.at(pulls.size() - 1));

                this->prev_coords = Point_2(cur_coords);
                this->cur_coords = Point_2(this->cur_coords + this->learning_rate * this->momentum);
                this->reflex_vertex = false;
            }
                
        }

    private:
        // save the current coords and the reflex vertex in whose reflex area the guard is in
        Point_2 cur_coords, reflex_area_vertex, prev_coords;
        Arrangement_2 visibility_region;
        // set gamma = 1 for not using momentum
        // set pull_attraction = 0 for not using the pull
        double area, learning_rate{0.5}, gamma{0.8}, pull_attraction{1};
        Vector_2 momentum{0, 0};
        bool reflex_vertex = false, reflex_area = false, pull_onto = true, pull_capping = true;
};