#include "utils.hpp"
#include <CGAL/Vector_2.h>

typedef Kernel::Vector_2                                                            Vector_2;

class Gradient {
    public:
        Gradient(std::vector<Vector_2> g, std::vector<Vector_2> p, std::vector<Point_2> r) : gradients(g), pulls(p), reflex_vertices(r) {};

        std::vector<Vector_2> get_gradients() {
            return this->gradients;
        }

        std::vector<Vector_2> get_pulls() {
            return this->pulls;
        }

        std::vector<Point_2> get_reflex_vertices() {
            return this->reflex_vertices;
        }

    private:
        std::vector<Vector_2> gradients, pulls;
        std::vector<Point_2> reflex_vertices;
};