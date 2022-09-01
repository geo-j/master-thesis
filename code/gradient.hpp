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

        void scale_gradient(float factor) {
            // std::cout << "gradient before " << this->gradients[this->gradients.size() - 1] << std::endl;
            this->gradients[this->gradients.size() - 1] *= factor;
            // std::cout << "gradient after " << this->gradients[this->gradients.size() - 1] << std::endl;

        }

    private:
    
        /* vector<Vector_2> gradients:                         vector of all the gradients of the guard around n reflex vertices it sees: 
        *                                                           indices 0, ..., n contain the gradient as computed for each reflex vertex
        *                                                           index n + 1 contains the sum of the gradients
        * 
        * vector<Vector_2> pulls:                              vector of all the pulls of the guard toward the n reflex vertices it sees:
        *                                                           indices 0, ..., n contain the gradient as computed for each reflex vertex
        *                                                           index n + 1 contains the sum of the gradients
        * 
        * vector<Point_2> reflex_vertices:                     vector of all the reflex vertices the guard sees
        */
        std::vector<Vector_2> gradients, pulls;
        std::vector<Point_2> reflex_vertices;
};