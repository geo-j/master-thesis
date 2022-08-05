#include <vector>
#include <string>

#include <CGAL/Polygon_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel                   Kernel;
typedef CGAL::Polygon_2<Kernel>                                             Polygon_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Segment_2                                                   Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                                  Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                       Arrangement_2;

typedef Kernel::FT                                                          FT;
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
typedef CGAL::Cartesian<NT> K;



template<typename type>
void push_back_unique(std::vector<type> &v, type element) {
    auto it = std::find(v.begin(), v.end(), element);

    if (it == v.end())
        v.push_back(element);
}

/* distance function
* :in param Point_2 p1:		source point
* :in param Point_2 p2:		destination point
* :return double:			square distance between the two input points
*
* This function computes the square distance p1^2 + p2^2 between two points p1 and p2.
*/
inline double distance(Point_2 p1, Point_2 p2) {
	return CGAL::to_double(((p1.x() - p2.x()) * (p1.x() - p2.x()) + (p1.y() - p2.y()) * (p1.y() - p2.y())).exact());
}

/* min_dist_reflex_vertices method
* :in param std::vector<Point_2> reflex_vertices:       the vector of all reflex vertices seen by a guard
* :return double:                                       the minimum distance between all pairs of reflex vertices seen by a guard
* 
* This method computes the minimum distance between all pairs of reflex vertices seen by a guard
*/
inline double min_dist_reflex_vertices(std::vector<Point_2> reflex_vertices) {
	double D = -1;
	for (auto q : reflex_vertices)
		for (auto r : reflex_vertices)
			if (q != r) {
				// std::cout << "more vertices\n";
				if (D == -1)
					D = distance(q, r);
				else {
					double d = distance(q, r);
					if (d < D)
						D = d;
				}
			}
	
	return D;
}

/* get_number function
* :in param string s:	string that needs to be converted to a double
* :return double:		input string converted to a double
*
* This function checks whether the input string is a number or a fraction. If it's a fraction, it takes the numerator and denominator on each side of the fraction and returns their fraction value.
*/
inline double get_number(std::string s) {
	double x, y;
	size_t frac = s.find('/');

	// check if input string is fraction, and if so, take the numerator and denominator on both sides of the fraction and return their result.
	if (frac != std::string::npos) {
		x = stod(s.substr(0, frac));
		y = stod(s.substr(frac + 1, -1));
		
		return x / y;
	} else
		return stod(s);
}


// TODO: should probably take into account polygons with holes
/* arrangement to polygon function, as adapted from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArrangementFunctions.cpp
* :in param    Arrangement_2 arrangement: input arrangement to be converted to a polygon
* :return   Polygon_2     polygon:        output polygon converted from the given arrangement
*/
inline Polygon_2 arrangement_to_polygon(Arrangement_2 &arrangement) {
	std::vector<Point_2> vertices;

	// loop around the arrangement
	auto eit = *arrangement.unbounded_face()->inner_ccbs_begin();
	do {
		vertices.push_back(eit->source()->point());
	} while (++ eit != *arrangement.unbounded_face()->inner_ccbs_begin());

	// TODO: could a polygon that actually has 3 collinear boundary points be a problem?
	// if 3 points are collinear, it means that the middle point is extra from the joined arrangement, so we can delete it
	// so that the polygon comparison works
	auto i = 0;
	while (i < vertices.size()) {
		// check whether the last 2 points in the vertex vector are collinear with the first one
		if (i == vertices.size() - 2) {
			if (CGAL::collinear(vertices.at(i), vertices.at(i + 1), vertices.at(0))) {
				vertices.erase(vertices.begin() + i + 1);
			} else
				i ++;
		// check whether the last point in the vertex vector is collinear with the first 2 points
		} else if (i == vertices.size() - 1) {
			if (CGAL::collinear(vertices.at(i), vertices.at(0), vertices.at(1))) {
				vertices.erase(vertices.begin());
			} else
				i ++;
		// check if any of the other points are collinear
		} else {
			if (CGAL::collinear(vertices.at(i), vertices.at(i + 1), vertices.at(i + 2))) {
				vertices.erase(vertices.begin() + i + 1);
			} else
				i ++;
		}
	}

	return Polygon_2(vertices.begin(), vertices.end()); // clockwise order!
}


/* polygon to arrangement function, as adapted from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArrangementFunctions.cpp
* :param   Polygon_2     polygon:       input polygon to be converted to an arrangement
* :return    Arrangement_2 arrangement: output arrangement converted from a polygon
*/
inline Arrangement_2 polygon_to_arrangement(Polygon_2 polygon) {
	std::vector<Segment_2> edges;
    Arrangement_2 arrangement;

	//Add points to arrangement
	for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++ it)
		edges.push_back(Segment_2(*it, *((it + 1 == polygon.vertices_end()) ? polygon.vertices_begin() : it + 1)));

	CGAL::insert_non_intersecting_curves(arrangement, edges.begin(), edges.end());

    return arrangement;
}

/* get_reflex_vertices method, as adapted from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArtGallery.cpp
* :in param Arrangement_2 arrangement:	the arrangement whose reflex vertices need to be computed
*
*  This method adds all the reflex vertices of am arrangement in the reflex_vertices vector.
*/
inline std::vector<Point_2> get_reflex_vertices(Arrangement_2 arrangement) {
	std::vector<Point_2> reflex_vertices;

	// Identify reflex vertices
	// use do/while for circular loop
	// The polygon is a "hole" in the unbounded face of the arrangement, thus a clockwise inner_ccb
	auto eit = *arrangement.unbounded_face()->inner_ccbs_begin();

	do {
		//Left turn, because the boundary is clockwise...
		if (CGAL::orientation(eit->prev()->source()->point(), eit->prev()->target()->point(), eit->target()->point()) == CGAL::LEFT_TURN)
			reflex_vertices.push_back(eit->source()->point());

	} while (++ eit != *arrangement.unbounded_face()->inner_ccbs_begin());

	return reflex_vertices;
}

/* is_inside_arrangement method
* :in param Arrangement_2 arrangement:	the arrangement which needs to be checked whether it contains the point inside
* :in param Point_2 p:  				the point that needs to be checked whether it is inside the arrangement
* :return bool:         				true if r is inside the arrangement
*                           				false otherwise
* 
* This method checks whether a point p is inside the input arrangement
*/
inline bool is_inside_arrangement(Arrangement_2 arrangement, Point_2 p) {
	// find where in the visibility region the guard is placed
	CGAL::Arr_naive_point_location<Arrangement_2> pl(arrangement);
	auto obj = pl.locate(p);

	// The query point may be a face on the visibility region boundary
	auto *face = boost::get<Arrangement_2::Face_const_handle>(&obj);
	
	// if the query point is within a face, and the face is bounded, then it is inside the arrangement
	if (face) {
		if (!(*face)->is_unbounded())
			return true;
	}
	// if the query point is not within a face, but it is on the arrangement boundary, we consider it to still be inside the arrangement
	else {
		auto *edge = boost::get<Arrangement_2::Halfedge_const_handle>(&obj);
		
		if (edge)
			return true;
	}

	return false;
}

/* compute_area function
* :in param Arrangement_2 arrangement:	the arrangement whose area needs to be computed
* :return double:						the area of the arrangement
*
* This function computes the area of an arrangement without holes. It first converts it to a simple polygon, then uses the built-in CGAL function to compute its area and return it.
*/
inline double compute_area(Arrangement_2 &arrangement) {
	auto visibility_polygon = arrangement_to_polygon(arrangement);
	return -CGAL::to_double(visibility_polygon.area());
}
