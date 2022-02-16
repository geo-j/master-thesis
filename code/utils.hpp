#include <vector>

#include <CGAL/Polygon_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>


// TODO: think about how the kernel would need to be changed
typedef CGAL::Exact_predicates_exact_constructions_kernel                   Kernel;
typedef CGAL::Polygon_2<Kernel>                                             Polygon_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Segment_2                                                   Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                                  Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                       Arrangement_2;


template<typename type>
void push_back_unique(std::vector<type> &v, type element) {
    auto it = std::find(v.begin(), v.end(), element);

    if (it == v.end())
        v.push_back(element);
}

/* arrangement to polygon function, as adapted from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArrangementFunctions.cpp
* :param    Arrangement_2 arrangement: input arrangement to be converted to a polygon
* :return   Polygon_2     polygon:     output polygon converted from the given arrangement
*/
Polygon_2 arrangement_to_polygon(Arrangement_2 &arrangement) {
	std::vector<Point_2> vertices;

	// loop around the arrangement
	auto eit = *arrangement.unbounded_face()->inner_ccbs_begin();
	do {
		vertices.push_back(eit->source()->point());
	} while (++ eit != *arrangement.unbounded_face()->inner_ccbs_begin());

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
Arrangement_2 polygon_to_arrangement(Polygon_2 polygon) {
	std::vector<Segment_2> edges;
    Arrangement_2 arrangement;

	//Add points to arrangement
	for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++ it)
		edges.push_back(Segment_2(*it, *((it + 1 == polygon.vertices_end()) ? polygon.vertices_begin() : it + 1)));

	CGAL::insert_non_intersecting_curves(arrangement, edges.begin(), edges.end());

    return arrangement;
}