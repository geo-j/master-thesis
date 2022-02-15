#include <vector>


// TODO: think about how the kernel would need to be changed
typedef CGAL::Exact_predicates_exact_constructions_kernel                   Kernel;
typedef CGAL::Polygon_2<Kernel>                                             Polygon_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Segment_2                                                   Segment_2;


template<typename type>
void push_back_unique(std::vector<type> &v, type element) {
    auto it = std::find(v.begin(), v.end(), element);

    if (it == v.end())
        v.push_back(element);
}

/* arrangement to polygon function, as taken from Simon's implementation https://github.com/simonheng/AGPIterative/blob/main/ArtGalleryCore/ArrangementFunctions.cpp
* :param    Arrangement_2 arrangement: input arrangement to be converted to a polygon
* :return   Polygon_2     polygon:     output polygon converted from the given arrangement
*/
Polygon_2 arrangement_to_polygon(Arrangement_2 arrangement) {
	std::vector<Point_2> vertices;

	auto eit = arrangement.unbounded_face()->inner_ccbs_begin();
	do {
		vertices.push_back(eit->source()->point());
	} while (++ eit != arrangement.unbounded_face()->inner_ccbs_begin());

	return Polygon_2(vertices.begin(), vertices.end()); //clockwise order!
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