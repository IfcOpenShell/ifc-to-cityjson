#ifndef GLOBAL_EXECUTION_CONTEXT_H
#define GLOBAL_EXECUTION_CONTEXT_H

#include "processing.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// State that is relevant accross the different radii for which the
// program is executed. Used to map back the semantics from the original 
// IFC input to the newly constructed polyhedra.
template <typename TreeKernel>
struct global_execution_context : public execution_context {
	typedef CGAL::Polyhedron_3<TreeKernel> TreeShapeType;
	typedef CGAL::AABB_face_graph_triangle_primitive<TreeShapeType, CGAL::Default, CGAL::Tag_false> Primitive;
	typedef CGAL::AABB_traits<TreeKernel, Primitive> AAbbTraits;
	typedef typename AAbbTraits::Bounding_box Bounding_box;
	typedef CGAL::AABB_tree<AAbbTraits> AAbbTree;
	typedef typename AAbbTree::Primitive_id Primitive_id;
	typedef std::vector<std::pair<
		ifcopenshell::geometry::taxonomy::style*,
		std::list<cgal_shape_t::Facet_handle>>> segmentation_return_type;

	AAbbTree tree;

	std::list<ifcopenshell::geometry::taxonomy::style> styles;
	std::map<std::pair<double, std::pair<double, double>>, typename decltype(styles)::iterator> diffuse_to_style;

	// A reference is kept to the original shapes in a std::list.
	// Later an aabb tree is used map eroded triangle centroids
	// back to the original elements to preserve semantics.
	std::list<TreeShapeType> triangulated_shape_memory;
	std::map<typename TreeShapeType::Facet_handle, typename decltype(styles)::const_iterator> facet_to_style;

#ifdef GEOBIM_DEBUG
	simple_obj_writer obj_;
#endif

	global_execution_context();

	void operator()(shape_callback_item& item);

	void finalize();

	segmentation_return_type segment(const cgal_shape_t& input);
};

#endif