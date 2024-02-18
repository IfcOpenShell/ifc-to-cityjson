#ifndef WRITER_H
#define WRITER_H

#include "processing.h"

#include <ifcgeom/kernels/cgal/CgalKernel.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/make_shared.hpp>

#include <nlohmann/json.hpp>

#include <array>
#include <fstream>

template <typename P>
struct vertex_cache {
	std::map<P, size_t> vertex_indices;
	std::vector<P> vertex_points;
};

struct print_and_clear_point_cache {
	std::ostream& obj_;

	print_and_clear_point_cache(std::ostream& obj)
	: obj_(obj) 
	{}

	void operator()(boost::blank&) {}
	
	template <typename T>
	void operator()(T& c) {
		for (auto& p : c.vertex_points) {
			obj_ << "v "
				<< p.cartesian(0) << " "
				<< p.cartesian(1) << " "
				<< p.cartesian(2) << "\n";
		}
		c.vertex_points.clear();
	}
};

// Abstract writer class that takes triangular facets.
struct abstract_writer {
	typedef CGAL::Simple_cartesian<double>::Point_3 P3;
	std::vector<P3>* point_lookup;

	boost::variant<boost::blank, vertex_cache<cgal_point_t>, vertex_cache<P3>> cache;

	std::array<P3, 3> points_from_facet(std::vector<std::vector<size_t>>::const_iterator f) {
		return {
			(*point_lookup)[(*f)[0]],
			(*point_lookup)[(*f)[1]],
			(*point_lookup)[(*f)[2]]
		};
	}

	std::array<P3, 3> points_from_facet(std::list<std::vector<std::vector<size_t>>::const_iterator>::const_iterator f) {
		return points_from_facet(*f);
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(cgal_shape_t::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<P3, 3> points_from_facet(CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>>::Facet_handle f) {
		return {
				f->facet_begin()->vertex()->point(),
				f->facet_begin()->next()->vertex()->point(),
				f->facet_begin()->next()->next()->vertex()->point()
		};
	}

	std::array<size_t, 3> point_indices_from_facet(std::vector<std::vector<size_t>>::const_iterator f) {
		return {
			(*f)[0],
			(*f)[1],
			(*f)[2]
		};
	}

	std::array<size_t, 3> point_indices_from_facet(std::list<std::vector<std::vector<size_t>>::const_iterator>::const_iterator f) {
		return point_indices_from_facet(*f);
	}

	template <typename T>
	std::array<size_t, 3> point_indices_from_facet(T t) {
		auto arr = points_from_facet(t);
		typedef typename decltype(arr)::value_type U;
		if (cache.which() == 0) {
			cache = vertex_cache<U>();
		}
		auto& C = boost::get<vertex_cache<U>>(cache);
		std::array<size_t, 3> idxs;
		for (int i = 0; i < 3; ++i) {
			auto it = C.vertex_indices.find(arr[i]);
			if (it == C.vertex_indices.end()) {
				auto n = C.vertex_indices.size();
				idxs[i] = C.vertex_indices[arr[i]] = n;
				C.vertex_points.push_back(arr[i]);
			}
			else {
				idxs[i] = it->second;
			}
		}
		return idxs;
	}

	std::array<Kernel_::Point_3, 3> points_from_facet(std::list<cgal_shape_t::Facet_handle>::iterator f) {
		return points_from_facet(*f);
	}

	bool has_finalized = false;
	virtual void do_finalize() = 0;
	virtual ~abstract_writer() {
	}
	void finalize() {
		if (!has_finalized) {
			has_finalized = true;
			do_finalize();
		}
	}
};

// OBJ writer for CGAL facets paired with a style
struct simple_obj_writer : public abstract_writer {
	int group_id = 1;
	int vertex_count = 1;
	std::ofstream obj, mtl;
	rgb GRAY;
	std::string fn_prefix_;

	simple_obj_writer(const std::string& fn_prefix)
		: fn_prefix_(fn_prefix)
		, GRAY(0.6, 0.6, 0.6)
	{}

	void begin() {
		obj.open((fn_prefix_).c_str());
		mtl.open((fn_prefix_.substr(0, fn_prefix_.size() - 4) + ".mtl").c_str());
		obj << "mtllib " << fn_prefix_ << ".mtl\n";
	}

	// @todo make the Kernel or Point_3 type discoverable from this template
	template <typename It>
	void operator()(const item_info* info, It begin, It end) {
		auto diffuse = info && info->style && info->style->diffuse ? rgb(info->style->diffuse.ccomponents()) : GRAY;

		obj << "g " << (info ? info->guid : "unknown") << "\n";
		obj << "usemtl m" << group_id << "\n";
		mtl << "newmtl m" << group_id << "\n";
		mtl << "kd " << diffuse.r() << " " << diffuse.g() << " " << diffuse.b() << "\n";

		group_id++;

		/*
		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			for (int i = 0; i < 3; ++i) {
				obj << "v "
					<< points[i].cartesian(0) << " "
					<< points[i].cartesian(1) << " "
					<< points[i].cartesian(2) << "\n";
			}
			obj << "f "
				<< (vertex_count + 0) << " "
				<< (vertex_count + 1) << " "
				<< (vertex_count + 2) << "\n";
			vertex_count += 3;
		}
		*/

		std::vector<std::array<size_t, 3>> fs;
		for (auto it = begin; it != end; ++it) {
			fs.push_back(point_indices_from_facet(it));
		}

		boost::apply_visitor(print_and_clear_point_cache(obj), cache);		

		for (auto& f : fs) {
			obj << "f "
				<< f[0] + 1 << " "
				<< f[1] + 1 << " "
				<< f[2] + 1 << "\n";
		}
	}

	void do_finalize() {}
};

namespace {

	struct predicate_always {
		bool operator()(const Eigen::Vector3d&) const {
			return true;
		}
	};

	struct predicate_is_up {
		bool operator()(const Eigen::Vector3d& norm) const {
			return norm(2) > 0.;
		}
	};

	std::string map_semantics(const std::string& ifc, const Eigen::Vector3d& norm) {
		static predicate_always always;
		static predicate_is_up is_up;

		static std::vector<std::pair<std::pair<std::string, std::function<bool(const Eigen::Vector3d&)>>, std::string>> mappings {
			{{"IfcSlab", is_up}, "RoofSurface"},
			{{"IfcSlab", always}, "GroundSurface"},
			{{"IfcWall", always}, "WallSurface"},
			{{"IfcWindow", always}, "Window"},
			{{"IfcDoor", always}, "Door"},
		};

		for (auto& m : mappings) {
			if (m.first.first == ifc && m.first.second(norm)) {
				return m.second;
			}
		}

		return "ClosureSurface";
	};
}

struct city_json_writer : public abstract_writer {
	rgb GRAY;

	using json = nlohmann::json;

	std::string filename;

	std::vector<std::array<double, 3>> vertices;
	std::vector<std::vector<std::vector<std::vector<int>>>> boundaries;
	std::vector<std::vector<int>> boundary_materials;
	std::vector<json> boundary_semantics;
	std::vector<int> boundary_semantics_values;

	json materials;

	city_json_writer(const std::string& fn_prefix)
		: filename(fn_prefix)
		, materials(json::array())
		, GRAY(0.6, 0.6, 0.6)
	{		
		boundaries.emplace_back();
		boundary_materials.emplace_back();
	}

	template <typename It>
	void operator()(const item_info* info, It begin, It end) {
		auto diffuse = info && info->style && info->style->diffuse ? rgb(info->style->diffuse.ccomponents()) : GRAY;

		json material = json::object();
		material["name"] = "material-" + boost::lexical_cast<std::string>(materials.size());
		material["diffuseColor"] = std::array<double, 3>{diffuse.r(), diffuse.g(), diffuse.b()};
		material["specularColor"] = std::array<double, 3>{0., 0., 0.};
		material["shininess"] = 0.;
		material["isSmooth"] = false;
		materials.push_back(material);

		for (auto it = begin; it != end; ++it) {
			auto points = points_from_facet(it);
			std::vector<int> faces;
			for (int i = 0; i < 3; ++i) {
				faces.push_back(vertices.size());
				vertices.push_back({ {
					CGAL::to_double(points[i].cartesian(0)),
					CGAL::to_double(points[i].cartesian(1)),
					CGAL::to_double(points[i].cartesian(2))
				} });
			}

			Eigen::Vector3d a, b, c;
			a << vertices[faces[0]][0], vertices[faces[0]][1], vertices[faces[0]][2];
			b << vertices[faces[1]][0], vertices[faces[1]][1], vertices[faces[1]][2];
			c << vertices[faces[2]][0], vertices[faces[2]][1], vertices[faces[2]][2];
			Eigen::Vector3d norm = (b - a).cross(c - a);

			boundaries.back().push_back({ faces });
			boundary_materials.back().push_back(materials.size() - 1);
			json json_type = json::object();
			json_type["type"] = map_semantics(info ? info->entity_type : "x", norm);
			boundary_semantics.push_back(json_type);  
			boundary_semantics_values.push_back(boundary_semantics_values.size());
		}
	}

	void do_finalize() {
		json city;

		city["type"] = "CityJSON";
		city["version"] = "1.0";
		city["extensions"] = json::object();
		city["metadata"]["referenceSystem"] = "urn:ogc:def:crs:EPSG::2355";
		city["vertices"] = vertices;
		city["appearance"]["materials"] = materials;

		auto& building1 = city["CityObjects"]["id-1"];
		building1["type"] = "Building";
		building1["geographicalExtent"] = std::array<double, 6>{0, 0, 0, 1, 1, 1};

		json geom = json::object();
		geom["type"] = "Solid";
		geom["lod"] = 2;
		geom["boundaries"] = boundaries;
		geom["semantics"]["values"][0] = boundary_semantics_values;
		geom["semantics"]["surfaces"] = boundary_semantics;
		geom["material"]["diffuse"]["values"] = boundary_materials;
		building1["geometry"].push_back(geom);

		std::ofstream(filename.c_str()) << city;
	}
};

struct external_element_collector : public abstract_writer {
	using json = nlohmann::json;

	std::string filename;
	const std::list<item_info*>& all_infos;
	std::set<const item_info*> part_of_exterior;


	external_element_collector(const std::string& fn_prefix, const std::list<item_info*>& infos)
		: filename(fn_prefix)
		, all_infos(infos)
	{}

	template <typename It>
	void operator()(const item_info* info, It begin, It end) {
		if (info) {
			part_of_exterior.insert(info);
		}
	}

	void do_finalize() {
		json data = json::array();

		for (auto& info : all_infos) {
			json object = json::object();
			object["guid"] = info->guid;
			object["is_external"] = part_of_exterior.find(info) != part_of_exterior.end();

			data.push_back(object);
		}

		std::ofstream(filename.c_str()) << data;
	}
};

struct polyhedron_collector : public abstract_writer {
	std::list<IfcGeom::Element*> elems;

	template <typename It>
	void operator()(const item_info* info, It begin, It end) {
		std::vector<std::array<size_t, 3>> fs;
		for (auto it = begin; it != end; ++it) {
			fs.push_back(point_indices_from_facet(it));
		}
		if (cache.which() != 1) {
			return;
			// @nb make sure exact segmentation is on.
			// @todo serialize CgalShapeSimple() otherwise
		}

		cgal_shape_t P;
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(boost::get<vertex_cache<cgal_point_t>>(cache).vertex_points, fs, P);
		P.normalize_border();

		IfcGeom::ConversionResults shapes;
		shapes.push_back(IfcGeom::ConversionResult(0, new ifcopenshell::geometry::CgalShape(P), info ? info->style : ifcopenshell::geometry::taxonomy::style::ptr(nullptr)));

		ifcopenshell::geometry::Settings settings;
		settings.get<ifcopenshell::geometry::settings::WeldVertices>().value = false;

		IfcGeom::BRepElement brep(
			info ? info->id : 0,
			info ? info->parent_id : 0,
			info ? info->name : std::string(""),
			info ? info->entity_type : std::string(""),
			info ? info->guid : std::string(""),
			"exterior", // context
			// @todo should we have an option to use local coordinates? (i.e multiple with placement inverse here?)
			ifcopenshell::geometry::taxonomy::make<ifcopenshell::geometry::taxonomy::matrix4>(),
			boost::make_shared<IfcGeom::Representation::BRep>(settings, info ? info->entity_type : std::string(""), std::to_string(info ? info->id : 0) + "-" + std::to_string(elems.size()) + "-exterior", shapes), // boost::shared_ptr<IfcGeom::Representation::BRep>& geometry
			// @todo can this remain nullptr safely?
			nullptr // IfcUtil::IfcBaseEntity* product
		);

		// @todo based on settings
		auto tri = new IfcGeom::TriangulationElement(brep);
		elems.push_back(tri);
	}

	void do_finalize() {}
};

#endif
