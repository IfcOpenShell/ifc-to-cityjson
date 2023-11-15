#ifndef SETTINGS_H
#define SETTINGS_H

// @tfk these are excluded on purpose, because when referenced from ifcconvert relative
// paths are needed. when referenced from the converter installed paths are assumed.

// #include <ifcgeom/IteratorSettings.h>
// #include <ifcparse/IfcFile.h>

#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include <string>
#include <vector>

// Global settings that are derived from the command line invocation
// parsed by Boost.ProgramOptions
struct geobim_settings {
	std::vector<std::string> input_filenames;

	std::string output_filename;	
	std::string cityjson_output_filename;
	std::string obj_output_filename;
	std::string json_output_filename;
	
	std::vector<std::string> radii;
	bool apply_openings, apply_openings_posthoc, debug, exact_segmentation, minkowski_triangles, no_erosion, spherical_padding;
	ifcopenshell::geometry::Settings settings;
	boost::optional<std::set<std::string>> entity_names;
	bool entity_names_included;
	std::vector<IfcParse::IfcFile*> file;
	boost::optional<size_t> threads;
};

// Parse the command line settings and do basic initialization
int parse_command_line(geobim_settings& settings, int argc, char** argv);

#endif
