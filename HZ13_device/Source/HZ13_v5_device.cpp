/*
  Hornof-Zhang 2013 variation on Williams 67 device created by Dave Kieras, 2004.
  Modified by Kieras, 2013, 2014, 2015
*/


#include "HZ13_v5_device.h"
#include "EPICLib/Geometry.h"
#include "EPICLib/Output_tee_globals.h"
#include "EPICLib/Numeric_utilities.h"
#include "EPICLib/Random_utilities.h"
#include "EPICLib/Symbol_utilities.h"
#include "EPICLib/Assert_throw.h"
#include "EPICLib/Device_exception.h"
#include "EPICLib/Standard_Symbols.h"
#include "EPICLib/Statistics.h"

#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
//#include <cassert>
#include <cmath>
#include <fstream>

namespace GU = Geometry_Utilities;
//using GU::Point;
using namespace std;

/*
Visual field is 39 X 30 degrees, with 75 objects distributed at random,
using 96 different search displays in file search_fields.txt.

The actual trial sequences are in the file all_trials.txt, and represent
22 subjects who viewed each each search display under a particular cue
condition, for a total of 2112 trials, each condition presented 264 times.

Version 3 Trial specs are always used. If All is specified, then all conditions 
in the trial spec are run for a total of 2112 trials. 

Number, or some combination of {Color, Size, Shape} can be specified and the
corresponding single condition is run. This single condition consists of the matching
condition in the trial specs, so that an individual condition is a subset of the 
complete experiment stimuli in the same condition.

The number of replications is the number of passes made over the trial specs,
so three replications of All would be 2112 * 3 = 6336 trials.

Stimuli are in three sizes: 
2.8, 1.6, 0.8 degrees, called large, medium, and small as "encoded size"
actual visual size differs depending on shape.

5 colors:
blue, green, yellow, red, purple

5 shapes:
circle, semicircle, triangle, square, cross

each stimulus has a two digit number 1-75, that was 0.26 DVA in height (.52 wide X .26 high we'll say)

probe is mouse target XX, then color, size, shape, number



trial start
probe cue presented
after study, user points, clicks on XX
cue disappears,
search screen appears
point to and click on target.

We have the actual stimuli and trial sequences.
There are 96 search fields, numbered 24-119, each with a target.
There are 22 subjects, each saw the search fields once, with a 
random assignment of condition to search field, and in random order.

One repetition through this set of trials will yield 2112 trials, each
search field presented 22 times, each condition presented 264 times.

Since a repetition with trials data contains all of the conditions,
statistics on results has to be accumulated for each condition.

Eight conditions: Color/~Color X Size/~Size X Shape/~Shape
Show proportion matches to target (regardless of whether cued)
Color, Size, Shape
Show proportion matches as a function of object size (small, medium, large).
Show saccade distance for each property.
Color [small, medium, large] proportion, distance
Size ditto
Shape ditto

Data accumulation: fixation proportions and means computed in each trial (denominator = number of fixations)
with some variables, like RT being a single data point for that condition,
then all variables averaged over trials for each condition (denominator = number of trials in that condition).

Use Q99 criterion from data reduction to set excess fixations criterion in each condition; 
do not include trials with excess fixations in the means.
*/




// handy internal constants
const Symbol Fixation_point_c("Fixation_point");
const Symbol Display_c("Display");
const Symbol probe_mouse_target_c("probe_mouse_target");
const Symbol probe_mouse_target_text_c("XX");
const Symbol probe_label_c("probe_label");
const Symbol probe_shape_c("probe_shape");
const Symbol probe_color_c("probe_color");
const Symbol probe_size_c("probe_size");
const Symbol unspecified_cue_c("---");

const Symbol Small_c("Small");
const Symbol Medium_c("Medium");
const Symbol Large_c("Large");

// from the HZ13 display properties
// const double degrees_per_pixel_c = 0.02586;
//const double display_center_x_c = 800.;
//const double display_center_y_c = 600.;

constexpr double display_hor_pixel_size_c = 1600.;
constexpr double display_ver_pixel_size_c = 1200;
constexpr double display_hor_pixel_center_c = display_hor_pixel_size_c / 2.;
constexpr double display_ver_pixel_center_c = display_ver_pixel_size_c / 2.;
constexpr double display_hor_mm_size_c = 432;
constexpr double display_ver_mm_size_c = 324.;
//constexpr double display_hor_mm_per_pixel_c = display_hor_mm_size_c /display_hor_pixel_size_c;
//constexpr double display_ver_mm_per_pixel_c = display_ver_mm_size_c /display_ver_pixel_size_c;
constexpr double display_viewing_distance_mm_c = 600.;
constexpr double display_viewing_distance_pixels_h_c = display_viewing_distance_mm_c * display_hor_pixel_size_c / display_hor_mm_size_c;
constexpr double display_viewing_distance_pixels_v_c = display_viewing_distance_mm_c * display_ver_pixel_size_c / display_ver_mm_size_c;
const GU::Point display_center_pixels_c(display_hor_pixel_center_c, display_ver_pixel_center_c);
const GU::Point initial_pixel_eye_position_c(800., 675);  // assumed initial eye position in pixels - see Zhang email 9/15/15

GU::Point pixel_Point_to_DVA_Point(GU::Point pixel_point)
{
		// using provided display parameters
		return GU::Point(
            GU::degrees_subtended(pixel_point.x - display_hor_pixel_center_c, display_viewing_distance_pixels_h_c),
            GU::degrees_subtended(pixel_point.y - display_ver_pixel_center_c, display_viewing_distance_pixels_v_c)
            );
}
const GU::Point initial_eye_position_c(pixel_Point_to_DVA_Point(initial_pixel_eye_position_c));  // assumed initial eye position in pixels


/* x, y gaussian noise parameters to be applied to eye location sent to device
values from Yunfeng 12/14/15
Mean of the residual error (DoV): x=-0.018879890470142 y=-0.0650516458420223
Standard deviation of the residual error (DoV): x=0.296691737023752 y=0.350642391504673
*/
const double x_noise_mean_c = -0.018879890470142;
const double x_noise_sd_c = 0.296691737023752;
//const double x_noise_sd_c = 1.0; // test
//const double x_noise_sd_c = 0.75; // test
//const double x_noise_sd_c = 0.1; // test
const double y_noise_mean_c = -0.0650516458420223;
const double y_noise_sd_c = 0.350642391504673;
//const double y_noise_sd_c = 1.0; // test
//const double y_noise_sd_c = 0.75; // test
//const double y_noise_sd_c = 0.1; // test

GU::Point add_measurement_error_noise(GU::Point measured_location)
{
    double noise_x = normal_random_variable(x_noise_mean_c, x_noise_sd_c);
    double noise_y = normal_random_variable(y_noise_mean_c, y_noise_sd_c);
    GU::Point noisy_fixation_location(measured_location.x + noise_x, measured_location.y + noise_y);
    return noisy_fixation_location;
}

#define CHECK_DENOMINATOR_OUTPUT
#undef CHECK_DENOMINATOR_ASSERT
void VisSearch_device::Trial_accumulator::check_denominators(Output_tee& device_out) const
{
    int nfix = n_trial_fixations;
    int n_inbetween = inbetween_fixations.get_count();
    int n_AOI = AOI_fixations.get_count();
	
	
#ifdef CHECK_DENOMINATOR_OUTPUT
    
    if(
    (nfix != (n_inbetween + n_AOI)) ||
    (abs(n_inbetween - inbetween_dwell.get_n()) > 1) ||
    (nfix != inbetween_fixations.get_n() || nfix != AOI_fixations.get_n()) ||
    (nfix != color_matches.get_n() || nfix != size_matches.get_n() || nfix != shape_matches.get_n() ) ||
    (nfix != revisits.get_n())
    ) {
    device_out << "Trial counts: nfix = " << nfix << ", n_inbetween = " << n_inbetween << ", n_AOI = " << n_AOI << endl;
    device_out << "Trial denominators: nfix = " << nfix << ", inbetween = " << inbetween_fixations.get_n() << ", AOI = " << AOI_fixations.get_n() << endl;
    device_out << "inbetween_dwell denominator =  " << inbetween_dwell.get_n() << endl;
	// something odd if unequal, but not by 1
	if((n_inbetween != inbetween_dwell.get_n()) && (n_inbetween - 1 != inbetween_dwell.get_n()))
		device_out << "n_inbetween " << n_inbetween << " inbetween_dwell denominator = " << inbetween_dwell.get_n() << endl;
    if(nfix != AOI_fixations.get_n())
		device_out << "Trial denominators: nfix = " << nfix << ", inbetween = " << inbetween_fixations.get_n() << ", AOI = " << AOI_fixations.get_n() << endl;
//    device_out << "dwell denominator =  " << dwell.get_n() << endl;

    if(n_AOI != dwell.get_n() && (n_AOI - 1 != dwell.get_n()) && (n_AOI - 2 != dwell.get_n()))
		device_out << "n_AOI " << n_AOI << " dwell denominator = " << dwell.get_n() << endl;
    device_out << "color_matches count =  " << color_matches.get_count() << ", denominator =  " << color_matches.get_n() << endl;
    device_out << "color_match_distance denominator =  " << color_match_distance.get_n() << endl;
    device_out << "color_mismatch_distance denominator =  " << color_mismatch_distance.get_n() << endl;
    device_out << "size_matches count =  " << size_matches.get_count() << ", denominator =  " << size_matches.get_n() << endl;
    device_out << "shape_matches count =  " << shape_matches.get_count() << ", denominator =  " << shape_matches.get_n() << endl;
    device_out << "size accumulators:" << endl;
    for(auto& accumulator : size_accumulators) {
        device_out << "size all fixations count = " << accumulator.all_fixations.get_count()
            << ", denominator =  " << accumulator.all_fixations.get_n() << endl;
        device_out << "size saccade_distance denominator = " << accumulator.saccade_distance.get_n() << endl;
        device_out << "size color_matches count = " << accumulator.color_matches.get_count()
            << ", denominator =  " << accumulator.color_matches.get_n() << endl;
        }
    device_out << "revisits count =  " << revisits.get_count() << ", denominator =  " << revisits.get_n() << endl;
    device_out << "unique_objects_fixated =  " << unique_objects_fixated << ", mean_fixation_frequency  =  " << mean_fixation_frequency << endl;
    }
#endif

#ifdef CHECK_DENOMINATOR_ASSERT
    Assert(nfix == inbetween_fixations.get_n());

    Assert(n_inbetween == inbetween_dwell.get_n()|| (n_inbetween - 1 == inbetween_dwell.get_n()) || (n_inbetween - 2 == inbetween_dwell.get_n())); // one of these could be last
    Assert(n_inbetween == inbetween_saccade_distance.get_n());
    Assert(n_inbetween == inbetween_closest_distance.get_n());
    Assert(nfix == AOI_fixations.get_n());
    
    // last fixation on click target is not included in dwell calculation
    Assert(n_AOI == dwell.get_n() || (n_AOI - 1 == dwell.get_n()) || (n_AOI - 2 == dwell.get_n()));  // last fixation on click target is not included in dwell calculation
    Assert(n_AOI == saccade_distance.get_n());
    Assert(n_AOI == fixation_closest_distance.get_n());
    Assert(nfix == revisits.get_n());
    Assert(nfix == color_matches.get_n());
    Assert(nfix == size_matches.get_n());
    Assert(nfix == shape_matches.get_n());
    for(auto& accumulator : size_accumulators) {
        accumulator.check_denominators(nfix);
        }
    Assert(nfix == immediate_repeats.get_n()); // a proportion of total
    Assert(n_AOI >= unique_objects_fixated);
    Assert(unique_objects_fixated == mean_fixation_frequency.get_n());
    Assert(n_AOI >= repeat_lags.get_n());  // denominator of repeat_lags depends on how many repeat lags there are

            Assert(color_matches.get_count() == color_match_distance.get_n());  // all of these denominators depend on number of matches/mismatches
            // 6/9/16 now counting color match/mismatch and color_matches only if color is cued, so no relation to number of AOI fixations
//            Assert(n_AOI - color_matches.get_count() == color_mismatch_distance.get_n());
            Assert(size_matches.get_count() == size_match_distance.get_n());
//            Assert(n_AOI - size_matches.get_count() == size_mismatch_distance.get_n());
            Assert(shape_matches.get_count() == shape_match_distance.get_n());
//            Assert(n_AOI - shape_matches.get_count() == shape_mismatch_distance.get_n());
 #endif
	
}


string make_condition_label_string(int cue_color, int cue_size, int cue_shape);

	
VisSearch_device::Search_object::Search_object(int id_, GU::Point location_, GU::Size size_, const Symbol& size_name_, int isize_, const Symbol& color_, int icolor_,
	const Symbol& shape_, int ishape_) :
	id(id_), location(location_), size(size_), size_name(size_name_), isize(isize_), color(color_), icolor(icolor_), shape(shape_), ishape(ishape_)
{	
//	name = concatenate_to_Symbol("Search", id);
	label = int_to_Symbol(id);
}

VisSearch_device::VisSearch_device(const std::string& device_name, Output_tee& ot) :
		Device_base(device_name, ot), 
		state(START), condition_string("1 300 AOI All"), do_all_conditions(true)
{
	// initialize the various display parameters
	parse_condition_string();
	
	// map string in search_fields file to array index
	color_index_map = {
		{"blue", 0},
		{"green", 1},
		{"yellow", 2},
		{"red", 3},
		{"purple", 4}
		};
	// map array index to Symbols used in simulation; search to get index from Symbol
    //          0       1       2           3       4
	colors = {Blue_c, Green_c, Yellow_c, Pink_c, Teal_c};

	shape_index_map = {
		{"circle", 0},
		{"square", 1},
		{"cross", 2},
		{"triangle", 3},
		{"semi-circle", 4}
		};
//              0          1        2          3        4
	shapes = {Circle_c, Square_c, Cross_c, Triangle_c, Top_Semicircle_c};

	size_index_map = {
		{"small", 0},
		{"medium", 1},
		{"large", 2}
		};
//                  0       1           2
	size_names = {Small_c, Medium_c, Large_c};
	size_values = {0.8, 1.6, 2.8}; // diameters
    AOI_radius_values = {1.0, 1.4, 2.0}; // object radius + 0.6
    
    	shape_size_factors = {
		{shape_index_map.at("circle"), GU::Size(1.0, 1.0)},
		{shape_index_map.at("square"), GU::Size(0.8862269, 0.8862269)},
		{shape_index_map.at("cross"), GU::Size(1.079906, 1.079906)},
		{shape_index_map.at("triangle"), GU::Size(1.346774, 1.16634)},
		{shape_index_map.at("semi-circle"), GU::Size(1.414214, 0.707107)}
		};

/* from Yunfeng 12/12/15 Q99 for nfixations in conditions
Number	295 NOnly
Shape	253
Color	73
Color+Shape	78
Size	190
Size+Shape	123
Color+Size	66
All	58 ClrSizShp

color 0, 1
	size 0, 1
		shape 0, 1
*/
	max_fixations_criterion[0][0][0] = 295; // NOnly
	max_fixations_criterion[0][0][1] = 253; // Shape
	max_fixations_criterion[0][1][0] = 190; // Size
	max_fixations_criterion[0][1][1] = 123; // Size+Shape
	max_fixations_criterion[1][0][0] = 73; // Color
	max_fixations_criterion[1][0][1] = 78; // Color+Shape
	max_fixations_criterion[1][1][0] = 66; // Color+Size
	max_fixations_criterion[1][1][1] = 58; // Color+Size+Shape
	
	cursor_location = GU::Point(0., 0.);
	
	initialize();

}

void VisSearch_device::parse_condition_string()
{
	// build an error message string in case we need it
	string error_msg(condition_string);
	error_msg += "\n Should be: number of repetitions (int > 0), max number of fixations (int > 0), either Cls or AOI for fixation criterion, and either All, Number, or any combination of {Size, Color, Shape} for cue condition";
	istringstream iss(condition_string);

	// Parse the number of repetitions and max n fixatons
	int nr, nf;
	iss >> nr >> nf;
	// do all error checks first
	if(!iss)
		throw Device_exception(this, string("Incorrect condition string: ") + error_msg);
	if(nr <= 0)
		throw Device_exception(this, string("Number of repetitions must be positive: ") + error_msg);
	if(nf <= 0)
		throw Device_exception(this, string("Max number of fixations must be positive: ") + error_msg);
	// now go with it
	n_repetitions = nr;
	max_n_fixations	= nf;
    
    // look for the fixation criterion
    string crit;
    iss >> crit;
    if(crit == "Cls")
        using_AOI = false;
    else if(crit == "AOI")
        using_AOI = true;
    else
        throw Device_exception(this, string("Missing fixation criterion: ") + error_msg);

	// look for the condition words
	current_trial_spec = Trial_spec();
	single_condition_spec = Trial_spec();
	do_all_conditions = false;
    bool number_only = false;
	string word;
	while(iss >> word) {
        if(word == "All") {
            do_all_conditions = true;
            }
       else if(word == "Number") {
            if(do_all_conditions)
                throw Device_exception(this, string("All conditions was specified: ") + error_msg);
            number_only = true;
            single_condition_spec.cue_size = 0;
            single_condition_spec.cue_color = 0;
            single_condition_spec.cue_shape = 0;
            }
        else if(word == "Size") {
           if(do_all_conditions || number_only)
                throw Device_exception(this, string("All or Number was specified: ") + error_msg);
            single_condition_spec.cue_size = 1;
            }
        else if(word == "Color") {
           if(do_all_conditions || number_only)
                throw Device_exception(this, string("All or Number was specified: ") + error_msg);
           single_condition_spec.cue_color = 1;
            }
        else if(word == "Shape") {
           if(do_all_conditions || number_only)
                throw Device_exception(this, string("All or Number was specified: ") + error_msg);
             single_condition_spec.cue_shape = 1;
            }
        else
			throw Device_exception(this, string("Unrecognized condition or cue specification: ") + error_msg);
		}
}

void VisSearch_device::set_parameter_string(const string& condition_string_)
{
	condition_string = condition_string_;
	parse_condition_string();
}

string VisSearch_device::get_parameter_string() const
{
	return condition_string;
}

void VisSearch_device::initialize()
{
	state = START;
	object_counter = 0;
    // check containers
    Assert(color_index_map.size() == n_colors);
    Assert(colors.size() == n_colors);
    Assert(size_index_map.size() == n_sizes);
    Assert(size_names.size() == n_sizes);
    Assert(size_values.size() == n_sizes);
    Assert(shape_index_map.size() == n_shapes);
    Assert(shapes.size() == n_shapes);

}

void VisSearch_device::initialize_data_collection()
{
	current_trial_spec_index = 0;
	current_repetition = 0;
    current_search_field = 0; // starting subscript
	
	// doing all conditions, start at first trial spec
	if(do_all_conditions) {
		current_trial_spec = trial_specs[current_trial_spec_index];
		}
	// if a single condition, advance to the next (first) trial spec that matches the condition
	else {
		bool first_single_condition_spec_found = find_next_single_condition_trial_spec();
		Assert(first_single_condition_spec_found);
		}
	
    n_trials = 0;
	n_probe_fixations = 0;
    n_terminated = 0;

	fixated_objects.clear();
	fixated_object_frequencies.clear();
	last_fixation_numbers.clear();

	last_fixation_location = GU::Point();
	last_fixation_target_name = Nil_c;
	last_fixated_Search_object = Search_object();
    trial_accumulators.reset();
    reset_condition_accumulators();
	
}

// this returns the string because we need it in more than one place
string make_condition_label_string(int cue_color, int cue_size, int cue_shape)
{
    string result;
    if(!cue_color && !cue_size && !cue_shape){
        result = "NOnly";
        return result;
        }
    if(cue_color)
        result += "Clr";
    if(cue_size)
        result += "Siz";
    if(cue_shape)
        result += "Shp";
    return result;
}


void VisSearch_device::reset_condition_accumulators()
{
    for(int cue_color = 0; cue_color < 2; cue_color++)
        for(int cue_size = 0; cue_size < 2; cue_size++)
            for(int cue_shape = 0; cue_shape < 2; cue_shape++)
                condition_accumulators[cue_color][cue_size][cue_shape].reset();
}

void VisSearch_device::handle_Start_event()
{
	if(device_out) {
		device_out << processor_info() << "Eye location noise parameters: X: N("
        << x_noise_mean_c << ", " << x_noise_sd_c << "), Y: N("
        << y_noise_mean_c << ", " << y_noise_sd_c << ")" << endl;
		device_out << "Max number of fixations per trial in conditions:\n";
		for(int cue_color = 0; cue_color < 2; cue_color++)
			for(int cue_size = 0; cue_size < 2; cue_size++)
				for(int cue_shape = 0; cue_shape < 2; cue_shape++) {
					device_out << make_condition_label_string(cue_color, cue_size, cue_shape) << ": "
						<< max_fixations_criterion[cue_color][cue_size][cue_shape] << endl;
					}
		}
//	if(device_out)
//		device_out << processor_info() << "received Start_event" << endl;
	if(out_file.is_open())  // if it is open from previous run, close it.
		out_file.close();
    out_file.open("HZ13_output_file.txt");
    if(!out_file.is_open())
        throw Device_exception(this, "Could not open output file!");
    out_file << get_human_prs_filename() << endl;
	load_trial_data_files();
    initialize();
	initialize_data_collection();

 	schedule_delay_event(500);
}

void VisSearch_device::handle_Stop_event()
{
//	if(device_out)
//		device_out << processor_info() << "received Stop_event" << endl;
//	output_statistics();
		
}

void VisSearch_device::load_trial_data_files()
{
	// open and load search fields file
	search_fields_file.open("search_fields.txt");
	if(!search_fields_file.good())
		throw Device_exception(this, "could not open search_fields.txt");

	search_fields.clear();
	search_fields.resize(n_search_fields);
	for(int i = 0; i < n_search_fields; i++)
		load_search_field();
//	device_out << n_search_fields << " search fields read" << endl;

	search_fields_file.close();
		

	trial_specs_file.open("all_trials.txt");
	if(!trial_specs_file.good())
		throw Device_exception(this, "could not open all_trials.txt");
	
	trial_specs.clear();
	for(int i = 0; i < n_trials_per_repetition; i++)
		load_trial_spec();
	device_out << trial_specs.size() << " trial specs read" << endl;
	for(int icolor = 0; icolor < 2; icolor++)
		for(int isize = 0; isize < 2; isize++)
			for(int ishape = 0; ishape < 2; ishape++)
				device_out << make_condition_label_string(icolor, isize, ishape) << ' ' << trials_per_condition[icolor][isize][ishape] << endl;
	
	Assert(trial_specs.size() == n_trials_per_repetition);

	trial_specs_file.close();
}

void VisSearch_device::load_search_field()
{
	int in_search_field_num, in_target_ID;
	
	search_fields_file >> in_search_field_num >> in_target_ID;
	if(!search_fields_file)
		throw Device_exception(this, "could not read search field number or target ID from search fields file");
	
	int search_field_index = in_search_field_num - n_practice_trials;
	Assert(search_field_index >= 0 && search_field_index < n_search_fields);
	search_fields[search_field_index].target_index = in_target_ID - 1;
	search_fields[search_field_index].objects.resize(n_objects);
	for(int i = 0; i < n_objects; i++) {
		int in_ID;
		string in_size, in_color, in_shape;
		double in_x, in_y;
		search_fields_file >> in_ID >> in_size >> in_color >> in_shape >> in_x >> in_y;
		if(!search_fields_file)
			throw Device_exception(this, "could not read object information from search fields file");
		int icolor = lookup_index(color_index_map, in_color);
        Assert(icolor >= 0 && icolor < colors.size());
		int isize = lookup_index(size_index_map, in_size);
        Assert(isize >= 0 && isize < size_values.size());
		int ishape = lookup_index(shape_index_map, in_shape);
        Assert(ishape >= 0 && ishape < shapes.size());
		// get the actual object size - for semicircles, object center might be off ..
		GU::Size object_simple_size{size_values[isize], size_values[isize]};
		GU::Size shape_size_factor = shape_size_factors[ishape];
		GU::Size object_actual_size{object_simple_size.h * shape_size_factor.h, object_simple_size.v * shape_size_factor.v};
		//device_out << in_size << ' ' << object_simple_size << ' ' << in_shape << ' ' << object_actual_size << endl;
		// using provided display parameters
        // use more accurate angle subtended information
        GU::Point object_center_pixels(in_x,in_y);
		GU::Point object_center_DVA = pixel_Point_to_DVA_Point(object_center_pixels);
//		object_center.x = (in_x - display_center_x_c) * degrees_per_pixel_c;
//		object_center.y = (in_y - display_center_y_c) * degrees_per_pixel_c;
		int ID = in_ID - 1;
		Assert(ID >= 0 && ID < n_objects);
		search_fields[search_field_index].objects[ID] =
//			Search_object(in_ID, object_center, object_simple_size,
			Search_object(in_ID, object_center_DVA, object_actual_size,
            size_names[isize], isize, colors[icolor], icolor, shapes[ishape], ishape);
		}
}

void VisSearch_device::load_trial_spec()
{
	int in_search_field_ID;
	trial_specs_file >> in_search_field_ID;
	if(!trial_specs_file)
		throw Device_exception(this, "could not read search field ID in trial spec");
	bool color_cue, size_cue, shape_cue;
	trial_specs_file >> color_cue >> size_cue >> shape_cue;
	if(!trial_specs_file)
		throw Device_exception(this, "could not read cue information in trial spec");
	int search_field_ID = in_search_field_ID - n_practice_trials;
	trial_specs.push_back(Trial_spec(search_field_ID, color_cue, size_cue, shape_cue));
	trials_per_condition[int(color_cue)][int(size_cue)][int(shape_cue)]++;
}

// look for the symbol in the index map, return the index or throw an exception
int VisSearch_device::lookup_index(const map<string, int>& themap, const string& key)
{
	auto it = themap.find(key);
	if(it == themap.end()) 
		throw Device_exception(this, string("could not find ") + key + string(" in map"));
	return it->second;
}



void VisSearch_device::handle_Delay_event(const Symbol& type, const Symbol& datum, 
		const Symbol& object_name, const Symbol& property_name, const Symbol& property_value)
{	
	switch(state) {
		case START:
			state = PRESENT_CURSOR;
			schedule_delay_event(300);
			break;
		case PRESENT_CURSOR:
			present_cursor();
			state = START_TRIAL;
			schedule_delay_event(200);
            break;
		case START_TRIAL:
            first_mouse_move_happened = false; // reset to wait for first mouse move to happen
			// reset excess fixations termination flag
			trial_terminated = false;
			state = PRESENT_PROBE;
			schedule_delay_event(100);
			break;
		case STOPPING_TRIAL:
			remove_stop_signal();
			do_trial_end();
//			state = START_TRIAL;
//			schedule_delay_event(1000);
			break;
		case PRESENT_PROBE:
//            device_out << "---------- starting trial ----------" << endl;
			prepare_trial();
			present_probe();
            state = WAITING_FOR_PRECUE_RESPONSE;
			break;
		case PRESENT_SEARCH_OBJECTS:
            cursor_in_motion = false;
			present_search_objects();
			state = WAITING_FOR_RESPONSE;
			break;
		case SHUTDOWN:
            out_file.close();
			stop_simulation();
			break;
		case WAITING_FOR_RESPONSE:
		default:
			throw Device_exception(this, "Device delay event in unknown or improper device state");
			break;
		}
}

void VisSearch_device::present_cursor()
{
	make_visual_object_appear(Cursor_name_c, cursor_location, GU::Size(1., 1.));
	set_visual_object_property(Cursor_name_c, Color_c, Black_c);
	set_visual_object_property(Cursor_name_c, Shape_c, Cursor_Arrow_c);
}

void VisSearch_device::prepare_trial()
{
	last_fixation_location = GU::Point();
	last_fixation_target_name = Nil_c;
	last_fixated_Search_object = Search_object();
	fixated_objects.clear();
	fixated_object_frequencies.clear();
	last_fixation_numbers.clear();

	search_objects.clear();
    trial_accumulators.reset();
    current_search_field = current_trial_spec.search_field_index;

//	if(!(current_search_field < search_fields.size()))
//        device_out << current_search_field << ' ' << search_fields.size() << ' ' << current_repetition << endl;

    Assert(current_search_field < search_fields.size());
	search_objects = search_fields[current_search_field].objects;
	int target_index = search_fields[current_search_field].target_index;
	
	// save a copy for convenience
	target = search_objects[target_index];
	
//	device_out << "No of locations:" << all_locations.size() << endl;
//	device_out << "No of search objects:" << search_objects.size() << endl;

}

/* size, color, shape, number, from top to bottom */
void VisSearch_device::present_probe()
{
    // set the current condition_label_string color size shape
    condition_label_string = make_condition_label_string(current_trial_spec.cue_color, current_trial_spec.cue_size, current_trial_spec.cue_shape);

	// give each object a unique name
	object_counter++;	// increment for each trial
    probe_mouse_target_name = concatenate_to_Symbol(probe_mouse_target_c, object_counter);
	probe_size_object_name = concatenate_to_Symbol(probe_size_c, object_counter);
	probe_color_object_name = concatenate_to_Symbol(probe_color_c, object_counter);
	probe_shape_object_name = concatenate_to_Symbol(probe_shape_c, object_counter);
	probe_label_object_name = concatenate_to_Symbol(probe_label_c, object_counter);

	// present the probe as a set of objects in the center of the display
	make_visual_object_appear(probe_mouse_target_name, GU::Point(0.0, 1.0), GU::Size(1.0, 0.3));
	make_visual_object_appear(probe_size_object_name, GU::Point(0.0, 0.6), GU::Size(1.0, 0.3));
	make_visual_object_appear(probe_color_object_name, GU::Point(0.0, 0.2), GU::Size(1.0, 0.3));
	make_visual_object_appear(probe_shape_object_name, GU::Point(0.0, -0.2), GU::Size(1.0, 0.3));
	make_visual_object_appear(probe_label_object_name, GU::Point(0.0, -0.6), GU::Size(0.6, 0.3));
	
	set_visual_object_property(probe_mouse_target_name, Text_c, probe_mouse_target_text_c);
	set_visual_object_property(probe_mouse_target_name, Position_c, Top_c);
	set_visual_object_property(probe_mouse_target_name, Above_c, probe_size_object_name);

	// if the cue flag for each cue type is true, provide the text for the cue; otherwise set to "---"
	if(current_trial_spec.cue_size)
		set_visual_object_property(probe_size_object_name, Text_c, target.size_name);
	else
		set_visual_object_property(probe_size_object_name, Text_c, unspecified_cue_c);
	set_visual_object_property(probe_size_object_name, Below_c, probe_mouse_target_name);
	set_visual_object_property(probe_size_object_name, Above_c, probe_color_object_name);

	if(current_trial_spec.cue_color)
		set_visual_object_property(probe_color_object_name, Text_c, target.color);
	else
		set_visual_object_property(probe_color_object_name, Text_c, unspecified_cue_c);
	set_visual_object_property(probe_color_object_name, Below_c, probe_size_object_name);
	set_visual_object_property(probe_color_object_name, Above_c, probe_shape_object_name);

	if(current_trial_spec.cue_shape)
		set_visual_object_property(probe_shape_object_name, Text_c, target.shape);
	else
		set_visual_object_property(probe_shape_object_name, Text_c, unspecified_cue_c);
	set_visual_object_property(probe_shape_object_name, Below_c, probe_color_object_name);
	set_visual_object_property(probe_shape_object_name, Above_c, probe_label_object_name);

	set_visual_object_property(probe_label_object_name, Text_c, target.label);
	set_visual_object_property(probe_label_object_name, Position_c, Bottom_c);
	set_visual_object_property(probe_label_object_name, Below_c, probe_shape_object_name);
}

void VisSearch_device::present_search_objects()
{
	for(int i = 0; i < n_objects; i++) {
		Search_object& obj = search_objects[i];
		ostringstream oss;
		// object_counter is incremented for each trial, ensuring a unique name across trials
		oss << "S" << obj.id << "_" << object_counter;
        string str_label = oss.str();
		obj.name = Symbol(str_label);
        str_label += "L";
        obj.label_name = Symbol(str_label);
        
        // present the object label as a separate object so it has its own size
		make_visual_object_appear(obj.name, obj.location, obj.size);
		set_visual_object_property(obj.name, Encoded_size_c, obj.size_name);
		set_visual_object_property(obj.name, Color_c, obj.color);
		set_visual_object_property(obj.name, Shape_c, obj.shape);
        make_visual_object_appear(obj.label_name, obj.location, GU::Size(0.52, 0.26));
		set_visual_object_property(obj.label_name, Inside_c, obj.name);
		set_visual_object_property(obj.label_name, Text_c, obj.label);
        
        
//		device_out << obj.name << ' ' << obj.location << ' ' << obj.size << ' ' << obj.size_name << ' '
//			<< obj.color << ' ' << obj.shape << ' ' << obj.label << endl;
		}
		
	stimulus_onset_time = get_time();
    first_mouse_move_happened = false; // reset to wait for first mouse move to happen
	n_trial_fixations = 0;
	last_fixation_dwell = 0;
	last_fixation_start_time = -1; // check for valid setting when used

}

void VisSearch_device::remove_probe()
{
	make_visual_object_disappear(probe_mouse_target_name);
	make_visual_object_disappear(probe_size_object_name);
	make_visual_object_disappear(probe_color_object_name);
	make_visual_object_disappear(probe_shape_object_name);
	make_visual_object_disappear(probe_label_object_name);
}

void VisSearch_device::remove_stimulus()
{	
	for(int i = 0; i < n_objects; ++i) {
		make_visual_object_disappear(search_objects[i].name);
		make_visual_object_disappear(search_objects[i].label_name);
		}
}

// make the stop signal a sound cue so that it can be perceived regardless
// of eye location
void VisSearch_device::present_stop_signal()
{
	stop_signal_object_name = concatenate_to_Symbol("StopSignal", n_trials);
	make_auditory_sound_event(stop_signal_object_name, Signal_c,
		GU::Point(0, 0), Symbol("Beep"), 500, 500);
//	make_visual_object_appear(stop_signal_object_name, GU::Point(0.5,0.5), GU::Size(30, 20));
//	set_visual_object_property(stop_signal_object_name, Color_c, Red_c);
//	set_visual_object_property(stop_signal_object_name, Shape_c, Rectangle_c);
}

void VisSearch_device::remove_stop_signal()
{
//	make_visual_object_disappear(stop_signal_object_name);
	stop_signal_object_name = Nil_c;
}


// here if an eyemovement start event
void VisSearch_device::handle_Eyemovement_Start_event(const Symbol& target_name, GU::Point new_location)
{
/*		device_out<< get_time() << ' ' << "Eye movement start received "
             << target_name << ' ' << new_location
            << " device state, mouse move, terminated " << state << ' ' << first_mouse_move_happened << ' ' << trial_terminated << endl;
*/
/*    if(first_mouse_move_happened  || trial_terminated)
		device_out<< get_time() << ' ' << "Eye movement start received after first mouse move or trial terminated "
             << target_name << ' ' << new_location
            << " device state, mouse move, terminated " << state << ' ' << first_mouse_move_happened << ' ' << trial_terminated << endl;
*/
	// last fixation in the trial is ignored because there is not another eye movement started until we are done
	if(state != WAITING_FOR_RESPONSE) {
//		device_out << "skip when not yet waiting for response" << endl;
		return;
		}

	// ignore if trial has been terminated due to excess fixations
	if(trial_terminated)
		return;
    
    // should we ignore if response has been started?
	
/* 021916 using first and last fixations, but not using dwell in last fixation
	// first fixation is ignored
	if(!n_trial_fixations) {
//		device_out << "skip first fixation" << endl;
		return;
		}
*/
	//	device_out << "em start at " << get_time() << endl;
	// here, aggregate the saved information from the previous fixation
	// but only if we actually have a bonafide start time for the fixation
	if(last_fixation_start_time >= 0) {
		last_fixation_dwell = get_time() - last_fixation_start_time;
		do_aggregate_fixation();
		}
}

// here if an eyemovement end event; identify what is fixated, but don't aggregate the stats
// except for spacial case in which we are terminating the trial - currently only partial aggregation
void VisSearch_device::handle_Eyemovement_End_event(const Symbol& fixation_target_name, GU::Point org_fixation_location)
{
/*		device_out<< get_time() << " Eye movement end received "
            << fixation_target_name << ' ' << org_fixation_location
            << " device state, mouse move, terminated " << state << ' ' << first_mouse_move_happened << ' ' << trial_terminated << endl;
*/
    // apply a noise function to the x and y coordinates of the fixation_location
    GU::Point fixation_location(add_measurement_error_noise(org_fixation_location));
    

	if(state != WAITING_FOR_RESPONSE) {
        // save when this fixation started in case actual search fixations are about to start
        last_fixation_start_time = get_time();
        // do nothing further
		return;
        }
	
	// if trial is already terminated, do nothing further
    if(trial_terminated)
		return;
	
	// save the information about this fixation to be accumulated when next one starts
	last_fixation_start_time = get_time();
	
//	device_out << "em end at " << get_time() << endl;
	
	// note saccade length
	last_fixation_target_name = fixation_target_name;
	last_saccade_distance = cartesian_distance(last_fixation_location, fixation_location);
	last_fixation_location = fixation_location;

	/*** fixation counted at eye movement end - start of fixation ***/
	n_trial_fixations++;

	// if maximum number of fixations has been exceeded, the trial is terminated
	// nothing further should be processed, but device must continue responding to events
	// because simulated human will not stop behaving immediately.

//	if(n_trial_fixations >= max_n_fixations) {
//	if(n_trial_fixations > max_n_fixations) {
	if(n_trial_fixations > max_fixations_criterion[current_trial_spec.cue_color][current_trial_spec.cue_size][current_trial_spec.cue_shape]) {
/*		if(device_out)
			device_out << processor_info() << "Excess fixations: "
				<< condition_label_string << ' ' << n_trials << ' '
				<< n_trial_fixations << endl;
*/
		// last fixation is not going to be counted, so decrement n_trial_fixations so
		// that accumulated statistics denominators will agree with it. - 032616
		n_trial_fixations--;
		
		trial_terminated = true; // trial is terminated due to excess fixations
//		trial_accumulators.excess_fixations = trial_terminated;
//        device_out << "*** TRIAL TERMINATED DUE TO EXCESS FIXATIONS: " << n_trial_fixations << " ***" << endl;
		n_terminated++;
		present_stop_signal();
        rt = get_time() - stimulus_onset_time; // record RT here, since trial is being terminated
		state = STOPPING_TRIAL;
		schedule_delay_event(1000);
		return;  // doing nothing further
		}
	
	// here if fixation should be processed

    // check in case the production rules attempted to look at part of the probe ...
	// device_out << "eyemovement_end at " << org_fixation_location << ' ' << fixation_target_name << endl;
	// count probe fixations
	if(match_probe_object(fixation_target_name)) {
//        device_out << "**** Fixation on probe object! ****" << endl;
		n_probe_fixations++;
		return;
		}
		
    // check to make sure object name is known = throws if not
//    find_object(eyemovement_target_name);

    // because of saccade noise, instead of taking the event's word for it, we compute the actual closest object ...

    // a function object class for comparing two objects in terms of the distance from the fixation_location
    struct Dist_from_fix {
        Dist_from_fix (GU::Point fix_loc_) : fix_loc(fix_loc_) {}
        bool operator() (const Search_object& obj1, const Search_object& obj2) const
            {
                double d1 = cartesian_distance(obj1.location, fix_loc);
                double d2 = cartesian_distance(obj2.location, fix_loc);
                return d1 < d2;
            }
        GU::Point fix_loc;
        };
 
    // designate the current fixated object using the current fixation location
    // always reset this so that if not using_AOI, data accumulation code is still valid
    last_fixation_is_inbetween = false;
    last_fixation_is_inAOI = false;
	// now find the closest object to the fixation location
	auto closest_it = min_element(search_objects.begin(), search_objects.end(), Dist_from_fix(fixation_location));
	Assert(closest_it != search_objects.end());
	last_closest_Search_object = *closest_it;
	last_closest_Search_object_distance = cartesian_distance(last_closest_Search_object.location, fixation_location);
    if(!using_AOI) {
        last_fixated_Search_object = last_closest_Search_object;
        }
    else { // using_AOI
        vector<Search_object> within_AOI_objs;
        copy_if(search_objects.begin(), search_objects.end(), back_inserter(within_AOI_objs),
            [this, fixation_location](const Search_object& obj) {return cartesian_distance(obj.location, fixation_location) <= AOI_radius_values[obj.isize];} );
        if(within_AOI_objs.empty()) {
//            device_out << "fixation is in-between objects" << endl;
            last_fixated_Search_object = Search_object();   // default ctor'd search object;
            last_fixation_is_inbetween = true;
            }
        else if(within_AOI_objs.size() == 1) {
            last_fixated_Search_object = *within_AOI_objs.begin();
			last_fixation_is_inAOI = true;
			}
        else if(within_AOI_objs.size() > 1) {
            device_out << "fixation is within multiple object AOIs, using closest" << endl;
            auto closest_AOIit = min_element(within_AOI_objs.begin(), within_AOI_objs.end(), Dist_from_fix(fixation_location));
            last_fixated_Search_object = *closest_AOIit;
			last_fixation_is_inAOI = true;
            }
        }

//	device_out << get_time() << " eyemovement_end at " << org_fixation_location << ' ' << fixation_target_name << ' ' << last_fixated_Search_object.name << endl;

}

// aggregate statistics about last fixation; called when eye starts to move to next fixation ...
void VisSearch_device::do_aggregate_fixation()
{
//	device_out << "aggregate fixation on closest object " << last_fixated_Search_object.name << endl;
	// pick out the set of accumulators for this condition
    // count how many times the fixation was in-between,or within AOI

	Assert(using_AOI && (last_fixation_is_inbetween && !last_fixation_is_inAOI || !last_fixation_is_inbetween && last_fixation_is_inAOI));
    trial_accumulators.inbetween_fixations.update(last_fixation_is_inbetween);  // accumulate proportions in this trial
    trial_accumulators.AOI_fixations.update(last_fixation_is_inAOI);

    // 021916 - do not accumulate dwell if first_mouse_move_happened - we must be aggregating the last fixation
    
    // 032316 - Policy:
    // All Trial match fixation tallies use ALL FIXATIONS in the denominator, so all proportions reported are a % of ALL FIXATIONS
    // All Trial averages like saccade distance use the number of fixations of that type (e.g. AOI) in their denominator
    // The dwell on the last fixation is not used in the average dwell
	
    // if the fixation was in-between, accumulate additional IBO statistics
    if(last_fixation_is_inbetween) {
        if(!first_mouse_move_happened) {
            trial_accumulators.inbetween_dwell.update(last_fixation_dwell);
//			device_out << "dwell on IBO " << last_fixated_Search_object.name << " = " << last_fixation_dwell << endl;
			}
        trial_accumulators.inbetween_saccade_distance.update(last_saccade_distance);
		trial_accumulators.inbetween_closest_distance.update(last_closest_Search_object_distance);
		}
	
    // following aggregations are conditional on inAOI, following above policy
    // count how many times the last_fixated_Search_object inAOI was not the eye movement target
	trial_accumulators.closest_not_target.update(last_fixation_is_inAOI && (last_fixated_Search_object.name != last_fixation_target_name));

    // calculate averages of AOI fixations - only in AOI counted here
	bool immediate_repeat = false;
    bool revisit_repeat = false;
    int repeat_lag = 0;
    
    if(last_fixation_is_inAOI) {
        if(!first_mouse_move_happened) {
            trial_accumulators.dwell.update(last_fixation_dwell);
//			device_out << "dwell on AOI " << last_fixated_Search_object.name << " = " << last_fixation_dwell << endl;
            }
        trial_accumulators.saccade_distance.update(last_saccade_distance);
        double distance_from_closest = cartesian_distance(last_fixated_Search_object.location, last_fixation_location);
        trial_accumulators.fixation_closest_distance.update(distance_from_closest);
    
        bool is_repeated = fixated_objects.find(last_fixated_Search_object.label) != fixated_objects.end();
        if(!is_repeated)
            fixated_objects.insert(last_fixated_Search_object.label);

        // classify a repeat fixation as either an immediate repeat or a revisit (at least one fixation between the two visits)
        fixated_object_frequencies[last_fixated_Search_object.label]++;
        if(fixated_object_frequencies[last_fixated_Search_object.label] > 1) {
            repeat_lag = n_trial_fixations - last_fixation_numbers[last_fixated_Search_object.label];
//			average repeat lag, inAOI objects only
			trial_accumulators.repeat_lags.update(repeat_lag);
			if(repeat_lag == 1) {
				immediate_repeat = true;
                revisit_repeat = false;
                }
            else {
                immediate_repeat = false;
                revisit_repeat = true;
                }
            }
        last_fixation_numbers[last_fixated_Search_object.label] = n_trial_fixations;
 /*
        device_out << "Fixation on " << last_fixated_Search_object.name << ' ' << last_fixated_Search_object.label << " is " << ((is_repeated) ? "repeated" : "first")
            << ((immediate_repeat) ? " immediate " : " ") << ((revisit_repeat) ? " revisit " : " ") << "lag: " << repeat_lag << endl;
*/
        // following should be true but only if an object was identified
        Assert(last_fixated_Search_object.isize >= 0 && last_fixated_Search_object.isize < n_sizes);

        // compute match/mismatch, fixation distances only if on an AOI object and if the dimension is cued; N is number of AOI fixations
        if(current_trial_spec.cue_color) {
            if(last_fixated_Search_object.icolor == target.icolor)
                trial_accumulators.color_match_distance.update(last_saccade_distance);
            else
                trial_accumulators.color_mismatch_distance.update(last_saccade_distance);
            }

        if(current_trial_spec.cue_size) {
            if(last_fixated_Search_object.isize == target.isize)
                trial_accumulators.size_match_distance.update(last_saccade_distance);
            else
                trial_accumulators.size_mismatch_distance.update(last_saccade_distance);
            }

        if(current_trial_spec.cue_shape) {
            if(last_fixated_Search_object.ishape == target.ishape)
                trial_accumulators.shape_match_distance.update(last_saccade_distance);
            else
                trial_accumulators.shape_mismatch_distance.update(last_saccade_distance);
            }
        
        // accumulate fixation distances for the three sizes
        trial_accumulators.size_accumulators[last_fixated_Search_object.isize].saccade_distance.update(last_saccade_distance);
        }
	
    // counting proportions, with ALL fixations in the denominator, which is why these are not inside the above conditional
	
    // count proportion of fixations that are an immediate repeat
    trial_accumulators.immediate_repeats.update(last_fixation_is_inAOI && immediate_repeat);
    // count proportion of fixations that are an revisit repeat
    trial_accumulators.revisits.update(last_fixation_is_inAOI && revisit_repeat);
	
	// count matches, distances on each dimension, only if cued
    trial_accumulators.color_matches.update(current_trial_spec.cue_color && last_fixation_is_inAOI && last_fixated_Search_object.icolor == target.icolor);
    trial_accumulators.size_matches.update(current_trial_spec.cue_size && last_fixation_is_inAOI && last_fixated_Search_object.isize == target.isize);
    trial_accumulators.shape_matches.update(current_trial_spec.cue_shape && last_fixation_is_inAOI && last_fixated_Search_object.ishape == target.ishape);
	// count matches for the three different sizes
	for(int i = 0; i < n_sizes; i++) {
		trial_accumulators.size_accumulators[i].all_fixations.update(
			last_fixation_is_inAOI && i == last_fixated_Search_object.isize);
		trial_accumulators.size_accumulators[i].color_matches.update(current_trial_spec.cue_color &&
			last_fixation_is_inAOI && last_fixated_Search_object.icolor == target.icolor && i == last_fixated_Search_object.isize);
		trial_accumulators.size_accumulators[i].size_matches.update(current_trial_spec.cue_size &&
			last_fixation_is_inAOI && last_fixated_Search_object.isize == target.isize && i == last_fixated_Search_object.isize);
		trial_accumulators.size_accumulators[i].shape_matches.update(current_trial_spec.cue_shape &&
			last_fixation_is_inAOI && last_fixated_Search_object.ishape == target.ishape && i == last_fixated_Search_object.isize);
		}
}


// here if a ply event is received
void VisSearch_device::handle_Ply_event(const Symbol& cursor_name, const Symbol& target_name,
		GU::Point new_location, GU::Polar_vector)
{
	if(trial_terminated) // do nothing if trial has been terminated
		return;
	if((state != WAITING_FOR_PRECUE_RESPONSE) && (state != WAITING_FOR_RESPONSE))
		throw Device_exception(this, "Ply received while not waiting for a response");
    
	// update the cursor position
	// sanity check
	Assert(Cursor_name_c == cursor_name);
	cursor_location = new_location;
	set_visual_object_location(Cursor_name_c, cursor_location);
	current_pointed_to_object = target_name;
	if(get_trace() && Trace_out && current_pointed_to_object != "nil")
		Trace_out << processor_info() << "Ply to: " << current_pointed_to_object << endl;
	else if(get_trace() && Trace_out) // shows intermediate points
		Trace_out << processor_info() << "Ply to: " << current_pointed_to_object << endl;
    
    // calculate rt from when mouse first starts to move
    if(state == WAITING_FOR_RESPONSE && !first_mouse_move_happened) {
        rt = get_time() - stimulus_onset_time;
        // this flag means that a response is in progress
        first_mouse_move_happened = true;
//dwell on AO		device_out << get_time()  << " first_mouse_move_happened" << endl;
        
        // 021916 - This is the signal to aggregate the data from the last fixation, but not the dwell ...
        do_aggregate_fixation();
        }
    
}


// here if a keystroke event is received - mouse buttons in this case
// here sequence to next trial state, or setup the next trial
void VisSearch_device::handle_Keystroke_event(const Symbol& key_name)
{
//    device_out << processor_info() << "Keystroke " << key_name << " received" <<
//            << " device state, mouse move, terminated " << state << ' ' << first_mouse_move_happened << ' ' << trial_terminated << endl;
	if(state == WAITING_FOR_PRECUE_RESPONSE) {
		remove_probe();
		state = PRESENT_SEARCH_OBJECTS;
		schedule_delay_event(1000); 
		return;
		}
    //
	else if(state != WAITING_FOR_RESPONSE) {
        // provide some information about what is going on here
        device_out << get_time() << " Keystroke " << key_name << " received while not waiting for a response " <<
        "device state, mouse move, terminated " << state << ' ' << first_mouse_move_happened << ' ' << trial_terminated << endl;
        // do nothing if the trial is being terminated; throw if not
       if(!trial_terminated)
            throw Device_exception(this, "Keystroke received while not waiting for a response");
        return;
        }
    // now that response has been received, reset the first mouse move happened flag
    first_mouse_move_happened = false;
    
     // identify the pointed-to object - should be the target
	const Search_object& pointed_to_obj = find_object(current_pointed_to_object);
	if(pointed_to_obj.label != target.label)
		throw Device_exception(this, "keystroke received for incorrect current pointed-to object");

   do_trial_end();
    }

void VisSearch_device::do_trial_end()
{
	// make the stimulus disappear
	remove_stimulus();
    
	
	// check the last object we looked at
//	const Search_object& fixated_obj = find_object(eyemovement_target_name);
//	if(fixated_obj.label != target.label)
//        throw Device_exception(this, "keystroke received for incorrect last fixated object");

    // no longer, because HZ stopped RT clock when mouse started to move
    // at this point, the RT is determined
//	rt = get_time() - stimulus_onset_time;
 		// identify the condition and accumulate the stats
		// pick out the set of accumulators for this condition
		Assert(current_trial_spec.cue_color >= 0 && current_trial_spec.cue_color < 2);
		Assert(current_trial_spec.cue_size >= 0 && current_trial_spec.cue_size < 2);
		Assert(current_trial_spec.cue_shape >= 0 && current_trial_spec.cue_shape < 2);
		Condition_accumulator& cond_accum =
			condition_accumulators[current_trial_spec.cue_color][current_trial_spec.cue_size][current_trial_spec.cue_shape];
	
	cond_accum.excess_fixations.update(trial_terminated);  // accumulate proportion of terminated trials
    
    // compute final statistics for this trial
	// if trial was not terminated, accumulate its data for the condition:
	if(!trial_terminated) {
//		trial_accumulators.excess_fixations = trial_terminated;
		trial_accumulators.rt = rt;
	//	cond_accum.n_trial_fixations.update(n_trial_fixations - 1);  // last fixation is counted in n_trial_fixations, so subtract one
		trial_accumulators.n_trial_fixations = n_trial_fixations;  // ??? last fixation is counted in n_trial_fixations, so subtract one
		trial_accumulators.unique_objects_fixated = fixated_objects.size();
		Mean_accumulator mean_fixation_frequency;
		for(auto pair : fixated_object_frequencies)
			mean_fixation_frequency.update(pair.second);
	//	cond_accum.mean_fixation_frequency.update(mean_fixation_frequency.get_mean());
		trial_accumulators.mean_fixation_frequency = mean_fixation_frequency.get_mean();
    
	/*    device_out << trial_accumulators.n_trial_fixations << ' ' << trial_accumulators.inbetween_fixations.get_n() << endl;
		device_out << trial_accumulators.n_trial_fixations << ' ' << trial_accumulators.inbetween_dwell.get_n() << endl;
		device_out << trial_accumulators.n_trial_fixations << ' ' << trial_accumulators.inbetween_dwell.get_n() << endl;

		device_out << "revisits " << trial_accumulators.revisits.get_count() << ' ' << trial_accumulators.revisits.get_n() << endl;
		device_out << "imm repeats " << trial_accumulators.immediate_repeats.get_count() << ' ' << trial_accumulators.immediate_repeats.get_n() << endl;
	*/

		trial_accumulators.check_denominators(device_out);
    

		// accumulate the trial data into the data for the condition
		cond_accum.update(trial_accumulators);
		cond_accum.check_denominators();

		} // end of acccumulate data if trial was not terminated
    
/*
	device_out << "trial_end " << n_trial_fixations << ' ' << cond_accum.n_trial_fixations.get_total() << ' '
		<< cond_accum.inbetween_fixation.get_n() << ' ' << cond_accum.AOI_fixation.get_n() << ' '
		<< cond_accum.n_trial_fixations.get_n() << ' '<< ' ' << cond_accum.rt.get_n()
		<< endl;
*/
	
/*	if(device_out)
		device_out << processor_info() << "Trial: " << n_trials << ", RT: " << rt 
			<< ", fixations: " << n_trial_fixations << endl;
*/

	if(device_out && (n_trials % 500) == 0)
		device_out << processor_info() << "Repetition: " << current_repetition << " Trial within repetition: " << current_trial_spec_index << " Overall trial: " << n_trials << endl;


    if(setup_next_trial()) {	// returns true if the next condition should be run
        // schedule a new trial
        state = START_TRIAL;
        schedule_delay_event(5000);
        }
    else {
        // stop the simulation after a delay
        state = SHUTDOWN;
        schedule_delay_event(500);
        }
    // do nothing further
}

// a trial has been completed - do what needs to be done
// return false if time to stop run; set up next trial if not.
bool VisSearch_device::setup_next_trial()
{
	n_trials++;
    
    // decide what to do for the next trial, if any
	// we run the current trial through the trial specs, where we fetch out a search field.
    // current_trial_spec_index is used as a subscript into the trial_specs,
	// so max value is n_trials_per_repetition - 1
	current_trial_spec_index++;
 //   device_out << current_trial_spec_index << ' ' << n_trials_per_repetition << endl;
	if(do_all_conditions) {
		if(current_trial_spec_index < n_trials_per_repetition) {
			Assert(current_trial_spec_index >= 0 && current_trial_spec_index < trial_specs.size());
			current_trial_spec = trial_specs[current_trial_spec_index];
			return true;	// keep going
			}
		// here if finished with all trials in this repetition; roll over into next repetition
		current_repetition++;
		if(current_repetition < n_repetitions) {
			current_trial_spec_index = 0;
			current_trial_spec = trial_specs[current_trial_spec_index];
			return true;	// keep going
			}
		}
	else {
	// here for single condition case, single_condition_spec set when run options are selected.
	// set the current_trial_spec_index to the next one matching the condition
		bool result = find_next_single_condition_trial_spec();
		// if we now have another single condition trial spec, go ahead and use it
		if(result) {
			return true;	// keep going
			}
		// here if finished with all trials in this repetition; roll over into next repetition
		else {
			current_repetition++;
			if(current_repetition < n_repetitions) {
				current_trial_spec_index = 0;
				bool next_first_single_condition_spec_trial_spec_found = find_next_single_condition_trial_spec();
				Assert(next_first_single_condition_spec_trial_spec_found);
				return true;	// keep going
				}
			}
		}
		
	// here if finished all repetitions
	output_statistics();    // dump the data
	return false;   // time to stop
}

// this function assumes that current_trial_spec_index has already been incremented after the current trial completes
bool VisSearch_device::find_next_single_condition_trial_spec()
{
	// if we are doing a single condition, scan through the trials specs looking for a match to the current condition
	// it will be the first trial; current_trial_spec_index is the index into trial_specs, not the actual trial in the run.
	for(; current_trial_spec_index < n_trials_per_repetition; current_trial_spec_index++) {
		current_trial_spec = trial_specs[current_trial_spec_index];
		if(current_trial_spec == single_condition_spec)
			return true;
		}
	return false;
}



void VisSearch_device::update_statistics()
{



}

bool VisSearch_device::match_probe_object(const Symbol& name) const
{
	return (name == probe_size_object_name || name == probe_color_object_name || name == probe_shape_object_name || name == probe_label_object_name);
}

const VisSearch_device::Search_object& VisSearch_device::find_object(const Symbol& name) const
{
	auto it = search_objects.begin();
	for(; it != search_objects.end(); ++it)
		if(it->name == name)
			return *it;
	if(it == search_objects.end())
		throw Device_exception(this, string("Looked at an unknown search object: ") + name.str());
	return *it;
}

void VisSearch_device::output_statistics()
{

	if(!out_file)
		return;
    out_file << fixed << setprecision(2);
    
//	device_out << "*** Condition: " << condition_label_string << endl;
	out_file << "N trials = " << '\t' << n_trials << endl;
//	device_out << setw(12) << " " << '\t' << setw(8) << "Mean" << '\t' << setw(8) << "tot N" << endl;
//	device_out << setw(12) << " " << setw(8) << "Mean" << setw(8) << "N Cor" << setw(8) << "Acc" << setw(8) << "tot N" << endl;

//	device_out << setw(12) << " " << '\t' << setw(8) << "Prop." << '\t' << setw(8) << "count" << '\t' << setw(8) << "tot N" << endl;

    // output differently depending on whether using trials specs or not
    // if we are:
	/* output lines as follows:
	Condition description using "---" to mean unspecified, Color, Size, Shape appearing in that order
	Overall  statistics: 
	mean rt, number of trials, mean saccade distance, mean number of saccades/trial, proportion of saccades that are repeat fixations,
	proportion of matches on color, size, shape to target object regardless of the cue type
	Then for each fixated object size, small, medium, large,
	saccade distance, number of saccades / tria. , prop. color matches, prop. size matches, prop. shape matches.
	Note: the proportion of matches corresponds to matches to the target object regardless of the cue type,
	so they must be interpreted with the cue type in mind.
	
	
	RT, N_fixations/trial, Fixations: Prop S, M, L, Color Match Prop: All, S, M, L, Size Match Prop: All, S, M, L,
		Shape Match Prop: All, S, M, L.
	
	*/
	
	output_condition(1, 0, 0);	// Color
	output_condition(1, 1, 0);	// Color+Size
	output_condition(1, 0, 1);	// Color+Shape
	output_condition(1, 1, 1);	// Color+Size+Shape
	output_condition(0, 1, 0);	// Size
	output_condition(0, 1, 1);	// Size+Shape
	output_condition(0, 0, 1);	// Shape
	output_condition(0, 0, 0);	// Number only
}

void VisSearch_device::output_condition(int cue_color, int cue_size, int cue_shape)
{
	/* output lines as follows:
	Condition description using "---" to mean unspecified, Color, Size, Shape appearing in that order
	Overall  statistics: 
	mean rt, number of trials, mean saccade distance, mean number of saccades/trial, proportion of saccades that are repeat fixations,
	proportion of matches on color, size, shape to target object regardless of the cue type
	Then for each fixated object size, small, medium, large,
	saccade distance, number of saccades / tria. , prop. color matches, prop. size matches, prop. shape matches.
	Note: the proportion of matches corresponds to matches to the target object regardless of the cue type,
	so they must be interpreted with the cue type in mind. ???
	
	
	Prop Excess fixations, Prop closest-not-target or IBO if using AOI, RT, N_fixations/trial, saccade-distance, closest-distance, Prop repeats;
	Fixations: Prop S, M, L, Color Match Prop: All, S, M, L, Size Match Prop: All, S, M, L, Shape Match Prop: All, S, M, L.
	proportion of immediate repeats, unique objects fixated, mean fixation frequency, average lag,
    color match/mismatch distance, ditto size, ditto shape
	
	
	2/13/16 saccade_distance is counting only those that are AOI fixations
	*/
    
	const Condition_accumulator& cond_accum = condition_accumulators[cue_color][cue_size][cue_shape];
    cond_accum.check_denominators();
    string condition_label = make_condition_label_string(cue_color, cue_size, cue_shape);

 //   Assert(cond_accum.inbetween_fixation.get_n() == cond_accum.AOI_fixation.get_n());
//	Assert(cond_accum.n_trial_fixations.get_total() == cond_accum.inbetween_fixation.get_n() + cond_accum.AOI_fixation.get_n());
	device_out << n_repetitions << ' ' << condition_label << ' ' << trials_per_condition[cue_color][cue_size][cue_shape] << ' '
		<< n_probe_fixations << ' ' << cond_accum.n_trial_fixations.get_total() << ' '
		<< cond_accum.n_trial_fixations.get_n() << ' ' << cond_accum.rt.get_n() << endl;

	out_file << setw(6) << condition_label << '\t';
	
//	int stat_n_trials = cond_accum.rt.get_n(); // get the number of trials in this condition
//	int stat_n_fixations = cond_accum.n_trial_fixations.get_total(); // get the total number of fixations counted in this condition
	// output overall statistics
    // output proportion of trials with excess fixations
    out_file << cond_accum.excess_fixations.get_proportion() << '\t';
	
    if(using_AOI) {
        // if using AOI, output the proportion of fixations that are in-between
        out_file << setw(6) << cond_accum.inbetween_fixations.get_mean() << '\t';
        }
    else {
        // if not using AOI, output proportion of fixations where the closest object is not the saccade target
        out_file << setw(6) << cond_accum.closest_not_target.get_mean() << '\t';
        }
	// output RT, mean number of fixations per trial, mean saccade distance
	out_file << setw(6) << cond_accum.rt.get_mean() / 1000. << '\t'
		<< setw(6) << cond_accum.n_trial_fixations.get_mean() << '\t'
        << setw(6) << cond_accum.inbetween_fixations.get_mean() << '\t'
        << setw(6) << cond_accum.inbetween_dwell.get_mean() << '\t'
        << setw(6) << cond_accum.inbetween_saccade_distance.get_mean() << '\t'
        << setw(6) << cond_accum.inbetween_closest_distance.get_mean() << '\t'
		<< setw(6) << cond_accum.AOI_fixations.get_mean() << '\t'
		<< setw(6) << cond_accum.dwell.get_mean() << '\t'
		<< setw(6) << cond_accum.saccade_distance.get_mean() << '\t';
	// output the distance to the closest, proportion of repeats
	out_file << setw(6) << cond_accum.fixation_closest_distance.get_mean() << '\t';
	out_file << setw(6) << cond_accum.revisits.get_mean() << '\t';

	// output the overall proportion of fixations on small, medium, large objects
/*	for(int isize = 0; isize < n_sizes; isize++) {
        double pfix_size = (stat_n_fixations > 0) ? cond_accum.size_accumulators[isize].saccade_distance.get_n() / double(stat_n_fixations) : 0.;
        out_file << setw(6) << pfix_size << '\t';
		}
*/
	// mean proportion fixations on each size
	for(int isize = 0; isize < n_sizes; isize++) {
		out_file << setw(6) << cond_accum.size_accumulators[isize].all_fixations.get_mean() << '\t';
		}
	// mean saccade distance to each size
	for(int isize = 0; isize < n_sizes; isize++) {
		out_file << setw(6) << cond_accum.size_accumulators[isize].saccade_distance.get_mean() << '\t';
		}

	// output proportion of fixations on matching color all sizes, then S, M, L
    out_file << setw(6) << cond_accum.color_matches.get_mean() << '\t';
	for(int isize = 0; isize < n_sizes; isize++) {
		out_file << setw(6) << cond_accum.size_accumulators[isize].color_matches.get_mean() << '\t';
		}
	// output proportion of fixations on matching size all sizes, then S, M, L
	out_file << setw(6) << cond_accum.size_matches.get_mean() << '\t';
	for(int isize = 0; isize < n_sizes; isize++) {
		out_file << setw(6) << cond_accum.size_accumulators[isize].size_matches.get_mean() << '\t';
		}
	// output proportion of fixations on matching shape all sizes, then S, M, L
	out_file << setw(6) << cond_accum.shape_matches.get_mean() << '\t';
	for(int isize = 0; isize < n_sizes; isize++) {
		out_file << setw(6) << cond_accum.size_accumulators[isize].shape_matches.get_mean() << '\t';
		}
	
	// output the proportion of immediate repeats, unique objects fixated, mean fixation frequency, average lag
	out_file << setw(6) << cond_accum.immediate_repeats.get_mean() << '\t'
		<< setw(6) << cond_accum.unique_objects_fixated.get_mean() << '\t'
		<< setw(6) << cond_accum.mean_fixation_frequency.get_mean() << '\t'
		<< setw(6) << cond_accum.repeat_lags.get_mean() << '\t';


    // for backwards compatibility here output the distances for color match/mismatch distance, ditto size, ditto shape
    out_file << setw(6) << cond_accum.color_match_distance.get_mean() << '\t' << setw(6) << cond_accum.color_mismatch_distance.get_mean() << '\t';
    out_file << setw(6) << cond_accum.size_match_distance.get_mean() << '\t' << setw(6) << cond_accum.size_mismatch_distance.get_mean() << '\t';
    out_file << setw(6) << cond_accum.shape_match_distance.get_mean() << '\t' << setw(6) << cond_accum.shape_mismatch_distance.get_mean();

	out_file << endl;

/*               device_out << setw(6) << cond_accum.rt.get_mean() / 1000. << '\t'
		<< setw(6) << stat_n_trials << '\t'
                    << setw(6) << cond_accum.saccade_distance.get_mean() << '\t'
		<< setw(6) << stat_n_fixations / double(stat_n_trials) << '\t'
                    << setw(6) << cond_accum.repeats.get_proportion() << '\t'
                    << setw(6) << cond_accum.color_matches.get_proportion() << '\t'
                    << setw(6) << cond_accum.size_matches.get_proportion() << '\t'
                    << setw(6) << cond_accum.shape_matches.get_proportion() << endl;
		
                // output fixation results for each object size
                for(int isize = 0; isize < n_sizes; isize++) {
                    const Fixation_accumulator& fix_accum = cond_accum.size_accumulators[isize];
                    device_out
	//		<< left << setw(6) <<size_names[isize] << '\t'
			<< setw(6) << fix_accum.saccade_distance.get_mean() << '\t'
			<< setw(6) << fix_accum.saccade_distance.get_n() / stat_n_fixations << '\t'
			<< setw(6) << fix_accum.color_matches.get_proportion() << '\t'
			<< setw(6) << fix_accum.size_matches.get_proportion() << '\t'
			<< setw(6) << fix_accum.shape_matches.get_proportion() << '\t'
			<< endl;
                    }
*/
}


