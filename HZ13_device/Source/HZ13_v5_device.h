#ifndef VISSEARCH_DEVICE_H
#define VISSEARCH_DEVICE_H

/*
Visual field is 39 X 30 degrees, divided up into a 13 X 13 ??? grid of 3 degree square.
Objects occupy 75 of the squares, chosen at random, except the center is 
reserved for the target description

stimuli are in three sizes: 
2.8, 1.6, 0.8 degrees, called large, medium, and small

5 colors:
blue, green, yellow, orange, pink ???

5 shapes:
circle, semicircle, triangle, square, cross

each stimulus has a two digit number 1-75

probe is mouse target XX, then color, size, shape, number


trial start
probe cue presented
after study, user points, clicks on XX
cue disappears,
search screen appears
point to and click on target.

time RT from stimulus appear to mouse starts to move

do not record first and last fixation in the trial (like HZ13 reduction)

*/


#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cassert>

#include "EPICLib/Device_base.h"
#include "EPICLib/Symbol.h"
#include "EPICLib/Statistics.h"
#include "EPICLib/Geometry.h"
#include "EPICLib/Assert_throw.h"

namespace GU = Geometry_Utilities;

// see implementation file for conversion functions and constants 
// pertaining to screen layout

class VisSearch_device : public Device_base {
public:
//	VisSearch_device(const std::string& id, Output_tee& ot, int n_total_trials_);
	VisSearch_device(const std::string& id, Output_tee& ot);
			
	virtual void initialize();
	virtual void set_parameter_string(const std::string&);
	virtual std::string get_parameter_string() const;

	virtual void handle_Start_event();
	virtual void handle_Stop_event();
	virtual void handle_Delay_event(const Symbol& type, const Symbol& datum, 
		const Symbol& object_name, const Symbol& property_name, const Symbol& property_value);
    virtual void handle_Ply_event(const Symbol& cursor_name, const Symbol& target_name,
		GU::Point new_location, GU::Polar_vector);
	virtual void handle_Keystroke_event(const Symbol& key_name);
	virtual void handle_Eyemovement_Start_event(const Symbol& target_name, GU::Point new_location);
	virtual void handle_Eyemovement_End_event(const Symbol& target_name, GU::Point new_location);
		
private:
	struct Search_object {
		Symbol name;
		Symbol label_name;
		int id;	// an ID number - turned into object name and label
		GU::Point location;
		GU::Size size;
		Symbol size_name;
		int isize;
		Symbol color;
		int icolor;
		Symbol shape;
		int ishape;
		Symbol label;
		
		// default ctor needed for container
		Search_object() : id(-1) {}
		// set name, label from id,
		Search_object(int id_, GU::Point location_, GU::Size size_,
			const Symbol& size_name_, int isize_, const Symbol& color_, int icolor_,
			const Symbol& shape_, int ishape_);
		};

	struct Search_field {
		int target_index;
		std::vector<Search_object> objects;
		};
		
    // holds a trial search field and condition cues;
    // cues coded 0 - not cued, 1 - cued for use as data accumulation subscripts
    struct Trial_spec {
		Trial_spec (int search_field_index_ = 0, int cue_color_ = 0, int cue_size_ = 0, int cue_shape_ = 0) :
			search_field_index(search_field_index_), cue_color(cue_color_), cue_size(cue_size_), cue_shape(cue_shape_) {}
		int search_field_index;
		int cue_color;
		int cue_size;
		int cue_shape;
		// equality means cue conditions are the same - for implementing single-condition scan of trial specs
		bool operator== (const Trial_spec& rhs) const
			{return cue_color == rhs.cue_color && cue_size == rhs.cue_size && cue_shape == rhs.cue_shape;}
	};
    
	enum State_e {START, PRESENT_CURSOR, START_TRIAL, STOPPING_TRIAL,
		PRESENT_PROBE, WAITING_FOR_PRECUE_RESPONSE, PRESENT_SEARCH_OBJECTS,
		WAITING_FOR_RESPONSE,SHUTDOWN};
	
	State_e state;
	
	// experiment parameters
	const int n_search_fields = 96;
	const int n_practice_trials = 24;	// have to deduct from input search field number
	const int n_trials_per_repetition = 22 * n_search_fields;
    const int n_conditions = 8;

    // experimental conditions control
	std::string condition_string;
    bool do_all_conditions;             // do all conditions
	int n_repetitions = 1;			// number of repeats of trial specifications

	std::vector<Search_field> search_fields; // the search fields to cycle through
	std::vector<Trial_spec> trial_specs;	 // the specs for each trial in a repetition

	Trial_spec single_condition_spec;		// if only one condition is used
	Trial_spec current_trial_spec;		// the current condition being used for this trial
    int current_repetition;
    int current_trial_spec_index;		// this is subscript into trial_specs
    int current_search_field;
   

	// stimulus parameters
	static const int n_objects = 75;
	static const int n_colors = 5;
	static const int n_sizes = 3;
	static const int n_shapes = 5;
	
	// stimulus generation
	// search field file
	std::ifstream search_fields_file;
	std::ifstream trial_specs_file;
    std::ofstream out_file;
	void load_trial_data_files();
	void load_search_field();
	void load_trial_spec();
	
	// maps and vectors to translate search_field entries into indices and Symbols
	std::map<std::string, int> color_index_map;
	std::vector<Symbol> colors;
	std::map<std::string, int> size_index_map;
	std::vector<Symbol> size_names;
	std::vector<double> size_values;
	std::vector<double> AOI_radius_values;
	std::map<std::string, int> shape_index_map;
	std::vector<Symbol> shapes;
	std::map<int, Geometry_Utilities::Size> shape_size_factors; // shape index to size factor
	typedef std::vector<Symbol> labels_t;
	labels_t labels;

	// display state
	long object_counter;
	Symbol probe_mouse_target_name;
	Symbol probe_size_object_name;
	Symbol probe_color_object_name;
	Symbol probe_shape_object_name;
	Symbol probe_label_object_name;

	// this holds the set of search objects put on the display
	// cell [0] has object with label of 1
	std::vector<Search_object> search_objects;
				
	int target_num;
	Search_object target;
	GU::Point cursor_location;
	Symbol current_pointed_to_object;
    Symbol stop_signal_object_name;
	
	// data accumulation
    bool using_AOI = false;  // true if fixations are resolved using AOI and last_fixation_is_inbetween, false if closest used
    bool first_mouse_move_happened = false; // record RT when first mouse move happens
	long stimulus_onset_time;
    long rt;
	int n_trial_fixations;
	int n_trials;
	int n_probe_fixations;
	int max_n_fixations = 300;	// stop trial if this many fixations or more
	bool trial_terminated = false;	// flag so show that trial has been terminated
    bool cursor_in_motion = false; // flag that cursor is in motion for a response
	int n_terminated = 0;		// counts number of terminated trials

	Symbol last_fixation_target_name;  // from eye movement end event
	GU::Point last_fixation_location;	// from eye movement end event
	double last_saccade_distance;
	long last_fixation_start_time = 0;
	long last_fixation_dwell = 0;
	Search_object last_closest_Search_object;
	Search_object last_fixated_Search_object;
	double last_closest_Search_object_distance;
    bool last_fixation_is_inbetween = false;
	bool last_fixation_is_inAOI = false;
//	GU::Point old_location;
	std::set<Symbol> fixated_objects;
	std::map<Symbol, int> fixated_object_frequencies; // object name, #times fixated
	std::map<Symbol, int> last_fixation_numbers;	// object name, fixation number last fixated
	// lag will be this fixation - last_fixation_number

	std::string condition_label_string;
    
    // Data aggregation concept for equal trial weight (non-homogenous trials)
    // During each trial accumulate proportions and means where denominator = number of fixations in that trial.
    // At end of each trial accumulate mean proportion and means where denominator = number of trials
    
    // accumulate these statistics during each trial, then accumulate for each condition
    struct Trial_accumulator {
//        bool excess_fixations;  // false or true
        Proportion_accumulator closest_not_target;
        double rt;
        int n_trial_fixations;  // accumulates number of fixations in a trial
        Proportion_accumulator inbetween_fixations;  // inbetween accumulated only if using AOI
        Mean_accumulator inbetween_dwell;
        Mean_accumulator inbetween_saccade_distance;
        Mean_accumulator inbetween_closest_distance;
        Proportion_accumulator AOI_fixations;		// accumulated only if using AOI
        Mean_accumulator dwell;
        Mean_accumulator saccade_distance;
        Mean_accumulator fixation_closest_distance;
        Proportion_accumulator revisits;
        
		Proportion_accumulator color_matches;
		Proportion_accumulator size_matches;
		Proportion_accumulator shape_matches;
        
        // This is used to accumulate statistics for each object size
        // denominator = number of fixations on objects of a given size in a trial
        struct Trial_Fixation_accumulator {
			Proportion_accumulator all_fixations;	// counts proportion of fixations for this object size, regardless of matches
            Mean_accumulator saccade_distance;
            Proportion_accumulator color_matches;
            Proportion_accumulator size_matches;
            Proportion_accumulator shape_matches;
            void reset() {
                all_fixations.reset();
                saccade_distance.reset();
                color_matches.reset();
                size_matches.reset();
                shape_matches.reset();
                }
            void check_denominators(int n) const {
                Assert(n == all_fixations.get_n());
                Assert(saccade_distance.get_n() == all_fixations.get_count());
                Assert(n == color_matches.get_n());
                Assert(n == size_matches.get_n());
                Assert(n == shape_matches.get_n());
                }
            };
    
		// an array of accumulators, one for each size
		Trial_Fixation_accumulator size_accumulators[n_sizes];

        Proportion_accumulator immediate_repeats;
		int unique_objects_fixated;
		double mean_fixation_frequency;
        Mean_accumulator repeat_lags;
        
        Mean_accumulator color_match_distance;
        Mean_accumulator color_mismatch_distance;
        Mean_accumulator size_match_distance;
        Mean_accumulator size_mismatch_distance;
        Mean_accumulator shape_match_distance;
        Mean_accumulator shape_mismatch_distance;
        
       void reset() {
//			excess_fixations = false;
            closest_not_target.reset();
            rt = 0;;
            n_trial_fixations = 0;
            inbetween_fixations.reset();
            inbetween_dwell.reset();
            inbetween_saccade_distance.reset();
            inbetween_closest_distance.reset();
            AOI_fixations.reset();
            dwell.reset();
            saccade_distance.reset();
            fixation_closest_distance.reset();
            revisits.reset();

            color_matches.reset();
            size_matches.reset();
            shape_matches.reset();
           
            for(auto& accumulator : size_accumulators)
                accumulator.reset();

            immediate_repeats.reset();
 			unique_objects_fixated = 0;
			mean_fixation_frequency = 0;
			repeat_lags.reset();
           
            color_match_distance.reset();
            color_mismatch_distance.reset();
            size_match_distance.reset();
            size_mismatch_distance.reset();
            shape_match_distance.reset();
            shape_mismatch_distance.reset();
            }
        
        void check_denominators(Output_tee& ot) const;
        };
    
    // denominator = number of trials in that condition  // not finished editing yet!
    // mean proportions are accumulated, not pooled proportions
    struct Condition_accumulator {
        Proportion_accumulator excess_fixations; // accumulates proportion of terminated trials
        Mean_accumulator closest_not_target; // accumulates mean proportion
        Mean_accumulator rt;
        Mean_accumulator n_trial_fixations;  // accumulates number of fixations in a trial
        Mean_accumulator inbetween_fixations;  // inbetween accumulated only if using AOI
        Mean_accumulator inbetween_dwell;
        Mean_accumulator inbetween_saccade_distance;
        Mean_accumulator inbetween_closest_distance;
        Mean_accumulator AOI_fixations;		// accumulated only if using AOI
        Mean_accumulator dwell;
        Mean_accumulator saccade_distance;
        Mean_accumulator fixation_closest_distance;
        Mean_accumulator revisits;
        
		Mean_accumulator color_matches;
		Mean_accumulator size_matches;
		Mean_accumulator shape_matches;
        
        // This is used to accumulate statistics for each object size
        // denominator = number of fixations on objects of a given size in a trial
        struct Condition_Fixation_accumulator {
            Mean_accumulator all_fixations;
            Mean_accumulator saccade_distance;
            Mean_accumulator color_matches;
            Mean_accumulator size_matches;
            Mean_accumulator shape_matches;
            void reset() {
                all_fixations.reset();
                saccade_distance.reset();
                color_matches.reset();
                size_matches.reset();
                shape_matches.reset();
                }
            void update(const Trial_accumulator::Trial_Fixation_accumulator& trial_size_accumulator) {
                all_fixations.update(trial_size_accumulator.all_fixations.get_proportion());
                saccade_distance.update(trial_size_accumulator.saccade_distance.get_mean());
                color_matches.update(trial_size_accumulator.color_matches.get_proportion());
                size_matches.update(trial_size_accumulator.size_matches.get_proportion());
                shape_matches.update(trial_size_accumulator.shape_matches.get_proportion());
                }
            void check_denominators(int n) const {
                Assert(n == all_fixations.get_n());
                Assert(n == saccade_distance.get_n());
                Assert(n == color_matches.get_n());
                Assert(n == size_matches.get_n());
                Assert(n == shape_matches.get_n());
                }
            };
    
		// an array of accumulators, one for each size
		Condition_Fixation_accumulator size_accumulators[n_sizes];

        Mean_accumulator immediate_repeats;
		Mean_accumulator unique_objects_fixated;
		Mean_accumulator mean_fixation_frequency;
        Mean_accumulator repeat_lags;
        
        Mean_accumulator color_match_distance;
        Mean_accumulator color_mismatch_distance;
        Mean_accumulator size_match_distance;
        Mean_accumulator size_mismatch_distance;
        Mean_accumulator shape_match_distance;
        Mean_accumulator shape_mismatch_distance;
        
       void reset() {
			excess_fixations.reset();
            closest_not_target.reset();
            rt.reset();
            n_trial_fixations.reset();
            inbetween_fixations.reset();
            inbetween_dwell.reset();
            inbetween_saccade_distance.reset();
            inbetween_closest_distance.reset();
            AOI_fixations.reset();
            dwell.reset();
            saccade_distance.reset();
            fixation_closest_distance.reset();
            revisits.reset();

            color_matches.reset();
            size_matches.reset();
            shape_matches.reset();
           
            for(auto& accumulator : size_accumulators)
                accumulator.reset();

            immediate_repeats.reset();
 			unique_objects_fixated.reset();
			mean_fixation_frequency.reset();
			repeat_lags.reset();
           
            color_match_distance.reset();
            color_mismatch_distance.reset();
            size_match_distance.reset();
            size_mismatch_distance.reset();
            shape_match_distance.reset();
            shape_mismatch_distance.reset();
            }
        
        // update from a trial accumulator
        void update(const Trial_accumulator& trial) {
//			excess_fixations.update(trial.excess_fixations);  // updated in trial_end processing
            closest_not_target.update(trial.closest_not_target.get_proportion());
            rt.update(trial.rt);
            n_trial_fixations.update(trial.n_trial_fixations);
            inbetween_fixations.update(trial.inbetween_fixations.get_proportion());
            inbetween_dwell.update(trial.inbetween_dwell.get_mean());
            inbetween_saccade_distance.update(trial.inbetween_saccade_distance.get_mean());
            inbetween_closest_distance.update(trial.inbetween_closest_distance.get_mean());
            AOI_fixations.update(trial.AOI_fixations.get_proportion());
            dwell.update(trial.dwell.get_mean());
            saccade_distance.update(trial.saccade_distance.get_mean());
            fixation_closest_distance.update(trial.fixation_closest_distance.get_mean());
            revisits.update(trial.revisits.get_proportion());

            color_matches.update(trial.color_matches.get_proportion());
            size_matches.update(trial.size_matches.get_proportion());
            shape_matches.update(trial.shape_matches.get_proportion());
           
            for(int i = 0; i < 3; i++)
                size_accumulators[i].update(trial.size_accumulators[i]);

            immediate_repeats.update(trial.immediate_repeats.get_proportion());;
 			unique_objects_fixated.update(trial.unique_objects_fixated);
			mean_fixation_frequency.update(trial.mean_fixation_frequency);;
			repeat_lags.update(trial.repeat_lags.get_mean());
           
            color_match_distance.update(trial.color_match_distance.get_mean());
            color_mismatch_distance.update(trial.color_mismatch_distance.get_mean());
            size_match_distance.update(trial.size_match_distance.get_mean());
            size_mismatch_distance.update(trial.size_mismatch_distance.get_mean());
            shape_match_distance.update(trial.shape_match_distance.get_mean());
            shape_mismatch_distance.update(trial.shape_mismatch_distance.get_mean());
            }
        
        void check_denominators() const {
            int rt_n = rt.get_n(); // should be number of trials in this condition
            Assert(rt_n == n_trial_fixations.get_n());
            Assert(rt_n == inbetween_fixations.get_n());
            Assert(rt_n == inbetween_dwell.get_n());
            Assert(rt_n == inbetween_saccade_distance.get_n());
            Assert(rt_n == inbetween_closest_distance.get_n());
            Assert(rt_n == AOI_fixations.get_n());
            Assert(rt_n == dwell.get_n());
            Assert(rt_n == saccade_distance.get_n());
            Assert(rt_n == fixation_closest_distance.get_n());
            Assert(rt_n == revisits.get_n());
            Assert(rt_n == color_matches.get_n());
            Assert(rt_n == size_matches.get_n());
            Assert(rt_n == shape_matches.get_n());
            for(auto& accumulator : size_accumulators)
                accumulator.check_denominators(rt_n);
            Assert(rt_n == immediate_repeats.get_n());
            Assert(rt_n == unique_objects_fixated.get_n());
            Assert(rt_n == mean_fixation_frequency.get_n());
            Assert(rt_n == repeat_lags.get_n());
            Assert(rt_n == color_match_distance.get_n());
            Assert(rt_n == color_mismatch_distance.get_n());
            Assert(rt_n == size_match_distance.get_n());
            Assert(rt_n == size_mismatch_distance.get_n());
            Assert(rt_n == shape_match_distance.get_n());
            Assert(rt_n == shape_mismatch_distance.get_n());
            }
        
        
        };
    
    Trial_accumulator trial_accumulators;    // just one is needed
    
    Condition_accumulator condition_accumulators[2][2][2];

	/* from Yunfeng 12/12/15
	Number	295 NOnly 0
Shape	253  Shp 1
Color	73   Clr 2
Color+Shape	78 ClrShp 3
Size	190 Siz 4
Size+Shape	123 SizShp 5
Color+Size	66 ClrSiz 6
All	58 ClrSizShp 7
*/
//	array<int, 8> excess_fixations_in_trial_criteria = {295, 253, 73, 78, 190, 123, 66, 58};

	int max_fixations_criterion[2][2][2] = {};
	int trials_per_condition[2][2][2] = {};
		
	// helpers
	void parse_condition_string();
	int lookup_index(const std::map<std::string, int>& themap, const std::string& key);
	bool find_next_single_condition_trial_spec();

	bool setup_next_trial();
	void initialize_data_collection();
    void reset_condition_accumulators();
    void present_cursor();
	void prepare_trial();
	void present_probe();
	void present_search_objects();
	void present_search_object(Search_object& search_object, const char * name_prefix);
	void remove_probe();
	void remove_stimulus();
	void do_aggregate_fixation();
    void do_trial_end();
    void present_stop_signal();
    void remove_stop_signal();
	bool match_probe_object(const Symbol& name) const;
	const Search_object& find_object(const Symbol& name) const;

    void update_statistics();
	void output_condition(int cue_color, int cue_size, int cue_shape);
	void output_statistics();

	// rule out default ctor, copy, assignment
	VisSearch_device(const VisSearch_device&);
	VisSearch_device& operator= (const VisSearch_device&);
};

#endif


