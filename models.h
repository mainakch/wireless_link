#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define C 299792458 // speed of light in SI units
#define PI 3.141592654
#define MAX_TRANSMITTERS 200 // maximum number of transmitters, real or virtual
#define MAX_RECEIVERS 2000 // maximum number of receivers
#define MAX_TIME 2000 // maximum number of receivers
#define LINE_LENGTH 256


typedef struct 
{
  double frequency;
  double wavelength;
  double delta_time;
  double max_limit;
  double min_limit;
  double boundary_tolerance;
  int total_time;

  bool generate_path;
  bool print_measurements;

} Simulation;

typedef struct 
{
  Simulation *sim;
  double velocity[3];
  double position[3];
  double acceleration_factor;
  double distance_from_actual_source;
  bool is_static;
} SpatialMotionModel;

typedef struct
{
  Simulation *sim;
  double power_in_dBm;
  double start_time;
  double end_time;
  double initial_phase;
  double doppler_offset;
} TransmissionModel;

typedef struct
{
  Simulation *sim;
  double distance;
} PropagationModel;

typedef struct
{
  SpatialMotionModel smm;
  TransmissionModel tm;
  int id;
} GeneralNode;

typedef struct
{
  GeneralNode gn;
  bool is_real_transmitter;
}  Transmitter;

typedef struct
{
  GeneralNode gn;
  double recv_noise_power;
}  Receiver;

typedef struct
{
  Simulation * sim;
  int num_receivers;
  int num_transmitters;
  int num_virtual_transmitters;
  Receiver receivers_array[MAX_RECEIVERS];
  Transmitter transmitters_array[MAX_TRANSMITTERS];
  int time;
  GeneralNode *node_array[MAX_RECEIVERS+MAX_TRANSMITTERS];
  int total_times[MAX_RECEIVERS+MAX_TRANSMITTERS];
  int total_gaussian_samples;
  double gaussiansamples[10*MAX_TIME];
  double **_nodepath;
} Environment;

int id();
void init_simulation(Simulation * sim_not_null);
void init_spatial_motion_model(SpatialMotionModel * smm, Simulation * sim_not_null);
void init_transmission_model(TransmissionModel * tm, Simulation * sim_not_null);
void init_propagation_model(PropagationModel * pm, Simulation * sim_not_null);
void init_general_node(GeneralNode * gn, Simulation * sim_not_null);
void init_transmitter(Transmitter * tn, Simulation * sim_not_null);
void init_receiver(Receiver * rc, Simulation * sim_not_null);
void init_environment(Environment * env, Simulation * sim_not_null);
void destroy_environment(Environment * env);
double _gaussrand(); // http://c-faq.com/lib/gaussian.html
void add_transmitter_with_position(Environment * env, double * pos);
void add_static_scatterer_with_position(Environment * env, double * pos);
void add_static_receiver(Environment * env, double * pos);
void _update_spatial_parameters(SpatialMotionModel * smm, GeneralNode * gn);
void update_all_locations(Environment * env);
void update_virtual_transmitters(Environment * env, int use_only_real);
double distance(GeneralNode * gn1, GeneralNode * gn2);
void readout_receiver_array(Environment * env);
void compute_shift(PropagationModel * pm, double * power_attn, double * phase_shift, Environment * env);
void populate_from_file(Environment * env, const char * filename);
void handle_request(Environment * env, FILE * fp, const char * req_type);
void print_to_file(Environment * env, const char * filename, bool print_state);
