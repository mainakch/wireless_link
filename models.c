#include "models.h"

int id()
{
  static int id=0;
  return (id++);
}

void init_simulation(Simulation * sim_not_null)
{
  sim_not_null->frequency = 60e9; // 60 Ghz 
  sim_not_null->wavelength = C/sim_not_null->frequency;
  sim_not_null->delta_time = 0.1;
  sim_not_null->max_limit = 1000;
  sim_not_null->min_limit = 0;
  sim_not_null->boundary_tolerance = 1;
  sim_not_null->total_time = 200;
}

void init_spatial_motion_model(SpatialMotionModel * smm)
{
  int ctr;
  for (ctr=0; ctr<3; ctr++)
  {
    smm->velocity[ctr] = 0;
    smm->position[ctr] = 5;
  }
  smm->distance_from_actual_source = 0;
  smm->acceleration_factor = 20;
  smm->is_static = false;
}

void init_transmission_model(TransmissionModel *tm)
{
  tm->power_in_dBm = -200;
  tm->start_time = 0;
  tm->end_time = 1000;
  tm->initial_phase = 0;
  tm->doppler_offset = 0;
}

void init_propagation_model(PropagationModel * pm)
{
  pm->distance = 0;
}

void init_general_node(GeneralNode * gn)
{
  init_spatial_motion_model(&(gn->smm));
  init_transmission_model(&(gn->tm));
  gn->id = id();
}

void init_transmitter(Transmitter * tn)
{
  init_general_node(&(tn->gn));
  tn->is_real_transmitter = false;
}

void init_receiver(Receiver * rc)
{
  init_general_node(&(rc->gn));
  rc->recv_noise_power = pow(10, -16);
}

void init_environment(Environment * env)
{
  if (!(env->num_transmitters > 0 && env->num_receivers > 0)
    || env->total_time < 1)
  {
    fprintf(stderr, "Please set non zero tx, rx and total time!\n");
    exit(EXIT_FAILURE);
  }

  env->_num_transmitters_ctr = 0;
  env->_num_receivers_ctr = 0;
  env->num_virtual_transmitters = 0;
  env->time = 0;
  init_environment_malloc(env);
}

void init_environment_malloc(Environment * env)
{
  fprintf(stderr, "Inside malloc for environment\n");

  env->receivers_array = calloc(env->num_receivers, sizeof(Receiver));
  env->transmitters_array = calloc(env->num_transmitters, sizeof(Transmitter));
  env->node_array = calloc((env->num_receivers + env->num_transmitters), sizeof(GeneralNode *));
  env->unit_power_gaussian_noise = calloc((env->num_receivers), sizeof(complex double));
}

void init_filereader(Filereader * fr)
{
  printf("Opening %s for writing \n", fr->output_filename);
  printf("Opening %s for reading \n", fr->input_filename); 
  
  fr->infile = fopen(fr->input_filename, "r");
  fr->outfile = fopen(fr->output_filename, "w");  
}

void destroy_environment(Environment * env)
{
  free(env->receivers_array);
  free(env->transmitters_array);
  free(env->node_array);
  free(env->unit_power_gaussian_noise);
}

void destroy_filereader(Filereader * fr)
{
  if (fr->infile != NULL)  fclose(fr->infile);
  if (fr->outfile != NULL)  fclose(fr->outfile);
}

double _gaussrand() // http://c-faq.com/lib/gaussian.html
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}

/* void add_node_to_environment_array(Environment * env, GeneralNode * gn) */
/* { */
/*   env->node_array[gn->id] = gn; */
/* } */

/* void print_node_path(Environment * env) */
/* { */
/*   printf("Printing node paths\n"); */
/*   int ctr, ctr1; */
/*   ctr1=0; */
/*   double * db = env->_nodepath[ctr1]; */
/*   while (db != NULL) */
/*   { */
/*     for (ctr=0; ctr<60; ctr++) */
/*     { */
/*       printf("%lf ", *(db+ctr)); */
/*     } */
/*     printf("\n"); */
/*     db = env->_nodepath[++ctr1]; */
/*   } */
/*   printf("Total number of nodes is %d \n", ctr1); */
/* } */

/* void print_current_locations(Environment * env) */
/* { */
/*   int ctr=0, ctr1; */
/*   GeneralNode * gn = env->node_array[ctr]; */
/*   while(gn != NULL) */
/*   { */
/*     printf("Node %d \n", ctr); */
/*     for(ctr1=0; ctr1<3; ctr1++) */
/*     { */
/*       printf("%lf ", gn->smm.position[ctr1]); */
/*     } */
/*     printf("\n"); */
/*     gn = env->node_array[++ctr]; */
/*   } */
/* } */

double gaussenv(Environment * env)
{
  static int numcalls=0;
  //return env->gaussiansamples[numcalls++];
  return _gaussrand();
}

void parse_input(int argc, char * argv[], Filereader * fr)
{
  strcpy(fr->input_filename, "/tmp/input.txt");
  strcpy(fr->output_filename, "/tmp/output.txt");

  struct option long_options[] = 
  {
    {"input_filename", required_argument, 0, 'i'},
    {"output_filename", required_argument, 0, 'o'},
    {0, 0, 0, 0}
  };

  int c=0, option_index;
  c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
  while (c != -1)
  {
    switch (c)
    {
      case 0:
        break;
      case 'i':
        strcpy(fr->input_filename, optarg);
        break;
      case 'o':
        strcpy(fr->output_filename, optarg);
	break;
      default:
        fprintf(stderr, "Error parsing options\n");
	exit(EXIT_FAILURE);
    }
    c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
  }
}
