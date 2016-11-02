#include "models.h"

int id()
{
  static int id=0;
  return (id++);
}

void init_simulation(Simulation * sim_not_null)
{
  sim_not_null->generate_path = true;
  sim_not_null->print_measurements = true;
  sim_not_null->frequency = 60e9; // 60 Ghz 
  sim_not_null->wavelength = C/sim_not_null->frequency;
  sim_not_null->delta_time = 0.1;
  sim_not_null->max_limit = 1000;
  sim_not_null->min_limit = 0;
  sim_not_null->boundary_tolerance = 1;
  sim_not_null->total_time = 200;
}

void init_spatial_motion_model(SpatialMotionModel * smm, Simulation * sim_not_null)
{
  smm->sim = sim_not_null;
  int ctr;
  for (ctr=0; ctr<3; ctr++)
  {
    smm->velocity[ctr] = 0;
    smm->position[ctr] = smm->sim->max_limit/2;
  }
  smm->distance_from_actual_source = 0;
  smm->acceleration_factor = 20;
  smm->is_static = false;
}

void init_transmission_model(TransmissionModel *tm, Simulation * sim_not_null)
{
  tm->power_in_dBm = -200;
  tm->start_time = 0;
  tm->end_time = 1000;
  tm->initial_phase = 0;
  tm->doppler_offset = 0;
}

void init_propagation_model(PropagationModel * pm, Simulation * sim_not_null)
{
  pm->sim = sim_not_null;
  pm->distance = 0;
}

void init_general_node(GeneralNode * gn, Simulation * sim_not_null)
{
  init_spatial_motion_model(&(gn->smm), sim_not_null);
  init_transmission_model(&(gn->tm), sim_not_null);
  gn->id = id();
}

void init_transmitter(Transmitter * tn, Simulation * sim_not_null)
{
  init_general_node(&(tn->gn), sim_not_null);
  tn->is_real_transmitter = false;
}

void init_receiver(Receiver * rc, Simulation * sim_not_null)
{
  init_general_node(&(rc->gn), sim_not_null);
  rc->recv_noise_power = pow(10, -16);
}

void init_environment(Environment * env, Simulation * sim_not_null)
{
  env->num_receivers = 0;
  env->num_transmitters = 0;
  env->num_virtual_transmitters = 0;
  env->sim = sim_not_null;
  env->time = 0;
  env->_nodepath = malloc((MAX_TRANSMITTERS+MAX_RECEIVERS)*sizeof(double *));
  env->total_gaussian_samples = 0;

  int ctr;
  for (ctr=0; ctr<MAX_TRANSMITTERS+MAX_RECEIVERS; ctr++)
  {
    env->total_times[ctr] = 0;
  }
  if (!sim_not_null->generate_path)
  {
    populate_from_file(env, "/tmp/input.txt");
  }
}

void destroy_environment(Environment * env)
{
  int ctr=0;
  while(env->_nodepath[ctr] != NULL)
  {
    free(env->_nodepath[ctr]);
    ctr++;
  }
  free(env->_nodepath);
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

void add_transmitter_with_position(Environment * env, double * pos)
{
  Transmitter tn;
  init_transmitter(&tn, env->sim);
  int ctr;
  for (ctr=0; ctr<3; ctr++)
  {
    tn.gn.smm.position[ctr] = *(pos+ctr);
  }
  tn.gn.tm.power_in_dBm = 0; // set power to 0 dBm instead of -200
  tn.is_real_transmitter = true;
  env->transmitters_array[env->num_transmitters] = tn;
  GeneralNode * gnp =   &(env->transmitters_array[env->num_transmitters].gn);
  add_node_to_environment_array(env, gnp);
  (env->num_transmitters)++;
}

void add_static_scatterer_with_position(Environment * env, double * pos)
{
  add_transmitter_with_position(env, pos);
  Transmitter * tn = &(env->transmitters_array[env->num_transmitters]);
  tn->is_real_transmitter = false;
  tn->gn.smm.acceleration_factor = 0;
  tn->gn.smm.is_static = true;
  GeneralNode * gnp = &tn->gn;
  add_node_to_environment_array(env, gnp);
  (env->num_transmitters)++;  
}

void add_static_receiver(Environment * env, double * pos)
{
  Receiver rx;
  init_receiver(&rx, env->sim);
  rx.gn.smm.acceleration_factor = 0;
  rx.gn.smm.is_static = true;
  env->receivers_array[env->num_receivers] = rx;
  int ctr;
  for (ctr=0; ctr<3; ctr++)
  {
    env->receivers_array[env->num_receivers].gn.smm.position[ctr] = *(pos+ctr);
  }
  GeneralNode * gnp = &(env->receivers_array[env->num_receivers].gn);
  add_node_to_environment_array(env, gnp);
  (env->num_receivers)++;
}

void _update_spatial_parameters(SpatialMotionModel * smm, GeneralNode * gn)
{
  int ctr;

  for (ctr=0; ctr<2; ctr++) // no updates for the third dimension
  {
    if (smm->position[ctr]> smm->sim->max_limit - smm->sim->boundary_tolerance)
    {
      smm->position[ctr] = smm->sim->max_limit - 2*smm->sim->boundary_tolerance;
      smm->velocity[ctr] *= -1;
    }
	
    if (smm->position[ctr]< smm->sim->min_limit + smm->sim->boundary_tolerance)
    {
      smm->position[ctr] = smm->sim->min_limit + 2*smm->sim->boundary_tolerance;
      smm->velocity[ctr] *= -1;
    }    
    
    smm->position[ctr] += smm->velocity[ctr]*smm->sim->delta_time;
    smm->velocity[ctr] += smm->acceleration_factor*_gaussrand()*smm->sim->delta_time;

  }
  /* printf("Nodepath %d 1 ", gn->id); */
  /* for (ctr=0; ctr<3; ctr++) */
  /* { */
  /*   printf("%f %f ", smm->position[ctr], smm->velocity[ctr]); */
  /* } */
  /* printf("\n"); */
}

void update_all_locations(Environment * env)
{
  int ctr;
  if (env->sim->generate_path)
  {
    for (ctr=0; ctr<env->num_transmitters; ctr++)
    {
      Transmitter * tn = &(env->transmitters_array[ctr]);
      _update_spatial_parameters(&(tn->gn.smm), &(tn->gn));
    }
    for (ctr=0; ctr<env->num_receivers; ctr++) 
    {
      Receiver * rcv = &(env->receivers_array[ctr]);
      _update_spatial_parameters(&(rcv->gn.smm), &(rcv->gn));
    }
  }
  else
  {
    ctr = 0;
    int ctr1;
    while(env->_nodepath[ctr] != NULL && env->time<env->sim->total_time)
    {
      GeneralNode * gn = env->node_array[ctr];
      double * db = env->_nodepath[ctr];
      db += 6*env->time;
      if (!gn->smm.is_static)
      {
	for(ctr1=0; ctr1<3; ctr1++)
	{
	  // printf("%lf %lf\n", *db, *(db+1));
	  gn->smm.position[ctr1] = *db;
	  gn->smm.velocity[ctr1] = *(db+1);
	  db += 2;
	}
      }
      ctr++;
    }
  }
  env->time++;
}

void update_virtual_transmitters(Environment * env, int use_only_real)
{
  if(use_only_real)
  {
    // ignore the virtual transmitters in env->transmitters_array
    int ctr1, ctr2;
    env->num_virtual_transmitters = 0;
    for (ctr1=0; ctr1<env->num_transmitters; ctr1++)
    {
      for (ctr2=0; ctr2<env->num_transmitters; ctr2++)
      {
	Transmitter * tx1 = env->transmitters_array + ctr1;
	Transmitter * tx2 = env->transmitters_array + ctr2;
	if((tx1->is_real_transmitter) &&
	   (!(tx2->gn.tm.power_in_dBm<-199)))
	{
	  // ctr1 points at a transmitter
	  // ctr2 points at a scatterer

	  Transmitter tn = env->transmitters_array[env->num_transmitters +
						   env->num_virtual_transmitters];
	  init_transmitter(&tn, env->sim);
	  double dist = distance(&(tx1->gn), &(tx2->gn));
	  tn.gn.smm.distance_from_actual_source += dist;
	  tn.gn.tm.start_time += dist/C;
	  tn.gn.tm.end_time += dist/C;
	  env->num_virtual_transmitters++;
	}
      }
    }
  }
  else
  {
    // TODO: update code for virtual transmitters
  }
}

double distance(GeneralNode * gn1, GeneralNode *gn2)
{
  double * pos1 = gn1->smm.position;
  double * pos2 = gn2->smm.position;
  double res = 0;
  int ctr;
  for(ctr=0; ctr<3; ctr++)
  {
    res += pow(*(pos1+ctr) - *(pos2+ctr), 2);
  }

  return pow(res, 0.5);
}

void readout_receiver_array(Environment * env)
{
  double complex *rx_output;
  rx_output = malloc(env->num_receivers*sizeof(double complex));
  Receiver * rx;
  Transmitter * tx;
  int tot_tx = env->num_transmitters + env->num_virtual_transmitters;
  int ctr=0;
  for (rx=env->receivers_array; rx<env->receivers_array + env->num_receivers; rx++)
  {
    rx_output[ctr] = 0;
    for (tx=env->transmitters_array; tx<env->transmitters_array + tot_tx; tx++)
    {
      double dist = distance(&(tx->gn), &(rx->gn));
      PropagationModel pm;
      init_propagation_model(&pm, env->sim);
      pm.distance = dist;
      double pow_attn, phase_shift;
      compute_shift(&pm, &pow_attn, &phase_shift, env);
      double recv_power_dBm = tx->gn.tm.power_in_dBm - pow_attn;
      rx_output[ctr] += pow(10, 0.5*recv_power_dBm/10)*cexp(I*phase_shift);
    }
    rx_output[ctr] += pow(rx->recv_noise_power/2, 0.5)*(gaussenv(env) + I*gaussenv(env));
    if (env->sim->print_measurements)
    {
      printf("Receiver %d %d %10.6e %10.4e \n", rx->gn.id, env->time, creal(rx_output[ctr]), cimag(rx_output[ctr]));
    }
    ctr++;
  }
  free(rx_output);
}

void compute_shift(PropagationModel * pm, double * power_attn, double * phase_shift, Environment * env)
{
  double dist = pm->distance/env->sim->wavelength;
  
  * power_attn = 20*log10(4*PI*dist);
  * phase_shift = 2*PI*dist;
}

void populate_from_file(Environment * env, const char * filename)
{
  FILE * fp;
  char mode[] = "r";
  char fmt[] = "%s";
  char buff[LINE_LENGTH];
  fp = fopen(filename, mode);
  while (fscanf(fp, fmt, buff) != EOF)
  {
    handle_request(env, fp, buff);
  }
  fclose(fp);
}

void handle_request(Environment * env, FILE * fp, const char * req_type)
{
  int ctr;

  if (!strcmp(req_type, "Timedelta"))
  {
    fscanf(fp, "%lf", &env->sim->delta_time);
  }
  else if (!strcmp(req_type, "Totaltime"))
  {
    fscanf(fp, "%d", &env->sim->total_time);
  }
  else if (!strcmp(req_type, "Transmitter"))
  {
    Transmitter * tx = &env->transmitters_array[env->num_transmitters];
    init_transmitter(tx, env->sim);
    fscanf(fp, "%d", &tx->gn.id);
    fscanf(fp, "%lf", &tx->gn.tm.power_in_dBm);
    fscanf(fp, "%lf", &tx->gn.tm.start_time);
    fscanf(fp, "%lf", &tx->gn.tm.end_time);
    add_node_to_environment_array(env, &tx->gn);
    env->num_transmitters++;
  }
  else if (!strcmp(req_type, "Receiver"))
  {
    Receiver * rx = &env->receivers_array[env->num_receivers];
    init_receiver(rx, env->sim);
    fscanf(fp, "%d", &rx->gn.id);
    double np;
    fscanf(fp, "%lf", &np);
    rx->recv_noise_power = pow(10, np/10);
    add_node_to_environment_array(env, &rx->gn);
    env->num_receivers++;
  }
  else if (!strcmp(req_type, "Scatterer"))
  {
    Transmitter * sc = &env->transmitters_array[env->num_transmitters];
    init_transmitter(sc, env->sim);
    fscanf(fp, "%d", &sc->gn.id);
    add_node_to_environment_array(env, &sc->gn);
    env->num_transmitters++;
  }
  else if (!strcmp(req_type, "Nodepath"))
  {
    int id, ctr;
    char ch;
    fscanf(fp, "%d", &id);
    fscanf(fp, " %c", &ch);
    // printf("ID and CH %d, %c\n", id, ch);
    if (ch=='s')
    {
      // node is static
      env->node_array[id]->smm.is_static = true;
      for (ctr=0; ctr<3; ctr++)
      {
	fscanf(fp, "%lf", &env->node_array[id]->smm.position[ctr]);
      }      
    }
    else
    {
      int num_inst;
      fscanf(fp, "%d", &num_inst);
      //printf("Number of numbers %d \n", num_inst);
      // node is dynamic
      if (env->_nodepath[id] == NULL)
      {
	env->_nodepath[id] = malloc(6*MAX_TIME*sizeof(double));
      }
      env->node_array[id]->smm.is_static = false;

      for (ctr=6*env->total_times[id]; ctr<6*(env->total_times[id]+num_inst); ctr++)
      {
	fscanf(fp, "%lf", &env->_nodepath[id][ctr]);	
	//printf("Read in %lf into %d, %d\n", env->_nodepath[id][ctr], id, ctr);
      }
      env->total_times[id] += num_inst;
    }
  }
  else if (!strcmp(req_type, "Gaussianrand"))
  {
    int num_elem, ctr;
    fscanf(fp, "%d", &num_elem);
   
    for (ctr=env->total_gaussian_samples; ctr<num_elem+env->total_gaussian_samples; ctr++)
    {
      fscanf(fp, "%lf", &env->gaussiansamples[ctr]);
    }
    env->total_gaussian_samples += num_elem;
  }
}

void print_to_file(Environment * env, const char * filename, bool print_state)
{
  FILE * fp;
  fp = fopen(filename, "a");

  if (fp != NULL)
  {
    if (print_state)
    {
    
      fprintf(fp, "Timedelta %lf\n", env->sim->delta_time);
      fprintf(fp, "Totaltime %d\n", env->sim->total_time);
      Transmitter * tx;
      for (tx=env->transmitters_array; tx!= env->transmitters_array
  	   + env->num_transmitters; tx++)
      {
        if (tx->gn.tm.power_in_dBm > -190)
        {
  	fprintf(fp, "Transmitter %d %lf %lf %lf\n",
  		tx->gn.id, tx->gn.tm.power_in_dBm,
  		tx->gn.tm.start_time,
  		tx->gn.tm.end_time);
        }
        else
        {
  	fprintf(fp, "Scatterer %d\n", tx->gn.id);
        }
      }
      Receiver * rx;
      for (rx=env->receivers_array; rx!=env->receivers_array
  	   + env->num_receivers; rx++)
      {
        fprintf(fp, "Receiver %d %lf\n",
		rx->gn.id, 10*log10(rx->recv_noise_power));
      }

      //print static node positions
      int ctr=0;
      GeneralNode * gn = env->node_array[ctr];
      while (gn != NULL)
      {
        if (gn->smm.is_static)
        {
	  fprintf(fp, "Nodepath %d s ", gn->id);
	  int ctr1;
	  for (ctr1=0; ctr1<3; ctr1++)
	  {
	    fprintf(fp, "%lf ", gn->smm.position[ctr1]);
	  }
	  fprintf(fp, "\n");
        }
	gn = env->node_array[++ctr];
      }

      int num_gauss = 5*MAX_TIME;
      num_gauss = 1;
      fprintf(fp, "Gaussianrand %d ", num_gauss);
      for (ctr=0; ctr<num_gauss; ctr++)
      {
        fprintf(fp, "%lf ", _gaussrand());
      }
      fprintf(fp, "\n");
    }
    else
    {
      int ctr = 0;
      GeneralNode * gn = env->node_array[ctr];
      while (gn != NULL)
      {
        if (!gn->smm.is_static)
        {
	  fprintf(fp, "Nodepath %d d 1 ", gn->id);
	  int ctr1;
	  for (ctr1=0; ctr1<3; ctr1++)
	  {
	    fprintf(fp, "%lf ", gn->smm.position[ctr1]);
	    fprintf(fp, "%lf ", gn->smm.velocity[ctr1]);	  	  
	  }
	  fprintf(fp, "\n");
        }
	gn = env->node_array[++ctr];
      }
    }
    fclose(fp);
  }
}

void add_node_to_environment_array(Environment * env, GeneralNode * gn)
{
  env->node_array[gn->id] = gn;
}

void print_node_path(Environment * env)
{
  int ctr, ctr1;
  ctr1=0;
  double * db = env->_nodepath[ctr1];
  while (db != NULL)
  {
    for (ctr=0; ctr<60; ctr++)
    {
      printf("%lf ", *(db+ctr));
    }
    printf("\n");
    db = env->_nodepath[++ctr1];
  }
  printf("Total number of nodes is %d \n", ctr1);
}

void print_current_locations(Environment * env)
{
  int ctr=0, ctr1;
  GeneralNode * gn = env->node_array[ctr];
  while(gn != NULL)
  {
    printf("Node %d \n", ctr);
    for(ctr1=0; ctr1<3; ctr1++)
    {
      printf("%lf ", gn->smm.position[ctr1]);
    }
    printf("\n");
    gn = env->node_array[++ctr];
  }
}

double gaussenv(Environment * env)
{
  static int numcalls=0;
  return env->gaussiansamples[numcalls++];
}
