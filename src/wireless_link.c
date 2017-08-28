#include "raytracing.h"

int main(int argc, char *argv[])
{
	struct simulation *sim = init_simulation(argc, argv);

	if (sim == 0) {
	        exit(1);
	}

        update_environment_from_file_sim(sim);
	int ctr = 0;
	while (!(update_environment_from_file_sim(sim))) {
	        fprintf(stderr, "%d\n", ctr);
		populate_tx_paths(sim->env);
		populate_receiver_ray_ribbons(sim->env);
		update_all_receiver_ray_ribbons(sim->env);
	        readout_all_signals_buffer(sim->env);
	        printout_all_signals_buffer(sim->env, sim->fr->outfile);
		printout_path_nariman(sim->env, sim->fr->outfile);
		++ctr;
	}
        destroy_simulation(sim);
}
