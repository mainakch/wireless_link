#include "raytracing.h"

void *custom_malloc(size_t size)
{
        void *tmp = get_or_free_memory(size, 0, 0);
        return tmp;
}

void *custom_calloc(size_t num, size_t size)
{
        void *tmp = custom_malloc(num * size);
        memset(tmp, 0, num * size);
        return tmp;
}

void custom_free(void *ptr)
{
        get_or_free_memory(0, ptr, 1);
}

unsigned long hash(unsigned char *str)
{
        // http://www.cse.yorku.ca/~oz/hash.html
        unsigned long hash = 5381;
        int c;

        while (c = *str++)
                hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

        return hash;
}

unsigned long hash_ptr(void *ptr)
{
        char str[16];
        sprintf(str, "%p", ptr);
        return hash(str)%4096;
}

unsigned long hash_int(size_t num)
{
        char str[16];
        sprintf(str, "%zx", num);
        return hash(str)%256;
}

struct memory_register *init_memory_register()
{
        struct memory_register *mem_reg = calloc(1,
                                                 sizeof(struct memory_register));
        mem_reg->malloc_array = calloc(256, sizeof(struct memory_block *));
        mem_reg->free_array = calloc(4096, sizeof(struct memory_block *));
        return mem_reg;
}

void destroy_memory_register(struct memory_register *mem_reg)
{
        int ctr;
        for (ctr = 0; ctr < 256; ++ctr) {
                if (*(mem_reg->malloc_array + ctr) != 0) {
                        destroy_memory_block(*(mem_reg->malloc_array + ctr));
                }
        }
        for (ctr = 0; ctr < 4096; ++ctr) {
                if (*(mem_reg->free_array + ctr) != 0) {
                        destroy_memory_block(*(mem_reg->free_array + ctr));
                }
        }
        free(mem_reg->malloc_array);
        free(mem_reg->free_array);
        free(mem_reg);
}

struct memory_block *init_memory_block(size_t size)
{
        struct memory_block *mem_blk = malloc(sizeof(struct memory_block));
        mem_blk->size = size;
        mem_blk->ptr = malloc(size);
        mem_blk->sizehash = (int) hash_int(size);
        mem_blk->ptrhash = (int) hash_ptr(mem_blk->ptr);
        mem_blk->next = 0;
        return mem_blk;
}

void destroy_memory_block(struct memory_block *mem_blk)
{
        if (mem_blk->next != 0) {
                destroy_memory_block(mem_blk->next);
        }
        free(mem_blk->ptr);
        free(mem_blk);
}

void *get_or_free_memory(size_t size, void *ptr, int zero_for_get)
{
        // negative zero_for_get to destroy all
        static struct memory_register *mem_reg = 0;
        void *ptr_n = 0;

        if (mem_reg == 0) {
                mem_reg = init_memory_register();
        }

        // free allocated memory if zero_for_get is positive
        if (zero_for_get > 0) {
                free_memory(ptr, mem_reg);
        }

        if (zero_for_get == 0) {
                ptr_n = malloc_memory(size, mem_reg);
        }

        if (mem_reg != 0 && zero_for_get < 0) {
                destroy_memory_register(mem_reg);
                mem_reg = 0;
        }

        return ptr_n;
}

static void free_memory(void *ptr, struct memory_register *mem_reg)
{
        // find block
        int hash_p = (int) hash_ptr(ptr);
        struct memory_block *mem_blk = *(mem_reg->free_array + hash_p);
        struct memory_block *mem_prv = 0;

        while (mem_blk->ptr != ptr) {
                mem_prv = mem_blk;
                mem_blk = mem_blk->next;
        }

        // remove block from free_array, zero memory
        if (mem_prv != 0) {
                mem_prv->next = mem_blk->next;
        } else {
                *(mem_reg->free_array + hash_p) = mem_blk->next;
                // handles case when next is null also
        }
        mem_blk->next = 0;

        // move to appropriate location in malloc_array
        struct memory_block *m_mem_blk = *(mem_reg->malloc_array +
                                           mem_blk->sizehash);
        // add to the beginning of malloc_array linked list
        if (m_mem_blk != 0) mem_blk->next = m_mem_blk;
        *(mem_reg->malloc_array + mem_blk->sizehash) = mem_blk;
}

static void *malloc_memory(size_t size, struct memory_register *mem_reg)
{
        // find block
        int hash_i = (int) hash_int(size);
        struct memory_block *mem_blk = *(mem_reg->malloc_array + hash_i);
        struct memory_block *mem_prv = 0;

        // find memory block of appropriate size since
        // multiple sizes may map to the same hash
        while (mem_blk != 0 && mem_blk->size != size) {
                mem_prv = mem_blk;
                mem_blk = mem_blk->next;
        }

        if (mem_blk == 0) {
                mem_blk = init_memory_block(size);
        } else {
                // disconnect mem_blk from malloc_array
                if (mem_prv == 0) {
                        *(mem_reg->malloc_array + hash_i) = mem_blk->next;
                } else {
                        mem_prv->next = mem_blk->next;
                }
                mem_blk->next = 0;
        }
        assert(mem_blk->next == 0);

        // add mem_blk to the beginning of appropriate link of free_array
        struct memory_block *mem_blk_f = *(mem_reg->free_array
                                           + mem_blk->ptrhash);
        if (mem_blk_f == 0) {
                *(mem_reg->free_array + mem_blk->ptrhash)
                        = mem_blk;
        } else {
                mem_blk->next = mem_blk_f->next;
                mem_blk_f->next = mem_blk;
        }

        return mem_blk->ptr;
}

void release_all_blocks()
{
        get_or_free_memory(0, 0, -1);
}

int counter()
{
        static int ctr = 0;
        return ctr++;
}

struct receiver_ray_ribbon_ll_node *init_receiver_ray_ribbon_ll_node() {
        struct receiver_ray_ribbon_ll_node *rnll =
                custom_calloc(1, sizeof(struct receiver_ray_ribbon_ll_node));
        return rnll;
}

struct receiver_ray_ribbon *init_receiver_ray_ribbon(
        struct receiver *rx,
        const struct ray_ribbon *ribbon,
        struct environment *env)
{
        struct receiver_ray_ribbon *rrbn_new =
                custom_calloc(1, sizeof(struct receiver_ray_ribbon));
        const struct perfect_reflector **prconst =
                (const struct perfect_reflector **) env->prarray;
        rrbn_new->ribbon = refine_ray_ribbon_image(ribbon->start_tx,
                                                   ribbon,
                                                   rx,
                                                   prconst);
        rrbn_new->start_tx = ribbon->start_tx;

        if (rrbn_new->ribbon == 0) {
                destroy_receiver_ray_ribbon(rrbn_new);
                return 0;
        }
        return rrbn_new;
}

bool populate_if_ray_ribbon_doesnt_exist
(const struct ray_ribbon *rb1, struct receiver *rx, struct environment *env) {
        // if type is the same declare to be the same
        struct receiver_ray_ribbon_ll_node *rnll = rx->rlln;
        struct receiver_ray_ribbon_ll_node *rnprev = 0;

        while (rnll != 0) {
                assert(rnll->rrbn != 0);
                assert(rnll->rrbn->ribbon != 0);
                if (check_same_type(rb1, rnll->rrbn->ribbon)) return false;
                rnprev = rnll;
                rnll = rnll->next;
        }

        struct receiver_ray_ribbon *rrbn_new =
                init_receiver_ray_ribbon(rx, rb1, env);
        if (!rrbn_new) return false;
        // add rb1 to rx
        if (rnprev == 0) {
                rx->rlln = init_receiver_ray_ribbon_ll_node();
                rx->rlln->rrbn = rrbn_new;
        } else {
                rnprev->next = init_receiver_ray_ribbon_ll_node();
                rnprev->next->rrbn = rrbn_new;
        }

        return true;
}

void populate_receiver_ray_ribbons(struct environment *env)
{
        if (env->tx_paths_updated_rx_paths_updated) {
                return;
        }
        fprintf(stderr, "Updating receiver ribbons\n");

        for (int ctrrx = 0; ctrrx < env->num_receivers; ++ctrrx) {
                struct receiver *rx = *(env->receivers_array + ctrrx);
                int ctrprx = 0;
                int ctrtx = 0;
                struct ray_ribbon_array *rbnarr = *(env->tx_paths + ctrtx);
                while (rbnarr != 0) {
                        int ctrp = 0;
                        struct ray_ribbon *rbn = *(rbnarr->ribbons);
                        while (rbn != 0) {
                                rbn->start_tx = *(env->transmitters_array
                                                  + ctrtx);
                                bool added = populate_if_ray_ribbon_doesnt_exist(
                                        rbn, rx, env);
                                if (added) fprintf(stderr, "Added ribbon\n");
                                ++ctrp;
                                rbn = *(rbnarr->ribbons + ctrp);
                        }

                        ++ctrtx;
                        rbnarr = *(env->tx_paths + ctrtx);
                }
        }
        env->tx_paths_updated_rx_paths_updated = true;
}

void update_all_receiver_ray_ribbons(struct environment *env)
{
        for (int ctrrx = 0; ctrrx < env->num_receivers; ++ctrrx) {
                struct receiver *rx = *(env->receivers_array + ctrrx);
                update_receiver_ray_ribbons(rx, env);
                update_receiver_ray_ribbons_signal_buffer(rx, env);
        }
}

bool update_receiver_ray_ribbons(struct receiver *rx,
                                struct environment *env)
{
        // this function updates the ray ribbon spatial parameters
        // based on the spatial locations
        // of the transmitter and the receiver, returns true

        struct receiver_ray_ribbon_ll_node *rlln = rx->rlln;
        const struct perfect_reflector **prconst =
                (const struct perfect_reflector **) env->prarray;
        while (rlln != 0) {
                struct ray_ribbon *refined_ribbon = 0;
                if (rlln->rrbn->ribbon) {
                        refined_ribbon = refine_ray_ribbon_image(
                                rlln->rrbn->ribbon->start_tx,
                                rlln->rrbn->ribbon,
                                rx,
                                prconst);
                        destroy_ray_ribbon(rlln->rrbn->ribbon);
                }
                if (refined_ribbon != 0) {
                        refined_ribbon->start_tx = rlln->rrbn->ribbon->start_tx;
                        refined_ribbon->end_rx = rx;
                        rlln->rrbn->ribbon = refined_ribbon;
                        update_receiver_ribbon_delay_dopplers(rlln->rrbn, env);
                } else {
                        // tx path update is triggered whenever ray is lost
                        // env->tx_paths_updated = false;
                        rlln->rrbn->ribbon = 0;
                }

                rlln = rlln->next;
        }
        return true;
}

void update_receiver_ray_ribbons_signal_buffer(struct receiver *rx,
                                              struct environment *env)
{
        // this function updates the buffer just before readout and removes
        // all outdated signals

        // code to update buffers
        struct receiver_ray_ribbon_ll_node *rlln = rx->rlln;
        struct receiver_ray_ribbon_ll_node *rllnprev = 0;
        while (rlln != 0) {
                /* if (rlln->rrbn->ribbon == 0 && rlln->rrbn->signal->next == 0 && */
                /*     rlln->rrbn->signal->receiver_read && */
                /*     rlln->rrbn->signal->transmit_time */
                /*        + rlln->rrbn->signal->delay + */
                /*     env->delta_time < env->time) { */
                if (rlln->rrbn->ribbon == 0) {
                        if (rllnprev != 0) {
                                rllnprev->next = rlln->next;
                                rlln->next = 0;
                                destroy_receiver_ray_ribbon_ll_node(rlln);
                                rlln = rllnprev->next;
                        } else {
                                rx->rlln = rlln->next;
                                rlln->next = 0;
                                destroy_receiver_ray_ribbon_ll_node(rlln);
                                rlln = rx->rlln;
                        }

                } else {
                        // add current signal first
                        struct signal_buffer *signal = rlln->rrbn->signal;
                        if (signal == 0) {
                                signal = init_signal_buffer();
                                assert(rlln->rrbn->start_tx != 0);
                                signal->signal = rlln->rrbn->start_tx->baseband_signal;
                                signal->transmit_time = env->time;
                                signal->delay = rlln->rrbn->delay;
                                rlln->rrbn->signal = signal;
                        } else {
                                while(signal->next != 0) {
                                        signal = signal->next;
                                }

                                signal->next = init_signal_buffer();
                                signal->next->delay = rlln->rrbn->delay;
                                signal->next->signal =
                                        rlln->rrbn->start_tx->baseband_signal;
                                signal->next->transmit_time = env->time;
                        }

                        // remove signal if there are newer signals in the buffer
                        while (rlln->rrbn->signal != 0 && rlln->rrbn->signal->next != 0
                               && rlln->rrbn->signal->next->transmit_time
                               + rlln->rrbn->signal->next->delay < env->time) {
                                rlln->rrbn->signal =
                                        destroy_signal_buffer_first(rlln->rrbn->signal);
                        }

                        rllnprev = rlln;
                        rlln = rlln->next;
                }
        }
}

void destroy_receiver_ray_ribbon_ll_node(struct receiver_ray_ribbon_ll_node *p)
{
        destroy_receiver_ray_ribbon(p->rrbn);
        custom_free(p);
}

void destroy_receiver_ray_ribbon(struct receiver_ray_ribbon *rrbn)
{
        destroy_signal_buffer(rrbn->signal);
        if (rrbn->ribbon != 0) {
                destroy_ray_ribbon(rrbn->ribbon);
        }
        custom_free(rrbn);
}

struct signal_buffer *init_signal_buffer()
{
        struct signal_buffer *sgn =
                custom_calloc(1, sizeof(struct signal_buffer));
        return sgn;
}

void destroy_signal_buffer(struct signal_buffer *sgn)
{
        if (sgn == 0) return;
        if (sgn->next != 0) {
                destroy_signal_buffer(sgn->next);
        }
        custom_free(sgn);
}

struct signal_buffer *destroy_signal_buffer_first(struct signal_buffer *sgn)
{
        struct signal_buffer *sgnnext = sgn->next;
        sgn->next = 0;
        destroy_signal_buffer(sgn);
        return sgnnext;
}

struct ray_ribbon_array *init_ray_ribbon_array(int number)
{
        struct ray_ribbon_array *rarr =
                custom_calloc(1, sizeof(struct ray_ribbon_array));
        rarr->ribbons = custom_calloc(number, sizeof(struct ray_ribbon *));
        rarr->max_len = number;
        return rarr;
}

void populate_ray_ribbon_array(struct transmitter *tx,
                               const struct perfect_reflector **ref_arr,
                               struct ray_ribbon_array *rarr,
                               int num_divs, int num_ref, bool single_type)
{
        populate_ray_ribbon_array_long(tx, ref_arr, rarr, num_ref,
                                       -PI / 2, PI / 2, PI / num_divs,
                                       0, 2 * PI, 2 * PI / num_divs,
                                       single_type);
}

void populate_ray_ribbon_array_long(struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    struct ray_ribbon_array *rarr,
                                    int num_ref,
                                    const double phi_start,
                                    const double phi_end,
                                    const double phi_delta,
                                    const double thet_start,
                                    const double thet_end,
                                    const double thet_delt,
                                    bool single_type)
{
        int num_points = (1 + floor((phi_end - phi_start) / phi_delta)) *
                (1 + floor((thet_end - thet_start) / thet_delt));
        if (_RAYTRACING_DEBUG) fprintf(stderr, "Number is %d\n", num_points);
        double complex *angles = custom_malloc(num_points * sizeof(double complex));

        int ctr = 0;

        double phi, theta;
        for (phi = phi_start; phi < phi_end; phi += phi_delta) {
                for (theta = thet_start; theta < thet_end; theta += thet_delt) {
                        *(angles + ctr) = phi + I * theta;
                        ++ctr;
                }
        }
        populate_ray_ribbon_array_full_copy(tx, ref_arr, num_ref, ctr,
                                       angles, rarr, single_type);
        custom_free(angles);
}

void populate_ray_ribbon_array_full_copy(const struct transmitter *tx,
                                         const struct perfect_reflector **ref_arr,
                                         int num_ref, int num_points,
                                         const double complex *angles,
                                         struct ray_ribbon_array *rarr,
                                         bool single_type)
{
        int ctr = 0;
        double phi, theta;

        struct ribbon_node *rn = init_chain_of_ribbon_nodes(6);
        struct ray_ribbon *rb = init_ray_ribbon(rn);
        rb->start_tx = tx;

        for (ctr = 0; ctr < num_points; ++ctr) {
                phi = creal(*(angles + ctr));
                theta = cimag(*(angles + ctr));
                cblas_dcopy(3, tx->gn->smm->position, 1,
                            rn->current->point, 1);
                double direction[3] = {cos(phi),
                                       sin(phi) * cos(theta),
                                       sin(phi) * sin(theta)};
                cblas_dcopy(3, direction, 1,
                            rn->current->unit_direction, 1);

                bool hit_des = process_vertical_chain_nomalloc(rn,
                                                               ref_arr,
                                                               num_ref);
                if (hit_des) {
                        add_ray_ribbon_copy(rarr, rb, single_type);
                }
        }
        destroy_ray_ribbon(rb);
}

struct ray_ribbon *init_ray_ribbon(struct ribbon_node *rn)
{
        struct ray_ribbon *rb = custom_calloc(1, sizeof(struct ray_ribbon));
        rb->head = rn;
        return rb;
}

struct ribbon_node *init_ribbon_node()
{
        struct ribbon_node *rn = custom_calloc(1, sizeof(struct ribbon_node));
        rn->current = custom_calloc(1, sizeof(struct half_infinite_ray));
        rn->ctr = counter();
        rn->surface_index = -1;
        return rn;
}

struct ribbon_node *init_ribbon_node_from_copy(const struct ribbon_node *rn)
{
        struct ribbon_node *rn_tmp = init_ribbon_node();
        rn_tmp->hit_destination_patch = rn->hit_destination_patch;
        rn_tmp->num_reflections = rn->num_reflections;
        rn_tmp->ctr = rn->ctr;
        rn_tmp->surface_index = rn->surface_index;
        *(rn_tmp->current) = *(rn->current);
        return rn_tmp;
}

struct ribbon_node *init_ribbon_node_from_points(const double *pt1,
                                                 const double *pt2)
{
        struct ribbon_node *rn_tmp = init_ribbon_node();
        cblas_dcopy(3, pt1, 1, rn_tmp->current->point, 1);
        cblas_dcopy(3, pt2, 1, rn_tmp->current->end_pt, 1);

        double diff_v[3];
        diff(pt1, pt2, diff_v);
        rn_tmp->current->length = normalize_unit_vector(diff_v);
        cblas_dcopy(3, diff_v, 1, rn_tmp->current->unit_direction, 1);
        return rn_tmp;
}

struct ribbon_node *init_chain_of_ribbon_nodes(int length_of_node)
{
        int ctr = 0;

        struct ribbon_node *rn = 0;
        struct ribbon_node *rn_tmp = 0;
        struct ribbon_node *rn_prev = 0;

        while (ctr < length_of_node) {
                rn_tmp = init_ribbon_node();
                rn_tmp->num_reflections = ctr;

                if (rn == 0) {
                        rn = rn_tmp;
                } else {
                        rn_prev->down = rn_tmp;
                }

                rn_prev = rn_tmp;
                ++ctr;
        }

        return rn;
}

struct ribbon_node *copy_chain_of_ribbon_nodes(const struct ribbon_node *rn)
{
        struct ribbon_node *rn_new = 0;
        struct ribbon_node *rn_tmp = 0;
        struct ribbon_node *rn_tmp_prev = 0;
        while (rn != 0) {
                rn_tmp = init_ribbon_node();
                rn_tmp->hit_destination_patch = rn->hit_destination_patch;
                rn_tmp->num_reflections = rn->num_reflections;
                rn_tmp->ctr = rn->ctr;
                rn_tmp->surface_index = rn->surface_index;
                *(rn_tmp->current) = *(rn->current);

                if (rn_new == 0) {
                        rn_new = rn_tmp;
                } else {
                        rn_tmp_prev->down = rn_tmp;
                }
                rn = rn->down;
                rn_tmp_prev = rn_tmp;
        }

        return rn_new;
}

struct ribbon_node *copy_chain_of_ribbon_nodes_till_dest(const struct ribbon_node *rn)
{
        struct ribbon_node *rn_new = 0;
        struct ribbon_node *rn_tmp = 0;
        struct ribbon_node *rn_tmp_prev = 0;
        bool isfinal = false;
        while (rn != 0 && !isfinal) {
                rn_tmp = init_ribbon_node();
                rn_tmp->hit_destination_patch = rn->hit_destination_patch;
                rn_tmp->num_reflections = rn->num_reflections;
                rn_tmp->ctr = rn->ctr;
                rn_tmp->surface_index = rn->surface_index;
                *(rn_tmp->current) = *(rn->current);

                if (rn_new == 0) {
                        rn_new = rn_tmp;
                } else {
                        rn_tmp_prev->down = rn_tmp;
                }
                isfinal = rn->hit_destination_patch;
                rn = rn->down;
                rn_tmp_prev = rn_tmp;
        }

        return rn_new;
}

struct ray_ribbon *copy_ray_ribbon(const struct ray_ribbon *rb, bool till_dest)
{
        struct ray_ribbon *rbnew = custom_malloc(sizeof(struct ray_ribbon));
        if (till_dest) {
                rbnew->head = copy_chain_of_ribbon_nodes_till_dest(rb->head);
        } else {
                rbnew->head = copy_chain_of_ribbon_nodes(rb->head);
        }
        rbnew->start_tx = rb->start_tx;
        rbnew->end_rx = rb->end_rx;
        rbnew->delay = rb->delay;
        rbnew->doppler = rb->doppler;
        rbnew->integrated_doppler_phase = rb->integrated_doppler_phase;
        rbnew->gain = rb->gain;
        rbnew->reflection_phase = rb->reflection_phase;
        return rbnew;
}

void destroy_ray_ribbon(struct ray_ribbon *rb)
{
        if (rb == 0) return;
        destroy_ray_ribbon_nodes(rb);
        custom_free(rb);
}

void destroy_ray_ribbon_nodes(struct ray_ribbon *rb)
{
        if (rb == 0) return;
        struct ribbon_node *rn = rb->head;
        destroy_chain_of_ribbon_nodes(rn);
        rb->head = 0;
}

void destroy_ray_ribbon_array(struct ray_ribbon_array *array)
{
        if (array == 0) return;
        int ctr = 0;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon_nodes(*(array->ribbons + ctr));
                custom_free(*(array->ribbons + ctr));
                ++ctr;
        }
        custom_free(array->ribbons);
        custom_free(array);
}

void destroy_ray_ribbon_array_all_but_first(struct ray_ribbon_array *array)
{
        if (array == 0) return;
        int ctr = 1;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon_nodes(*(array->ribbons + ctr));
                custom_free(*(array->ribbons + ctr));
                ++ctr;
        }
        custom_free(array->ribbons);
        custom_free(array);
}

void destroy_ray_ribbon_array_ribbons(struct ray_ribbon_array *array)
{
        int ctr = 0;
        if (array == 0) return;
        while (*(array->ribbons + ctr) != NULL) {
                destroy_ray_ribbon(*(array->ribbons + ctr));
                ++ctr;
        }
        custom_free(array->ribbons);
        array->max_len = 0;
        array->current_len = 0;
        array->ribbons = 0;
}

void destroy_chain_of_ribbon_nodes(struct ribbon_node *rn)
{
        if (rn == NULL) return;
        if (rn->current != NULL) custom_free(rn->current);
        struct ribbon_node *rndown = rn->down;
        custom_free(rn);
        destroy_chain_of_ribbon_nodes(rndown);
}

void destroy_ribbon_node(struct ribbon_node *rn) {
        if (rn == NULL) return;
        if (rn->current != NULL) custom_free(rn->current);
        custom_free(rn);
}

bool check_same_type(const struct ray_ribbon *ray_rb1,
                     const struct ray_ribbon *ray_rb2)
{
        struct ribbon_node *rn1 = ray_rb1->head;
        struct ribbon_node *rn2 = ray_rb2->head;
        // same type if same transmitter and same reflectors
        if (ray_rb1->start_tx != ray_rb2->start_tx) return false;

        bool hit1, hit2;

        while(rn1 != NULL && rn2 != NULL) {
                hit1 = rn1->hit_destination_patch;
                hit2 = rn2->hit_destination_patch;

                if (rn1->surface_index != rn2->surface_index) return false;
                rn1 = rn1->down;
                rn2 = rn2->down;
        }
        assert(hit1 || hit2);

        return (hit1 == hit2);
        //return (rn1 == NULL && rn2 == NULL);
}

bool add_ray_ribbon(struct ray_ribbon_array *array, struct ray_ribbon *rb,
                    bool single_type)
{
        if (rb == NULL) return false;

        if (single_type) {
                int ctr = 0;
                // check for type here
                while (*(array->ribbons + ctr) != 0) {
                        if(check_same_type(*(array->ribbons + ctr), rb)) {
                                return false;
                        }
                        ++ctr;
                }
        }

        *(array->ribbons + array->current_len) = rb;
        array->current_len++;
        // this enforces null termination
        *(array->ribbons + array->current_len) = 0;
        return true;
}

bool add_ray_ribbon_copy(struct ray_ribbon_array *array,
                         const struct ray_ribbon *rb,
                         bool single_type)
{
        if (rb == NULL) return false;
        if (single_type) {
                int ctr = 0;
                // check for type here
                while (*(array->ribbons + ctr) != 0) {
                        if(check_same_type(*(array->ribbons + ctr), rb)) {
                                return false;
                        }
                        ++ctr;
                }
        }

        *(array->ribbons + array->current_len) = copy_ray_ribbon(rb, true);
        array->current_len++;
        // this enforces null termination
        *(array->ribbons + array->current_len) = 0;
        return true;
}

double complex compute_intersection(const struct half_infinite_ray *hr,
                                    const struct perfect_reflector *pr)
{
        double t, sgn;
        double diff[3];
        cblas_dcopy(3, hr->point, 1, diff, 1);
        cblas_daxpy(3, -1, pr->center_point, 1, diff, 1);
        t = -cblas_ddot(3, diff, 1, pr->unit_normal, 1)
                /cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
        //if (_RAYTRACING_DEBUG) fprintf(stderr, "t obtained is %lf\n", t);

        // check if t lies within the bounds of the patch

        if (t < INFINITY) {
                // Verify the signs
                cblas_daxpy(3, t, hr->unit_direction, 1, diff, 1);
                double lengtht = cblas_ddot(3, diff, 1,
                                            pr->unit_length_normal, 1);
                double widtht = cblas_ddot(3, diff, 1,
                                           pr->unit_width_normal, 1);

                if (fabs(lengtht) > pr->length / 2 ||
                    fabs(widtht) > pr->width / 2 ||
                    t < 0) {
                        t = INFINITY;
                }
        }

        sgn = cblas_ddot(3, hr->unit_direction, 1, pr->unit_normal, 1);
        return t + I * sgn;
}

bool process_vertical_chain_nomalloc(struct ribbon_node *rn,
                                     const struct perfect_reflector **pr,
                                     int num_reflections)
{
        // this function computes whether a ray can hit the
        // destination after a max num_reflections
        int ctr = 0, ctrindex = -1, num_reflectors = 0;
        double tmin = INFINITY, sgn = -1;

        const struct perfect_reflector *prsurf = *(pr + ctr);
        while(prsurf != NULL) {
                double complex dbl = compute_intersection(rn->current, prsurf);
                rn->current->length = creal(dbl);
                if (creal(dbl) < tmin && ctr != rn->surface_index) {
                        tmin = creal(dbl);
                        sgn = cimag(dbl);
                        ctrindex = ctr;
                }
                ++ctr;
                prsurf = *(pr + ctr);
        }
        num_reflectors = ctr;

        if (sgn>0 || tmin>1e5) return false;
        if (ctrindex == num_reflectors - 1) {
                rn->hit_destination_patch = true;
                return true;
        } else {
                rn->hit_destination_patch = false;
        }

        if (rn->num_reflections > num_reflections) return false;

        // only case remaining is if there is intersection with
        // reflector and number of reflections is small

        // update starting point
        struct ribbon_node *rn_next = rn->down;

        cblas_dcopy(3, rn->current->point, 1,
                    rn_next->current->point, 1);
        cblas_daxpy(3, tmin, rn->current->unit_direction, 1,
                    rn_next->current->point, 1);
        rn_next->surface_index = ctrindex;

        // update ending point of previous ray
        cblas_dcopy(3, rn_next->current->point, 1, rn->current->end_pt, 1);

        // next update direction
        cblas_dcopy(3, rn->current->unit_direction, 1,
                    rn_next->current->unit_direction, 1);
        const struct perfect_reflector *prsurface = pr[ctrindex];
        double factor = -2*cblas_ddot(3, rn->current->unit_direction,
                                      1, prsurface->unit_normal, 1);
        cblas_daxpy(3, factor, prsurface->unit_normal, 1,
                    rn_next->current->unit_direction, 1);

        return process_vertical_chain_nomalloc(rn_next, pr, num_reflections);
}

void print_ray_ribbon(const struct ray_ribbon *rb)
{
        fprintf(stderr, "Printing rayribbon: \n\n");
        struct ribbon_node *rn = rb->head;
        int ctr = 0;
        while(rn != NULL) {
                fprintf(stderr, "Level %d:\n", ctr);
                print_ribbon_node(rn);
                rn = rn->down;
                ++ctr;
        }
        fprintf(stderr, "Delay doppler for rayribbon are: %lf ns; %lf;\n",
                (10e9) * rb->delay, rb->doppler);

}

void print_receiver_ray_ribbon(const struct receiver_ray_ribbon *rb)
{
        fprintf(stderr, "Printing receiver rayribbon: \n\n");
        assert(rb->ribbon);
        struct ribbon_node *rn = rb->ribbon->head;
        int ctr = 0;
        while(rn != NULL) {
                fprintf(stderr, "Level %d:\n", ctr);
                print_ribbon_node(rn);
                rn = rn->down;
                ++ctr;
        }
        fprintf(stderr, "Delay doppler for receiver rayribbon are:",
                "%lf ns; %lf;\n",
                (10e9) * rb->delay, rb->doppler);

}

void print_ray_ribbon_flattened(const struct ray_ribbon *rb)
{
        struct ribbon_node *rn = rb->head;
        int ctr = 0;
        while (rn != NULL) {

                fprintf(stderr, "(%lf, %lf, %lf) -- (%lf, %lf, %lf) ",
                        rn->current->point[0], rn->current->point[1],
                        rn->current->point[2], rn->current->end_pt[0],
                        rn->current->end_pt[1], rn->current->end_pt[2]);
                fprintf(stderr, "Surface index: %d; ", rn->surface_index);
                rn = rn->down;
                ++ctr;
        }
        fprintf(stderr, "\n");
}

void print_ray_ribbon_array(const struct ray_ribbon_array *rarr)
{
        struct ray_ribbon * rb;
        rb = *(rarr->ribbons);
        int ctr = 0;
        while (rb != NULL) {
                fprintf(stderr, "Printing ribbon %d ", ctr);
                print_ray_ribbon_flattened(rb);
                ++ctr;
                rb = *(rarr->ribbons + ctr);
        }
}

void print_vertical_strip(const struct ribbon_node *rn)
{
        int ctr = 0;
        while (rn != NULL) {
                fprintf(stderr, "Level %d\n", ctr++);
                print_ribbon_node(rn);
                rn = rn->down;
        }
}

void print_ribbon_node(const struct ribbon_node *rn)
{
        if (rn == NULL) return;
        fprintf(stderr, "Starting point: ");
        int ctr = 0;
        for(ctr = 0; ctr < 3; ++ctr) {
                fprintf(stderr, "%lf ", rn->current->point[ctr]);
        }

        fprintf(stderr, "Unit direction: ");
        for(ctr = 0; ctr < 3; ++ctr) {
                fprintf(stderr, "%lf ",
                        rn->current->unit_direction[ctr]);
        }

        fprintf(stderr, "Ending point: ");
        for(ctr = 0; ctr < 3; ++ctr) {
                fprintf(stderr, "%lf ", rn->current->end_pt[ctr]);
        }

        fprintf(stderr,
                "Hit dest: %d, num reflec: %d, Surf index: %d",
                rn->hit_destination_patch,
                rn->num_reflections, rn->surface_index);
        fprintf(stderr, "\n");
}

int count_segments(const struct ribbon_node *rn)
{
        int ctr = 0;
        while (rn != NULL) {
                rn = rn->down;
                ++ctr;
        }
        return ctr;
}

void invert_spherical_angles(const double *unit_vector, double *phi,
                             double *thet)
{
        *thet = atan(unit_vector[2]/unit_vector[1]);
        *phi = acos(unit_vector[0]);
        if (sin(*phi) * sin(*thet) / unit_vector[2] < 0 )
                *phi = 2 * PI - (*phi);
}

struct ray_ribbon *refine_ray_ribbon_image(const struct transmitter *tx,
                                           const struct ray_ribbon *rb,
                                           const struct receiver *rx,
                                           const struct perfect_reflector **pr)
{
        if (rb == 0) return 0;
        // first count number of reflectors
        int cnt = count_ribbon_nodes(rb->head);

        // allocate addresses for virtual points
        double **virtual_points = custom_malloc(cnt * sizeof(double *));
        double *zero_pt = custom_calloc(3, sizeof(double));
        struct ribbon_node **rnnodes = custom_calloc(cnt,
                                              sizeof(struct ribbon_node *));
        const struct perfect_reflector **prref = custom_calloc(
                cnt - 1, sizeof(struct perfect_reflector *));

        int ctr = 0;
        *virtual_points = custom_malloc(3 * sizeof(double));
        cblas_dcopy(3, tx->gn->smm->position, 1, *(virtual_points), 1);

        // compute reflected points
        struct ribbon_node *rn = rb->head;

        for (ctr = 0; ctr < cnt - 1; ++ctr) {
                *(rnnodes + ctr) = rn;
                *(virtual_points + ctr + 1) = custom_malloc(3 * sizeof(double));
                cblas_dcopy(3, *(virtual_points + ctr), 1,
                            *(virtual_points + ctr + 1), 1);
                int st_in = rn->down->surface_index;
                const struct perfect_reflector *prr = *(pr + st_in);
                *(prref + ctr) = prr;
                reflect(prr->center_point, prr->unit_normal, zero_pt,
                        *(virtual_points + ctr + 1));
                rn = rn->down;
        }
        *(rnnodes + ctr) = rn;

        // now work backwards to get points of intersection
        double *ptprev = rx->gn->smm->position;

        bool validrayribbon = true;
        for (ctr = cnt - 1; ctr > 0; --ctr) {
                struct ribbon_node *rn = init_ribbon_node_from_points(
                        *(virtual_points + ctr), ptprev);
                rn->surface_index = (*(rnnodes + ctr))->surface_index;
                rn->num_reflections = (*(rnnodes + ctr))->num_reflections;

                double complex tsgn = compute_intersection(
                        rn->current, *(prref + ctr - 1));
                rn->current->length = creal(tsgn);
                cblas_dcopy(3, ptprev, 1, rn->current->end_pt, 1);

                if (creal(tsgn) > 1e5) {
                        validrayribbon = false;
                        ctr = 0;
                } else {
                        cblas_daxpy(3, tsgn, rn->current->unit_direction,
                                    1, rn->current->point, 1);
                        *(rnnodes + ctr) = rn;
                        ptprev = rn->current->point;
                }
        }

        // construct final ray ribbon
        struct ray_ribbon *rbfinal = 0;
        if (validrayribbon) {
                // construct ray ribbon
                struct ribbon_node *rn =
                        init_ribbon_node_from_points(tx->gn->smm->position,
                                                     ptprev); // init_ribbon_node();
                struct ribbon_node *rninit = rn;
                rn->surface_index = -1;
                rn->num_reflections = 0;

                for (ctr = 1; ctr < cnt; ++ctr) {
                        rn->down = *(rnnodes + ctr);
                        rn = rn->down;
                        rn->num_reflections = ctr;
                }
                rn->hit_destination_patch = true;
                rbfinal = init_ray_ribbon(rninit);
                rbfinal->start_tx = tx;
                rbfinal->end_rx = rx;
        }

        // destroy temp
        for (ctr = 0; ctr < cnt; ++ctr) {
                custom_free(*(virtual_points + ctr));
        }
        custom_free(virtual_points);
        custom_free(prref);
        custom_free(rnnodes);
        custom_free(zero_pt);
        return rbfinal;
}

int count_ribbon_nodes(const struct ribbon_node *rn)
{
        int cnt = 0;
        bool isfinal = false;
        while (!isfinal) {
                cnt++;
                isfinal = rn->hit_destination_patch;
                rn = rn->down;
        }
        return cnt;
}

long type_ray_ribbon(const struct ray_ribbon *rb)
{
        int ctr=0;
        struct ribbon_node *rn = rb->head;
        while (rn != NULL) {
                ctr = MAX_SURFACES * ctr + rn->surface_index + 1;
                rn = rn->down;
        }
        return ctr;
}

struct ray_ribbon_array *throw_three_dim_ray_ribbon(struct transmitter *tn,
                                                    const struct
                                                    perfect_reflector **p,
                                                    int num_ref,
                                                    const double phi_start,
                                                    const double phi_end,
                                                    const double phi_incr,
                                                    const double thet_start,
                                                    const double thet_end,
                                                    const double thet_incr)
{
        struct ray_ribbon_array *rarr =
                custom_malloc(sizeof(struct ray_ribbon_array));
        populate_ray_ribbon_array_long(tn, p, rarr,
                                       num_ref, phi_start, phi_end,
                                       phi_incr, thet_start, thet_end,
                                       thet_incr, true);
        return rarr;
}

struct ribbon_node *get_last_ribbon_node(const struct ray_ribbon *rb)
{
        if (rb == 0) {
                fprintf(stderr, "Unexpected error! Should be non null!\n");
        }
        struct ribbon_node *rn = rb->head;
        while (rn->down != NULL) {
                rn = rn->down;
        }
        return rn;
}

void populate_tx_paths(struct environment *env)
{
	if (env->time_index % env->refresh_time == 0) {
                fprintf(stderr, "Refreshing rays at time index: %d\n",
                        env->time_index);
                env->tx_paths_updated = false;
        }

        if (env->tx_paths_updated) return;
        clear_tx_paths(env);
        add_receiver_patch(env, 10);
        const struct perfect_reflector **prconst =
                (const struct perfect_reflector **) env->prarray;

        // now populate individual paths

        struct ray_ribbon_array *rb_arr;
        for (int ctr = 0; ctr < env->num_transmitters; ++ctr) {
                struct transmitter *tx = *(env->transmitters_array + ctr);
                rb_arr = init_ray_ribbon_array(30);
                //populate_ray_ribbon_array(tx, prconst, rb_arr, 600, 3, true);

                populate_ray_ribbon_array_long(tx, prconst,
                                    rb_arr,
                                    3,
                                               -PI,
                                    PI,
                                    0.01,
                                    0,
                                    2 * PI,
                                    0.01,
                                    true);
                *(env->tx_paths + ctr) = rb_arr;
        }

        destroy_last_reflector(env);
        env->tx_paths_updated = true;
        env->tx_paths_updated_rx_paths_updated = false;
}

void update_receiver_ribbon_delay_dopplers(struct receiver_ray_ribbon *rb,
                                  const struct environment *env) {
        if (rb == 0) return;
        // compute delay
        double dist = 0;
        rb->reflection_phase = -1;
        struct ribbon_node *rn = rb->ribbon->head;
        while (rn != NULL) {
                dist += length_ribbon_node(rn);
                rn = rn->down;
                rb->reflection_phase++;
        }
        rb->delay = dist/C;
        // free space path loss
        rb->gain = 1 / (4 * PI * dist / env->wavelength);

        // compute doppler
        rb->doppler = compute_doppler(rb->ribbon, env);
}

double compute_doppler(const struct ray_ribbon *rb,
                       const struct environment *env) {
        struct perfect_reflector *pr;
        struct ribbon_node *rn = rb->head->down;
        double src_pos[3];
        double src_vel[3];

        cblas_dcopy(3, rb->start_tx->gn->smm->position, 1, src_pos, 1);
        cblas_dcopy(3, rb->start_tx->gn->smm->velocity, 1, src_vel, 1);

        while (rn != NULL) {
                pr = *(env->prarray + rn->surface_index);
                reflect(pr->center_point, pr->unit_normal, src_vel, src_pos);
                rn = rn->down;
        }
        double rel_pos[3];
        double rel_vel[3];
        diff(src_pos,
             rb->end_rx->gn->smm->position, rel_pos);
        diff(src_vel,
             rb->end_rx->gn->smm->velocity, rel_vel);
        double tmp = cblas_dnrm2(3, rel_pos, 1);
        return -cblas_ddot(3, rel_vel, 1, rel_pos, 1) / tmp / env->wavelength;
}

double length_ribbon_node(const struct ribbon_node *rn) {
        double diff_vector[3];
        diff(rn->current->point, rn->current->end_pt, diff_vector);
        return cblas_dnrm2(3, diff_vector, 1);
}

void reflect(const double *pos1, const double *n1, double *vel, double *pos) {
        double diff_vec[3];
        double r1[3];
        double v1[3];

        diff(pos1, pos, diff_vec);
        reflection_operation(diff_vec, n1, r1);
        reflection_operation(vel, n1, v1);

        // now update vel and pos
        cblas_dcopy(3, pos1, 1, pos, 1);
        cblas_daxpy(3, 1, r1, 1, pos, 1);
        cblas_dcopy(3, v1, 1, vel, 1);
}

void reflection_operation(const double *v1, const double *n1, double *vref) {
        cblas_dcopy(3, v1, 1, vref, 1);
        double tmp = cblas_ddot(3, n1, 1, vref, 1);
        cblas_daxpy(3, -2 * tmp, n1, 1, vref, 1);
}

void readout_all_signals(struct environment *env, FILE *fpout) {
        static bool first_call = true;
        if (first_call) {
                first_call = false;
                if (fpout != NULL) {
                        fprintf(fpout, "time\treceiver\t"
                                "real\timag\n");
                }
        }
        double complex signal;
        struct ray_ribbon_array *rba;
        int ctr = 0;
        rba = *(env->env_paths + ctr);
        struct receiver *rx = (*(env->receivers_array + ctr));
        while (rba != 0) {
                signal = 0;
                int ctr1 = 0;
                struct ray_ribbon *rb = *(rba->ribbons + ctr1);
                while (rb != 0) {

                        rb->integrated_doppler_phase = fmod(
                                (rb->integrated_doppler_phase +
                                 rb->doppler * env->delta_time), 1);

                        double phase = 0;
                        phase += rb->integrated_doppler_phase
                                + rb->reflection_phase -
                                (env->frequency + rb->doppler) * rb->delay;
                        signal += rb->gain * cexp(2 * PI * phase * I)
                                * pow(10, rb->start_tx->gn->tm->power_in_dBm/10)
                                * rb->start_tx->baseband_signal;
                        ctr1++;
                        rb = *(rba->ribbons + ctr1);
                }

                double rx_noise_std = pow(rx->recv_noise_power, 0.5);
                signal += rx_noise_std *
                        (*(env->unit_power_gaussian_noise + ctr));
                double real_sig = creal(signal);
                double imag_sig = cimag(signal);
                if (fpout != NULL) {
                        fprintf(fpout, "%lf\t%d\t"
                                "%e\t%e\n",
                                env->time, ctr, real_sig, imag_sig);
                }
                ++ctr;
                rba = *(env->env_paths + ctr);
                rx = (*(env->receivers_array + ctr));
        }
}

void readout_all_signals_buffer(struct environment *env, FILE *fpout) {
        static bool first_call = true;
        env->time_index++;
        if (first_call) {
                first_call = false;
                if (fpout != NULL) {
                        fprintf(fpout, "time\treceiver\t"
                                "real\timag\n");
                }
        }

        double complex signal;
        for (int ctr = 0; ctr < env->num_receivers; ++ctr) {
                struct receiver *rx = (*(env->receivers_array + ctr));
                signal = 0;
                int ctr1 = 0;
                struct receiver_ray_ribbon_ll_node *rlln = rx->rlln;
                while (rlln != 0) {
                        struct receiver_ray_ribbon *rrbn = rlln->rrbn;
                        rrbn->integrated_doppler_phase = fmod(
                                (rrbn->integrated_doppler_phase +
                                 rrbn->doppler * env->delta_time), 1);

                        double phase = 0;
                        assert(rrbn->signal);
                        phase += rrbn->integrated_doppler_phase
                                + rrbn->reflection_phase -
                                (env->frequency + rrbn->doppler)
                                * rrbn->signal->delay;

                        if (rrbn->signal->delay +
                            rrbn->signal->transmit_time
                            < env->time) {
                                double txpower =
                                        rrbn->start_tx->gn->tm->power_in_dBm/10;
                                signal += rrbn->gain
                                        * cexp(2 * PI * phase * I)
                                        * pow(10, txpower)
                                        * rrbn->signal->signal;
                                rrbn->signal->receiver_read = true;
                        }

                        rlln = rlln->next;
                }

                double rx_noise_std = pow(rx->recv_noise_power, 0.5);
                signal += rx_noise_std *
                        (*(env->unit_power_gaussian_noise + ctr));
                double real_sig = creal(signal);
                double imag_sig = cimag(signal);
                if (fpout != NULL) {
                        fprintf(fpout, "%lf\t%d"
                                "\t%e\t%e\n",
                                env->time, ctr, real_sig, imag_sig);
                }
        }
}

void clear_tx_paths(struct environment *env) {
        // clear existing paths
        int ctr = 0;
        while (*(env->tx_paths + ctr) != 0) {
                destroy_ray_ribbon_array(*(env->tx_paths + ctr));
                *(env->tx_paths + ctr) = 0;
                ++ctr;
        }
}

void clear_env_paths(struct environment *env) {
        // clear existing paths
        int ctr = 0;
        while (*(env->env_paths + ctr) != 0) {
                destroy_ray_ribbon_array(*(env->env_paths + ctr));
                *(env->env_paths + ctr) = 0;
                ++ctr;
        }
}

// Deprecated functions

// deprecated
void populate_ray_ribbon_array_full_malloc(const struct transmitter *tx,
                                    const struct perfect_reflector **ref_arr,
                                    int num_ref, int num_points,
                                    const double complex *angles,
                                    struct ray_ribbon_array *rarr,
                                    bool single_type)
{
        int ctr = 0;
        double phi, theta;
        for (ctr = 0; ctr < num_points; ++ctr) {
                phi = creal(*(angles + ctr));
                theta = cimag(*(angles + ctr));
                struct ribbon_node *rn = init_ribbon_node();
                cblas_dcopy(3, tx->gn->smm->position, 1,
                            rn->current->point, 1);
                double direction[3] = {cos(phi),
                                       sin(phi) * cos(theta),
                                       sin(phi) * sin(theta)};
                cblas_dcopy(3, direction, 1,
                            rn->current->unit_direction, 1);

                bool hit_des = process_vertical_chain(rn,
                                                      ref_arr,
                                                      num_ref);
                if (hit_des) {
                        struct ray_ribbon *rb = init_ray_ribbon(rn);
                        rb->start_tx = tx;
                        bool ribbon_added = add_ray_ribbon(rarr, rb,
                                                           single_type);
                        if (!ribbon_added) destroy_ray_ribbon(rb);

                } else {
                        destroy_chain_of_ribbon_nodes(rn);
                }
        }
}

// deprecated
bool process_vertical_chain(struct ribbon_node *rn,
                            const struct perfect_reflector **pr,
                            int num_reflections)
{
        // this function computes whether a ray can hit the
        // destination after a max num_reflections
        int ctr=0, ctrindex=-1, num_reflectors=0;
        double tmin = INFINITY, sgn = -1;

        const struct perfect_reflector *prsurf = *(pr + ctr);
        while(prsurf != NULL) {
                double complex dbl = compute_intersection(rn->current, prsurf);
                rn->current->length = creal(dbl);
                if (creal(dbl) < tmin && ctr != rn->surface_index) {
                        tmin = creal(dbl);
                        sgn = cimag(dbl);
                        ctrindex = ctr;
                }
                ++ctr;
                prsurf = *(pr + ctr);
        }
        num_reflectors = ctr;

        if (sgn>0 || tmin>1e5) return false;
        if (ctrindex == num_reflectors - 1) {
                rn->hit_destination_patch = true;
                return true;
        }

        if (rn->num_reflections > num_reflections) return false;

        // only case remaining is if there is intersection with
        // reflector and number of reflections is small

        // update starting point
        struct ribbon_node *rn_next = custom_malloc(sizeof(struct ribbon_node));
        rn_next->down = 0;
        rn_next->current = custom_calloc(1, sizeof(struct half_infinite_ray));
        rn_next->hit_destination_patch = false;
        rn_next->num_reflections = rn->num_reflections + 1;

        cblas_dcopy(3, rn->current->point, 1,
                    rn_next->current->point, 1);
        cblas_daxpy(3, tmin, rn->current->unit_direction, 1,
                    rn_next->current->point, 1);
        rn_next->surface_index = ctrindex;

        // update ending point of previous ray
        cblas_dcopy(3, rn_next->current->point, 1, rn->current->end_pt, 1);

        // next update direction
        cblas_dcopy(3, rn->current->unit_direction, 1,
                    rn_next->current->unit_direction, 1);
        const struct perfect_reflector *prsurface = pr[ctrindex];
        double factor = -2*cblas_ddot(3, rn->current->unit_direction,
                                      1, prsurface->unit_normal, 1);
        cblas_daxpy(3, factor, prsurface->unit_normal, 1,
                    rn_next->current->unit_direction, 1);

        // update pointers
        rn->down = rn_next;
        return process_vertical_chain(rn_next, pr, num_reflections);
}

// deprecated
struct ray_ribbon *refine_ray_image_ribbon(const struct transmitter *tx,
                                     const struct ray_ribbon *rb,
                                     const struct receiver *rx,
                                     const struct perfect_reflector **pr)
{
        struct ray_ribbon_array *node_array_mod = generate_nearby_ribbons(
                tx, pr, 3, rb);

        // return if insufficient number of ribbons
        int len_array = 0;
        while (*(node_array_mod->ribbons + len_array) != 0) {
                len_array++;
        }
        if (len_array < 3) {
                destroy_ray_ribbon_array(node_array_mod);
                fprintf(stderr, "Insufficient ribbons generated\n");
                return 0;
        }

        const double *point = rx->gn->smm->position;
        if (_RAYTRACING_DEBUG) {
                fprintf(stderr, "Position\n");
                print_vector(point);
        }
        int ctr_typ_so_far = 0, ctr = 0;
        double weights[3];
        struct ray_ribbon *rbn = *(node_array_mod->ribbons);
        struct ribbon_node *rn = 0;

        while ((!is_close_ribbon(rbn, point)) && ctr < 100) {
                compute_averaging_coefficients(point, node_array_mod, weights);
                rn = init_ribbon_node();
                compute_average_ribbon_node(rn, node_array_mod, weights);
                destroy_ray_ribbon_nodes(rbn);
                bool has_hit = process_vertical_chain(rn, pr, 3);
                if (!has_hit) {
                        fprintf(stderr, "Unexpected error. Destroying rn.\n");
                        destroy_chain_of_ribbon_nodes(rn);
                        ctr = 101;
                        break;
                }
                (*(node_array_mod->ribbons))->head = rn;
                rbn =  *(node_array_mod->ribbons);
                ++ctr;
        }

        destroy_ray_ribbon_array_all_but_first(node_array_mod);
        if (ctr < 100) {
                rbn->start_tx = tx;
                rbn->end_rx = rx;
                return rbn;
        } else {
                fprintf(stderr, "Ray did not converge!\n");
                destroy_ray_ribbon(rbn);
                return 0;
        }
        // return 0 if ray ribbon does not converge
}

// deprecated
bool is_close_ribbon(const struct ray_ribbon *rb, const double *point)
{
        return isclose(rb->head, point);
}

// deprecated
bool isclose(const struct ribbon_node *rn, const double *point)
{
        double diff[3];
        while (rn->down != NULL) {
                rn = rn->down;
        }
        cblas_dcopy(3, rn->current->end_pt, 1, diff, 1);
        cblas_daxpy(3, -1, point, 1, diff, 1);
        if (cblas_dnrm2(3, diff, 1) < 1e-6)  return true;
        return false;
}

// deprecated
void compute_average_ribbon_node(struct ribbon_node *rn,
                                const struct ray_ribbon_array *rba,
                                double *weights)
{
        double thet, phi;
        double thetaav = 0, phiav = 0;
        int ctr = 0;
        struct ribbon_node *node;
        for (; ctr < 3; ++ctr) {
                node = (*(rba->ribbons + ctr))->head;
                invert_spherical_angles(node->current->unit_direction,
                                        &phi, &thet);
                thetaav += *(weights+ctr) *thet;
                phiav += *(weights+ctr) *phi;
        }
        rn->current->unit_direction[0] = cos(phiav);
        rn->current->unit_direction[1] = sin(phiav)*cos(thetaav);
        rn->current->unit_direction[2] = sin(phiav)*sin(thetaav);
        rn->hit_destination_patch = 0;
        rn->num_reflections = 0;
        rn->ctr = 1;

        cblas_dcopy(3, node->current->point, 1,
                    rn->current->point, 1);
}

// deprecated
void compute_averaging_coefficients(const double *point,
                                    const struct ray_ribbon_array *rba,
                                    double *weights)
{
        gsl_matrix *mat = gsl_matrix_alloc(3, 2);
        gsl_permutation *perm = gsl_permutation_alloc(3);
        gsl_vector *x = gsl_vector_alloc(2);
        gsl_vector *b = gsl_vector_alloc(3);
        gsl_vector *tau = gsl_vector_alloc(2);
        gsl_vector *residual = gsl_vector_alloc(3);

        int c0, c1;
        for (c0=0; c0<3; c0++) {
                /* fprintf(stderr, ANSI_COLOR_GREEN); */
                /* print_ray_ribbon_array(rba); */
                /* fprintf(stderr, ANSI_COLOR_RESET); */
                const struct ribbon_node *node =
                        get_last_ribbon_node(*(rba->ribbons));
                gsl_vector_set(b, c0, *(point + c0)
                               - node->current->end_pt[c0]);
                for (c1=0; c1<2; c1++) {
                        const struct ribbon_node *nodeset =
                                get_last_ribbon_node(*(rba->ribbons + c1 + 1));
                        gsl_matrix_set(mat, c0, c1,
                                       nodeset->current->end_pt[c0]
                                       - node->current->end_pt[c0]);
                }
        }

        gsl_linalg_QR_decomp(mat, tau);
        gsl_linalg_QR_lssolve(mat, tau, b, x, residual);

        *(weights) = 1 - gsl_vector_get(x, 0) - gsl_vector_get(x, 1);
        *(weights + 1) = gsl_vector_get(x, 0);
        *(weights + 2) = gsl_vector_get(x, 1);

        gsl_matrix_free(mat);
        gsl_permutation_free(perm);
        gsl_vector_free(x);
        gsl_vector_free(b);
        gsl_vector_free(tau);
        gsl_vector_free(residual);
}

// deprecated
void remove_ribbon_node_duplicates(struct ribbon_node *rn) {
        struct ribbon_node *rn_orig = rn;
        struct ribbon_node *rn_tmp;
        while (rn != NULL) {
                if ((rn->down != NULL) &&
                    (length_ribbon_node(rn->down) < 1e-4)) {
                        rn_tmp = rn->down;
                        rn->down = rn->down->down;
                        destroy_ribbon_node(rn_tmp);
                }
                rn = rn->down;
        }
        // update num_reflections
        rn = rn_orig;
        int ctr = 0;
        while (rn != NULL) {
                rn->num_reflections = ctr;
                ++ctr;
                rn = rn->down;
        }
}

// deprecated
void update_ribbon_delay_dopplers(struct ray_ribbon *rb,
                                  const struct environment *env) {
        if (rb == 0) return;
        // compute delay
        double dist = 0;
        rb->reflection_phase = -1;
        struct ribbon_node *rn = rb->head;
        while (rn != NULL) {
                dist += length_ribbon_node(rn);
                rn = rn->down;
                rb->reflection_phase++;
        }
        rb->delay = dist/C;
        // free space path loss
        rb->gain = 1 / (4 * PI * dist / env->wavelength);

        // compute doppler
        rb->doppler = compute_doppler(rb, env);
}

// deprecated
void populate_env_paths(struct environment *env)
{
        add_receiver_patch(env, 20);

        // first sound the channel
        const struct perfect_reflector **prconst =
                (const struct perfect_reflector **) env->prarray;

        // now generate rayribbons for each receiver
        // clear existing ray ribbon arrays
        clear_env_paths(env);

        // ctr loops over receivers
        int ctr = 0;
        struct receiver *rx = *(env->receivers_array);
        while (ctr < env->num_receivers) {
                // ctrtx loops over transmitters
                int ctrtx = 0;
                struct ray_ribbon_array *rb_arr =
                        init_ray_ribbon_array(30); // if too low, can cause bugs


                struct ray_ribbon_array *rba = *(env->tx_paths);
                struct transmitter *tx = *(env->transmitters_array);
                while (ctrtx < env->num_transmitters) {
                        // ctr2 loops over rays in tx path ray ribbon from a
                        // part. tx

                        int ctr2 = 0;
                        struct ray_ribbon *rb;
                        rb = *(rba->ribbons + ctr2);
                        while(rb != 0) {
                                struct ray_ribbon *tmprb =
                                        refine_ray_ribbon_image(tx,
                                                                rb, rx, prconst);
                                bool stat;
                                if (tmprb == 0) {
                                        fprintf(stderr, "Null ribbon for ray %d"
                                                " of tx %d at rx %d\n",
                                                ctr2, ctrtx, ctr);
                                } else {
                                        stat = add_ray_ribbon(rb_arr,
                                                              tmprb, true);
                                }
                                if (!stat && tmprb != 0) {
                                        fprintf(stderr, "Unexpected error! "
                                                "stat should always be true!\n");
                                }
                                ctr2++;
                                rb = *(rba->ribbons + ctr2);
                        }
                        ctrtx++;
                        rba = *(env->tx_paths + ctrtx);
                        tx = *(env->transmitters_array + ctrtx);
                }
                *(env->env_paths + ctr) = rb_arr;
                ++ctr;
                *(env->env_paths + ctr) = 0;
                rx = *(env->receivers_array + ctr);
        }
        destroy_last_reflector(env);
}

// deprecated

struct ray_ribbon_array *generate_nearby_ribbons(const struct transmitter *tx,
                                                 const struct
                                                 perfect_reflector **ref_arr,
                                                 int num_ref,
                                                 const struct ray_ribbon *rb)
{
        double phi, theta;
        invert_spherical_angles(rb->head->current->unit_direction,
                                &phi, &theta);
        struct ray_ribbon_array *rarr = init_ray_ribbon_array(4);

        double complex *angles = custom_malloc(3 * sizeof(double complex));
        *angles = (phi - 0.00001) + I * (theta + 0.00002);
        *(angles + 1) = (phi + 0.000039) + I * (theta + 0.000029);
        *(angles + 2) = (phi - 0.000041) + I * (theta - 0.00007);
        /* double fact = (1e-3)/RAND_MAX; */
        /* *angles = (phi - fact * rand()) + I * (theta + fact * rand()); */
        /* *(angles + 1) = phi + fact * rand() + I * (theta + fact * rand()); */
        /* *(angles + 2) = phi + fact * rand() + I * (theta + fact * rand()); */

        populate_ray_ribbon_array_full_copy(tx, ref_arr, num_ref, 3, angles, rarr,
                                       false);
        custom_free(angles);
        return rarr;
}

// deprecated

void update_env_paths_delay_dopplers(struct environment *env) {
        int ctr = 0;
        struct ray_ribbon_array *rba = *(env->env_paths + ctr);
        while (rba != 0) {
                int ctr1 = 0;
                struct ray_ribbon *rb = *(rba->ribbons + ctr1);
                while (rb != NULL) {
                        update_ribbon_delay_dopplers(rb, env);
                        ctr1++;
                        rb = *(rba->ribbons + ctr1);
                }
                ++ctr;
                rba = *(env->env_paths + ctr);
        }
}
