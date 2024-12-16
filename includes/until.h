//
// Created by ra on 03.12.24.
//

#ifndef UNTIL_H
#define UNTIL_H
#include "vec.h"

//
// Created by ra on 03.12.24.
//



void runge_kutta_step(vec *r, vec *v, vec B, vec E, long double q, long double m, long double dt);

vec getE(array_array_double *phi_data, long double y);

bool random_chance(int percent);

void save_data_to_file(char* fname, vec* result, size_t size);

#endif //UNTIL_H
