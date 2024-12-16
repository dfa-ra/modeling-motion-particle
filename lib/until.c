//
// Created by ra on 03.12.24.
//

#include <float.h>
#include <stdbool.h>
#include <time.h>

#include "../includes/vec.h"
#include "../includes/arr.h"


void runge_kutta_step(vec *r, vec *v, vec B, vec E, long double q, long double m, long double dt) {
    vec k1_v, k2_v, k3_v, k4_v;
    vec k1_r, k2_r, k3_r, k4_r;

    vec a = vec_mul(vec_sum(E, vec_cross(*v, B)), q / m);

    k1_v = vec_mul(a, dt);
    k1_r = vec_mul(*v, dt);

    a = vec_mul(vec_sum(E, vec_cross(vec_sum(*v, vec_mul(k1_v, 0.5)), B)), q / m);
    k2_v = vec_mul(a, dt);
    k2_r = vec_mul(vec_sum(*v, vec_mul(k1_v, 0.5)), dt);

    a = vec_mul(vec_sum(E, vec_cross(vec_sum(*v, vec_mul(k2_v, 0.5)), B)), q / m);
    k3_v = vec_mul(a, dt);
    k3_r = vec_mul(vec_sum(*v, vec_mul(k2_v, 0.5)), dt);

    a = vec_mul(vec_sum(E, vec_cross(vec_sum(*v, k3_v), B)), q / m);
    k4_v = vec_mul(a, dt);
    k4_r = vec_mul(vec_sum(*v, k3_v), dt);

    *v = vec_sum(*v, vec_mul(vec_sum(vec_sum(k1_v, vec_mul(k2_v, 2.0)), vec_sum(vec_mul(k3_v, 2.0), k4_v)), 1.0 / 6.0));
    *r = vec_sum(*r, vec_mul(vec_sum(vec_sum(k1_r, vec_mul(k2_r, 2.0)), vec_sum(vec_mul(k3_r, 2.0), k4_r)), 1.0 / 6.0));
}

#include <math.h>

vec getE(array_array_double *phi_data, long double y) {

    const long double EPSILON = 1e-5;


    for (size_t i = 0; i < phi_data->cols; i++) {
        long double temp = phi_data->data[0][i];
        // Проверка нахождения точки
        if (fabsl(y - temp) < EPSILON) {
            long double delta_phi;
            long double delta_d;

            printf("---------------\n");
            if (i == 0) {
                printf("%Lg, %Lg :1: \n", phi_data->data[1][i], phi_data->data[1][i + 1]);
                printf("%Lg, %Lg :1: \n", phi_data->data[0][i], phi_data->data[0][i + 1]);
                delta_phi = fabsl(phi_data->data[1][i] - phi_data->data[1][i + 1]);
                delta_d = fabsl(phi_data->data[0][i] - phi_data->data[0][i + 1]);
            } else if (i == phi_data->cols - 1) {

                printf("%Lg, %Lg :2: \n", phi_data->data[1][i], phi_data->data[1][i - 1]);
                printf("%Lg, %Lg :2: \n", phi_data->data[0][i], phi_data->data[0][i - 1]);
                delta_phi = fabsl(phi_data->data[1][i] - phi_data->data[1][i - 1]);
                delta_d = fabsl(phi_data->data[0][i] - phi_data->data[0][i - 1]);
            } else {

                printf("%Lg, %Lg :3: \n", phi_data->data[0][i - 1], phi_data->data[0][i + 1]);
                printf("%Lg, %Lg :3: \n", phi_data->data[1][i - 1], phi_data->data[1][i + 1]);
                delta_phi = fabsl(phi_data->data[1][i - 1] - phi_data->data[1][i + 1]);
                delta_d = fabsl(phi_data->data[0][i - 1] - phi_data->data[0][i + 1]);
            }

            if (delta_phi / delta_d > LDBL_MAX) {
                fprintf(stderr, "Переполнение при вычислении E\n");
                exit(1);
            }
            printf("%Lg, %Lg : \n", delta_phi, delta_d);
            return (vec){0, delta_phi / delta_d, 0};
        }
    }
    printf("ZERO");
    return (vec){0, 0, 0};
}


bool random_chance(int percent) {

    if (percent < 0 || percent > 100) {
        fprintf(stderr, "Процент должен быть в диапазоне от 0 до 100\n");
        return false;
    }
    int random_number = rand() % 100; // Случайное число от 0 до 99
    return random_number < percent;
}



//
void save_data_to_file(char* fname, vec* result, size_t real_steps) {

    FILE* file = fopen(fname, "w");

    if (file == NULL) {
        perror("Ошибка открытия файла 0_0");
        return;
    }

    for (size_t i = 0; i < real_steps; i++) {
        fprintf(file, "%.13Le %.13Le %.13Le\n", result[i].x, result[i].y, result[i].z);
    }
    fclose(file);

};