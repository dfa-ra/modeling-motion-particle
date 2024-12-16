#include <float.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#include "includes/arr.h"
#include "includes/vec.h"
#include "includes/phisics.h"
#include "includes/until.h"


void save_data_to_file(char* fname, vec* result, size_t size);


// Основная программа
int main() {

    int flag = 1;
    // Исходные данные
    long double Ek = 5 * e;         // Кинетическая энергия (Дж)
    long double alpha = M_PI / 4;  // Угол вылета (радианы)
    vec v = {                      // Начальная скорость
        sqrtl(2 * Ek / me) * cosl(alpha),
        sqrtl(2 * Ek / me) * sinl(alpha),
        0
    };
    vec r = {0, 0, 0};             // Начальная позиция
    vec B = (vec){0, 0, 0};            // Магнитное поле
    vec E = (vec){0, 0, 0};

    long double dt = 1e-14;        // Шаг времени (с)
    long double T = 1e-10;         // Общее время (с)
    size_t steps = T / dt;         // Количество шагов
    size_t real_steps = 0;
    vec result[steps];


    FILE* file = fopen("/media/ra/_work/ra/ITMO/PHISICS/projectC/sourse/data.txt", "r");

    array_array_double phi_data = read_array_array_double_from_file(file);

    if (phi_data.data == NULL ) {
        return 1;
    }

    // Численный метод Эйлера
    for (size_t i = 0; i < steps; i++) {
        if (i == 230 && flag == 1) {
            flag = 0;
            B = (vec){0, 0, 6};
        }
        else if (flag == 0) {
            E = vec_mul(getE(&phi_data, r.y * 100), 100);
        }

        runge_kutta_step(&r, &v, B, E, e, me, dt);

        result[i] = r;
        real_steps++;
        if (r.y <= 0) {
            result[i] = r;
            real_steps++;
            flag = -1;
            B = (vec){0, 0, 0};
            E = (vec){0, 0, 0};
            break;
        };
    }

    save_data_to_file("/media/ra/_work/ra/ITMO/PHISICS/projectC/drawer/data.txt", result, real_steps - 1);
    free_array_array_double(phi_data);
    return 0;
}
