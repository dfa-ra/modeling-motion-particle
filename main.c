#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "includes/vec.h"

#define e 1.602e-19   // Заряд электрона (Кл)
#define me 9.109e-31  // Масса электрона (кг)

void runge_kutta_step(vec *r, vec *v, vec B, long double q, long double m, long double dt) {
    vec k1_v, k2_v, k3_v, k4_v;
    vec k1_r, k2_r, k3_r, k4_r;

    vec a = vec_mul(vec_cross(*v, B), q / m);

    k1_v = vec_mul(a, dt);
    k1_r = vec_mul(*v, dt);

    a = vec_mul(vec_cross(vec_sum(*v, vec_mul(k1_v, 0.5)), B), q / m);
    k2_v = vec_mul(a, dt);
    k2_r = vec_mul(vec_sum(*v, vec_mul(k1_v, 0.5)), dt);

    a = vec_mul(vec_cross(vec_sum(*v, vec_mul(k2_v, 0.5)), B), q / m);
    k3_v = vec_mul(a, dt);
    k3_r = vec_mul(vec_sum(*v, vec_mul(k2_v, 0.5)), dt);

    a = vec_mul(vec_cross(vec_sum(*v, k3_v), B), q / m);
    k4_v = vec_mul(a, dt);
    k4_r = vec_mul(vec_sum(*v, k3_v), dt);

    *v = vec_sum(*v, vec_mul(vec_sum(vec_sum(k1_v, vec_mul(k2_v, 2.0)), vec_sum(vec_mul(k3_v, 2.0), k4_v)), 1.0 / 6.0));
    *r = vec_sum(*r, vec_mul(vec_sum(vec_sum(k1_r, vec_mul(k2_r, 2.0)), vec_sum(vec_mul(k3_r, 2.0), k4_r)), 1.0 / 6.0));
}



// Основная программа
int main() {
    // Исходные данные
    long double Ek = 5 * e;         // Кинетическая энергия (Дж)
    long double alpha = M_PI / 4;  // Угол вылета (радианы)
    vec v = {                      // Начальная скорость
        sqrtl(2 * Ek / me) * cosl(alpha),
        sqrtl(2 * Ek / me) * sinl(alpha),
        0
    };
    vec r = {0, 0, 0};             // Начальная позиция
    vec B = {0, 0, 6};             // Магнитное поле
    long double dt = 1e-13;        // Шаг времени (с)
    long double T = 1e-10;         // Общее время (с)
    size_t steps = T / dt;         // Количество шагов

    // Открытие файла для записи
    FILE* file = fopen("/media/ra/_work/ra/ITMO/PHISICS/projectC/drawer/data.txt", "w");
    if (file == NULL) {
        perror("Ошибка открытия файла");
        return EXIT_FAILURE;
    }

    // Численный метод Эйлера
    for (size_t i = 0; i < steps; i++) {
        runge_kutta_step(&r, &v, B, e, me, dt);
        fprintf(file, "%.13Le %.13Le %.13Le\n", r.x, r.y, r.z);
    }


    fclose(file);
    printf("Траектория записана в файл trajectory.txt\n");
    return 0;
}
