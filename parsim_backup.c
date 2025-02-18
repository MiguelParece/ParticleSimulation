#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1


typedef struct {
    double x, y;    // Position
    double vx, vy;  // Velocity
    double m;       // Mass
} particle_t;

unsigned int seed;
void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

double rnd_uniform01()
{
    int seed_in = seed;
    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);
    return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}

double rnd_normal01()
{
    double u1, u2, z, result;
    do
    {
        u1 = rnd_uniform01();
        u2 = rnd_uniform01();
        z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
    } while (result < 0 || result >= 1);
    return result;
}

void init_particles(long seed, double side, long ncside, long long n_part, particle_t* par)
{
    double (*rnd01)() = rnd_uniform01;
    long long i;

    if (seed < 0)
    {
        rnd01 = rnd_normal01;
        seed = -seed;
    }

    init_r4uni(seed);

    for (i = 0; i < n_part; i++)
    {
        par[i].x = rnd01() * side;
        par[i].y = rnd01() * side;
        par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
        par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;

        par[i].m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <seed> <side_length> <grid_size> <n_particles> <n_timesteps>\n", argv[0]);
        return 1;
    }

    // Parse command line arguments
    int seed = atol(argv[1]);
    double side_length = atof(argv[2]);
    long grid_size = atol(argv[3]);
    long long n_particles = atoll(argv[4]);
    long n_timesteps = atol(argv[5]);

    // Allocate memory for particles
    particle_t* particles = malloc(n_particles * sizeof(particle_t));
    if (!particles) {
        fprintf(stderr, "Error allocating particles array\n");
        return 1;
    }

    // Initialize particles
    init_particles(seed, side_length, grid_size, n_particles, particles);

    // Measure execution time
    // TODO: Uncomment for the second delivery
    //double exec_time;
    //exec_time = -omp_get_wtime();

    //simulation();  // You'll implement this function

    //exec_time += omp_get_wtime();
    //fprintf(stderr, "%.1fs\n", exec_time);

    //print_result();  // You'll implement this function

    // Clean up
    free(particles);
    return 0;
}