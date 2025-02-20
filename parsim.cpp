#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>
#include <string>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

class RandomGenerator
{
private:
    unsigned int seed;

public:
    RandomGenerator(int input_seed) : seed(input_seed + 987654321) {}

    double uniform01()
    {
        int seed_in = seed;
        seed ^= (seed << 13);
        seed ^= (seed >> 17);
        seed ^= (seed << 5);
        return 0.5 + 0.2328306e-09 * (seed_in + static_cast<int>(seed));
    }

    double normal01()
    {
        double u1, u2, z, result;
        do
        {
            u1 = uniform01();
            u2 = uniform01();
            z = std::sqrt(-2 * std::log(u1)) * std::cos(2 * M_PI * u2);
            result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
        } while (result < 0 || result >= 1);
        return result;
    }
};

class Particle
{

public:
    double x, y;   // Position
    double vx, vy; // Velocity
    double m;      // Mass
    bool alive;    // has Collided

    Particle() : x(0), y(0), vx(0), vy(0), m(0), alive(true) {}

    // Add methods for particle behavior here
    void updatePosition(double dt)
    {
        x += vx * dt;
        y += vy * dt;
    }

    void updatePositionAndVelocity(double dt, double ax, double ay)
    {
        x += vx * dt + 0.5 * ax * dt * dt;
        y += vy * dt + 0.5 * ay * dt * dt;

        vx += ax * dt; // Optionally update velocity
        vy += ay * dt;
    }

    std::pair<double, double> calculateForceBetweenParticles(const Particle &p2)
    {
        double dx = p2.x - x;
        double dy = p2.y - y;
        // TODO: Why epsilon?
        double distSquared = dx * dx + dy * dy + EPSILON2;
        double dist = std::sqrt(distSquared);
        double forceMagnitude = (G * m * p2.m) / distSquared;

        // Normalize direction
        double fx = forceMagnitude * (dx / dist);
        double fy = forceMagnitude * (dy / dist);

        return {fx, fy};
    }

    std::pair<double, double> calculateForceWithCell(const Cell &c)
    {
        double dx = c.x - x;
        double dy = c.y - y;
        // TODO: Why epsilon?
        double distSquared = dx * dx + dy * dy + EPSILON2;
        double dist = std::sqrt(distSquared);
        double forceMagnitude = (G * m * c.m) / distSquared;

        // Normalize direction
        double fx = forceMagnitude * (dx / dist);
        double fy = forceMagnitude * (dy / dist);

        return {fx, fy};
    }

    void applyForce(double fx, double fy, double dt)
    {
        double ax = fx / m;
        double ay = fy / m;

        updatePositionAndVelocity(ax, ay, dt);
    }
};

class Cell
{
public:
    unsigned int unique_id; // Unique identifier
    double x, y;            // Center of Mass
    double m;               // Mass

    Cell(unsigned int id = 0) : unique_id(id), x(0), y(0), m(0) {}

    void addParticle(Particle* p) {
        m += p->m;
        x += p->m * p->x;
        y += p->m * p->y;
    }

    std::pair<double, double> getCOM(){
        return {x/m, y/m};
    }

    void resetCOM(){
        m = 0;
        x = 0;
        y = 0;
    }
    
};

struct CellBounds
{
    long start, end;
};

class ParticleSimulation
{
private:
    // std::vector<Particle> particles;         // particle data
    // std::vector<Particle*> cell_sorted;        // array de ponteiros para a data que servirá como middle man para o sort
    //                                            // assim a actual data esta sempre no mesmo sitio. (mais eficiente ) sera worth tho ?
    // std::vector<CellBounds> cell_boundaries; // cell boundaries estes bounds sao o inicio e o fim (indices) do array cell_sorted assim se 
    //                                          // fizermos cell_boundaries[0].start / .end sabemos quais particulas estao na cell 0 rapidamente 

    // RandomGenerator rng;
    // double side_length;
    // long grid_size;

    std::vector<Particle> particles;
    std::vector<std::vector<Particle *>> cellParticles;
    std::vector<Cell> cells;
    RandomGenerator rng;
    double side_length;
    long grid_size;

public:
    ParticleSimulation(long seed, double side, long ncside, long long n_part)
        : particles(n_part), rng(seed), side_length(side), grid_size(ncside)
    {
        initializeSimulation();
    }

    //TODO: tem que se chamar initializeParticles.
    void initializeSimulation()
    {
        auto rnd01 = [this]()
        { return rng.uniform01(); };

        cellParticles.resize(grid_size * grid_size);
        cells.resize(grid_size * grid_size);

        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].x = rnd01() * side_length;
            particles[i].y = rnd01() * side_length;
            particles[i].vx = (rnd01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].vy = (rnd01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].m = rnd01() * 0.01 * (grid_size * grid_size) /
                             particles.size() / G * EPSILON2;

            int cell_x = static_cast<int>(particles[i].x / (side_length / grid_size));
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));

            // Convert 2D cell index to 1D index
            int cell_index = cell_y * grid_size + cell_x;

            cellParticles[cell_index].push_back(&particles[i]);
            cells[cell_index].addParticle(&particles[i]);
            
        }
    }
    
    void simulate(long n_time_steps)
    {
        for (long i = 0; i < n_time_steps; i++)  //TODO main loop with the right steps
        {
        //Calculate Cell center of mass

            //Calculate force for particles
            //Update position and velocity
            //Check collisons
        }
    }

    const std::vector<Particle> &getParticles() const
    {
        return particles;
    }
};

int main(int argc, char *argv[])
{
    try
    {
        if (argc != 6)
        {
            throw std::runtime_error("Usage: " + std::string(argv[0]) +
                                     " <seed> <side_length> <grid_size> <n_particles> <n_timesteps>");
        }

        int seed = std::stol(argv[1]);
        double side_length = std::stod(argv[2]);
        long grid_size = std::stol(argv[3]);
        long long n_particles = std::stoll(argv[4]);
        long n_timesteps = std::stol(argv[5]);

        ParticleSimulation simulation(seed, side_length, grid_size, n_particles);
        simulation.simulate(n_timesteps);

        // TODO: Add timing and result printing

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
