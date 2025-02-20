#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>

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
private:
    void updatePositionAndVelocity(double ax, double ay)
    {
        x += vx * DELTAT + 0.5 * ax * DELTAT * DELTAT;
        y += vy * DELTAT + 0.5 * ay * DELTAT * DELTAT;

        vx += ax * DELTAT;
        vy += ay * DELTAT;
    }
public:
    double x, y;   // Position
    double vx, vy; // Velocity
    double m;      // Mass
    double fx, fy;  // Force
    bool alive;    // has Collided

    Particle() : x(0), y(0), vx(0), vy(0), m(0), alive(true) {}

    void calculateForceBetweenParticles(Particle *p2) {
        double dx = p2->x - x;
        double dy = p2->y - y;
    
        // Small constant to avoid division by zero issues
        double distSquared = dx * dx + dy * dy;
        double dist = std::sqrt(distSquared);
        if (dist == 0) {
            return; // Avoid self-force
        }
    
        double forceMagnitude = (G * m * p2->m) / distSquared;
    
        // Normalize and apply force
        double fx_add = forceMagnitude * (dx / dist);
        double fy_add = forceMagnitude * (dy / dist);
    
        fx += fx_add;
        fy += fy_add;
        p2->fx -= fx_add;
        p2->fy -= fy_add;
    }

    void calculateForceWithCell(const Cell *c)
    {
        double dx = c->x - x;
        double dy = c->y - y;
    
        // Small constant to avoid division by zero issues
        double distSquared = dx * dx + dy * dy;
        double dist = std::sqrt(distSquared);
        if (dist == 0) {
            return; // Avoid self-force
        }
    
        double forceMagnitude = (G * m * c->m) / distSquared;
    
        // Normalize and apply force
        fx += forceMagnitude * (dx / dist);
        fy += forceMagnitude * (dy / dist);
    }


    void applyForce()
    {
        double ax = fx / m;
        double ay = fy / m;

        updatePositionAndVelocity(ax, ay);
    }

    double getDistance(Particle* p) {
        return sqrt(pow((x - p->x), 2) + pow((y - p->y), 2));
    }
};

class Cell
{
public:
    double mx, my;          // Center of Mass
    double m;               // Mass
    int x, y;               // Cell position

    Cell(unsigned int id = 0) : mx(0), my(0), m(0) , x(0), y(0) {}

    void addParticle(Particle* p) {
        m += p->m;
        mx += p->m * p->x;
        my += p->m * p->y;
    }

    void addParticle(Particle* p) {
        if (m == 0) { 
            mx = p->x;
            my = p->y;
        } else {
            mx = (mx * m + p->m * p->x) / (m + p->m);
            my = (my * m + p->m * p->y) / (m + p->m);
        }
        m += p->m;
    }
    
};

class ParticleSimulation
{
private:
    std::vector<Particle> particles;
    std::vector<std::vector<Particle *>> cellParticles;
    std::vector<Cell> cells;
    RandomGenerator rng;
    double side_length;
    long grid_size;
    long collisions;

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

        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].x = rnd01() * side_length;
            particles[i].y = rnd01() * side_length;
            particles[i].vx = (rnd01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].vy = (rnd01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].m = rnd01() * 0.01 * (grid_size * grid_size) /
                             particles.size() / G * EPSILON2;
            
        }
    }

    void updateCOM(){
        cellParticles.assign(grid_size * grid_size, std::vector<Particle*>{});
        cells.assign(grid_size * grid_size, Cell{});  
        for (size_t i = 0; i < particles.size(); i++)
        {
            // Calculate cell index
            int cell_x = static_cast<int>(particles[i].x / (side_length / grid_size));
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));

            // Convert 2D cell index to 1D index
            int cell_index = cell_y * grid_size + cell_x;

            cellParticles[cell_index].push_back(&particles[i]);
            cells[cell_index].addParticle(&particles[i]);
            cells[cell_index].x = cell_x;
            cells[cell_index].y = cell_y;
        }
    }
    

    void updateForces(){
        for(int i = 0 ; i < cellParticles.size();i++){ // percorrer todas as cells
            for(int j = 0;j<cellParticles[i].size();i++){ // por cada por todas as particulas de cada cell
                   for(int k = j+1; j<cellParticles[i].size()-1;k++){ // ver todas as outras particulas dentro da mesma cell
                        cellParticles[i][j]->calculateForceBetweenParticles(cellParticles[i][k]);
                   }
                // ver as cells que estao arround desta particula

                for(int dx=-1;dx<1;dx++){
                    for(int dy=-1; dy<1;dy++){
                        if(cells[i].x==0 && cells[i].y==0){continue;} // salta a cell do meio 

                        int neighbour_cell_x = cells[i].x + dx;
                        int neighbour_cell_y = cells[i].y + dy;

                        //necessario loopback
                        neighbour_cell_x = (neighbour_cell_x + grid_size) % grid_size; // calcula o mirror caso seja preciso
                        neighbour_cell_y = (neighbour_cell_y + grid_size) % grid_size;
                        int neighbour_cell_index = neighbour_cell_y + neighbour_cell_y * grid_size;
                        cellParticles[i][j]->calculateForceWithCell(&cells[neighbour_cell_index]);
                    }
                }

            }
        }
    }

    void updatePositionAndVelocity(){
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].applyForce();
        }
    }

    void checkCollisions(){
        for (int i = 0; i < cellParticles.size(); i++) {
            std::unordered_set<Particle*> collisionSet;
            for (int j = 0; j < cellParticles[i].size(); j++) {
                for (int k = j + 1; k < cellParticles[i].size(); k++) {
                    // check distance between particles
                    if (cellParticles[i][j]->getDistance(cellParticles[i][k]) < EPSILON2) {
                        // if particles not in set, increment collision counter
                        if (collisionSet.count(cellParticles[i][j]) == 0 && collisionSet.count(cellParticles[i][j]) == 0) {
                            collisions++;
                        }
                        // if distance < epsilon, add particles to set
                        collisionSet.insert(cellParticles[i][j]);
                        collisionSet.insert(cellParticles[i][k]);
                    }
                    // if particles not in set, increment collision counter
                }
            }
            for (const auto& elem : collisionSet) {
                elem->alive = false;
            }
        }
    }
    
    void simulate(long n_time_steps)
    {
        for (long i = 0; i < n_time_steps; i++)  //TODO main loop with the right steps
        {
        //Calculate Cell center of mass
            updateCOM();
            //Calculate force for particles
            updateForces();
            //Update position and velocity
            updatePositionAndVelocity();
            //Check collisons
            checkCollisions();
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
