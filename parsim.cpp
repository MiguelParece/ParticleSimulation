#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <iomanip>
#include <algorithm>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

class RandomGenerator {
    private:
        unsigned int seed;
        bool useNormal;
    
    public:
        RandomGenerator(int input_seed) : seed(abs(input_seed) + 987654321), useNormal(input_seed < 0) {}
    
        double uniform01() {
            int seed_in = seed;
            seed ^= (seed << 13);
            seed ^= (seed >> 17);
            seed ^= (seed << 5);
            return 0.5 + 0.2328306e-09 * (seed_in + static_cast<int>(seed));
        }
    
        double normal01() {
            double u1, u2, z, result;
            do {
                u1 = uniform01();
                u2 = uniform01();
                z = std::sqrt(-2 * std::log(u1)) * std::cos(2 * M_PI * u2);
                result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
            } while (result < 0 || result >= 1);
            return result;
        }
    
        double getRandom01() {
            return useNormal ? normal01() : uniform01();
        }
    };

class Cell;

class Particle {
private:
    void updatePositionAndVelocity(double ax, double ay){
        x += vx * DELTAT + 0.5 * ax * DELTAT * DELTAT;
        y += vy * DELTAT + 0.5 * ay * DELTAT * DELTAT;

        vx += ax * DELTAT;
        vy += ay * DELTAT;
    }

public:
    double x, y;   // Position
    double vx, vy; // Velocity
    double m;      // Mass
    double fx, fy; // Force
    bool alive;    // has Collided

    Particle() : x(0), y(0), vx(0), vy(0), m(0), fx(0), fy(0), alive(true) {}

    void calculateForceBetweenParticles(Particle *p2);
    void calculateForceWithCell(const Cell *c);
    void applyForce(double sidelen);

    double getDistance(Particle *p)
    {
        return sqrt(pow((x - p->x), 2) + pow((y - p->y), 2));
    }
};

class Cell
{
public:
    double mx, my; // Center of Mass
    double m;      // Mass
    int x, y;      // Cell position

    Cell() : mx(0), my(0), m(0), x(0), y(0) {}

    void addParticle(Particle *p)
    {
        if (m == 0)
        {
            mx = p->x;
            my = p->y;
        }
        else
        {
            mx = (mx * m + p->m * p->x) / (m + p->m);
            my = (my * m + p->m * p->y) / (m + p->m);
        }
        m += p->m;
    }
};

void Particle::calculateForceWithCell(const Cell *c) 
{
    double dx = c->mx - x; // Changed from c->x to c->mx
    double dy = c->my - y; // Changed from c->y to c->my

    double distSquared = dx * dx + dy * dy;
    double dist = std::sqrt(distSquared);
    if (dist == 0)
    {
        return;
    }

    double forceMagnitude = (G * m * c->m) / distSquared;

    fx += forceMagnitude * (dx / dist);
    fy += forceMagnitude * (dy / dist);
}

void Particle::calculateForceBetweenParticles(Particle *p2)
{
    double dx = p2->x - x;
    double dy = p2->y - y;

    double distSquared = dx * dx + dy * dy;
    double dist = std::sqrt(distSquared);
    if (dist == 0)
    {
        return;
    }

    double forceMagnitude = (G * m * p2->m) / distSquared;

    double fx_add = forceMagnitude * (dx / dist);
    double fy_add = forceMagnitude * (dy / dist);

    fx += fx_add;
    fy += fy_add;
    p2->fx -= fx_add;
    p2->fy -= fy_add;
}

void Particle::applyForce(double sidelen){
    if (m == 0) {
        fx = 0;
        fy = 0;
        return;
    }
    
    double ax = fx / m;
    double ay = fy / m;

    updatePositionAndVelocity(ax, ay);
    
    x = fmod((x + sidelen), sidelen);
    y = fmod((y + sidelen), sidelen);

    // Reset forces for next iteration
    fx = 0;
    fy = 0;
}

class ParticleSimulation
{
private:
    std::vector<Particle> particles;
    std::vector<std::vector<Particle *>> cellParticles;
    std::vector<Cell> cells;
    RandomGenerator rng;
    double side_length;
    long grid_size = 0;
    long collisions = 0;

public:
    ParticleSimulation(long seed, double side, long ncside, long long n_part)
    : particles(n_part), 
    cellParticles(ncside * ncside),
    cells(ncside * ncside),
    rng(seed), 
    side_length(side), 
    grid_size(ncside)
    {
    init_particles();
    }

    void init_particles()
    {

        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].x = rng.getRandom01() * side_length;
            particles[i].y = rng.getRandom01() * side_length;
            particles[i].vx = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].vy = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].m = rng.getRandom01() * 0.01 * (grid_size * grid_size) /
                           particles.size() / G * EPSILON2;
        }
    }
    
    void updateCellParticles()
    {
        cellParticles.assign(grid_size * grid_size, std::vector<Particle *>{});
        cells.assign(grid_size * grid_size, Cell{}); //TODO: Precisamos disto?
        for (size_t i = 0; i < particles.size(); i++)
        {
            // Calculate cell index
            int cell_x = static_cast<int>(particles[i].x / (side_length / grid_size));
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));
            
            //TODO: REMOVE prints
            if (cell_x < 0 || cell_x >= grid_size || cell_y < 0 || cell_y >= grid_size) {
                std::cout << "[PANIC3] Cell out of bounds" << std::endl;
                continue;
            }

            // Convert 2D cell index to 1D index
            int cell_index = cell_y * grid_size + cell_x;

            cellParticles[cell_index].push_back(&particles[i]);
            cells[cell_index].x = cell_x;
            cells[cell_index].y = cell_y;
        }
    }

    void updateCOM()
    {
        for (int i = 0; i < cellParticles.size(); i++)
        {
            for (int j = 0; j < cellParticles[i].size(); j++)
            {
                cells[i].addParticle(cellParticles[i][j]);
            }
        }
    }

    void updateForces()
    {
        for (int i = 0; i < cellParticles.size(); i++)
        { // percorrer todas as cells

            std::vector<Cell> temp_cells; // Remove initial size
            temp_cells.reserve(8); // Reserve space for efficiency

            // calculate the neighbout cells
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    if (dx == 0 && dy == 0) continue; //Salta a cell do meio

                    int neighbour_cell_x = cells[i].x + dx;
                    int neighbour_cell_y = cells[i].y + dy;

                    Cell temp_cell = Cell();

                    // loop around mirror math
                    // perceber se temos de somar ou subtrair o grid lenght
                    if (neighbour_cell_x >= grid_size)
                    {
                        temp_cell.mx += side_length;
                    }
                    else if (neighbour_cell_x < 0)
                    {
                        temp_cell.mx -= side_length;
                    }
                    if (neighbour_cell_y >= grid_size)
                    {
                        temp_cell.my += side_length;
                    }
                    else if (neighbour_cell_y < 0)
                    {
                        temp_cell.my -= side_length;
                    }

                    // garantir que obtemos a mirror cell
                    neighbour_cell_x = (neighbour_cell_x + grid_size) % grid_size; // calcula o mirror caso seja preciso
                    neighbour_cell_y = (neighbour_cell_y + grid_size) % grid_size;

                    if (neighbour_cell_x < 0 || neighbour_cell_x >= grid_size || neighbour_cell_y < 0 || neighbour_cell_y >= grid_size) {
                        std::cout << "[PANIC1] Cell out of bounds" << std::endl;
                        continue;
                    }

                    // obter o index  cell
                    int cell_index = neighbour_cell_x + neighbour_cell_y * grid_size;
                    

                    // somar as coordenadas

                    temp_cell.mx += cells[cell_index].mx;
                    temp_cell.my += cells[cell_index].my;

                    temp_cell.m = cells[cell_index].m;
                    // adicionar a cell ao array

                    temp_cells.push_back(temp_cell);
                }
            }

            for (int j = 0; j < cellParticles[i].size(); j++)
            { // por cada por todas as particulas de cada cell
                if (j != cellParticles[i].size() - 1)
                { // evitar calculos duplicados
                    for (int k = j + 1; k < cellParticles[i].size(); k++)
                    { // ver todas as outras particulas dentro da mesma cell
                            if (cellParticles[i][j]->alive == true && cellParticles[i][k]->alive == true){
                                cellParticles[i][j]->calculateForceBetweenParticles(cellParticles[i][k]);
                            }
                    }
                }

                // calcular as forcas
                for (int l = 0; l < temp_cells.size(); l++)
                {
                    if (cellParticles[i][j]->alive == true){
                        cellParticles[i][j]->calculateForceWithCell(&temp_cells[l]);
                    }
                }
            }
        }
    }

    void updatePositionAndVelocity(double sidelen)
    {
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].applyForce(sidelen);
        }
        updateCellParticles();
    }

    void checkCollisions()
    {
        for (int i = 0; i < cellParticles.size(); i++)
        {   
            std::unordered_set<Particle *> collisionSet; //set used to temporarily store collided particles in a cell
            for (int j = 0; j < cellParticles[i].size(); j++)
            {
                if (cellParticles[i][j]->alive == true) { // Only check particles that are alive
                    for (int k = j + 1; k < cellParticles[i].size(); k++)
                    {
                        // if both particles are alive, check if distance between them is smaller than EPSILON
                        if (cellParticles[i][k]->alive == true &&
                            cellParticles[i][j]->getDistance(cellParticles[i][k]) < sqrt(EPSILON2))
                        {
                            // if particles not in set, new collision detected
                            //std::cout << std::fixed << std::setprecision(6) << "Collision of dist: " << cellParticles[i][j]->getDistance(cellParticles[i][k]) << std::endl;
                            if (collisionSet.count(cellParticles[i][j]) == 0 && collisionSet.count(cellParticles[i][k]) == 0)
                                collisions++;

                            // add particles to set since new they have collided
                            collisionSet.insert(cellParticles[i][j]);
                            collisionSet.insert(cellParticles[i][k]);
                        }
                    }
                }
            }
            for (const auto &elem : collisionSet) // set all particles inside collisionSet as "dead"
            {
                elem->alive = false;
                elem->m = 0;
            }
        }
    }

    void simulate(long n_time_steps)
    {
        // for (size_t j = 0; j < particles.size(); j++)
        //     {
        //         // std::cout << std::fixed << std::setprecision(3) << "Particle " << j << ": mass=" << particles[j].m << " x=" << particles[j].x << " y=" << particles[j].y << " vx=" << particles[j].vx << " vy=" << particles[j].vy << std::endl;
        //     }
        //     for (size_t j = 0; j < cells.size(); j++)
        //     {
        //         // std::cout << std::fixed << std::setprecision(3) << "Cell " << j << " x: " << cells[j].mx << " y: " << cells[j].my << " m: " << cells[j].m << std::endl;
        //     }
        updateCellParticles();
        for (long i = 0; i < n_time_steps; i++) // TODO main loop with the right steps
        {   
            // Calculate Cell center of mass
            updateCOM();
            // Calculate force for particles
            updateForces();
            // Update position and velocity
            updatePositionAndVelocity(side_length);
            // Check collisons
            checkCollisions();
            // std::cout << "t=" << i << std::endl;
            // for (size_t j = 0; j < particles.size(); j++)
            // {
            //     std::cout << std::fixed << std::setprecision(3) << "Particle " << j << ": mass=" << particles[j].m << " x=" << particles[j].x << " y=" << particles[j].y << " vx=" << particles[j].vx << " vy=" << particles[j].vy << std::endl;
            // }
        }

        std::cout << std::fixed << std::setprecision(3) << particles[0].x << " " << particles[0].y << std::endl;
        std::cout << collisions << std::endl;
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

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
