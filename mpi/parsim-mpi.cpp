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
#include <omp.h>
#include <mpi.h>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define EPSILON 0.005
#define DELTAT 0.1

#define DEBUG_MPI 0 

#if DEBUG_MPI
#define DEBUG_PRINT(rank, msg, ...)                                                         \
    {                                                                                       \
        char debug_buf[256];                                                                \
        snprintf(debug_buf, sizeof(debug_buf), "[Rank %d] " msg "\n", rank, ##__VA_ARGS__); \
        fprintf(stderr, "%s", debug_buf);                                                   \
        fflush(stderr);                                                                     \
    }
#else
#define DEBUG_PRINT(rank, msg, ...) \
    {                               \
    }
#endif

// Structure for MPI transfer of particles
struct ParticleData
{
    int id;         // Particle ID
    double x, y;    // Position
    double vx, vy;  // Velocity
    double m;       // Mass
    bool alive;     // Has collided
    int cell_index; // Which cell it belongs to
};

class RandomGenerator
{
private:
    unsigned int seed;
    bool useNormal;

public:
    RandomGenerator(int input_seed) : seed(abs(input_seed) + 987654321), useNormal(input_seed < 0) {}

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
            result = 0.5 + 0.15 * z;
        } while (result < 0 || result >= 1);
        return result;
    }

    double getRandom01()
    {
        return useNormal ? normal01() : uniform01();
    }
};

class Cell;

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
    int id;         // Unique identifier for each particle
    double x, y;    // Position
    double vx, vy;  // Velocity
    double m;       // Mass
    double fx, fy;  // Force
    bool alive;     // has Collided
    int cell_index; // which cell it belongs to
    int proc_owner; // Process that owns this particle

    Particle() : id(-1), x(0), y(0), vx(0), vy(0), m(0), fx(0), fy(0), alive(true), proc_owner(-1) {}

    Particle(const ParticleData &data, int particle_id = -1) : id(data.id),
                                                               x(data.x), y(data.y),
                                                               vx(data.vx), vy(data.vy),
                                                               m(data.m), fx(0), fy(0),
                                                               alive(data.alive),
                                                               cell_index(data.cell_index) {}

    // Convert to MPI transfer structure
    ParticleData toParticleData() const
    {
        ParticleData data;
        data.id = id;
        data.x = x;
        data.y = y;
        data.vx = vx;
        data.vy = vy;
        data.m = m;
        data.alive = alive;
        data.cell_index = cell_index;
        return data;
    }

    void calculateForceBetweenParticles(Particle *p2);
    void calculateForceWithCell(const Cell *c);
    void applyForce(double sidelen, long grid_size, std::vector<Cell> &cells,
                    int row_start, int row_end, int proc_rows, int rank);

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
    bool change_flag;
    bool is_ghost; // Is this a ghost cell?

    Cell() : mx(0), my(0), m(0), x(0), y(0), change_flag(false), is_ghost(false) {}

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

// Structure to send center of mass data between processes
struct CellCOMData
{
    double mx, my;
    double m;
    int x, y;
};

void Particle::calculateForceWithCell(const Cell *c)
{
    double dx = c->mx - x;
    double dy = c->my - y;

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

void Particle::applyForce(double sidelen, long grid_size, std::vector<Cell> &cells,
                          int row_start, int row_end, int proc_rows, int rank)
{
    if (m == 0 || !alive)
    {
        fx = 0;
        fy = 0;
        return;
    }

    double ax = fx / m;
    double ay = fy / m;

    updatePositionAndVelocity(ax, ay);

    // Wrap around the simulation space properly
    x = fmod(x + sidelen, sidelen);
    y = fmod(y + sidelen, sidelen);

    // Calculate cell coordinates
    int cell_y = static_cast<int>(y / (sidelen / grid_size));
    int cell_x = static_cast<int>(x / (sidelen / grid_size));

    // Ensure cell coordinates are within grid bounds - convert types to match
    cell_x = std::max(0, std::min(cell_x, static_cast<int>(grid_size) - 1));
    cell_y = std::max(0, std::min(cell_y, static_cast<int>(grid_size) - 1));

    // Check if particle is in this process's domain
    if (cell_y >= row_start && cell_y < row_end)
    {
        int new_cell_index = (cell_y - row_start) * grid_size + cell_x;

        if (new_cell_index >= 0 && new_cell_index < static_cast<int>(cells.size()))
        {
            if (cell_index != new_cell_index && cell_index >= 0 &&
                cell_index < static_cast<int>(cells.size()))
            {
                cells[cell_index].change_flag = true;
            }
            cell_index = new_cell_index;
        }
        else
        {
            // This shouldn't happen with the bounds checking above
            cell_index = -1;
        }
    }
    else
    {
        proc_owner = cell_y / proc_rows;
        proc_owner = std::max(0, proc_owner); 
        cell_index = -1;                      
    }

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
    long grid_size;
    long local_collisions = 0;
    long global_collisions = 0;
    int num_procs;
    int rank;
    int rows_per_proc;
    int my_row_start;
    int my_row_end;
    MPI_Datatype mpi_particle_type;
    MPI_Datatype mpi_cell_com_type;

    void initMPITypes()
    {
        MPI_Datatype types[8] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL, MPI_INT};
        int blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1}; 
        MPI_Aint offsets[8];

        offsets[0] = offsetof(ParticleData, id);
        offsets[1] = offsetof(ParticleData, x);
        offsets[2] = offsetof(ParticleData, y);
        offsets[3] = offsetof(ParticleData, vx);
        offsets[4] = offsetof(ParticleData, vy);
        offsets[5] = offsetof(ParticleData, m);
        offsets[6] = offsetof(ParticleData, alive);
        offsets[7] = offsetof(ParticleData, cell_index);

        MPI_Type_create_struct(8, blocklengths, offsets, types, &mpi_particle_type);
        MPI_Type_commit(&mpi_particle_type);

        // Create MPI datatype for CellCOMData
        MPI_Datatype com_types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
        int com_blocklengths[5] = {1, 1, 1, 1, 1};
        MPI_Aint com_offsets[5];

        com_offsets[0] = offsetof(CellCOMData, mx);
        com_offsets[1] = offsetof(CellCOMData, my);
        com_offsets[2] = offsetof(CellCOMData, m);
        com_offsets[3] = offsetof(CellCOMData, x);
        com_offsets[4] = offsetof(CellCOMData, y);

        MPI_Type_create_struct(5, com_blocklengths, com_offsets, com_types, &mpi_cell_com_type);
        MPI_Type_commit(&mpi_cell_com_type);
    }

public:
    ParticleSimulation(long seed, double side, long ncside, long long n_part)
        : rng(seed), side_length(side), grid_size(ncside),
          mpi_particle_type(MPI_DATATYPE_NULL),
          mpi_cell_com_type(MPI_DATATYPE_NULL)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        rows_per_proc = ncside / num_procs;
        if (rows_per_proc == 0)
        {
            rows_per_proc = 1;
        }

        if (rank == 0)
        {
            // Only initialize particles on rank 0
            particles.resize(n_part);
            init_particles();
        }
    }

    ~ParticleSimulation()
    {
        if (mpi_particle_type != MPI_DATATYPE_NULL)
        {
            MPI_Type_free(&mpi_particle_type);
            mpi_particle_type = MPI_DATATYPE_NULL;
        }
        if (mpi_cell_com_type != MPI_DATATYPE_NULL)
        {
            MPI_Type_free(&mpi_cell_com_type);
            mpi_cell_com_type = MPI_DATATYPE_NULL;
        }
    }

    void setupMPI()
    {
        initMPITypes();

        // Calculate row divisions
        my_row_start = rank * rows_per_proc;
        my_row_end = (rank == num_procs - 1) ? grid_size : (rank + 1) * rows_per_proc;

        int local_grid_size = (my_row_end - my_row_start) * grid_size;
        cells.resize(local_grid_size);
        cellParticles.resize(local_grid_size);
    }

    void distributeInitialParticles(long long n_part)
    {
        distributeParticles(n_part);
    }

    void init_particles()
    {
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i].id = i; // Assign unique ID
            particles[i].x = rng.getRandom01() * side_length;
            particles[i].y = rng.getRandom01() * side_length;
            particles[i].vx = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].vy = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].m = rng.getRandom01() * 0.01 * (grid_size * grid_size) /
                             particles.size() / G * EPSILON2;

            // Determine which process should own this particle
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));
            particles[i].proc_owner = cell_y / rows_per_proc;
            if (particles[i].proc_owner >= num_procs)
            {
                particles[i].proc_owner = num_procs - 1;
            }
        }
    }

    void distributeParticles(long long n_part)
    {
        if (rank == 0)
        {
            std::vector<int> particles_per_proc(num_procs, 0);
            for (const auto &p : particles)
            {
                particles_per_proc[p.proc_owner]++;
            }

            MPI_Bcast(particles_per_proc.data(), num_procs, MPI_INT, 0, MPI_COMM_WORLD);

            for (int dest_rank = 1; dest_rank < num_procs; dest_rank++)
            {
                std::vector<ParticleData> particles_to_send;
                for (const auto &p : particles)
                {
                    if (p.proc_owner == dest_rank)
                    {
                        particles_to_send.push_back(p.toParticleData());
                    }
                }

                MPI_Send(particles_to_send.data(), particles_to_send.size(),
                         mpi_particle_type, dest_rank, 0, MPI_COMM_WORLD);
            }

            // Keep only particles for rank 0
            std::vector<Particle> my_particles;
            for (const auto &p : particles)
            {
                if (p.proc_owner == 0)
                {
                    my_particles.push_back(p);
                }
            }
            particles = std::move(my_particles);
        }
        else
        {
            // Receive counts from root
            std::vector<int> particles_per_proc(num_procs);
            MPI_Bcast(particles_per_proc.data(), num_procs, MPI_INT, 0, MPI_COMM_WORLD);

            // Receive particles for this process
            int my_particle_count = particles_per_proc[rank];
            std::vector<ParticleData> received_particles(my_particle_count);

            MPI_Status status;
            MPI_Recv(received_particles.data(), my_particle_count,
                     mpi_particle_type, 0, 0, MPI_COMM_WORLD, &status);

            // Convert to local particles
            particles.resize(my_particle_count);
            for (int i = 0; i < my_particle_count; i++)
            {
                particles[i] = Particle(received_particles[i]);
            }
        }
    }

    void updateCellParticles()
    {
        int local_grid_size = (my_row_end - my_row_start) * grid_size;
        DEBUG_PRINT(rank, "Starting updateCellParticles with %zu cells", cells.size());

        // Make sure cellParticles has the right size
        if (cellParticles.size() != local_grid_size)
        {
            DEBUG_PRINT(rank, "Resizing cellParticles from %zu to %d", cellParticles.size(), local_grid_size);
            cellParticles.resize(local_grid_size);
        }

        std::vector<std::vector<Particle *>> new_cellParticles(local_grid_size);

        for (int i = 0; i < local_grid_size; i++)
        {
            if (i < static_cast<int>(cells.size()))
            {
                cells[i].change_flag = false;
            }
        }

        for (size_t i = 0; i < particles.size(); i++)
        {
            if (!particles[i].alive)
            {
                continue;
            }

            int cell_index = particles[i].cell_index;

            // Skip particles with invalid cell indices
            if (cell_index < 0 || cell_index >= local_grid_size)
            {
                continue;
            }

            new_cellParticles[cell_index].push_back(&particles[i]);
        }

        cellParticles = std::move(new_cellParticles);

        DEBUG_PRINT(rank, "Finished updateCellParticles");
    }

    void exchangeParticles()
    {

        DEBUG_PRINT(rank, "Starting exchangeParticles with %zu particles", particles.size());

        std::vector<std::vector<ParticleData>> particles_to_send(num_procs);

        // Identify particles to send to other processes
        for (size_t i = 0; i < particles.size(); i++)
        {
            if (particles[i].cell_index == -1 && particles[i].proc_owner != rank)
            {
                int target_proc = std::max(0, std::min(particles[i].proc_owner, num_procs - 1));
                particles_to_send[target_proc].push_back(particles[i].toParticleData());
            }
        }

        std::vector<int> send_counts(num_procs, 0);
        for (int i = 0; i < num_procs; i++)
        {
            send_counts[i] = particles_to_send[i].size();
        }

        std::vector<int> recv_counts(num_procs, 0);
        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<MPI_Request> requests;
        for (int i = 0; i < num_procs; i++)
        {
            if (i != rank && send_counts[i] > 0)
            {
                MPI_Request req;
                MPI_Isend(particles_to_send[i].data(), send_counts[i],
                          mpi_particle_type, i, 0, MPI_COMM_WORLD, &req);
                requests.push_back(req);
            }
        }

        std::vector<ParticleData> received_particles;
        for (int i = 0; i < num_procs; i++)
        {
            if (i != rank && recv_counts[i] > 0)
            {
                std::vector<ParticleData> temp_received(recv_counts[i]);
                MPI_Status status;
                MPI_Recv(temp_received.data(), recv_counts[i],
                         mpi_particle_type, i, 0, MPI_COMM_WORLD, &status);

                received_particles.insert(received_particles.end(), temp_received.begin(), temp_received.end());
            }
        }

        if (!requests.empty())
        {
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        }

        // Remove particles that were sent to other processes
        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                           [this](const Particle &p)
                           {
                               return p.cell_index == -1 && p.proc_owner != rank;
                           }),
            particles.end());

        for (const auto &p_data : received_particles)
        {
            Particle new_particle(p_data);

            // Update cell index for new domain
            int cell_y = static_cast<int>(new_particle.y / (side_length / grid_size));
            int cell_x = static_cast<int>(new_particle.x / (side_length / grid_size));

            // Ensure coordinates are within bounds - use static_cast to match types
            cell_x = std::max(0, std::min(cell_x, static_cast<int>(grid_size) - 1));
            cell_y = std::max(0, std::min(cell_y, static_cast<int>(grid_size) - 1));

            // Only add if the particle belongs to this process's domain
            if (cell_y >= my_row_start && cell_y < my_row_end)
            {
                new_particle.cell_index = (cell_y - my_row_start) * grid_size + cell_x;
                new_particle.proc_owner = rank;
                particles.push_back(new_particle);
            }
        }

        DEBUG_PRINT(rank, "Finished exchangeParticles. Total particles: %zu", particles.size());
    }

    void updateCOM()
    {
        int local_grid_size = (my_row_end - my_row_start) * grid_size;

        // Resize only if needed
        if (cellParticles.size() != local_grid_size)
            cellParticles.resize(local_grid_size);
        if (cells.size() != local_grid_size)
            cells.resize(local_grid_size);

        for (auto &cell_particles : cellParticles)
        {
            cell_particles.clear();
        }

        std::fill(cells.begin(), cells.end(), Cell());

        DEBUG_PRINT(rank, "updateCOM: Processing %zu particles", particles.size());

        for (size_t i = 0; i < particles.size(); i++)
        {
            // Calculate cell index
            int cell_x = static_cast<int>(particles[i].x / (side_length / grid_size));
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));

            // Extra boundary checks
            if (cell_x < 0)
                cell_x = 0;
            else if (cell_x >= grid_size)
                cell_x = grid_size - 1;

            if (cell_y < 0)
                cell_y = 0;
            else if (cell_y >= grid_size)
                cell_y = grid_size - 1;

            // Skip if particle is outside my domain
            if (cell_y < my_row_start || cell_y >= my_row_end)
            {
                DEBUG_PRINT(rank, "Particle %d outside domain: cell_y=%d, range=[%d,%d)",
                            particles[i].id, cell_y, my_row_start, my_row_end);
                continue;
            }

            // Convert 2D cell index to local 1D index
            int local_cell_y = cell_y - my_row_start;
            int cell_index = local_cell_y * grid_size + cell_x;

            // Extra bounds check
            if (cell_index < 0 || cell_index >= local_grid_size)
            {
                DEBUG_PRINT(rank, "WARNING: Invalid cell index %d (max=%d) for particle %d at (%f,%f)",
                            cell_index, local_grid_size - 1, particles[i].id, particles[i].x, particles[i].y);
                continue;
            }

            particles[i].cell_index = cell_index; 
            cellParticles[cell_index].push_back(&particles[i]);

            cells[cell_index].addParticle(&particles[i]);
            cells[cell_index].x = cell_x;
            cells[cell_index].y = cell_y;
        }

        DEBUG_PRINT(rank, "updateCOM: Completed. Exchanging boundary info");
        exchangeBoundaryCellInfo();
    }

    void exchangeBoundaryCellInfo()
    {
        int local_grid_size = (my_row_end - my_row_start) * grid_size;
        cells.resize(local_grid_size); // Remove any ghost cells

        std::vector<CellCOMData> top_row_data(grid_size, {0});     
        std::vector<CellCOMData> bottom_row_data(grid_size, {0});  
        std::vector<CellCOMData> ghost_bottom_row(grid_size, {0}); 
        std::vector<CellCOMData> ghost_top_row(grid_size, {0});   

        if (rank > 0 && local_grid_size >= grid_size)
        {
            for (int x = 0; x < grid_size; x++)
            {
                int cell_index = x; 
                if (cell_index >= 0 && cell_index < local_grid_size)
                {
                    top_row_data[x].mx = cells[cell_index].mx;
                    top_row_data[x].my = cells[cell_index].my;
                    top_row_data[x].m = cells[cell_index].m;
                    top_row_data[x].x = cells[cell_index].x;
                    top_row_data[x].y = cells[cell_index].y;
                }
            }
        }

        if (rank < num_procs - 1 && local_grid_size >= grid_size)
        {
            int last_row_index = local_grid_size - grid_size;
            for (int x = 0; x < grid_size; x++)
            {
                int cell_index = last_row_index + x;
                if (cell_index >= 0 && cell_index < local_grid_size)
                {
                    bottom_row_data[x].mx = cells[cell_index].mx;
                    bottom_row_data[x].my = cells[cell_index].my;
                    bottom_row_data[x].m = cells[cell_index].m;
                    bottom_row_data[x].x = cells[cell_index].x;
                    bottom_row_data[x].y = cells[cell_index].y;
                }
            }
        }

        MPI_Status status;

        if (rank > 0)
        {
            MPI_Send(top_row_data.data(), grid_size, mpi_cell_com_type, rank - 1, 0, MPI_COMM_WORLD);

            MPI_Recv(ghost_bottom_row.data(), grid_size, mpi_cell_com_type, rank - 1, 1, MPI_COMM_WORLD, &status);
        }

        if (rank < num_procs - 1)
        {
            MPI_Send(bottom_row_data.data(), grid_size, mpi_cell_com_type, rank + 1, 1, MPI_COMM_WORLD);

            MPI_Recv(ghost_top_row.data(), grid_size, mpi_cell_com_type, rank + 1, 0, MPI_COMM_WORLD, &status);
        }

        if (rank > 0)
        {
            for (int x = 0; x < grid_size; x++)
            {
                Cell ghost_cell;
                ghost_cell.mx = ghost_bottom_row[x].mx;
                ghost_cell.my = ghost_bottom_row[x].my - side_length; // Adjust for wrapping
                ghost_cell.m = ghost_bottom_row[x].m;
                ghost_cell.x = ghost_bottom_row[x].x;
                ghost_cell.y = ghost_bottom_row[x].y;
                ghost_cell.is_ghost = true;

                cells.push_back(ghost_cell);
            }
        }

        // For top ghost row (from next rank)
        if (rank < num_procs - 1)
        {
            for (int x = 0; x < grid_size; x++)
            {
                Cell ghost_cell;
                ghost_cell.mx = ghost_top_row[x].mx;
                ghost_cell.my = ghost_top_row[x].my + side_length; // Adjust for wrapping
                ghost_cell.m = ghost_top_row[x].m;
                ghost_cell.x = ghost_top_row[x].x;
                ghost_cell.y = ghost_top_row[x].y;
                ghost_cell.is_ghost = true;

                cells.push_back(ghost_cell);
            }
        }
    }

    void updateForces()
    {
        DEBUG_PRINT(rank, "Starting updateForces with %zu cells", cells.size());

        int local_grid_size = (my_row_end - my_row_start) * grid_size;

        MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel for
        for (int i = 0; i < local_grid_size; i++)
        {
            if (i < 0 || i >= static_cast<int>(cells.size()) || i >= static_cast<int>(cellParticles.size()))
            {
                continue;
            }

            std::vector<Cell> temp_cells;
            temp_cells.reserve(8);

            int cell_x = cells[i].x;
            int cell_y = cells[i].y;

            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    if (dx == 0 && dy == 0)
                        continue;

                    int neighbor_x = (cell_x + dx + grid_size) % grid_size;
                    int neighbor_y = cell_y + dy;

                    if (neighbor_x < 0 || neighbor_x >= grid_size)
                        continue;

                    Cell temp_cell;
                    temp_cell.mx = 0.0;
                    temp_cell.my = 0.0;
                    temp_cell.m = 0.0;

                    bool found_cell = false;

                    if (neighbor_y < my_row_start)
                    {
                        // Look for ghost cells from previous process
                        int ghost_start = local_grid_size;
                        for (int g = ghost_start; g < static_cast<int>(cells.size()); g++)
                        {
                            if (cells[g].x == neighbor_x &&
                                (cells[g].y == neighbor_y || cells[g].y == neighbor_y + grid_size))
                            {
                                temp_cell.mx = cells[g].mx;
                                temp_cell.my = cells[g].my;
                                temp_cell.m = cells[g].m;
                                found_cell = true;
                                break;
                            }
                        }

                        if (!found_cell && rank == 0)
                        {
                            // For rank 0, wrap to bottom
                            int wrap_y = grid_size - 1;
                            int wrap_index = (wrap_y - my_row_start) * grid_size + neighbor_x;
                            if (wrap_index >= 0 && wrap_index < local_grid_size)
                            {
                                temp_cell.mx = cells[wrap_index].mx;
                                temp_cell.my = cells[wrap_index].my - side_length;
                                temp_cell.m = cells[wrap_index].m;
                                found_cell = true;
                            }
                        }
                    }
                    else if (neighbor_y >= my_row_end)
                    {
                        // Look for ghost cells from next process
                        int ghost_start = local_grid_size;
                        for (int g = ghost_start; g < static_cast<int>(cells.size()); g++)
                        {
                            if (cells[g].x == neighbor_x &&
                                (cells[g].y == neighbor_y || cells[g].y == neighbor_y - grid_size))
                            {
                                temp_cell.mx = cells[g].mx;
                                temp_cell.my = cells[g].my;
                                temp_cell.m = cells[g].m;
                                found_cell = true;
                                break;
                            }
                        }

                        if (!found_cell && rank == num_procs - 1)
                        {
                            // For last rank, wrap to top
                            int wrap_y = 0;
                            int wrap_index = (wrap_y - my_row_start) * grid_size + neighbor_x;
                            if (wrap_index >= 0 && wrap_index < local_grid_size)
                            {
                                temp_cell.mx = cells[wrap_index].mx;
                                temp_cell.my = cells[wrap_index].my + side_length;
                                temp_cell.m = cells[wrap_index].m;
                                found_cell = true;
                            }
                        }
                    }
                    else
                    {
                        // Regular cell within this process's domain
                        int neighbor_index = (neighbor_y - my_row_start) * grid_size + neighbor_x;
                        if (neighbor_index >= 0 && neighbor_index < local_grid_size)
                        {
                            temp_cell.mx = cells[neighbor_index].mx;
                            temp_cell.my = cells[neighbor_index].my;
                            temp_cell.m = cells[neighbor_index].m;
                            found_cell = true;
                        }
                    }

                    if (found_cell && temp_cell.m > 0)
                    {
                        temp_cells.push_back(temp_cell);
                    }
                }
            }

            // Calculate forces between particles in the same cell
            if (i < static_cast<int>(cellParticles.size()) && !cellParticles[i].empty())
            {
                // Check for null pointers and add extra safety
                for (size_t j = 0; j < cellParticles[i].size(); j++)
                {
                    if (!cellParticles[i][j] || !cellParticles[i][j]->alive)
                        continue;

                    for (size_t k = j + 1; k < cellParticles[i].size(); k++)
                    {
                        if (!cellParticles[i][k] || !cellParticles[i][k]->alive)
                            continue;

                        cellParticles[i][j]->calculateForceBetweenParticles(cellParticles[i][k]);
                    }

                    for (const auto &neighbor_cell : temp_cells)
                    {
                        cellParticles[i][j]->calculateForceWithCell(&neighbor_cell);
                    }
                }
            }
        }

        DEBUG_PRINT(rank, "Finished updateForces");

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void checkCollisions()
    {
        int local_grid_size = (my_row_end - my_row_start) * grid_size;
        DEBUG_PRINT(rank, "Starting checkCollisions with %d cells", local_grid_size);

        for (int i = 0; i < local_grid_size && i < static_cast<int>(cellParticles.size()); i++)
        {
            auto &cell_particles = cellParticles[i];

            // Create a local copy of particles that we're working with to avoid modification during iteration
            std::vector<Particle *> valid_particles;

            for (auto *p : cell_particles)
            {
                if (p && p->alive)
                {
                    valid_particles.push_back(p);
                }
            }

            // Check for collisions between valid particles
            for (size_t j = 0; j < valid_particles.size(); j++)
            {
                for (size_t k = j + 1; k < valid_particles.size(); k++)
                {
                    double dx = valid_particles[j]->x - valid_particles[k]->x;
                    double dy = valid_particles[j]->y - valid_particles[k]->y;
                    double dist_squared = dx * dx + dy * dy;

                    if (dist_squared <= EPSILON2)
                    {
                        valid_particles[j]->alive = false;
                        valid_particles[k]->alive = false;
                        local_collisions++;
                    }
                }
            }
        }

        size_t original_size = particles.size();
        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                           [](const Particle &p)
                           { return !p.alive; }),
            particles.end());

        DEBUG_PRINT(rank, "Removed %zu collided particles", original_size - particles.size());
        DEBUG_PRINT(rank, "Finished checkCollisions. Local collisions: %ld", local_collisions);
    }

    void updateParticlePositions()
    {
        DEBUG_PRINT(rank, "Updating particle positions");

        for (auto &p : particles)
        {
            if (!p.alive)
                continue;

            double orig_x = p.x;
            double orig_y = p.y;
            int orig_cell = p.cell_index;

            // Update position and cell index
            p.applyForce(side_length, grid_size, cells, my_row_start, my_row_end, rows_per_proc, rank);

            // Validate new state
            if (p.cell_index >= 0 && p.cell_index < static_cast<int>(cells.size()))
            {
                // Valid cell index
            }
            else if (p.cell_index == -1)
            {
                // Marked for migration - ensure proc_owner is valid
                p.proc_owner = std::max(0, std::min(p.proc_owner, num_procs - 1));
            }
            else
            {
                // Invalid cell index
                p.cell_index = -1;
                p.proc_owner = std::max(0, std::min(p.proc_owner, num_procs - 1));
            }
        }

        DEBUG_PRINT(rank, "Finished updating particle positions");
    }

    void run(int timesteps)
    {
        for (int step = 0; step < timesteps; step++)
        {
            updateCOM();
            MPI_Barrier(MPI_COMM_WORLD);

            updateForces();
            MPI_Barrier(MPI_COMM_WORLD);

            updateParticlePositions();
            MPI_Barrier(MPI_COMM_WORLD);

            exchangeParticles();
            MPI_Barrier(MPI_COMM_WORLD);

            updateCellParticles();
            MPI_Barrier(MPI_COMM_WORLD);

            checkCollisions();
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // Collect total collisions from all processes
        MPI_Reduce(&local_collisions, &global_collisions, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    std::vector<Particle> gatherAllParticles()
    {
        int local_particle_count = particles.size();

        // Create arrays to store counts and displacements for each process
        std::vector<int> particle_counts(num_procs, 0);
        std::vector<int> particle_displs(num_procs, 0);

        // Gather the number of particles from each process
        MPI_Allgather(&local_particle_count, 1, MPI_INT,
                      particle_counts.data(), 1, MPI_INT,
                      MPI_COMM_WORLD);

        // Calculate displacements and total particle count
        int total_particles = 0;
        for (int i = 0; i < num_procs; i++)
        {
            particle_displs[i] = total_particles;
            total_particles += particle_counts[i];
        }

        // Convert local particles to ParticleData for transmission
        std::vector<ParticleData> local_particle_data(local_particle_count);
        for (int i = 0; i < local_particle_count; i++)
        {
            local_particle_data[i] = particles[i].toParticleData();
        }

        // Prepare vectors only on root process to save memory
        std::vector<ParticleData> all_particle_data;
        std::vector<Particle> all_particles;

        if (rank == 0)
        {
            all_particle_data.resize(total_particles);
        }

        // Gather all particles from all processes
        MPI_Gatherv(
            local_particle_data.data(), local_particle_count, mpi_particle_type,
            rank == 0 ? all_particle_data.data() : nullptr,
            particle_counts.data(), particle_displs.data(),
            mpi_particle_type, 0, MPI_COMM_WORLD);

        // Convert back to Particles on root process
        if (rank == 0)
        {
            all_particles.reserve(total_particles);
            for (const auto &p_data : all_particle_data)
            {
                all_particles.push_back(Particle(p_data));
            }
        }

        return all_particles;
    }

    Particle getParticle(int index) const
    {
        if (rank == 0)
        {
            if (index < particles.size())
            {
                return particles[index];
            }

            // If not found locally, do a global search
            Particle result;
            bool found = false;

            // Broadcast the index to search for
            int local_found = 0;
            MPI_Bcast(&index, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Gather information about which process has the particle
            for (int i = 1; i < num_procs; i++)
            {
                int has_particle = 0;

                // Send request to other processes
                MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Recv(&has_particle, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (has_particle)
                {
                    ParticleData p_data;
                    MPI_Recv(&p_data, 1, mpi_particle_type, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    result = Particle(p_data);
                    found = true;
                    break;
                }
            }

            // If not found in any process, return a default particle
            if (!found)
            {
                result.x = 0.0;
                result.y = 0.0;
            }

            return result;
        }
        else //Non-root
        {
            int requested_index;
            MPI_Bcast(&requested_index, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Check if we have the requested particle
            int has_particle = 0;
            ParticleData p_data;
            for (const auto &particle : particles)
            {
                if (particle.cell_index == requested_index)
                {
                    has_particle = 1;
                    p_data = particle.toParticleData();
                    break;
                }
            }

            // Send response
            MPI_Send(&has_particle, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            if (has_particle)
            {
                MPI_Send(&p_data, 1, mpi_particle_type, 0, 2, MPI_COMM_WORLD);
            }

            return Particle();
        }
    }

    int getCollisionCount() const
    {
        return global_collisions;
    }
};

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (argc != 6)
        {
            if (rank == 0)
            {
                std::cerr << "Usage: " << argv[0] << " <seed> <side> <ncside> <n_part> <timesteps>\n";
            }
            MPI_Finalize();
            return 1;
        }

        long seed = std::stol(argv[1]);
        double side = std::stod(argv[2]);
        long ncside = std::stol(argv[3]);
        long long n_part = std::stoll(argv[4]);
        int timesteps = std::stoi(argv[5]);

        ParticleSimulation simulation(seed, side, ncside, n_part);

        double exec_time = -omp_get_wtime();

        // Setup MPI structures
        simulation.setupMPI();

        // Distribute initial particles to the appropriate processes
        simulation.distributeInitialParticles(n_part);

        simulation.run(timesteps);

        exec_time += omp_get_wtime();

        if (rank == 0)
        {
            std::vector<Particle> all_particles = simulation.gatherAllParticles();

            Particle particle0;

            for (const auto &p : all_particles)
            {
                if (p.id == 0)
                {
                    particle0 = p;
                    break;
                }
            }

            std::cout << std::fixed << std::setprecision(3)
                      << particle0.x << " " << particle0.y << std::endl;

            std::cout << simulation.getCollisionCount() << std::endl;

            std::cerr << std::fixed << std::setprecision(1) << exec_time << "s" << std::endl;
        }
        else
        {
            // Non-root processes just gather their particles
            simulation.gatherAllParticles();
        }
    } 

    MPI_Finalize();
    return 0;
}