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

// Structure for MPI transfer of particles
struct ParticleData {
    double x, y;   // Position
    double vx, vy; // Velocity
    double m;      // Mass
    bool alive;    // Has collided
    int cell_index; // Which cell it belongs to
};

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
    void updatePositionAndVelocity(double ax, double ay) {
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
    int cell_index; // which cell it belongs to
    int proc_owner; // Process that owns this particle
    
    Particle() : x(0), y(0), vx(0), vy(0), m(0), fx(0), fy(0), alive(true), proc_owner(-1) {}

    // Conversion from MPI transfer structure
    Particle(const ParticleData& data) : 
        x(data.x), y(data.y), 
        vx(data.vx), vy(data.vy), 
        m(data.m), fx(0), fy(0), 
        alive(data.alive), 
        cell_index(data.cell_index) {}

    // Convert to MPI transfer structure
    ParticleData toParticleData() const {
        ParticleData data;
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

    double getDistance(Particle *p) {
        return sqrt(pow((x - p->x), 2) + pow((y - p->y), 2));
    }
};

class Cell {
public:
    double mx, my; // Center of Mass
    double m;      // Mass
    int x, y;      // Cell position
    bool change_flag;
    bool is_ghost;  // Is this a ghost cell?

    Cell() : mx(0), my(0), m(0), x(0), y(0), change_flag(false), is_ghost(false) {}

    void addParticle(Particle *p) {
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

// Structure to send center of mass data between processes
struct CellCOMData {
    double mx, my;
    double m;
    int x, y;
};

void Particle::calculateForceWithCell(const Cell *c) {
    double dx = c->mx - x;
    double dy = c->my - y;

    double distSquared = dx * dx + dy * dy;
    double dist = std::sqrt(distSquared);
    if (dist == 0) {
        return;
    }

    double forceMagnitude = (G * m * c->m) / distSquared;

    fx += forceMagnitude * (dx / dist);
    fy += forceMagnitude * (dy / dist);
}

void Particle::calculateForceBetweenParticles(Particle *p2) {
    double dx = p2->x - x;
    double dy = p2->y - y;

    double distSquared = dx * dx + dy * dy;
    double dist = std::sqrt(distSquared);
    if (dist == 0) {
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

void Particle::applyForce(double sidelen, long grid_size, std::vector<Cell>& cells,
                         int row_start, int row_end, int proc_rows, int rank) {
    if (m == 0) {
        fx = 0;
        fy = 0;
        return;
    }
    
    double ax = fx / m;
    double ay = fy / m;

    int cell_y = static_cast<int>(y / (sidelen / grid_size));
    int cell_x = static_cast<int>(x / (sidelen / grid_size));

    cell_index = (cell_y - row_start) * grid_size + cell_x;
    if (cell_y < 0 || cell_y >= grid_size) {
        cell_index = -1; // Mark for migration
    }
        
    updatePositionAndVelocity(ax, ay);
    
    // Wrap around the simulation space
    x = fmod((x + sidelen), sidelen);
    y = fmod((y + sidelen), sidelen);

    cell_y = static_cast<int>(y / (sidelen / grid_size));
    cell_x = static_cast<int>(x / (sidelen / grid_size));

    // Check if particle is now in this process's domain
    if (cell_y >= row_start && cell_y < row_end) {
        int actual_cell_index = (cell_y - row_start) * grid_size + cell_x;

        if (cell_index != actual_cell_index && cell_index >= 0) {
            cells[cell_index].change_flag = true;
            cell_index = actual_cell_index;
        }
    } else {
        // Particle now belongs to another process
        proc_owner = cell_y / proc_rows;
        cell_index = -1; // Mark for migration
    }

    // Reset forces for next iteration
    fx = 0;
    fy = 0;
}

class ParticleSimulation {
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

    void initMPITypes() {
        // Create MPI datatype for ParticleData
        MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL, MPI_INT};
        int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
        MPI_Aint offsets[7];
        
        offsets[0] = offsetof(ParticleData, x);
        offsets[1] = offsetof(ParticleData, y);
        offsets[2] = offsetof(ParticleData, vx);
        offsets[3] = offsetof(ParticleData, vy);
        offsets[4] = offsetof(ParticleData, m);
        offsets[5] = offsetof(ParticleData, alive);
        offsets[6] = offsetof(ParticleData, cell_index);
        
        MPI_Type_create_struct(7, blocklengths, offsets, types, &mpi_particle_type);
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
        // Initialize MPI
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Create MPI datatypes
        initMPITypes();
        
        // Divide grid rows among processes
        rows_per_proc = ncside / num_procs;
        if (rows_per_proc == 0) {
            rows_per_proc = 1;
        }
        
        my_row_start = rank * rows_per_proc;
        my_row_end = (rank == num_procs - 1) ? ncside : (rank + 1) * rows_per_proc;
        
        // Allocate cells for my portion of the grid, plus ghost cells
        int local_grid_size = (my_row_end - my_row_start) * ncside;
        cells.resize(local_grid_size);
        cellParticles.resize(local_grid_size);
        
        // Initialize particles if root process
        if (rank == 0) {
            particles.resize(n_part);
            init_particles();
        }
        
        // Distribute particles to appropriate processes
        distributeParticles(n_part);
    }


    ~ParticleSimulation() {
        // Free MPI datatypes BEFORE MPI_Finalize is called
        if (mpi_particle_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&mpi_particle_type);
            mpi_particle_type = MPI_DATATYPE_NULL;
        }
        if (mpi_cell_com_type != MPI_DATATYPE_NULL) {
            MPI_Type_free(&mpi_cell_com_type);
            mpi_cell_com_type = MPI_DATATYPE_NULL;
        }
    }
    

    void init_particles() {
        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].x = rng.getRandom01() * side_length;
            particles[i].y = rng.getRandom01() * side_length;
            particles[i].vx = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].vy = (rng.getRandom01() - 0.5) * side_length / grid_size / 5.0;
            particles[i].m = rng.getRandom01() * 0.01 * (grid_size * grid_size) /
                           particles.size() / G * EPSILON2;
                           
            // Determine which process should own this particle
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));
            particles[i].proc_owner = cell_y / rows_per_proc;
            if (particles[i].proc_owner >= num_procs) {
                particles[i].proc_owner = num_procs - 1;
            }
        }
    }
    
    void distributeParticles(long long n_part) {
        if (rank == 0) {
            // Count particles for each process
            std::vector<int> particles_per_proc(num_procs, 0);
            for (const auto& p : particles) {
                particles_per_proc[p.proc_owner]++;
            }
            
            // Send counts to all processes
            MPI_Bcast(particles_per_proc.data(), num_procs, MPI_INT, 0, MPI_COMM_WORLD);
            
            // Send particles to appropriate processes
            for (int dest_rank = 1; dest_rank < num_procs; dest_rank++) {
                std::vector<ParticleData> particles_to_send;
                for (const auto& p : particles) {
                    if (p.proc_owner == dest_rank) {
                        particles_to_send.push_back(p.toParticleData());
                    }
                }
                
                MPI_Send(particles_to_send.data(), particles_to_send.size(), 
                         mpi_particle_type, dest_rank, 0, MPI_COMM_WORLD);
            }
            
            // Keep only particles for rank 0
            std::vector<Particle> my_particles;
            for (const auto& p : particles) {
                if (p.proc_owner == 0) {
                    my_particles.push_back(p);
                }
            }
            particles = std::move(my_particles);
        } else {
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
            for (int i = 0; i < my_particle_count; i++) {
                particles[i] = Particle(received_particles[i]);
            }
        }
    }
    
    void updateCellParticles() {
        for (int i = 0; i < cellParticles.size(); i++) {
            if (cells[i].change_flag) { 
                // Collect particles that need to move
                std::vector<std::pair<Particle*, int>> particlesToMove;
                
                for (int k = 0; k < cellParticles[i].size(); k++) {
                    int cell_index = cellParticles[i][k]->cell_index;
                    if (cell_index != i && cell_index >= 0) {
                        particlesToMove.push_back({cellParticles[i][k], cell_index});
                    }
                }
                
                // Move particles to new cells
                for (auto& pair : particlesToMove) {
                    cellParticles[pair.second].push_back(pair.first);
                }
                
                // Remove moved particles from original cell
                if (!particlesToMove.empty()) {
                    auto newEnd = std::remove_if(cellParticles[i].begin(), cellParticles[i].end(),
                        [i](Particle* p) { return p->cell_index != i && p->cell_index >= 0; });
                    cellParticles[i].erase(newEnd, cellParticles[i].end());
                }
                
                cells[i].change_flag = false;
            }
        }
    }

    void exchangeParticles() {
        // Particles that need to be sent to other processes
        std::vector<std::vector<ParticleData>> particles_to_send(num_procs);
        
        // Identify particles to send to other processes
        std::vector<Particle*> particles_to_remove;
        for (auto& p : particles) {
            if (p.cell_index == -1 && p.proc_owner != rank) {
                particles_to_send[p.proc_owner].push_back(p.toParticleData());
                particles_to_remove.push_back(&p);
            }
        }
        
        // Exchange particle counts
        std::vector<int> send_counts(num_procs);
        for (int i = 0; i < num_procs; i++) {
            send_counts[i] = particles_to_send[i].size();
        }
        
        std::vector<int> recv_counts(num_procs);
        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        
        // Exchange particles
        std::vector<MPI_Request> requests;
        for (int i = 0; i < num_procs; i++) {
            if (i != rank && send_counts[i] > 0) {
                MPI_Request req;
                MPI_Isend(particles_to_send[i].data(), send_counts[i], 
                         mpi_particle_type, i, 0, MPI_COMM_WORLD, &req);
                requests.push_back(req);
            }
        }
        
        // Receive particles from other processes
        std::vector<ParticleData> received_particles;
        for (int i = 0; i < num_procs; i++) {
            if (i != rank && recv_counts[i] > 0) {
                std::vector<ParticleData> temp_received(recv_counts[i]);
                MPI_Status status;
                MPI_Recv(temp_received.data(), recv_counts[i], 
                         mpi_particle_type, i, 0, MPI_COMM_WORLD, &status);
                
                received_particles.insert(received_particles.end(), temp_received.begin(), temp_received.end());
            }
        }
        
        // Wait for all sends to complete
        if (!requests.empty()) {
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        }
        
        // Remove particles that were sent to other processes
        particles.erase(
            std::remove_if(particles.begin(), particles.end(), 
                [&particles_to_remove](const Particle& p) {
                    return std::find(particles_to_remove.begin(), particles_to_remove.end(), &p) != particles_to_remove.end();
                }
            ),
            particles.end()
        );
        
        // Add received particles
        for (const auto& p_data : received_particles) {
            Particle new_particle(p_data);
            
            // Update cell index for new domain
            int cell_y = static_cast<int>(new_particle.y / (side_length / grid_size));
            int cell_x = static_cast<int>(new_particle.x / (side_length / grid_size));
            new_particle.cell_index = (cell_y - my_row_start) * grid_size + cell_x;
            new_particle.proc_owner = rank;
            
            particles.push_back(new_particle);
        }
    }

    void updateCOM() {
        // Reset cells and cellParticles
        cellParticles.assign((my_row_end - my_row_start) * grid_size, std::vector<Particle *>{});
        cells.assign((my_row_end - my_row_start) * grid_size, Cell{});
        
        for (size_t i = 0; i < particles.size(); i++) {
            // Calculate cell index
            int cell_x = static_cast<int>(particles[i].x / (side_length / grid_size));
            int cell_y = static_cast<int>(particles[i].y / (side_length / grid_size));
            
            // Skip if particle is outside my domain
            if (cell_y < my_row_start || cell_y >= my_row_end) {
                continue;
            }
            
            // Convert 2D cell index to local 1D index
            int local_cell_y = cell_y - my_row_start;
            int cell_index = local_cell_y * grid_size + cell_x;
            
            if (cell_x < 0 || cell_x >= grid_size || local_cell_y < 0 || local_cell_y >= (my_row_end - my_row_start)) {
                std::cout << "[PANIC2] Cell out of bounds at rank " << rank << std::endl;
                continue;
            }
            
            particles[i].cell_index = cell_index; // set particle cell index
            cellParticles[cell_index].push_back(&particles[i]);
            
            cells[cell_index].addParticle(&particles[i]);
            cells[cell_index].x = cell_x;
            cells[cell_index].y = cell_y;
        }
        
        // Exchange cell center of mass information for boundary cells
        exchangeBoundaryCellInfo();
    }
    
    void exchangeBoundaryCellInfo() {
        // Validate process boundaries
        if (rank < 0 || rank >= num_procs) {
            std::cerr << "Rank " << rank << " out of bounds" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    
        // Collect COM data for top and bottom rows
        std::vector<CellCOMData> top_row_data(grid_size), bottom_row_data(grid_size);
        
        // Prepare data for boundary exchange
        try {
            // If not the first process, prepare top row to send to previous process
            if (rank > 0) {
                for (int x = 0; x < grid_size; x++) {
                    int cell_index = x; // First row in my domain
                    if (cell_index >= 0 && cell_index < static_cast<int>(cells.size())) {
                        top_row_data[x].mx = cells[cell_index].mx;
                        top_row_data[x].my = cells[cell_index].my;
                        top_row_data[x].m = cells[cell_index].m;
                        top_row_data[x].x = cells[cell_index].x;
                        top_row_data[x].y = cells[cell_index].y;
                    }
                }
            }
            
            // If not the last process, prepare bottom row to send to next process
            if (rank < num_procs - 1) {
                int last_row = (my_row_end - my_row_start - 1) * grid_size;
                for (int x = 0; x < grid_size; x++) {
                    int cell_index = last_row + x;
                    if (cell_index >= 0 && cell_index < static_cast<int>(cells.size())) {
                        bottom_row_data[x].mx = cells[cell_index].mx;
                        bottom_row_data[x].my = cells[cell_index].my;
                        bottom_row_data[x].m = cells[cell_index].m;
                        bottom_row_data[x].x = cells[cell_index].x;
                        bottom_row_data[x].y = cells[cell_index].y;
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in data preparation: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    
        // Safer communication with explicit synchronization
        std::vector<MPI_Request> requests;
        std::vector<MPI_Status> statuses;
        
        // Receive buffers
        std::vector<CellCOMData> ghost_bottom_row(grid_size), ghost_top_row(grid_size);
        
        try {
            // Send/Receive operations with explicit error checking
            if (rank > 0) {
                MPI_Request send_req, recv_req;
                MPI_Status send_status, recv_status;
                
                // Non-blocking send to previous process
                int send_err = MPI_Isend(top_row_data.data(), grid_size, mpi_cell_com_type, 
                                         rank - 1, 0, MPI_COMM_WORLD, &send_req);
                
                // Non-blocking receive from previous process
                int recv_err = MPI_Irecv(ghost_bottom_row.data(), grid_size, mpi_cell_com_type, 
                                         rank - 1, 1, MPI_COMM_WORLD, &recv_req);
                
                if (send_err != MPI_SUCCESS || recv_err != MPI_SUCCESS) {
                    std::cerr << "MPI communication error at rank " << rank << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                // Wait for both send and receive
                MPI_Wait(&send_req, &send_status);
                MPI_Wait(&recv_req, &recv_status);
            }
            
            // Repeat for bottom row communication with next process
            if (rank < num_procs - 1) {
                MPI_Request send_req, recv_req;
                MPI_Status send_status, recv_status;
                
                // Non-blocking send to next process
                int send_err = MPI_Isend(bottom_row_data.data(), grid_size, mpi_cell_com_type, 
                                         rank + 1, 1, MPI_COMM_WORLD, &send_req);
                
                // Non-blocking receive from next process
                int recv_err = MPI_Irecv(ghost_top_row.data(), grid_size, mpi_cell_com_type, 
                                         rank + 1, 0, MPI_COMM_WORLD, &recv_req);
                
                if (send_err != MPI_SUCCESS || recv_err != MPI_SUCCESS) {
                    std::cerr << "MPI communication error at rank " << rank << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                // Wait for both send and receive
                MPI_Wait(&send_req, &send_status);
                MPI_Wait(&recv_req, &recv_status);
            }
        } catch (const std::exception& e) {
            std::cerr << "Communication error at rank " << rank << ": " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    
        // Add ghost cells to local grid with additional safety checks
        try {
            // For bottom ghost row (coming from previous process)
            if (rank > 0) {
                for (int x = 0; x < grid_size; x++) {
                    if (x < static_cast<int>(ghost_bottom_row.size())) {
                        Cell ghost_cell;
                        ghost_cell.mx = ghost_bottom_row[x].mx;
                        ghost_cell.my = ghost_bottom_row[x].my - side_length;
                        ghost_cell.m = ghost_bottom_row[x].m;
                        ghost_cell.x = ghost_bottom_row[x].x;
                        ghost_cell.y = ghost_bottom_row[x].y;
                        ghost_cell.is_ghost = true;
                        
                        cells.push_back(ghost_cell);
                    }
                }
            }
            
            // For top ghost row (coming from next process)
            if (rank < num_procs - 1) {
                for (int x = 0; x < grid_size; x++) {
                    if (x < static_cast<int>(ghost_top_row.size())) {
                        Cell ghost_cell;
                        ghost_cell.mx = ghost_top_row[x].mx;
                        ghost_cell.my = ghost_top_row[x].my + side_length;
                        ghost_cell.m = ghost_top_row[x].m;
                        ghost_cell.x = ghost_top_row[x].x;
                        ghost_cell.y = ghost_top_row[x].y;
                        ghost_cell.is_ghost = true;
                        
                        cells.push_back(ghost_cell);
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Ghost cell addition error at rank " << rank << ": " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    void updateForces() {
        #pragma omp parallel for
        for (int i = 0; i < (my_row_end - my_row_start) * grid_size; i++) {
            std::vector<Cell> temp_cells;
            temp_cells.reserve(8);

            int cell_x = cells[i].x;
            int cell_y = cells[i].y;

            // Calculate neighboring cells
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (dx == 0 && dy == 0) continue; // Skip the center cell

                    int neighbor_x = (cell_x + dx + grid_size) % grid_size;
                    int neighbor_y = cell_y + dy;
                    
                    Cell temp_cell;
                    
                    // Handle wrap-around for x
                    if (neighbor_x < 0) temp_cell.mx -= side_length;
                    else if (neighbor_x >= grid_size) temp_cell.mx += side_length;
                    
                    // Handle y-direction neighbors, which might be in ghost cells
                    if (neighbor_y < my_row_start) {
                        // Use ghost cells from previous process or wrap around
                        if (rank > 0) {
                            // Look in ghost cells
                            int ghost_index = (my_row_end - my_row_start) * grid_size + neighbor_x;
                            if (ghost_index < cells.size()) {
                                temp_cell.mx = cells[ghost_index].mx;
                                temp_cell.my = cells[ghost_index].my;
                                temp_cell.m = cells[ghost_index].m;
                                temp_cells.push_back(temp_cell);
                            }
                        } else {
                            // Wrap around to bottom of simulation
                            int wrap_y = grid_size - 1;
                            int wrap_index = (wrap_y - my_row_start) * grid_size + neighbor_x;
                            if (wrap_index >= 0 && wrap_index < (my_row_end - my_row_start) * grid_size) {
                                temp_cell.mx = cells[wrap_index].mx;
                                temp_cell.my = cells[wrap_index].my - side_length;
                                temp_cell.m = cells[wrap_index].m;
                                temp_cells.push_back(temp_cell);
                            }
                        }
                    } else if (neighbor_y >= my_row_end) {
                        // Use ghost cells from next process or wrap around
                        if (rank < num_procs - 1) {
                            // Look in ghost cells
                            int ghost_offset = (my_row_end - my_row_start) * grid_size + grid_size + neighbor_x;
                            if (ghost_offset < cells.size()) {
                                temp_cell.mx = cells[ghost_offset].mx;
                                temp_cell.my = cells[ghost_offset].my;
                                temp_cell.m = cells[ghost_offset].m;
                                temp_cells.push_back(temp_cell);
                            }
                        } else {
                            // Wrap around to top of simulation
                            int wrap_y = 0;
                            int wrap_index = (wrap_y - my_row_start) * grid_size + neighbor_x;
                            if (wrap_index >= 0 && wrap_index < (my_row_end - my_row_start) * grid_size) {
                                temp_cell.mx = cells[wrap_index].mx;
                                temp_cell.my = cells[wrap_index].my + side_length;
                                temp_cell.m = cells[wrap_index].m;
                                temp_cells.push_back(temp_cell);
                            }
                        }
                    } else {
                        // Regular neighbor within this process's domain
                        int neighbor_index = (neighbor_y - my_row_start) * grid_size + neighbor_x;
                        if (neighbor_index >= 0 && neighbor_index < (my_row_end - my_row_start) * grid_size) {
                            temp_cell.mx = cells[neighbor_index].mx;
                            temp_cell.my = cells[neighbor_index].my;
                            temp_cell.m = cells[neighbor_index].m;
                            temp_cells.push_back(temp_cell);
                        }
                    }
                }
            }

            // Calculate forces between particles in the same cell
            #pragma omp parallel for
            for (size_t j = 0; j < cellParticles[i].size(); j++) {
                for (size_t k = j + 1; k < cellParticles[i].size(); k++) {
                    if (!cellParticles[i][j]->alive || !cellParticles[i][k]->alive)
                        continue;
                        
                    cellParticles[i][j]->calculateForceBetweenParticles(cellParticles[i][k]);
                }
            }

            // Calculate forces with neighboring cells
            #pragma omp parallel for
            for (size_t j = 0; j < cellParticles[i].size(); j++) {
                if (!cellParticles[i][j]->alive)
                    continue;
                    
                for (const auto& neighbor_cell : temp_cells) {
                    if (neighbor_cell.m > 0) {
                        cellParticles[i][j]->calculateForceWithCell(&neighbor_cell);
                    }
                }
            }
        }
    }

    void checkCollisions() {
        for (int i = 0; i < (my_row_end - my_row_start) * grid_size; i++) {
            for (size_t j = 0; j < cellParticles[i].size(); j++) {
                if (!cellParticles[i][j]->alive)
                    continue;
                    
                for (size_t k = j + 1; k < cellParticles[i].size(); k++) {
                    if (!cellParticles[i][k]->alive)
                        continue;
                        
                    double distance = cellParticles[i][j]->getDistance(cellParticles[i][k]);
                    if (distance <= EPSILON) {
                        cellParticles[i][j]->alive = false;
                        cellParticles[i][k]->alive = false;
                        local_collisions++;
                    }
                }
            }
        }
        
        // Remove collided particles
        particles.erase(
            std::remove_if(particles.begin(), particles.end(), 
                [](const Particle& p) { return !p.alive; }),
            particles.end()
        );
    }

    void updateParticlePositions() {
        for (auto& p : particles) {
            if (!p.alive) continue;
            
            p.applyForce(side_length, grid_size, cells, my_row_start, my_row_end, rows_per_proc, rank);
        }
    }

    void run(int timesteps) {
        for (int step = 0; step < timesteps; step++) {
            // Step 1: Determine center of mass for each cell
            updateCOM();
            
            // Step 2: Calculate forces
            updateForces();
            
            // Step 3: Update particle positions and velocities
            updateParticlePositions();
            
            // Step 4: Exchange particles between processes
            exchangeParticles();
            
            // Step 5: Update cell particle lists after movement
            updateCellParticles();
            
            // Step 6: Check for collisions
            checkCollisions();
        }
        
        // Collect total collisions from all processes
        MPI_Reduce(&local_collisions, &global_collisions, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    Particle getParticle(int index) const {
        if (rank == 0) {
            // If the index is within local particles, return it directly
            if (index < particles.size()) {
                return particles[index];
            }
            
            // If not found locally, do a global search
            Particle result;
            bool found = false;
            
            // Broadcast the index to search for
            int local_found = 0;
            MPI_Bcast(&index, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            // Gather information about which process has the particle
            for (int i = 1; i < num_procs; i++) {
                int has_particle = 0;
                
                // Send request to other processes
                MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Recv(&has_particle, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                if (has_particle) {
                    ParticleData p_data;
                    MPI_Recv(&p_data, 1, mpi_particle_type, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    result = Particle(p_data);
                    found = true;
                    break;
                }
            }
            
            // If not found in any process, return a default particle
            if (!found) {
                result.x = 0.0;
                result.y = 0.0;
            }
            
            return result;
        } else {
            // Non-root processes handle particle search
            int requested_index;
            MPI_Bcast(&requested_index, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            // Check if we have the requested particle
            int has_particle = 0;
            ParticleData p_data;
            for (const auto& particle : particles) {
                if (particle.cell_index == requested_index) {
                    has_particle = 1;
                    p_data = particle.toParticleData();
                    break;
                }
            }
            
            // Send response
            MPI_Send(&has_particle, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            if (has_particle) {
                MPI_Send(&p_data, 1, mpi_particle_type, 0, 2, MPI_COMM_WORLD);
            }
            
            return Particle(); // Non-root processes return empty particle
        }
    }

    int getCollisionCount() const {
        return global_collisions;
    }
};

int main(int argc, char *argv[]) {
    // Initialize MPI first
    MPI_Init(&argc, &argv);
    
    {  // Add a scope block to ensure destructor is called before MPI_Finalize
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Check command line arguments
        if (argc != 6) {
            if (rank == 0) {
                std::cerr << "Usage: " << argv[0] << " <seed> <side> <ncside> <n_part> <timesteps>\n";
            }
            MPI_Finalize();
            return 1;
        }

        // Parse command line arguments
        long seed = std::stol(argv[1]);
        double side = std::stod(argv[2]);
        long ncside = std::stol(argv[3]);
        long long n_part = std::stoll(argv[4]);
        int timesteps = std::stoi(argv[5]);

        // Run simulation
        double exec_time;
        ParticleSimulation simulation(seed, side, ncside, n_part);
        
        exec_time = -omp_get_wtime();
        simulation.run(timesteps);
        exec_time += omp_get_wtime();
        
        // Output results
        if (rank == 0) {
            // Print all particles
            for (int i = 0; i < n_part; i++) {
                Particle p = simulation.getParticle(i);
                std::cout << std::fixed << std::setprecision(3) << p.x << " " << p.y << std::endl;
            }
        
            // Get particle 0 specifically
            Particle p0 = simulation.getParticle(0);
            
            // Print particle 0 coordinates
            std::cout << std::fixed << std::setprecision(3) << p0.x << " " << p0.y << std::endl;
            
            // Print collision count
            std::cout << simulation.getCollisionCount() << std::endl;
            
            // Print execution time to stderr
            std::cerr << std::fixed << std::setprecision(1) << exec_time << "s" << std::endl;
        }
    }  // Simulation object destructor called here

    // Finalize MPI after the simulation object is destroyed
    MPI_Finalize();
    return 0;
}