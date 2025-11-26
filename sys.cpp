#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>

struct Vector3D 
{
    double x {0};
    double y {0};
    double z {0};
   
    Vector3D(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} {}
    
    Vector3D() = default;
    
    Vector3D operator*(double scalar) const 
    {
        return {x * scalar, y * scalar, z * scalar};
    }
    
    Vector3D& operator+=(const Vector3D& vect)
    {
        x += vect.x;
        y += vect.y;
        z += vect.z; 
        return *this;  
    } 
    
    Vector3D& operator*=(double scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar; 
        return *this;  
    } 
};

Vector3D operator+(const Vector3D& v1, const Vector3D& v2) 
{   
    return {Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)};
}

class Particle
{
public:
    Vector3D r;           // Позиция 
    Vector3D v; // Скорость   
    Vector3D Grad; // градиент потенциала 
    double m = 1; // Масса частицы в а.е.м   
};

class System {
    
    const double box_size = 10;      //работаем в приведенных единицах
    const double dt = 0.001;        
    const unsigned int N = 512;
        
    const double epsilon = 1;       
    const double sigma = 1;  
    const double tau = 1;
    const double T_target = 1.7;
    std::ofstream velocity_coord_file;

    std::vector<Particle> particles {N};

public:    
    System()
    {
        generate_particles();
        velocity_coord_file.open("velocity_coord.txt");
        if (!velocity_coord_file.is_open()) 
        {
            std::cerr << "Ошибка открытия файла!" << std::endl;
            
        }
    }
    ~System()
    {
        velocity_coord_file.close();  
    }

private:    
    double random_double(double min, double max) 
    {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(min, max);
        return dist(gen);
    }

    void generate_particles()
    { 
        unsigned int number = 0;
        double x, y, z;
        double cube_size = box_size / 8;
        for(unsigned int i = 0; i < 8; i++) 
            for(unsigned int j = 0; j < 8; j++)
                for(unsigned int k = 0; k < 8; k++)
                {
                    x = cube_size / 2 + i * cube_size;
                    y = cube_size / 2 + j * cube_size;              //решетчатая инициализация начальных положений
                    z = cube_size / 2 + k * cube_size;
                    particles[number].r = Vector3D(x, y, z);
                    particles[number].v = Vector3D(0, 0, 0);
                    number++;
                }
    }

    double calculate_coordinate_of_substance(double x1, double x2)
    {
        double dx = x2 - x1;
        dx -= box_size * std::round(dx / box_size); 
        return dx;
    }
    
    Vector3D calculate_vect_of_substance(Particle p1, Particle p2)
    { // кладем эту разность в вектор
        Vector3D dr;
        dr.x = calculate_coordinate_of_substance(p1.r.x, p2.r.x);
        dr.y = calculate_coordinate_of_substance(p1.r.y, p2.r.y);
        dr.z = calculate_coordinate_of_substance(p1.r.z, p2.r.z);
        return dr;
    }

    Vector3D calculate_gradient(Particle p1, Particle p2)
    {
        Vector3D r = calculate_vect_of_substance(p2, p1);
        double r_sq = r.x * r.x + r.y * r.y + r.z * r.z;
        
        if(std::sqrt(r_sq) < 1e-10) 
            return Vector3D(0, 0, 0);  // защита от деления на 0
        if(r_sq > 9 * sigma * sigma) 
            return Vector3D(0, 0, 0);  // обрезание
        
        double r_abs = std::sqrt(r_sq);
        double inv_r = 1.0 / r_abs;
        double sr = sigma * inv_r;
        double sr6 = std::pow(sr, 6);  // (σ/r)⁶
        double sr12 = sr6 * sr6;    // (σ/r)¹²
        double force_mag = 24 * epsilon * (2 * sr12 - sr6) * inv_r;
        return r * (force_mag * inv_r);
    }
    
    
    double periodic_conditions_on_coordinats(double x)
    {   
        x = std::fmod(x + box_size / 2, box_size);
        if (x < 0)
            x += box_size;
        return x - box_size / 2;
    }
        
    Vector3D periodic_conditions_on_vect(Vector3D r)
    {   
        return Vector3D(periodic_conditions_on_coordinats(r.x), periodic_conditions_on_coordinats(r.y), periodic_conditions_on_coordinats(r.z));
    }
     
    double calculate_temperature()
    {   
        // Сначала удаляем движение центра масс
        remove_center_of_mass_motion();
        
        double E_kin = 0;
        for(unsigned int i = 0; i < N; i++)
        {   
            double sq_v_abs = particles[i].v.x * particles[i].v.x 
                            + particles[i].v.y * particles[i].v.y 
                            + particles[i].v.z * particles[i].v.z;
            E_kin += sq_v_abs * particles[i].m / 2;
        }
        
        // Число степеней свободы: 3N - 3 (после удаления v_cm)
        return 2 * E_kin / (3 * N - 3);
    }

    void remove_center_of_mass_motion()
    {
        Vector3D v_cm(0, 0, 0);
        double total_mass = 0;
    
        // Вычисляем скорость центра масс
        for(unsigned int i = 0; i < N; i++) {
            v_cm += particles[i].v * particles[i].m;
            total_mass += particles[i].m;
        }
        v_cm *= (1.0 / total_mass);
        
        // Вычитаем скорость центра масс
        for(unsigned int i = 0; i < N; i++) {
            particles[i].v += v_cm * (-1);
        }
    }

public: 
    void apply_berendsen_thermostat()
    {
        double T_current = calculate_temperature();
        
        // Защита от деления на ноль и отрицательных температур
        if(T_current <= 0) 
            return;
        
        double lambda = sqrt(1 + (dt / tau) * (T_target / T_current - 1));
        
        // Применяем ко всем частицам
        for(unsigned int i = 0; i < N; i++) 
        {
            particles[i].v *= lambda;
        }
    }

    void verlet_scheme_iteration()
    {   
        // Сохраняем старые силы
        std::vector<Vector3D> old_gradients(N);
        for(unsigned i = 0; i < N; i++) 
        {
            old_gradients[i] = particles[i].Grad;
        }
        
        for(unsigned i = 0; i < N; i++) 
        {
            particles[i].r = particles[i].r + particles[i].v * dt +
                             old_gradients[i] * (dt * dt / (2 * particles[i].m));
            particles[i].r = periodic_conditions_on_vect(particles[i].r);
        }
        
        std::vector<Vector3D> new_gradients(N, Vector3D(0, 0, 0));
        for(unsigned i = 0; i < N; i++) 
        {
            for(unsigned j = i + 1; j < N; j++) 
            {   
                // Оптимизация: j = i+1
                Vector3D grad_ij = calculate_gradient(particles[i], particles[j]);
                new_gradients[i] += grad_ij;
                new_gradients[j] += grad_ij * (-1); // 3-й закон Ньютона
            }
            particles[i].Grad = new_gradients[i];
        }
        
        for(unsigned i = 0; i < N; i++) 
        {
            particles[i].v += (old_gradients[i] + new_gradients[i]) *
                              (dt / (2 * particles[i].m));
        }
    }
    void write_velocity_coordinate_to_file()
    {
        if(!velocity_coord_file.is_open())
            return;
            
        for(unsigned int j = 0; j < 512 ; j++)
        {
            velocity_coord_file << particles[j].v.x  << " " << std::endl;
        }
    }
};
   


int main()
{   
    const double N_termostat = 1000;
    unsigned int N_iter = 10000;
    System system = System();
    
    for(unsigned int i = 0; i < N_iter; i++)
    {  
        system.verlet_scheme_iteration();

        if (i < N_termostat)
        {
            system.apply_berendsen_thermostat();
        }
        if((i % 5000 == 0) && (i > N_termostat))
        {
            std::cout << i << std::endl;
            system.write_velocity_coordinate_to_file();
        }

    }
    
 
    return 0;
}