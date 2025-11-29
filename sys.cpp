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
    Vector3D r0;            // позиция в нулевой момент времени(для среднеквадратичного смещения)
    Vector3D R;            // позиция без учета периодических граничных условий
    Vector3D r;           // Позиция 
    Vector3D v;          // Скорость   
    Vector3D Grad;      // градиент потенциала 
    double m = 1;      // Масса частицы в а.е.м   
};

class System {
    
    const double box_size = 10;      //работаем в приведенных единицах
    const double dt = 0.001;        
    const unsigned int N = 512;
        
    const double epsilon = 1;       
    const double sigma = 1;  
    const double tau = 0.1;
    const double T_target = 3;
    std::ofstream velocity_coord_file;
    std::ofstream velocity_module_file;
    std::ofstream energy_file;
    std::ofstream mean_square_offset_file;

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
        velocity_module_file.open("velocity_module.txt");
        if (!velocity_module_file.is_open()) 
        {
            std::cerr << "Ошибка открытия файла!" << std::endl;
            
        }
        energy_file.open("energy.txt");
        if (!energy_file.is_open()) 
        {
            std::cerr << "Ошибка открытия файла!" << std::endl;
            
        }
        mean_square_offset_file.open("mean_square_offset.txt");
        if (!mean_square_offset_file.is_open()) 
        {
            std::cerr << "Ошибка открытия файла!" << std::endl;
            
        }




    }
    ~System()
    {
        velocity_coord_file.close(); 
        velocity_module_file.close();
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
                    particles[number].r = particles[number].r0 = particles[number].R = Vector3D(x, y, z);
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
    std::vector<double> calculate_sr6(Vector3D r)
    {
        
        double r_sq = r.x * r.x + r.y * r.y + r.z * r.z;
        
        if(std::sqrt(r_sq) < 1e-10) 
            return {0,0};  // защита от деления на 0
        if(r_sq > 9 * sigma * sigma) 
            return {0,0};  // обрезание
        
        double r_abs = std::sqrt(r_sq);
        double inv_r = 1.0 / r_abs;
        double sr = sigma * inv_r;
        return {std::pow(sr, 6),inv_r};

    }

    Vector3D calculate_gradient(Particle p1, Particle p2)
    {   
        Vector3D r = calculate_vect_of_substance(p2, p1);
        std::vector<double> sr_6 = calculate_sr6(r);
        double sr6 = sr_6[0];  // (σ/r)⁶
        double sr12 = sr6 * sr6;    // (σ/r)¹²
        double force_mag = 24 * epsilon * (2 * sr12 - sr6) * sr_6[1];
        return r * (force_mag * sr_6[1]);
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
            particles[i].r += particles[i].v * dt +
                             old_gradients[i] * (dt * dt / (2 * particles[i].m)); 
            particles[i].R += particles[i].v * dt +
                              old_gradients[i] * (dt * dt / (2 * particles[i].m));
                                                                                   /*для R не применяем граничных условий
                                                                                    R нигде не используется, 
                                                                                    кроме поиска среднеквадратичного смещения*/
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
   
    double calculate_potential(Particle p1,Particle p2) //потенциал взаимодействия двух частиц
    {
        Vector3D r = calculate_vect_of_substance(p1, p2);
        return 4 * epsilon * (std::pow(calculate_sr6(r)[0],2) - calculate_sr6(r)[0]);
    }
    
    double calculate_potential_of_system()   //потенциал системы
    {
        double U = 0;
        for(unsigned int i = 0;i < 512; i++)
            for(unsigned int j = i + 1; j < 512; j++)
            {   
                U += calculate_potential(particles[i], particles[j]);
            }
        return U;
    }
    double calculate_E_kin(Particle p)  //кинетическая энергия одной частицы
    {
        double v_sq = p.v.x * p.v.x +
                      p.v.y * p.v.y + 
                      p.v.z * p.v.z;
        return p.m*v_sq/2;
    }
    double calculate_E_kin() //кинетическая энергия системы
    {
        double E_kin = 0;
        
        for(unsigned int i = 0; i < 512; i++)
        {   
            E_kin += calculate_E_kin(particles[i]);
        }
        return E_kin;
    }
    
    double calculate_energy() // эта функция вычисляет энергию системы
    {
        return calculate_E_kin() + calculate_potential_of_system();
    }
    double calculate_mean_square_offset()
    {
        double square_offset = 0;
        Vector3D delta_R;
        for(unsigned int i = 0; i < 512; i++)
        {
            delta_R = particles[i].r0 + particles[i].R*(-1);
            square_offset += std::pow(delta_R.x,2) + 
                                  std::pow(delta_R.y,2) +
                                  std::pow(delta_R.z,2);
        }
        return square_offset/512;

    }
        
    void write_energy_to_file()
    {
        if(!energy_file.is_open())
            return;
            
        energy_file << calculate_energy() << " " << std::endl;
        
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
    
    void write_velocity_module_to_file()
    {
        if(!velocity_module_file.is_open())
            return;
            
        for(unsigned int j = 0; j < 512 ; j++)
        {
            velocity_module_file << std::sqrt(particles[j].v.x * particles[j].v.x + 
                                   particles[j].v.y * particles[j].v.y + 
                                   particles[j].v.z * particles[j].v.z)  << " " << std::endl;
        }
    }

    
    void write_mean_square_offset_to_file()    
    {
        if(!mean_square_offset_file.is_open())
            return;
            
        mean_square_offset_file << calculate_mean_square_offset()  << " " << std::endl;
    }
    void write_everything_to_file()
    {
        write_energy_to_file();
        write_velocity_coordinate_to_file();
        write_velocity_module_to_file();
        write_mean_square_offset_to_file();
    }

};
   


int main()
{   
    const double N_termostat = 1000;
    unsigned int N_iter=10000;
    
    System system = System();

  for(unsigned int i = 0; i < N_iter; i++)
    {  
        system.verlet_scheme_iteration();

        if (i < N_termostat)
        {
            system.apply_berendsen_thermostat();
        }
        if((i % 100 == 0) && (i > N_termostat))
        {
            std::cout << i << std::endl;
            system.write_everything_to_file();
        }

    }
                
 
    return 0;
}