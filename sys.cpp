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
   
    Vector3D(double x_, double y_, double z_)
    { 
        x = x_;
        y = y_;
        z = z_;

    }
    
    Vector3D() = default;
    
    Vector3D operator*(double scalar) const 
    {
        return {x * scalar, y * scalar, z * scalar};
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
    const double dt = 0.01;        
    const unsigned int N = 100;
        
    const double epsilon = 1;       
    const double sigma = 1;  
    const double tau = 1;
    const double T_target = 1.7;


public:
    std::vector<Particle> particles {N};
    
    System()
    {
        generate_particles();
    };
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
        for(unsigned int i = 0; i < N; i++)
        {   
            double x = random_double(-box_size / 2, box_size / 2);
            double y = random_double(-box_size / 2, box_size / 2);
            double z = random_double(-box_size / 2, box_size / 2);
            Vector3D r_ = Vector3D(x, y, z);
            particles[i].r = r_;
            particles[i].v = Vector3D(0, 0, 0);
        }

    }
    double calculate_coordinate_of_substance(double x1, double x2)
    {  // считаем разность между координатами частицы с учетом периодических граничных условий

        //if (x2 >= x1)
        //{
           // if(x2 - x1 >= box_size / 2)
              //  return x2 - x1;
           // else
               // return -(box_size - x2 + x1);
      //  }
        //else
       // {
            //if(x2 - x1 <= -box_size / 2)
               // return box_size - x1 + x2;
            //else    
                //return x2 - x1;
        //}
        double dx = x2 - x1;
        dx -= box_size*std::round(dx/box_size); 
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
    
    Vector3D calculate_gradient(Particle p1, Particle p2) //считаем силу как градиент потенциала
    {
        Vector3D r = calculate_vect_of_substance(p2, p1); //именно в таком порядке
        double r_abs = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
        if((r_abs > 3 * sigma)||(r_abs < 1e-16)) return Vector3D(0,0,0);
        double sr = sigma/r_abs;
        double sr6 = sr * sr * sr * sr * sr * sr;  // (σ/r)^6
        double sr12 = sr6 * sr6;
        Vector3D r_norm = r*(1/r_abs);
        return r_norm*(24 * epsilon * (2 * sr12 - sr6) / r_abs) ;  //возвращает градиент потенциала от частицы 2 к частице 1



    }
    
    double periodic_conditions_on_coordinats(double x)
    {   
        x += box_size / 2;
        if((x  > box_size) || (x < 0))
        {
            x -= (floor(x / box_size)) * box_size;
        }
        return x - box_size / 2;
    }
    
    Vector3D periodic_conditions_on_vect(Vector3D r)
    {   
        return Vector3D(periodic_conditions_on_coordinats(r.x), periodic_conditions_on_coordinats(r.y), periodic_conditions_on_coordinats(r.z));
    }
    double calculate_temperature()
    {   
        double E_kin = 0;
        for(unsigned int i = 0;i < N;i++ )
        {   
            double sq_v_abs = particles[i].v.x * particles[i].v.x + particles[i].v.y * particles[i].v.y + particles[i].v.z * particles[i].v.z;
            E_kin += sq_v_abs*particles[i].m/2;
        }
        return 2 * E_kin / (3 * N - 3);
    }
    Vector3D berendsen_termostat_vector(Vector3D V)
    {
        double lambda = sqrt(1 + (dt/tau)*(T_target/calculate_temperature()-1));
        return V*lambda;
    }  

public: 
 
    void verlet_scheme_iteration()
    {   
        std::vector<Vector3D> Gradients{N}; //массив градиентов для вспех частиц
        for(unsigned int i = 0; i < N; i++)
        {   
            for(unsigned int j = 0; j < N; j++)
            { 
                if(i == j) continue;
                Gradients[i] = Gradients[i] + calculate_gradient(particles[i], particles[j]);
            }
            
            
        }
        for(unsigned int i = 0; i < N; i++)
        {   particles[i].Grad = Gradients[i];
            particles[i].r = particles[i].r + particles[i].v * dt + Gradients[i] * (-dt * dt / (2 * particles[i].m));
            particles[i].v = particles[i].v + Gradients[i] * (-dt/particles[i].m);
            particles[i].r = periodic_conditions_on_vect(particles[i].r);
        }
    }
    void berendsen_termostat()
    {
        for(unsigned i = 0; i < N; i++)
        {
            particles[i].v = berendsen_termostat_vector(particles[i].v);
        }
    }
};

int main()
{   
    const double N_termostat = 1000;
    unsigned int N_iter = 1000;
    System system = System();
    std::ofstream outfile;
    outfile.open("output.txt");
    std::cout << system.particles[0].r.x ;
    if (!outfile.is_open()) {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        return 1;
    }
    for(unsigned int i = 0; i < N_iter; i++)
    {  
        system.verlet_scheme_iteration();
        if (i==0)
        {
           
        }

        if (i<N_termostat)
        {
            system.berendsen_termostat();
        }
        
        std::cout << i << std::endl;

    }
    
    // Запись данных в файл
             // Запись данных в файл
    for(unsigned int j = 0; j <100 ; j++)
        {
            outfile << system.particles[j].v.x << " " << std::endl;
        }
    outfile.close();   
    

   
    return 0;
}