#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>

struct Vector3D 
{
    double x;
    double y;
    double z;

    Vector3D sum(Vector3D r1, Vector3D r2)
    {
        Vector3D result;
        result.x = r1.x + r2.x;
        result.y = r1.y + r2.y;
        result.z = r1.z + r2.z;
        return result;
    }
    Vector3D muliply(Vector3D r, double c)
    {
        Vector3D result;
        result.x = r.x * c;
        result.y = r.y * c;
        result.z = r.z * c;
        return result;
    }
    Vector3D(double x_, double y_, double z_)
    { 
        x = x_;
        y = y_;
        z = z_;

    }
    Vector3D() = default;
};


class Particle
{
public:
    Vector3D r;           // Позиция 
    Vector3D r_old; // Позиция на предыдущем шаге    
    Vector3D F; // Сила
    double m; // масса частицы в граммах   

};

class System {
    
    double box_size=1; //в мм
    //const double dt;      //длина шага по времени
public:
    std::vector<Particle> particles {100};
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

    void generate_particles(unsigned int N = 100)
    { 
        for(unsigned int i = 0; i < N; i++)
        {   
            double x = random_double(-box_size / 2, box_size / 2);
            double y = random_double(-box_size / 2, box_size / 2);
            double z = random_double(-box_size / 2, box_size / 2);
            Vector3D r_ = Vector3D(x, y, z);
            particles[i].r = particles[i].r_old = r_;
        }

    }
    
    void calculate_distance(Particle p1, Particle p2)
    {
        //между молекулой 1 и ближайшим изображением молекулы 2
    }
    
    void calculate_forces()
    {
        //to do
    }

public:   
    void verlet_scheme_iteration()
    {
        //to do
    }
    
    
};

int main()
{   
    System system = System();
    unsigned int N_iter = 100000;
    std::cout << system.particles[0].r.x ;
    for(unsigned int i = 1; i < N_iter; i++)
    {
        system.verlet_scheme_iteration();

    }
    
    return 0;
}