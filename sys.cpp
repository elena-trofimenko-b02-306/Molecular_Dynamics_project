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
    double m = 1; // Масса частицы в граммах   

};

class System {
    
    const double box_size = 1;      //в мм
    const double dt = 0.01;        //длина шага по времени в секундах
    const unsigned int N = 100;    // число частиц 
    double epsilon = 1;
    double sigma = 1;   

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

        if (x2 >= x1)
        {
            if(x2 - x1 >= box_size / 2)
                return x2 - x1;
            else
                return -(box_size - x2 + x1);
        }
        else
        {
            if(x2 - x1 <= -box_size / 2)
                return box_size - x1 + x2;
            else    
                return x2 - x1;
        }


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
        double r_abs = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
        return r * 24 * epsilon * (2 * pow(sigma / r_abs, 12) / r_abs + pow(sigma / r_abs, 6) / r_abs);  //возвращает градиент потенциала от частицы 2 к частице 1



    }
    
    double periodic_conditions_on_coordinats(double x)
    {
        if(x > box_size / 2) 
        {
            return x -box_size;
        }
        else if (x < -box_size / 2)
        {
            return x + box_size;
        }
        else return x;
    }
    
    Vector3D periodic_conditions_on_vect(Vector3D r)
    {   
        return Vector3D(periodic_conditions_on_coordinats(r.x), periodic_conditions_on_coordinats(r.y), periodic_conditions_on_coordinats(r.z));
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
        {
            particles[i].r = particles[i].r + particles[i].v * dt + Gradients[i] * (-dt * dt / 2 * particles[i].m);
            particles[i].v = particles[i].v + Gradients[i] * (-dt);
            particles[i].r = periodic_conditions_on_vect(particles[i].r);
        }
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
        std::cout << i << std::endl;

    }
    std::ofstream outfile;
    
    // Открытие файла для записи
    outfile.open("output.txt");
    
    // Проверка, успешно ли открылся файл
    if (!outfile.is_open()) {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        return 1;
    }
    
    // Запись данных в файл
    for(unsigned int i = 0; i < 100; i++)
    {
        outfile << system.particles[i].v.x << " " << std::endl;
    }
    // Закрытие файла
    outfile.close();
    return 0;
}