#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


const double Pi = acos(-1.);

class Fx : public Function
{
public:
    double operator()(double x) const
    {
        return sqrt(3) * sin(2 * Pi * x);
    }
};


class Fy : public Function
{
public:
    double operator()(double x) const
    {
        return 2.0 * sqrt(3) / 3.0 * cos(2 * Pi * x) + 2.0 / 3.0 * sqrt(sqrt(3) * fabs(sin(2 * Pi * x)));
    }
};




int main()
{
    std::vector<int> m={10,40,160};
    std::vector<Point> points((size_t)(m[0]+2), Point(0,0)); //  add m+1 consecutive marker points. Form the closed curve. The last one is the same as the first one
    Fx fx;
    Fy fy;

    std::cout << "Problem F" << std::endl;
    // m=10
    std::ofstream fout("./data/data_10.txt"); 
    if (!fout)
    { 
        std::cerr << "File open failed!" << std::endl;
        return 1;
    }
    
    // choose consecutive point
    for(int i=0; i<=m[0]; i++){
        double x = static_cast<double> (i)/ (m[0]+1);
        double y = static_cast<double> (i)/ (m[0]+1);
        
        points[i] = Point(x,y);
    }
    if(m[0]%2==0){
        points[m[0]/2]=Point(0.5,0.5); // change the point to the lowest consecutive point of the heart
    }


    BezierInterpolation Bezier_10(points,Fx(),Fy());
    std::vector<Point> Bezier_p_10=Bezier_10.curve(0.01);

    for(auto p:Bezier_p_10){
        fout << p.x << ", " << p.y << std::endl;
    }
    fout.close();

    // m=40
    points.resize((size_t)(m[1]+2),Point(0,0));
    std::ofstream fout2("./data/data_40.txt"); 
    if (!fout2)
    { 
        std::cerr << "File open failed!" << std::endl;
        return 1;
    }
    
    for(int i=0; i<=m[1]; i++){
        double x = static_cast<double> (i)/ (m[1]+1);
        double y = static_cast<double> (i)/ (m[1]+1);
        points[i] = Point(x,y);
    }
    if(m[1]%2==0){
        points[m[1]/2]=Point(0.5,0.5);
    }

    BezierInterpolation Bezier_40(points,Fx(),Fy());
    std::vector<Point> Bezier_p_40=Bezier_40.curve(0.01);
    

    for(auto p:Bezier_p_40){
        fout2 << p.x << ", " << p.y << std::endl;
    }
    fout2.close();

    // m=160
    points.resize((size_t)(m[2]+2),Point(0,0));
    std::ofstream fout3("./data/data_160.txt"); 
    if (!fout3)
    { 
        std::cerr << "File open failed!" << std::endl;
        return 1;
    }
    
    for(int i=0; i<=m[2]; i++){
        double x = static_cast<double> (i)/ (m[2]+1);
        double y = static_cast<double> (i)/ (m[2]+1);
        points[i] = Point(x,y);
    }
    if(m[2]%2==0){
        points[m[2]/2]=Point(0.5,0.5);
    }

    BezierInterpolation Bezier_160(points,Fx(),Fy());
    std::vector<Point> Bezier_p_160=Bezier_160.curve(0.01);


    for(auto p:Bezier_p_160){
        fout3 << p.x << ", " << p.y << std::endl;
    }
    fout3.close();

    return 0;
}