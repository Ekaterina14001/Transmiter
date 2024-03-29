﻿// find_transmission.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cassert>     
#include <string>

struct Point
{
    double x;
    double y;
    //offsets
    double AB;
    double AC;
    double BC;
};

class TestPoints
{
private:
    Point D;
    Point E;
    Point F;
    Point A;
    Point B;
    Point C;

    
public:
    TestPoints()
    {
        int start = 10;
        int end = 30;

        D.x = rand() % (end - start + 1) + start;
        D.y = rand() % (end - start + 1) + start;
        E.x = rand() % (end - start + 1) + start;
        E.y = rand() % (end - start + 1) + start;
        F.x = rand() % (end - start + 1) + start;
        F.y = rand() % (end - start + 1) + start;

        A.x = rand() % (end - start + 1) + start;
        A.y = rand() % (end - start + 1) + start;
        B.x = rand() % (end - start + 1) + start;
        B.y = rand() % (end - start + 1) + start;
        C.x = rand() % (end - start + 1) + start;
        C.y = rand() % (end - start + 1) + start;

        D.AB = sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2)) + sqrt(pow(D.x - B.x, 2) + pow(D.y - B.y, 2));
        D.AC = sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2)) + sqrt(pow(D.x - C.x, 2) + pow(D.y - C.y, 2));
        D.BC = sqrt(pow(D.x - B.x, 2) + pow(D.y - B.y, 2)) + sqrt(pow(D.x - C.x, 2) + pow(D.y - C.y, 2));
        E.AB = sqrt(pow(E.x - A.x, 2) + pow(E.y - A.y, 2)) + sqrt(pow(E.x - B.x, 2) + pow(E.y - B.y, 2));
        E.AC = sqrt(pow(E.x - A.x, 2) + pow(E.y - A.y, 2)) + sqrt(pow(E.x - C.x, 2) + pow(E.y - C.y, 2));
        E.BC = sqrt(pow(E.x - B.x, 2) + pow(E.y - B.y, 2)) + sqrt(pow(E.x - C.x, 2) + pow(E.y - C.y, 2));
        F.AB = sqrt(pow(F.x - A.x, 2) + pow(F.y - A.y, 2)) + sqrt(pow(F.x - B.x, 2) + pow(F.y - B.y, 2));
        F.AC = sqrt(pow(F.x - A.x, 2) + pow(F.y - A.y, 2)) + sqrt(pow(F.x - C.x, 2) + pow(F.y - C.y, 2));
        F.BC = sqrt(pow(F.x - B.x, 2) + pow(F.y - B.y, 2)) + sqrt(pow(F.x - C.x, 2) + pow(F.y - C.y, 2));
    }
    Point get_A()
    {
        return A;
    }
    Point get_B()
    {
        return B;
    }
    Point get_C()
    {
        return C;
    }
    Point get_D()
    {
        return D;
    }
    Point get_E()
    {
        return E;
    }
    Point get_F()
    {
        return F;
    }
};

class Transmission: public TestPoints
{
private:
    Point D;
    Point E;
    Point F;
    Point A;
    Point B;
    Point C;
    Point standartA;
    Point standartB;
    Point standartC;

    

    void set_D(double x, double y, double AB, double AC, double BC)
    {
        D.x = x;
        D.y = y;
        D.AB = AB;//AD-BD
        D.AC = AC; //AD-CD
        D.BC = BC; //BD-CD
    }
    void set_E(double x, double y, double AB, double AC, double BC)
    {
        E.x = x;
        E.y = y;
        E.AB = AB;//AE-BE
        E.AC = AC; //AE-CE
        E.BC = BC; //BE-CE
    }
    void set_F(double x, double y, double AB, double AC, double BC)
    {
        F.x = x;
        F.y = y;
        F.AB = AB;//AF-BF
        F.AC = AC; //AF-CF
        F.BC = BC; //BF-CF
    }
  
    //general equation is sqrt(pow(x-x0,2)+pow(y-y0,2))-sqrt(pow(x-x1,2)+pow(y-y1,2))-v(t0-t1)=0
    // D, F and E - (x,y), A, B, C - (x0,y0) or (x1, y1), so vt0 is distance (x,y) to (x0,y0), vt1 is distance (x,y) to (x1,y1)
    //system of 9 equations
    double equation(unsigned int i, Point &A, Point &B, Point &C) {
        assert(i >0 && i<10);
        
        switch (i) {//i-num of equation
        case 1:
            return sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2)) + sqrt(pow(D.x -B.x, 2) + pow(D.y - B.y, 2)) - D.AB;
        case 2:
            return sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2)) + sqrt(pow(D.x - C.x, 2) + pow(D.y - C.y, 2)) - D.AC;
        case 3:
            return sqrt(pow(D.x - B.x, 2) + pow(D.y - B.y, 2)) + sqrt(pow(D.x - C.x, 2) + pow(D.y - C.y, 2)) - D.BC;
        case 4:
            return sqrt(pow(E.x - A.x, 2) + pow(E.y - A.y, 2)) + sqrt(pow(E.x - B.x, 2) + pow(E.y - B.y, 2)) - E.AB;
        case 5:
            return sqrt(pow(E.x - A.x, 2) + pow(E.y - A.y, 2)) + sqrt(pow(E.x - C.x, 2) + pow(E.y - C.y, 2)) - E.AC;
        case 6:
            return sqrt(pow(E.x - B.x, 2) + pow(E.y - B.y, 2)) + sqrt(pow(E.x - C.x, 2) + pow(E.y - C.y, 2)) - E.BC;
        case 7:
            return sqrt(pow(F.x - A.x, 2) + pow(F.y - A.y, 2)) + sqrt(pow(F.x - B.x, 2) + pow(F.y - B.y, 2)) - F.AB;
        case 8:
            return sqrt(pow(F.x - A.x, 2) + pow(F.y - A.y, 2)) + sqrt(pow(F.x - C.x, 2) + pow(F.y - C.y, 2)) - F.AC;
        case 9:
            return sqrt(pow(F.x - B.x, 2) + pow(F.y - B.y, 2)) + sqrt(pow(F.x - C.x, 2) + pow(F.y - C.y, 2)) - F.BC;
        default:
            return 0.0;
        }
    }
    double derivative_x(int i) {
    switch (i) {//AB AC BC
    case 1:
        return -(D.x-A.x) / sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2));
    case 2:
        return -(D.x-B.x) / sqrt(pow(D.x - B.x, 2) + pow(D.x - B.y, 2));
    case 3:
        return -(D.x - C.x) / sqrt(pow(D.x - C.x, 2) + pow(D.x - C.y, 2));
  
    default:
        return 0.0;
    }
}

double derivative_y(int i) {
    switch (i)
    {
   
    case 1:
        return -(D.y - A.y) / sqrt(pow(D.x - A.x, 2) + pow(D.y - A.y, 2));
    case 2:
        return -(D.y - B.y) / sqrt(pow(D.x - B.x, 2) + pow(D.y - B.y, 2));
    case 3:
        return -(D.y - C.y) / sqrt(pow(D.x - C.x, 2) + pow(D.x - C.y, 2));
    
    default:
        return 0.0;
    }

}
    double goldenSectionSearch(int i, double grad) {
        double a = 0.0;
        double b = 1.0;
        double tau = 0.381966011;  // Значение (1 - tau) используется для разделения интервала
        double epsilon = 0.00001;   // Точность метода
        double x1 = a + (1 - tau) * (b - a);
        double x2 = a + tau * (b - a);

        while (abs(b - a) > epsilon) {
            Point A1;
            Point A2;
            Point B1;
            Point B2;
            Point C1;
            Point C2;
            A1.x=A.x + x1 * grad;
            A1.y=A.y + x1 * grad;
            B1.x=B.x + x1 * grad;
            B1.y=B.y + x1 * grad;
            C1.x=C.x + x1 * grad;
            C1.y=C.y + x1 * grad;
            A2.x = A.x + x2 * grad;
            A2.y = A.y + x2 * grad;
            B2.x = B.x + x2 * grad;
            B2.y = B.y + x2 * grad;
            C2.x = C.x + x2 * grad;
            C2.y = C.y + x2 * grad;
            if (equation(i, A1, B1, C1) < equation(i, A2, B2, C2)) {
                b = x2;
                x2 = x1;
                x1 = a + (1 - tau) * (b - a);
            }
            else {
                a = x1;
                x1 = x2;
                x2 = a + tau * (b - a);
            }
        }

        return (a + b) / 2;  // Возвращаем значение в середине интервала
    }
    void steepestDescent(int mIterations) {
        //set_first();
        int iteration = 0;
        while (iteration < mIterations) {
            // Вычисляем градиенты функций для переменных A.x, A.y, B.x, B.y, C.x, C.y
            Point grad_A;
            grad_A.x = derivative_x(1);
            grad_A.y = derivative_y(1);
            Point grad_B;
            grad_B.x = derivative_x(2);
            grad_B.y = derivative_y(2);
            Point grad_C;
            grad_C.x = derivative_x(3);
            grad_C.y = derivative_y(3);

            // Находим оптимальный шаг для каждого уравнения с помощью метода золотого сечения
            Point step_A;
            step_A.x = goldenSectionSearch(1, grad_A.x);
            step_A.y = goldenSectionSearch(2, grad_A.y);
            Point step_B;
            step_B.x = goldenSectionSearch(3, grad_B.x);
            step_B.y = goldenSectionSearch(4, grad_B.y);
            Point step_C;
            step_C.x = goldenSectionSearch(5, grad_C.x);
            step_C.y = goldenSectionSearch(6, grad_C.y);

            // Обновляем значения переменных с найденными оптимальными шагами
            A.x -= step_A.x * grad_A.x;
            A.y -= step_A.y * grad_A.y;
            B.x -= step_B.x * grad_B.x;
            B.y -= step_B.y * grad_B.y;
            C.x -= step_C.x * grad_C.x;
            C.y -= step_C.y * grad_C.y;

            iteration++;
        }
    }
    void print()
    {
        std::cout << "A, B, C, D, E, F: (" << A.x << "," << A.y << "), (" << B.x << "," << B.y << "), (" <<
            C.x << "," << C.y << "), (" << D.x << "," << D.y << "), (" << E.x << "," << E.y << "), ("
            << F.x << "," << F.y << ")" << std::endl;

        std::cout << "standartA, standartB, standartC: (" << standartA.x << "," << standartA.y << "), (" << standartB.x <<
            "," << standartB.y << "), (" << standartC.x << "," << standartC.y << ")"<< std::endl;
    }

    void print_differences()
    {
        std::cout << "Diferences: ||A-standartA||=" << sqrt(pow(A.x - standartA.x, 2) + pow(A.y - standartA.y, 2)) << ", ||B-standartB||="
            << sqrt(pow(B.x - standartB.x, 2) + pow(B.y - standartB.y, 2)) << ", ||C-standartC||="
            << sqrt(pow(C.x - standartC.x, 2) + pow(C.y - standartC.y, 2)) << std::endl;
    }
   
    public:
        Transmission(TestPoints& t)
        {
            D = t.get_D();
            E = t.get_E();
            F = t.get_F();
            standartA = t.get_A();
            standartB = t.get_B();
            standartC = t.get_C();
            A.x = 0;
            A.y = 0;
            B.x = 0;
            B.y = 0;
            C.x = 0;
            C.y = 0;
        }

    void compute()
    {
        print();
        steepestDescent(10);
        std::cout << "After calculation" << std::endl;
        print();
        print_differences();
      
    }
    

};

int main()
{
    TestPoints t1;
    Transmission T(t1);
    T.compute();
}


