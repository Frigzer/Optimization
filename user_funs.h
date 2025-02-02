#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// Lab_01

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
double GetFib(int n);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);

// Lab_02

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);

// Lab_03

matrix ff3T(matrix, matrix = NAN, matrix = NAN);

// Lab_04

matrix ff4T(matrix, matrix = NAN, matrix = NAN);

matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix hf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
//class Fun {
//public:
//    double x1;
//    double x12;
//    double x2;
//    double x22;
//    double c;
//    double x1x2;
//
//    // Konstruktor z parametrami
//    Fun(double a, double b, double f, double d, double e, double g) {
//        x1 = a;
//        x12 = b;
//        x2 = f;
//        x22 = d;
//        x1x2 = e;
//        c = g;
//    }
//
//    // Konstruktor domyœlny
//    Fun() {
//        x1 = 0;
//        x12 = 0;
//        x2 = 0;
//        x22 = 0;
//        c = 0;
//        x1x2 = 0; // Domyœlna wartoœæ dla nowego pola
//    }
//    Fun(const Fun& other) {
//        x1 = other.x1;
//        x12 = other.x12;
//        x2 = other.x2;
//        x22 = other.x22;
//        c = other.c;
//        x1x2 = other.x1x2;
//    }
//
//    // Prze³adowanie operatora +
//    Fun operator+(const Fun& other) const {
//        return Fun(
//            x1 + other.x1,
//            x12 + other.x12,
//            x2 + other.x2,
//            x22 + other.x22,
//            c + other.c,
//            x1x2 + other.x1x2 
//        );
//    }
//    void override(double a, double b, double f, double d, double e, double g) const {
//        x1 = a;
//        x12 = b;
//        x2 = f;
//        x22 = d;
//        x1x2 = e;
//        c = g;
//    }
//    
//    void divergence() const {
//        this->override(2 * x12 + x1x2, 0, x1x2 + 2 * x21, 0, 0, x1 + x2);
//    }
//
//};

// Lab_06

matrix ff6T(matrix x, matrix = NAN, matrix = NAN);
matrix df6(double t, matrix Y, matrix ud1, matrix = NAN);
matrix ff6R(matrix, matrix = NAN, matrix = NAN);