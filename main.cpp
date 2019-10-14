#include "stdio.h"
#include "math.h"
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream> 


const double G = 0.00000000006673;
// const double G = 1;

// const double timestep = 1.0;
const double timestep = 0.001;
const double kMASS = 1;
const double vi = 0.0;

// struct defining the mass, position, and velocity of each body
struct body {
    double mass;
    double pX;
    double pY;   
    double vX;
    double vY;
};

// places numBodies uniformly spaced along a line from 0 to 1
void generateBodies(std::vector<body> &bodies, int numBodies){
    double xcur = 0.0;
    double ycur = 0.0;
    double xStep = 1.0/numBodies;
    double yStep = 1.0/numBodies;
    for (auto i = bodies.begin(); i<bodies.end(); i++){
        i->mass = kMASS;
        i->vX = 0.0;
        i->vY = 0.0;
        i->pX = xcur;
        i->pY = ycur;
        xcur += xStep;
        ycur += yStep;
    }
}

// returns the rate (in millions of interactions / second)
double standardThreeBody(std::vector<body> &bodies, int numIterations){
    auto start_time = std::chrono::high_resolution_clock::now();

    for (auto i = bodies.begin(); i < bodies.end(); i++){
        for (auto j = bodies.begin(); j < bodies.end(); j++){
            if(i == j){
                continue;
            }
            double x_diff = (j->pX - i->pX);
            double y_diff = (j->pY - i->pY);
            double distance_squared = x_diff * x_diff + y_diff * y_diff;
            double inverse_distance = 1.0 / sqrt(distance_squared);

            // update by acceleration times time interval
            double product = G * j->mass * timestep / distance_squared;
            i->vX += x_diff * product * inverse_distance;
            i->vY += y_diff * product * inverse_distance;
        }
    }
    for (auto i = bodies.begin(); i < bodies.end(); i++){
        i->pX += i->vX * timestep;
        i->pY += i->vY * timestep;        
    }

    auto finish_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = (finish_time - start_time);
    double num_million_interactions = ((bodies.size() * (bodies.size() + 1)) / 2) / 1000000.0;
    return (num_million_interactions / time_span.count());
}

// similar to standardThreeBody, but uses newton's third to reduce the number 
// of calculations
double reducedThreeBody( std::vector<body> &bodies, int numIterations){
    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto i = bodies.begin(); i < bodies.end(); i++){
        for (auto j = i+1; j < bodies.end(); j++){
            double x_diff = (j->pX - i->pX);
            double y_diff = (j->pY - i->pY);
            double distance_squared = x_diff * x_diff + y_diff * y_diff;
            double inverse_distance = 1.0 / sqrt(distance_squared);
            // update by acceleration times time interval
            double product = G * timestep / distance_squared;
            double productj = -1.0 * i->mass * product;
            double producti = j->mass * product;
            i->vX += x_diff * producti * inverse_distance;
            i->vY += y_diff * producti * inverse_distance;
            j->vX += x_diff * productj * inverse_distance;
            j->vY += y_diff * productj * inverse_distance;
        }
    }
    for (auto i = bodies.begin(); i < bodies.end(); i++){
        i->pX += i->vX * timestep;
        i->pY += i->vY * timestep;
    }

    auto finish_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = (finish_time - start_time);
    double num_million_interactions = ((bodies.size() * (bodies.size() + 1)) / 2) / 1000000.0;
    return (num_million_interactions / time_span.count());
}

int main(){

    // std::vector<int> test_numbers = {15000};
    std::vector<int> test_numbers = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000};
    std::vector<body> bodies;
    std::vector<double> standard_performances(test_numbers.size());
    std::vector<double> reduced_performances(test_numbers.size());
    int ind = 0;
    std::vector<body> standard_bodies;
    std::vector<body> reduced_bodies;
    for (auto iter = test_numbers.begin(); iter < test_numbers.end(); iter++){
        bodies.clear();
        bodies.resize(*iter);
        generateBodies(bodies, *iter);
        standard_bodies = bodies;
        reduced_bodies = bodies;
        standard_performances[ind] = standardThreeBody(standard_bodies, test_numbers[ind]);
        reduced_performances[ind] = reducedThreeBody(reduced_bodies, test_numbers[ind]);
        ind += 1;
    }


    // write the performances to body_benchmarks.dat
    ////////////////////////////////////////////////////////////////////////////////////////
    // first column will be the number of bodies
    // second column will be standard performances
    // third column will be reduced performances
    std::ofstream out_benches;
    out_benches.open("body_benchmarks.dat"); // opens the file
    if( !out_benches ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    for (int i = 0; i < test_numbers.size(); i++){
        out_benches << test_numbers[i] << "\t\t" << standard_performances[i] << "\t\t"
                << reduced_performances[i] << std::endl;
    }
    out_benches.close();
    
    return 0;
}    