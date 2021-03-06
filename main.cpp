/*
 * main.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: d-w-h
 *
 *      This code solves the advection diffusion equation:
 *      Da/U*d2Ca/dz2 - dCa/dz + ra/U = 0
 *      Using the Gauss-Seidel iteration method.
 */

#include <math.h>
#include <stdio.h>
#include "solver.hpp"
#include "user_types.hpp"

double ra(double Ca) {
    /* Reaction rate law */
    double k = 1.0;
    return -k * Ca;
}

int main(int argc, char* argv[]) {
    p_params physical_parameters;
    g_params grid_parameters;
    s_data solver_data;

    /* Parameters */
    grid_parameters.num_nodes = 40;       //Number of nodes
    grid_parameters.L = 1.6;              //Length of domain
    physical_parameters.U = 1.0;          //Fluid velocity
    physical_parameters.Cao = 1.0;        //Inlet concentration
    physical_parameters.Da = 1.0;         //Diffusion coefficient

    /* Allocate data for solver results */
    solver_data.Ca = new double[grid_parameters.num_nodes];
    solver_data.z_c = new double[grid_parameters.num_nodes];

    /* Execute solver */
    solver(physical_parameters, grid_parameters, &solver_data);

    /* Print data */
    for(int i = 0; i < grid_parameters.num_nodes; ++i) {
        printf("i: %i, z: %f, Ca: %f\n", i, solver_data.z_c[i], solver_data.Ca[i]);
    }

    /* Deallocate data */
    delete [] solver_data.Ca;
    delete [] solver_data.z_c;

    return 0;
}

