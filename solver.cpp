/*
 * solver.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include "main.hpp"
#include "user_types.hpp"

double dra_dCa(double Ca) {
    double dCa = 0.0001;
    return (ra(Ca + dCa) - ra(Ca)) / dCa;
}

void solver(p_params physical_parameters, g_params grid_parameters, s_data* solver_data) {
    double L, U, Cin, Cao, Da, del_z;
    int num_nodes, i, outer_it, inner_it, max_outer_it, max_inner_it;

    /* Parameters */
    num_nodes = grid_parameters.num_nodes;
    L = grid_parameters.L;
    U = physical_parameters.U;
    Cao = physical_parameters.Cao;
    Da = physical_parameters.Da;

    max_outer_it = 300;
    max_inner_it = 300;

    /* Start simulation */
    del_z = L / num_nodes;
    double* Ca_prev_it = new double[num_nodes];

    /* Initialize Ca and z_c*/
    for(i = 0; i < num_nodes; ++i) {
        solver_data->Ca[i] = 0.0;
        solver_data->z_c[i] = i*del_z + 0.5*del_z;
    }

    /* Gauss Seidel iterations */
    outer_it = 0;
    while(outer_it < max_outer_it) {
        for(i = 0; i < num_nodes; ++i) {
            Ca_prev_it[i] = solver_data->Ca[i];
        }
        inner_it = 0;
        while(inner_it < max_inner_it) {
            //Inlet node
            Cin = (Cao + Da*solver_data->Ca[0]/(U*0.5*del_z)) / (1+Da/(U*0.5*del_z));
            // Left most node
            solver_data->Ca[0] = (Da*Cin/(0.5*del_z) + Da*solver_data->Ca[1]/del_z + U*Cin + ra(Ca_prev_it[0])*del_z - dra_dCa(Ca_prev_it[0])*Ca_prev_it[0]*del_z) / (Da/(0.5*del_z) + Da/del_z + U - dra_dCa(Ca_prev_it[0])*del_z);
            // Central nodes
            for(i = 1; i < num_nodes - 1; ++i) {
                solver_data->Ca[i] = (Da*solver_data->Ca[i-1]/del_z + Da*solver_data->Ca[i+1]/del_z + U*solver_data->Ca[i-1] + ra(Ca_prev_it[i])*del_z - dra_dCa(Ca_prev_it[i])*Ca_prev_it[i]*del_z) / (Da/del_z + Da/del_z + U - dra_dCa(Ca_prev_it[i])*del_z);
            }
            // Right most node
            solver_data->Ca[num_nodes-1] = (Da*solver_data->Ca[num_nodes-2]/del_z + U*solver_data->Ca[num_nodes-2] + ra(Ca_prev_it[num_nodes-1])*del_z - dra_dCa(Ca_prev_it[num_nodes-1])*Ca_prev_it[num_nodes-1]*del_z) / (Da/del_z + U - dra_dCa(Ca_prev_it[num_nodes-1])*del_z);

            ++inner_it;
        }

        ++outer_it;
    }
}
