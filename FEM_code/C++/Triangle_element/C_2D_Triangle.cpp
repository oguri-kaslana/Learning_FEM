#include "Triangle.h"
#include "Boundary.h"
#include "Calcaulate.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

int main()
{
    cout << scientific << setprecision(6);

    string element_type = "T6";

    std::filesystem::path output_root =
        LR"(D:\A_homework\有限元程序设计\Mid_homework)";

    std::filesystem::create_directories(output_root);

    auto double_to_name = [](double value) -> string
        {
            ostringstream oss;
            oss << fixed << setprecision(2) << value;

            string s = oss.str();

            for (int i = 0; i < s.size(); i++)
            {
                if (s[i] == '.')
                {
                    s[i] = 'p';
                }
            }

            return s;
        };

    auto EB_deflection = [](
        double x,
        double L,
        double H,
        double thickness,
        double E,
        double q) -> double
        {
            double P_total = q * H * thickness;
            double I = thickness * H * H * H / 12.0;

            double v = -P_total * x * x * (3.0 * L - x)
                / (6.0 * E * I);

            return v;
        };

    double H = 4.0;
    double thickness = 1.0;

    double E = 110000.0;
    double nu = 0.3;

    double q = 5.0;
    double qx = 0.0;
    double qy = -q;

    vector<double> L_list = {
        4.0,
        24.0,
        44.0,
        64.0,
        84.0
    };

    vector<double> seeds_list;

    if (element_type == "T3")
    {
        seeds_list = {
            2.0,
            1.0,
            0.5
        };
    }
    else if (element_type == "T6")
    {
        seeds_list = {
            2.0,
            1.0
        };
    }
    string result_prefix = element_type + "_";

    ofstream convergence_file(output_root / (result_prefix + "convergence.csv"));
    ofstream centerline_all_file(output_root / (result_prefix + "centerline_all.csv"));
    ofstream summary_file(output_root / (result_prefix + "summary.csv"));

    convergence_file
        << "element_type,L,seeds,number_nodes,number_elements,number_dof,"
        << "max_disp,max_node,max_ux,max_uy,relative_change\n";

    centerline_all_file
        << "element_type,node_id,L,seeds,x,y,ux_FEM,uy_FEM,uy_EB,abs_error,relative_error\n";

    summary_file
        << "element_type,L,seeds,number_nodes,number_elements,number_dof,"
        << "right_mid_node,right_mid_ux,right_mid_uy,right_mid_uy_EB,"
        << "right_mid_abs_error,right_mid_relative_error,"
        << "max_disp,max_node,max_ux,max_uy\n";

    for (int iL = 0; iL < L_list.size(); iL++)
    {
        double L = L_list[iL];

        double previous_max_disp = 0.0;

        double final_seeds = 0.0;
        int final_number_nodes = 0;
        int final_number_elements = 0;
        int final_number_dof = 0;

        int final_right_mid_node = -1;
        double final_right_mid_ux = 0.0;
        double final_right_mid_uy = 0.0;
        double final_right_mid_uy_EB = 0.0;
        double final_right_mid_abs_error = 0.0;
        double final_right_mid_relative_error = 0.0;

        double final_max_disp = 0.0;
        int final_max_node = 0;
        double final_max_ux = 0.0;
        double final_max_uy = 0.0;

        for (int iSeed = 0; iSeed < seeds_list.size(); iSeed++)
        {
            double seeds = seeds_list[iSeed];

            cout << "========================================" << endl;
            cout << "element_type = " << element_type << endl;
            cout << "L = " << L << endl;
            cout << "seeds = " << seeds << endl;
            cout << "========================================" << endl;

            string case_name =
                element_type + "_L_" + double_to_name(L)
                + "_seed_" + double_to_name(seeds);

            std::filesystem::path case_folder = output_root / case_name;

            std::filesystem::create_directories(case_folder);

            vector<vector<int>> nodes = create_nodes(L, H, seeds);
            vector<vector<double>> P3 = create_Points(L, H, seeds);
            vector<vector<int>> T3 = create_Topolocy(nodes);

            vector<vector<double>> P;
            vector<vector<int>> T;
            vector<vector<double>> D = create_D_plane_stress(E, nu);
            vector<vector<double>> K_global;
            vector<vector<double>> strain;

            if (element_type == "T3")
            {
                P = P3;
                T = T3;

                vector<vector<vector<double>>> J = create_Jacobi(P, T);
                vector<vector<double>> dN_local = create_dN_local();
                vector<vector<vector<double>>> B = create_B_triangle(dN_local, J);
                vector<vector<vector<double>>> K_local =
                    create_K_local(B, J, D, thickness);

                K_global = create_K_global(K_local, T, P);
            }
            else if (element_type == "T6")
            {
                pair<vector<vector<double>>, vector<vector<int>>> mesh6 =
                    create_T6_mesh(P3, T3);

                P = mesh6.first;
                T = mesh6.second;

                vector<vector<vector<double>>> K_local =
                    create_K_local_T6(P, T, D, thickness);

                K_global = create_K_global_general(K_local, T, P);
            }
            else
            {
                cout << "unknown element_type" << endl;
                return 0;
            }

            int number_nodes = P.size();
            int number_elements = T.size();
            int number_dof = number_nodes * 2;

            vector<double> F_global(number_dof, 0.0);

            if (element_type == "T3")
            {
                vector<vector<int>> node_q = find_right_boundary_edges(P, L);
                vector<vector<double>> values_q =
                    create_uniform_q_values(node_q, qx, qy);

                apply_q_boundary(F_global, P, node_q, values_q, thickness);
            }
            else if (element_type == "T6")
            {
                vector<vector<int>> node_q = find_right_boundary_edges_T6(P, T, L);
                vector<vector<double>> values_q =
                    create_uniform_q_values(node_q, qx, qy);

                apply_q_boundary_T6(F_global, P, node_q, values_q, thickness);
            }

            vector<vector<double>> K_origin = K_global;
            vector<double> F_origin = F_global;

            vector<vector<int>> node_u = find_left_fixed_nodes(P);
            vector<double> values_u = create_zero_u_values(node_u);

            apply_u_boundary(K_global, F_global, node_u, values_u);

            vector<double> U = solve_linear_system(K_global, F_global);

            if (element_type == "T3")
            {
                vector<vector<vector<double>>> J = create_Jacobi(P, T);
                vector<vector<double>> dN_local = create_dN_local();
                vector<vector<vector<double>>> B = create_B_triangle(dN_local, J);

                strain = create_strain(B, T, U);
            }
            else if (element_type == "T6")
            {
                strain = create_strain_T6(P, T, U);
            }

            vector<vector<double>> stress = create_stress(D, strain);
            vector<double> von_mises = create_von_mises_stress(stress);
            vector<double> R = create_reaction_force(K_origin, U, F_origin);

            double max_disp = 0.0;
            int max_node = 1;
            double max_ux = 0.0;
            double max_uy = 0.0;

            for (int i = 0; i < number_nodes; i++)
            {
                double ux = U[2 * i];
                double uy = U[2 * i + 1];

                double u_abs = sqrt(ux * ux + uy * uy);

                if (u_abs > max_disp)
                {
                    max_disp = u_abs;
                    max_node = i + 1;
                    max_ux = ux;
                    max_uy = uy;
                }
            }

            double relative_change = 0.0;

            if (iSeed > 0 && abs(max_disp) > 1.0e-12)
            {
                relative_change =
                    abs(max_disp - previous_max_disp) / abs(max_disp);
            }

            previous_max_disp = max_disp;

            double eps = 1.0e-8;

            int right_mid_node = -1;
            double right_mid_ux = 0.0;
            double right_mid_uy = 0.0;

            for (int i = 0; i < P.size(); i++)
            {
                double x = P[i][0];
                double y = P[i][1];

                if (abs(x - L) < eps && abs(y - H / 2.0) < eps)
                {
                    right_mid_node = i + 1;
                    right_mid_ux = U[2 * i];
                    right_mid_uy = U[2 * i + 1];
                    break;
                }
            }

            double right_mid_uy_EB =
                EB_deflection(L, L, H, thickness, E, q);

            double right_mid_abs_error =
                abs(right_mid_uy - right_mid_uy_EB);

            double right_mid_relative_error = 0.0;

            if (abs(right_mid_uy_EB) > 1.0e-12)
            {
                right_mid_relative_error =
                    right_mid_abs_error / abs(right_mid_uy_EB);
            }

            ofstream nodes_file(case_folder / "nodes.csv");

            nodes_file << "node_id,x,y\n";

            for (int i = 0; i < P.size(); i++)
            {
                nodes_file << i + 1 << ","
                    << P[i][0] << ","
                    << P[i][1] << "\n";
            }

            nodes_file.close();

            ofstream elements_file(case_folder / "elements.csv");

            if (element_type == "T3")
            {
                elements_file << "element_id,node1,node2,node3\n";

                for (int e = 0; e < T.size(); e++)
                {
                    elements_file << e + 1 << ","
                        << T[e][0] << ","
                        << T[e][1] << ","
                        << T[e][2] << "\n";
                }
            }
            else if (element_type == "T6")
            {
                elements_file << "element_id,node1,node2,node3,node4,node5,node6\n";

                for (int e = 0; e < T.size(); e++)
                {
                    elements_file << e + 1 << ","
                        << T[e][0] << ","
                        << T[e][1] << ","
                        << T[e][2] << ","
                        << T[e][3] << ","
                        << T[e][4] << ","
                        << T[e][5] << "\n";
                }
            }

            elements_file.close();

            ofstream displacement_file(case_folder / "displacement.csv");

            displacement_file << "node_id,x,y,ux,uy,u_abs\n";

            for (int i = 0; i < P.size(); i++)
            {
                double ux = U[2 * i];
                double uy = U[2 * i + 1];
                double u_abs = sqrt(ux * ux + uy * uy);

                displacement_file << i + 1 << ","
                    << P[i][0] << ","
                    << P[i][1] << ","
                    << ux << ","
                    << uy << ","
                    << u_abs << "\n";
            }

            displacement_file.close();

            ofstream strain_file(case_folder / "strain.csv");

            strain_file << "element_id,xc,yc,epsilon_x,epsilon_y,gamma_xy\n";

            for (int e = 0; e < T.size(); e++)
            {
                int n1 = T[e][0] - 1;
                int n2 = T[e][1] - 1;
                int n3 = T[e][2] - 1;

                double xc = (P[n1][0] + P[n2][0] + P[n3][0]) / 3.0;
                double yc = (P[n1][1] + P[n2][1] + P[n3][1]) / 3.0;

                strain_file << e + 1 << ","
                    << xc << ","
                    << yc << ","
                    << strain[e][0] << ","
                    << strain[e][1] << ","
                    << strain[e][2] << "\n";
            }

            strain_file.close();

            ofstream stress_file(case_folder / "stress.csv");

            stress_file << "element_id,xc,yc,sigma_x,sigma_y,tau_xy,von_mises\n";

            for (int e = 0; e < T.size(); e++)
            {
                int n1 = T[e][0] - 1;
                int n2 = T[e][1] - 1;
                int n3 = T[e][2] - 1;

                double xc = (P[n1][0] + P[n2][0] + P[n3][0]) / 3.0;
                double yc = (P[n1][1] + P[n2][1] + P[n3][1]) / 3.0;

                stress_file << e + 1 << ","
                    << xc << ","
                    << yc << ","
                    << stress[e][0] << ","
                    << stress[e][1] << ","
                    << stress[e][2] << ","
                    << von_mises[e] << "\n";
            }

            stress_file.close();

            ofstream reaction_file(case_folder / "reaction.csv");

            reaction_file << "node_id,direction,reaction\n";

            for (int k = 0; k < node_u.size(); k++)
            {
                int node_id = node_u[k][0];
                int direction = node_u[k][1];

                int dof = 2 * (node_id - 1) + (direction - 1);

                reaction_file << node_id << ","
                    << direction << ","
                    << R[dof] << "\n";
            }

            reaction_file.close();

            vector<int> centerline_nodes;

            for (int i = 0; i < P.size(); i++)
            {
                if (abs(P[i][1] - H / 2.0) < eps)
                {
                    centerline_nodes.push_back(i);
                }
            }

            sort(centerline_nodes.begin(), centerline_nodes.end(),
                [&](int a, int b)
                {
                    return P[a][0] < P[b][0];
                });

            ofstream centerline_file(case_folder / "centerline_EB_compare.csv");

            centerline_file
                << "node_id,L,seeds,x,y,ux_FEM,uy_FEM,uy_EB,abs_error,relative_error\n";

            for (int k = 0; k < centerline_nodes.size(); k++)
            {
                int id = centerline_nodes[k];

                double x = P[id][0];
                double y = P[id][1];

                double ux = U[2 * id];
                double uy = U[2 * id + 1];

                double uy_EB =
                    EB_deflection(x, L, H, thickness, E, q);

                double abs_error = abs(uy - uy_EB);

                double relative_error = 0.0;

                if (abs(uy_EB) > 1.0e-12)
                {
                    relative_error = abs_error / abs(uy_EB);
                }

                centerline_file
                    << id + 1 << ","
                    << L << ","
                    << seeds << ","
                    << x << ","
                    << y << ","
                    << ux << ","
                    << uy << ","
                    << uy_EB << ","
                    << abs_error << ","
                    << relative_error << "\n";

                centerline_all_file
                    << element_type << ","
                    << id + 1 << ","
                    << L << ","
                    << seeds << ","
                    << x << ","
                    << y << ","
                    << ux << ","
                    << uy << ","
                    << uy_EB << ","
                    << abs_error << ","
                    << relative_error << "\n";
            }

            centerline_file.close();

            convergence_file
                << element_type << ","
                << L << ","
                << seeds << ","
                << number_nodes << ","
                << number_elements << ","
                << number_dof << ","
                << max_disp << ","
                << max_node << ","
                << max_ux << ","
                << max_uy << ","
                << relative_change << "\n";

            final_seeds = seeds;
            final_number_nodes = number_nodes;
            final_number_elements = number_elements;
            final_number_dof = number_dof;

            final_right_mid_node = right_mid_node;
            final_right_mid_ux = right_mid_ux;
            final_right_mid_uy = right_mid_uy;
            final_right_mid_uy_EB = right_mid_uy_EB;
            final_right_mid_abs_error = right_mid_abs_error;
            final_right_mid_relative_error = right_mid_relative_error;

            final_max_disp = max_disp;
            final_max_node = max_node;
            final_max_ux = max_ux;
            final_max_uy = max_uy;

            cout << "number_nodes = " << number_nodes << endl;
            cout << "number_elements = " << number_elements << endl;
            cout << "number_dof = " << number_dof << endl;
            cout << "max_disp = " << max_disp << endl;
            cout << "relative_change = " << relative_change << endl;
            cout << "right_mid_uy = " << right_mid_uy << endl;
            cout << "right_mid_uy_EB = " << right_mid_uy_EB << endl;
            cout << "right_mid_relative_error = " << right_mid_relative_error << endl;
            cout << endl;
        }

        summary_file
            << element_type << ","
            << L << ","
            << final_seeds << ","
            << final_number_nodes << ","
            << final_number_elements << ","
            << final_number_dof << ","
            << final_right_mid_node << ","
            << final_right_mid_ux << ","
            << final_right_mid_uy << ","
            << final_right_mid_uy_EB << ","
            << final_right_mid_abs_error << ","
            << final_right_mid_relative_error << ","
            << final_max_disp << ","
            << final_max_node << ","
            << final_max_ux << ","
            << final_max_uy << "\n";
    }

    convergence_file.close();
    centerline_all_file.close();
    summary_file.close();

    cout << "All results have been written to:" << endl;
    cout << "D:\\A_homework\\有限元程序设计\\Mid_homework" << endl;

    return 0;
}