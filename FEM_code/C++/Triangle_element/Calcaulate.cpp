#include <vector>
#include <cmath>
using namespace std;

//矩阵乘法
vector<vector<double>> matrix_multiply(vector<vector<double>> A, vector<vector<double>> B)
{
    int rows_A = A.size();
    int cols_A = A[0].size();
    int rows_B = B.size();
    int cols_B = B[0].size();

    vector<vector<double>> C(rows_A, vector<double>(cols_B, 0.0));

    for (int i = 0; i < rows_A; i++)
    {
        for (int j = 0; j < cols_B; j++)
        {
            for (int k = 0; k < cols_A; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

//矩阵转置
vector<vector<double>> matrix_transpose(vector<vector<double>> A)
{
    int rows = A.size();
    int cols = A[0].size();

    vector<vector<double>> AT(cols, vector<double>(rows, 0.0));

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            AT[j][i] = A[i][j];
        }
    }

    return AT;
}

//2×2矩阵行列式
double determinant_2x2(vector<vector<double>> A)
{
    double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

    return detA;
}

//2×2矩阵求逆
vector<vector<double>> inverse_2x2(vector<vector<double>> A)
{
    vector<vector<double>> invA(2, vector<double>(2, 0.0));

    double detA = determinant_2x2(A);

    invA[0][0] = A[1][1] / detA;
    invA[0][1] = -A[0][1] / detA;
    invA[1][0] = -A[1][0] / detA;
    invA[1][1] = A[0][0] / detA;

    return invA;
}

//向量相乘
vector<double> matrix_vector_multiply(
    vector<vector<double>> A,
    vector<double> x)
{
    int rows = A.size();
    int cols = A[0].size();

    vector<double> b(rows, 0.0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            b[i] += A[i][j] * x[j];
        }
    }

    return b;
}

//求解位移线性方程组
vector<double> solve_linear_system(
    vector<vector<double>> K,
    vector<double> F)
{
    int n = K.size();

    for (int k = 0; k < n - 1; k++)
    {
        int pivot = k;

        for (int i = k + 1; i < n; i++)
        {
            if (abs(K[i][k]) > abs(K[pivot][k]))
            {
                pivot = i;
            }
        }

        if (pivot != k)
        {
            swap(K[k], K[pivot]);
            swap(F[k], F[pivot]);
        }

        for (int i = k + 1; i < n; i++)
        {
            double factor = K[i][k] / K[k][k];

            for (int j = k; j < n; j++)
            {
                K[i][j] -= factor * K[k][j];
            }

            F[i] -= factor * F[k];
        }
    }

    vector<double> U(n, 0.0);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++)
        {
            sum += K[i][j] * U[j];
        }

        U[i] = (F[i] - sum) / K[i][i];
    }

    return U;
}


//求解U
vector<double> get_element_displacement(
    vector<double> U,
    vector<int> element_nodes)
{
    vector<double> Ue(6, 0.0);

    int n1 = element_nodes[0] - 1;
    int n2 = element_nodes[1] - 1;
    int n3 = element_nodes[2] - 1;

    Ue[0] = U[2 * n1];
    Ue[1] = U[2 * n1 + 1];

    Ue[2] = U[2 * n2];
    Ue[3] = U[2 * n2 + 1];

    Ue[4] = U[2 * n3];
    Ue[5] = U[2 * n3 + 1];

    return Ue;
}

//求解应变
vector<vector<double>> create_strain(
    vector<vector<vector<double>>> B,
    vector<vector<int>> T,
    vector<double> U)
{
    int number_elements = T.size();

    vector<vector<double>> strain(number_elements, vector<double>(3, 0.0));

    for (int e = 0; e < number_elements; e++)
    {
        vector<double> Ue = get_element_displacement(U, T[e]);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                strain[e][i] += B[i][j][e] * Ue[j];
            }
        }
    }

    return strain;
}

//求解sigma
vector<vector<double>> create_stress(
    vector<vector<double>> D,
    vector<vector<double>> strain)
{
    int number_elements = strain.size();

    vector<vector<double>> stress(number_elements, vector<double>(3, 0.0));

    for (int e = 0; e < number_elements; e++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                stress[e][i] += D[i][j] * strain[e][j];
            }
        }
    }

    return stress;
}

//求解米塞斯应力
vector<double> create_von_mises_stress(
    vector<vector<double>> stress)
{
    int number_elements = stress.size();

    vector<double> von_mises(number_elements, 0.0);

    for (int e = 0; e < number_elements; e++)
    {
        double sx = stress[e][0];
        double sy = stress[e][1];
        double txy = stress[e][2];

        von_mises[e] = sqrt(sx * sx - sx * sy + sy * sy + 3.0 * txy * txy);
    }

    return von_mises;
}

//求解RF
vector<double> create_reaction_force(
    vector<vector<double>> K_origin,
    vector<double> U,
    vector<double> F_origin)
{
    vector<double> KU = matrix_vector_multiply(K_origin, U);

    int number_dof = U.size();

    vector<double> R(number_dof, 0.0);

    for (int i = 0; i < number_dof; i++)
    {
        R[i] = KU[i] - F_origin[i];
    }

    return R;
}