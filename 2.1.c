#include <stdio.h>
#include <math.h>

#define M_PI 3.141592
#define N 4
#define A -2
#define B 4

double original_function(double x)
{
    double y;
    y = pow(x, 3) - exp(x) + 1;
    return y;
}

double divided_differences(double list_x[], double list_y[], double dd[][N + 1])
{
    for (int i = 0; i <= N; i++)
    {
        dd[i][0] = list_y[i];
    }

    for (int j = 1; j <= N; j++)
    {
        for (int i = 0; i <= N - j; i++)
        {
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (list_x[i + j] - list_x[i]);
        }
    }
    return dd[0][N];
}

void print_divided_differences(double dd[][N + 1])
{
    printf("\nDivided_differences:\n");

    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j <= N - i; j++)
        {
            printf(" %10.5lf |", dd[i][j]);
        }
        printf("\n");
    }
}

double newton_backward_polynomial(double x, double list_x[], double dd[N + 1][N + 1])
{
    double result = dd[N][0];
    double term;
    for (int j = 1; j <= N; j++)
    {
        term = dd[N - j][j];
        for (int i = 0; i < j; i++)
        {
            term *= (x - list_x[N - i]);
        }
        result += term;
    }
    return result;
}

void print_newton_backward_polynomial(double list_x[], double dd[][N + 1])
{
    printf("P(x) = %.5lf", dd[N][0]);
    for (int j = 1; j <= N; j++)
    {
        printf(" + (%.5lf", dd[N - j][j]);
        for (int i = 0; i < j; i++)
        {
            printf(" * (x - %.5lf)", list_x[N - i]);
        }
        printf(")");
    }
    printf("\n");
}

void verify_polynomial(double list_x[], double list_y[], double dd[][N + 1])
{
    printf("\nTest\n");
    for (int i = 0; i <= N; i++)
    {
        double p_val = newton_backward_polynomial(list_x[i], list_x, dd);
        printf("x[%d] = %.6f, P(x) = %.6f, f(x) = %.6f, Error = %.6e\n", i, list_x[i], p_val, list_y[i], fabs(p_val - list_y[i]));
    }
}

void write_newton_backward_polynomial_to_file(const char* filename, double list_x[], double dd[][N + 1]) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Ошибка открытия файла");
        return;
    }

    // Записываем многочлен в файл
    fprintf(file, "P(x) = %.5f", dd[N][0]);
    for (int j = 1; j <= N; j++) {
        fprintf(file, " + (%.5f", dd[N - j][j]);
        for (int i = 0; i < j; i++) {
            fprintf(file, " * (x - %.5f)", list_x[N - i]);
        }
        fprintf(file, ")");
    }
    fprintf(file, "\n");

    fclose(file);
}

int main(void)
{
    double x, t;
    double list_x[N + 1], list_y[N + 1];
    double dd[N + 1][N + 1];

    //узлы
    for (int k = 0; k <= N; k++)
    {
        t = cos(M_PI * (2 * k + 1) / (2 * (N + 1)));
        x = (A + B) / 2 + ((B - A) / 2) * t;
        list_x[k] = x;
    }

    //значение в узлах
    for (int k = 0; k <= N; k++)
    {
        x = list_x[k];
        list_y[k] = original_function(x);
    }


    printf("nodes_x:\n");
    for (int k = 0; k <= N; k++)
    {
        printf("x[%d] = %.6f, y[%d] = %.6f\n", k, list_x[k], k, list_y[k]);
    }

    divided_differences(list_x, list_y, dd);
    print_divided_differences(dd);

    print_newton_backward_polynomial(list_x, dd);

    verify_polynomial(list_x, list_y, dd);

    write_newton_backward_polynomial_to_file("polynomial.txt", list_x, dd);

    return 0;
}