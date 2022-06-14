#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double dwalltime();

int main(int argc, char *argv[])
{

    /* Check command line parameters */
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((N % bs) != 0))
    {
        printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
        exit(1);
    }

    int *F, *fibonacci, check = 1, maxF, minF, countF, count;
    double *A, *B, *resAB, *C, *D, *resDF, *R, *fibF, val, timetick;

    /* Aloca memoria par las matrices */
    A = (double *)malloc(sizeof(double) * N * N);
    B = (double *)malloc(sizeof(double) * N * N);
    C = (double *)malloc(sizeof(double) * N * N);
    D = (double *)malloc(sizeof(double) * N * N);
    F = (int *)malloc(sizeof(int) * N * N);
    R = (double *)malloc(sizeof(double) * N * N);
    resAB = (double *)malloc(sizeof(double) * N * N);
    resDF = (double *)malloc(sizeof(double) * N * N);
    fibF = (double *)malloc(sizeof(double) * N * N);

    /* Aloca memoria para el arreglo*/
    fibonacci = (int *)malloc(sizeof(int) * 40);

    /* Inicializa las matrices*/
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i * N + j] = rand() % 41;
            B[j * N + i] = rand() % 41;
            C[j * N + i] = rand() % 41;
            D[i * N + j] = rand() % 41;
            F[j * N + i] = rand() % 41;
            resAB[i * N + j] = 0;
            resDF[i * N + j] = 0;
            fibF[j * N + i] = 0;
            R[i * N + j] = 0;
        }
    }

    /* Start timing */
    timetick = dwalltime();

    /* Calculo de los valores de fibonacci hasta 40, para recalcularlos cada vez */
    fibonacci[0] = 0;
    fibonacci[1] = 1;

    for (int i = 2; i < 40; i++)
    {
        fibonacci[i] = fibonacci[i - 1] + fibonacci[i - 2];
    }

    /* busca el max, min, y fib */
    minF = 41;
    maxF = -1;
    count = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (minF > F[j * N + i])
                minF = F[j * N + i];

            if (maxF < F[j * N + i])
                maxF = F[j * N + i];

            count += F[j * N + i];

            fibF[j * N + i] = (double)fibonacci[F[j * N + i]];
        }
    }

    /* calculo val */
    val = (maxF * minF) / (count / (N * N));

    /* Multiply Matrices */
    int i, j, k, p, q, r;
    int in, jn, ink, jnk, inkpn, jnkqn, inkpnq;

    for (i = 0; i < N; i += bs)
    {
        in = i * N;
        for (j = 0; j < N; j += bs)
        {
            jn = j * N;
            for (k = 0; k < N; k += bs)
            {
                ink = in + k;
                jnk = jn + k;
                for (p = 0; p < bs; p++)
                {
                    inkpn = ink + p * N;
                    for (q = 0; q < bs; q++)
                    {
                        inkpnq = inkpn + q;
                        jnkqn = jnk + q * N;
                        for (r = 0; r < bs; r++)
                        {
                            resAB[inkpnq] += A[inkpn + r] * B[jnkqn + r];
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < N; i += bs)
    {
        in = i * N;
        for (j = 0; j < N; j += bs)
        {
            jn = j * N;
            for (k = 0; k < N; k += bs)
            {
                ink = in + k;
                jnk = jn + k;
                for (p = 0; p < bs; p++)
                {
                    inkpn = ink + p * N;
                    for (q = 0; q < bs; q++)
                    {
                        inkpnq = inkpn + q;
                        jnkqn = jnk + q * N;
                        for (r = 0; r < bs; r++)
                        {
                            R[inkpnq] += resAB[inkpn + r] * C[jnkqn + r];
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < N; i += bs)
    {
        in = i * N;
        for (j = 0; j < N; j += bs)
        {
            jn = j * N;
            for (k = 0; k < N; k += bs)
            {
                ink = in + k;
                jnk = jn + k;
                for (p = 0; p < bs; p++)
                {
                    inkpn = ink + p * N;
                    for (q = 0; q < bs; q++)
                    {
                        inkpnq = inkpn + q;
                        jnkqn = jnk + q * N;
                        for (r = 0; r < bs; r++)
                        {
                            resDF[inkpnq] += D[inkpn + r] * fibF[jnkqn + r];
                        }
                    }
                }
            }
        }
    }

    /* Suma de matrices de resABC y resDF y multiplicacion con val */
    for (int i = 0; i < N; i++)
    {
        in = i * N;
        for (int j = 0; j < N; j++)
        {
            R[in + j] = (R[in + j] + resDF[in + j]) * val;
        }
    }

    /* Stop timeming */
    workTimeSequential = dwalltime() - timetick;
    printf("Tiempo en segundos secuencial: %f\n", workTimeSequential);

    /* Checking results */
    for (int i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            check = check && (R[i * N + j] == (N * N + N));
        }
    }

    if (check)
    {
        printf("\nResultado correcto de la ecuacion\n");
    }
    else
    {
        printf("\nResultado incorrecto de la ecuacion\n");
    }

    free(A);
    free(B);
    free(resAB);
    free(C);
    free(D);
    free(resDF);
    free(R);
    free(F);
    free(fibF);
    free(fibonacci);
}

double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}
