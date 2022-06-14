#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

// Variables compartidas
int thread_count, local_size, local_size_mul;
int N, bs;
int max_num = -1, min_num = 41, prom_num = 0;
double val, timetick, workTimeSequential;

// Arrays
double *A, *B, *resAB, *C, *D, *resDF, *R, *fibF;
int *F, *fibonacci, *ids;
int *prom_nums;

// Locks
pthread_mutex_t mutex;

// Barrier
pthread_barrier_t barrier;

double dwalltime();
void *equation(void *rank);
void sequential();

int main(int argc, char *argv[])
{
	int thread; // for

	/* Check command line parameters */
	if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((N % bs) != 0))
	{
		printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
		exit(1);
	}

	// Mutex
	pthread_mutex_init(&mutex, NULL);

	pthread_t *thread_handles;
	thread_count = atoi(argv[3]);
	local_size = N / thread_count;

	thread_handles = malloc(thread_count * sizeof(pthread_t));

	// Initilaize barrier
	pthread_barrier_init(&barrier, NULL, thread_count);

	/* Aloca memoria para los id de los thread*/
	ids = (int *)malloc(sizeof(int) * thread_count);

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

	/* Aloca memoria para el promedio*/
	prom_nums = (int *)malloc(sizeof(int) * thread_count);

	/* Aloca memoria para el arreglo de fibonacci*/
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

	/* Start timing paralels*/
	timetick = dwalltime();

	/* Calculo de los valores de fibonacci */
	fibonacci[0] = 0;
	fibonacci[1] = 1;

	for (int i = 2; i < 40; i++)
	{
		fibonacci[i] = fibonacci[i - 1] + fibonacci[i - 2];
	}

	/* pthred create */
	for (thread = 0; thread < thread_count; thread++)
	{
		ids[thread] = thread;
		pthread_create(&thread_handles[thread], NULL, equation, &ids[thread] /*(void *)thread*/);
	}

	for (thread = 0; thread < thread_count; thread++)
		pthread_join(thread_handles[thread], NULL);

	/* Stop timeming */
	double workTimeParallel = dwalltime() - timetick;
	printf("Tiempo en segundos paralelo: %f\n", workTimeParallel);

	/* Checking results */
	/*int check = 1;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
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
	}*/

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
	free(prom_nums);
	free(thread_handles);
	pthread_barrier_destroy(&barrier);
	pthread_mutex_destroy(&mutex);

	sequential();

	double speedup = workTimeSequential / workTimeParallel;
	double eficiencia = speedup / thread_count;
	printf("Speedup: %f\n", speedup);
	printf("Eficiencia: %f\n", eficiencia);
}

void *equation(void *rank)
{
	int *puntero, my_rank;
	puntero = (int *)rank;
	my_rank = *puntero;
	int my_first_row = my_rank * local_size;
	int my_last_row = (my_rank + 1) * local_size - 1;

	int mins, maxs, temp_thread;

	/*min, max y prom*/
	mins = 41;
	maxs = -1;
	temp_thread = 0;

	for (int i = my_first_row; i <= my_last_row; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (F[j * N + i] > maxs)
				maxs = F[i];

			if (F[j * N + i] < mins)
				mins = F[i];

			temp_thread += F[j * N + i];

			fibF[j * N + i] = (double)fibonacci[F[j * N + i]];
		}
	}

	// Average
	prom_nums[my_rank] = temp_thread;

	// Competimos una unica vez
	pthread_mutex_lock(&mutex);

	// Maximun
	if (maxs > max_num)
		max_num = maxs;

	// Minimun
	if (mins < min_num)
		min_num = mins;

	pthread_mutex_unlock(&mutex);

	pthread_barrier_wait(&barrier);

	// Calculamos el promedio despues de que todos los hilos hayan terminado, ademas del valor para la funcion
	// Lo calcula cada hilo

	prom_num = 0;
	for (int i = 0; i < thread_count; i++)
	{
		prom_num += prom_nums[i];
	}
	prom_num = (prom_num / (N * N));

	val = ((max_num * min_num) / (prom_num));

	/* Multiply Matrices */

	int i, j, k, p, q, r;
	int in, jn, ink, jnk, inkpn, jnkqn, inkpnq;

	// A mul B
	for (i = my_first_row; i <= my_last_row; i += bs)
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

	//  AB mul C
	for (i = my_first_row; i <= my_last_row; i += bs)
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

	//  D mul fibF
	for (i = my_first_row; i <= my_last_row; i += bs)
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

	/* Suma de matrices de resABC y resDF y multiplicacion con val */ // CORREJIR
	for (int i = my_first_row; i <= my_last_row; i++)
	{
		in = i * N;
		for (int j = 0; j < N; j++)
		{
			R[in + j] = (R[in + j] + resDF[in + j]) * val;
		}
	}

	pthread_exit(0);
}

void sequential()
{
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
