#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define COORDINATOR 0
#define MAX_THREADS 4

int main(int argc, char *argv[])
{
	int i, j, k, numProcs, rank, n, bs, stripSize, check = 1, threads, provide;
	double fnum;
	int localSum = 0, sum = 0, localMin = 41, min = 41, localMax = -1, max = -1;
	double *a, *b, *c, *d, *ab, *dff, *ff, *r;
	int *f, *fib;
	MPI_Status status;
	double commTimes[6], commTime, totalTime;

	/* Check command line parameters */
	if ((argc != 4) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
	{
		printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
		printf("cantidad de argumentos:%i", argc);
		exit(1);
	}

	threads = atoi(argv[3]);

	// MPI_Init(&argc, &argv);
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provide);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (n % numProcs != 0)
	{
		printf("El tama�o de la matriz debe ser multiplo del numero de procesos.\n");
		exit(1);
	}

	// calcular porcion de cada worker
	stripSize = n / numProcs;

	// Reservar memoria
	if (rank == COORDINATOR)
	{
		a = (double *)malloc(sizeof(double) * n * n);
		ab = (double *)calloc(sizeof(double), n * n);
		d = (double *)malloc(sizeof(double) * n * n);
		f = (int *)malloc(sizeof(int) * n * n);
		dff = (double *)calloc(sizeof(double), n * n);
		r = (double *)calloc(sizeof(double), n * n);
	}
	else
	{
		a = (double *)malloc(sizeof(double) * n * stripSize);
		ab = (double *)calloc(sizeof(double), n * stripSize);
		d = (double *)malloc(sizeof(double) * n * stripSize);
		f = (int *)malloc(sizeof(int) * n * stripSize);
		dff = (double *)calloc(sizeof(double), n * stripSize);
		r = (double *)calloc(sizeof(double), n * stripSize);
	}

	b = (double *)malloc(sizeof(double) * n * n);
	c = (double *)malloc(sizeof(double) * n * n);
	ff = (double *)malloc(sizeof(double) * n * n);
	fib = (int *)malloc(sizeof(int) * 40);

	// inicializar datos
	if (rank == COORDINATOR)
	{
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				a[i * n + j] = 1;

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				b[j * n + i] = 1;

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				c[j * n + i] = 1;

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				d[i * n + j] = 1;

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				f[j * n + i] = 1;
	}

	// inicializar datos de fib
	// lo hacen todos

	///////////////////////////////////
	////////////FIBONACCI//////////////
	///////////////////////////////////

	commTimes[0] = MPI_Wtime();

	fib[0] = 0;
	fib[1] = 1;

	for (int i = 2; i < 40; i++)
	{
		fib[i] = fib[i - 1] + fib[i - 2];
	}

	commTimes[1] = MPI_Wtime();

	///////////////////////////////////
	///////MAXIMO,MINIMO,PROMEDIO//////
	///////////////////////////////////

	/* distribuir datos*/
	MPI_Bcast(b, n * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
	MPI_Bcast(c, n * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	MPI_Scatter(a, stripSize * n, MPI_DOUBLE, a, stripSize * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
	MPI_Scatter(d, stripSize * n, MPI_DOUBLE, d, stripSize * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
	MPI_Scatter(r, stripSize * n, MPI_DOUBLE, r, stripSize * n, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	MPI_Scatter(f, n * stripSize, MPI_INT, f, n * stripSize, MPI_INT, COORDINATOR, MPI_COMM_WORLD);
	// MPI_Scatter(ff, n * stripSize, MPI_DOUBLE, ff, n * stripSize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

	/* ¿Como debo calcular los tiempos? */
	commTimes[2] = MPI_Wtime();

/* computar calculos parciales */
#pragma omp parallel
	{

#pragma omp parallel for reduction(+                                               \
								   : localSum) reduction(max                       \
														 : localMax) reduction(min \
																			   : localMin)
		for (i = 0; i < stripSize * n; i++)
		{
			localSum += f[i];

			if (f[i] > localMax)
				localMax = f[i];

			if (f[i] < localMin)
				localMin = f[i];

			ff[i + stripSize * n * rank] = (double)fib[f[i]];
		}

#pragma omp single
		{

			commTimes[3] = MPI_Wtime();

			// recolectar resultados parciales

			MPI_Allreduce(&localSum, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&localMin, &min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&localMax, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allgather(ff + stripSize * n * rank, n * stripSize, MPI_DOUBLE, ff, n * stripSize, MPI_DOUBLE, MPI_COMM_WORLD);

			// Calculo de fnum
			// lo hacen todos los procesos

			commTimes[4] = MPI_Wtime();
		}

		fnum = (min * max) / (sum / (n * n));

		///////////////////////////////////
		//////////MULTIPLICACION///////////
		///////////////////////////////////

		int p, q, s;
		int in, jn, ink, jnk, inkpn, jnkqn, inkpnq;

#pragma omp for private(j, k, p, q, s, in, jn, ink, jnk, inkpn, jnkqn, inkpnq) nowait
		for (i = 0; i < stripSize; i += bs)
		{
			in = i * n;
			for (j = 0; j < n; j += bs)
			{
				jn = j * n;
				for (k = 0; k < n; k += bs)
				{
					ink = in + k;
					jnk = jn + k;
					for (p = 0; p < bs; p++)
					{
						inkpn = ink + p * n;
						for (q = 0; q < bs; q++)
						{
							inkpnq = inkpn + q;
							jnkqn = jnk + q * n;
							for (s = 0; s < bs; s++)
							{
								ab[inkpnq] += a[inkpn + s] * b[jnkqn + s];
							}
						}
					}
				}
			}
		}

#pragma omp for private(j, k, p, q, s, in, jn, ink, jnk, inkpn, jnkqn, inkpnq) nowait
		for (i = 0; i < stripSize; i += bs)
		{
			in = i * n;
			for (j = 0; j < n; j += bs)
			{
				jn = j * n;
				for (k = 0; k < n; k += bs)
				{
					ink = in + k;
					jnk = jn + k;
					for (p = 0; p < bs; p++)
					{
						inkpn = ink + p * n;
						for (q = 0; q < bs; q++)
						{
							inkpnq = inkpn + q;
							jnkqn = jnk + q * n;
							for (s = 0; s < bs; s++)
							{
								r[inkpnq] += ab[inkpn + s] * c[jnkqn + s];
							}
						}
					}
				}
			}
		}

#pragma omp for private(j, k, p, q, s, in, jn, ink, jnk, inkpn, jnkqn, inkpnq) nowait
		for (i = 0; i < stripSize; i += bs)
		{
			in = i * n;
			for (j = 0; j < n; j += bs)
			{
				jn = j * n;
				for (k = 0; k < n; k += bs)
				{
					ink = in + k;
					jnk = jn + k;
					for (p = 0; p < bs; p++)
					{
						inkpn = ink + p * n;
						for (q = 0; q < bs; q++)
						{
							inkpnq = inkpn + q;
							jnkqn = jnk + q * n;
							for (s = 0; s < bs; s++)
							{
								dff[inkpnq] += d[inkpn + s] * ff[jnkqn + s];
							}
						}
					}
				}
			}
		}

#pragma omp for
		for (i = 0; i < stripSize * n; i++)
		{
			r[i] += (dff[i] * fnum);
		}

#pragma omp sigle
		{

			commTimes[5] = MPI_Wtime();

			MPI_Gather(r, n * stripSize, MPI_DOUBLE, r, n * stripSize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

			commTimes[6] = MPI_Wtime();
		}

	} // FIN pragma omp parallels

	MPI_Finalize();

	if (rank == COORDINATOR)
	{

		// Check results
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				check = check && (r[i * n + j] == n * n + n);
			}

		if (check)
		{
			printf("r resultado correcto\n");
		}
		else
		{
			printf("r resultado erroneo\n");
		}

		totalTime = commTimes[6] - commTimes[0];
		commTime = (commTimes[2] - commTimes[1]) + (commTimes[4] - commTimes[3]) + (commTimes[6] - commTimes[5]);
		printf("Multiplicacion de matrices (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n", n, totalTime, commTime);
	}

	free(a);
	free(b);
	free(c);
	free(d);
	free(ab);
	free(dff);
	free(ff);
	free(r);
	free(f);
	free(fib);

	return 0;
}
