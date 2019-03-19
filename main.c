#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define TEMP 0.001
#define MSQUARED 1.0
#define XSIZE 100
#define YSIZE 100

#define PI 3.14159265

// Box is periodic in x-direction, boundary is at y == 0, impose Dirichlet (phi=0) boundary condition at y == YSIZE - 1

double phi[XSIZE][YSIZE];

double unitrand() {
	return random() / (double) RAND_MAX;
}

void init() {
	int x, y;

	for (x = 0; x < XSIZE; x++)
		for (y = 0; y < YSIZE; y++)
			phi[x][y] = unitrand();
}

// Generate a random number from a normal distribution of standard deviation sigma
// This uses the Box-Muller method, maybe could be improved for speed

double normal(double sigma) {
	double u = unitrand();
	double v = unitrand();

	return(sigma * sqrt(-2.0 * log(u)) * cos(2.0 * PI * v));
}

void update(double temperature) {
	int x = (int) (unitrand() * XSIZE);
	int y = (int) (unitrand() * YSIZE);

	double c =  phi[x][y];
	double l, r, u, d;

	if (x != 0)
		l = phi[x-1][y];
	else
		l = phi[XSIZE - 1][y];
	if (x < XSIZE -1)
		r = phi[x+1][y];
	else
		r = phi[0][y];

	if (y < YSIZE - 1)
		u = phi[x][y + 1];
	else
		u = 0;
	if (y != 0)
		d = phi[x][y-1];
	else
		d = 2 * c - u;	// This makes the second derivative vanish

	double delta = normal(0.1 * (u - d + l - r));

	double deltaE = (l + r + u + d + (2.0 * MSQUARED + 8) * c) * delta - (-4 + MSQUARED) * delta * delta;
	if (unitrand() < exp(- deltaE / temperature)) {
		putchar('y');
		phi[x][y] += delta;
	} else {
		putchar('n');
	}
}

void show() {
	int x, y;
	FILE *h = fopen("/tmp/output.data", "w");
	if (!h) {
		printf("Cannot open pipe!\n");
		return;
	}

//	fprintf(h, "set xrange [1:%d]\nset yrange [1:%d]\nset grid\nset hidden3d\n$grid << EOD\n", XSIZE, YSIZE);
	for (x = 0; x < XSIZE; x++) {
		for (y = 0; y < YSIZE; y++)
			fprintf(h, "%f ", phi[x][y]);
		fprintf(h, "\n");
	}
//	fprintf(h, "EOD\nsplot '$grid' matrix with lines notitle\n");
	fclose(h);
}

void print_boundary() {
	int x;
	for (x = 0; x < XSIZE; x++)
		printf("%f\n", phi[x][0]);
}

int main()
{
	long t;
	init();
	for (t = 0L; t < 100000000L; t++)
		update(TEMP);
	show();
	print_boundary();
	return 0;
}
