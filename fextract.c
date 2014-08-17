/**
 * File:		fextract.c
 * Author:	Tobias Ramforth
 * Date:		09.10.2007
 *
 * Description:	feature extraction
 **/
#include <string.h>
#include <limits.h>
#include <math.h>

#include <mx/histogram.h>
#include <rs/memory.h>
#include <rs/messages.h>

#include <im/basics.h>
#include <im/errors.h>
#include <im/messages.h>
#include <im/data.h>
#include <im/filter.h>
#include <im/transform.h>
#include <im/statistics.h>

#define PEN_KERNEL
#include "basics.h"
#include "messages.h"
#include "errors.h"
#include "data.h"
#include "utils.h"
#include "segmentation.h"

#define BASELINE_DIST_CON 15
#define BASELINE_DIST 10

#define __trunc(x) x < 0 ? ceil(x) : floor(x) /* TODO function in mx? */
#define atand(x)   (atan(x) / M_PI * 180.0)

extern int fextract_verbose;


/* Berechnung einer kubischen Approximations-Spline */
/* aus Reinsch67:SBS                               */
void smooth_spline(
		mx_real_t *result,
		mx_real_t *xm,
		mx_real_t *ym,
		int num,
		mx_real_t s,
		mx_real_t p)
	{
	mx_real_t *x, *y, *dy, *a, *b, *c, *d;
	mx_real_t *r, *r1, *r2, *t, *t1, *u, *v;
	int i, n1, n2, m1, m2;
	mx_real_t e, f, f2, g, h, fxx;

	x = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	y = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	dy = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	a = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	b = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	c = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	d = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	r = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	r1 = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	r2 = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	t = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	t1 = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	u = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));
	v = (mx_real_t*)malloc(sizeof(mx_real_t) * (num + 2));

	for (i = 0; i < num; i++) {
		x[i + 1] = xm[i];
		y[i + 1] = ym[i];
		dy[i + 1] = 1.0;
	}

	n1 = 1;
	n2 = num;
	m1 = n1 - 1;
	m2 = n2 + 1;
	r[m1] = r[n1] = r1[n2] = r2[n2] = r2[m2] = u[m1] = u[n1] = u[n2] = u[m2] = 0;
	m1 = n1 + 1;
	m2 = n2 - 1;
	h = x[m1] - x[n1];
	f = (y[m1] - y[n1]) / h;

	for (i = m1; i <= m2; i++) {
		g = h;
		h = x[i + 1] - x[i];
		e = f;
		f = (y[i + 1] - y[i]) / h;
		a[i] = f-e;
		t[i] = 2 * (g + h) / 3;
		t1[i] = h / 3;
		r2[i] = dy[i - 1] / g;
		r[i] = dy[i + 1] / h;
		r1[i] = -dy[i] / g - dy[i] / h;
	}

	for (i = m1; i <= m2; i++) {
		b[i] = r[i] * r[i]+ r1[i] * r1[i] + r2[i] * r2[i];
		c[i] = r[i] * r1[i + 1] + r1[i] * r2[i + 1];
		d[i] = r[i] * r2[i + 2];
	}
	f2 = -s;

	while (1) {
		for (i = m1; i <= m2; i++) {
			r1[i - 1] = f * r[i - 1];
			r2[i - 2] = g * r[i - 2];
			r[i] = 1 / (p * b[i] + t[i] - f * r1[i - 1] - g * r2[i - 2]);
			u[i] = a[i] - r1[i - 1] * u[i - 1] - r2[i - 2] * u[i - 2];
			f = p * c[i] + t1[i] - h * r1[i - 1];
			g = h;
			h = d[i] * p;
		}

		for (i = m2; i >= m1; i--) {
			u[i] = r[i] * u[i] - r1[i] * u[i + 1] - r2[i] * u[i + 2];
		}
		e = h = 0;

		for (i = n1; i <= m2; i++) {
			g = h;
			h = (u[i + 1] - u[i]) / (x[i + 1] - x[i]);
			v[i] = (h - g) * dy[i] * dy[i];
			e = e + v[i] * (h - g);
		}
		g = v[n2] = -h * dy[n2] * dy[n2];
		e = e - g * h;
		g = f2;
		f2 = e * p * p;
		if (f2 >= s || f2 <= g) {
			break;
		}
		f = 0;
		h = (v[m1] - v[n1]) / (x[m1] - x[n1]);

		for (i = m1; i <= m2; i++) {
			g = h;
			h = (v[i + 1] - v[i]) / (x[i + 1] - x[i]);
			g = h - g - r1[i - 1] * r[i - 1] - r2[i - 2] * r[i - 2];
			f = f + g * r[i] * g;
			r[i] = g;
		}
		h = e - p * f;
		if (h <= 0) {
			break;
		}
		/*       p = p+(s-f2)/((sqrt(s/e)+p)*h); */
	}

	for (i = n1; i <= n2; i++) {
		a[i] = y[i] - p * v[i];
		c[i] = u[i];
	}

	for (i = n1; i <= m2; i++) {
		h = x[i + 1] - x[i];

		d[i] = (c[i + 1] - c[i]) / (3 * h);
		b[i] = (a[i + 1] - a[i]) / h - (h * d[i] + c[i]) * h;
	}

	for (i = n1; i <= m2; i++) {
		mx_real_t h;
		int xx;

		xx = x[i];

		for (xx = x[i]; xx <= x[i + 1]; xx++) {
			h = xx - x[i];
			fxx = ((d[i] * h + c[i]) * h + b[i]) * h + a[i];
			result[xx] = fxx;
		}
	}
	}

int compute_baselines_SPLINE(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) bin,
		im_std(plane,_t) plane,
		int *lower,
		int *upper)
	{
	im_std(plane,_t) dbg = NULL;
	int h, i, j, k, x, y;
	int min, max, down, up;
	int *lower_contour = NULL;
	int *upper_contour = NULL;
	int *maxima_x = NULL;
	int *maxima_y = NULL;
	int *minima_x = NULL;
	int *minima_y = NULL;
	mx_real_t *maxx = NULL;
	mx_real_t *maxy = NULL;
	mx_real_t *minx = NULL;
	mx_real_t *miny = NULL;
	int *x_contour_minima = NULL;
	int *y_contour_minima = NULL;
	int *x_contour_maxima = NULL;
	int *y_contour_maxima = NULL;
	int num_contour_minima = 0, num_min_outlier = 0;
	int num_contour_maxima = 0, num_max_outlier = 0, num_maxima = 0, num_minima = 0;
	mx_histogram_t *maxima_histogram = NULL;
	mx_real_t M1 = 0, M2 = 0, m, var, dist;
	mx_real_t t_s, t_qs, y_s, y_qs, u, beta = 0.0, zaehler, nenner;
	int norm, threshold, val;

	char name[80];
	char index[10];
	static int dbg_num = 0;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(bin);
	dimy = im_std(plane,_dimy)(bin);

	lower_contour = (int*)malloc(sizeof(int) * dimx);
	upper_contour = (int*)malloc(sizeof(int) * dimx);

	maxima_x = (int*)malloc(sizeof(int) * dimx);
	maxima_y = (int*)malloc(sizeof(int) * dimx);
	minima_x = (int*)malloc(sizeof(int) * dimx);
	minima_y = (int*)malloc(sizeof(int) * dimx);
	for (i = 0; i < dimx; i++) {
		maxima_x[i] = -1;
		maxima_y[i] = -1;
		minima_x[i] = -1;
		minima_y[i] = -1;
	}

	for (h = 1; h < dimy - 2; h++) {
		for (i = 0; i < dimx; i++) {
			max = -1;
			min = dimy;
			for (j = 0; j < 2; j++) {
				if (bin[h + j][i] == 0) {
					y = j;
					if (y > max) {
						max = y;
					}
					if (y < min) {
						min = y;
					}
				}
			}
			lower_contour[i] = max;
			upper_contour[i] = min;
		}
		for (i = 0, down = 0; i < dimx - 1; i++) {
			if (lower_contour[i] == -1 && lower_contour[i + 1] == 0) {
				down = 1;
			}	else if (down && lower_contour[i] == 0 && lower_contour[i + 1] == 0) {
				down = 1;
			} else if (down && lower_contour[i] == 0 && lower_contour[i + 1] == -1) {
				down = 0;
				num_contour_minima++;
				if (num_contour_minima == 1) {
					x_contour_minima = (int*)malloc(sizeof(int));
					y_contour_minima = (int*)malloc(sizeof(int));
				} else {
					x_contour_minima = (int*)realloc(x_contour_minima, sizeof(int)
							*num_contour_minima);
					y_contour_minima = (int*)realloc(y_contour_minima, sizeof(int)
							*num_contour_minima);
				}
				x_contour_minima[num_contour_minima - 1] = i;
				y_contour_minima[num_contour_minima - 1] = h + lower_contour[i];

				if (minima_x[i] != -1) {
					fprintf(stderr, "mehrere Minima in Spalte %d\n", i);
				}
				minima_x[i] = i;
				minima_y[i] = h + lower_contour[i];
			} else {
				down = 0;
			}
		}

		for (i = 0, up = 0; i < dimx - 1; i++) {
			if (upper_contour[i] == dimy && upper_contour[i + 1] == 1) {
				up = 1;
			} else if (up && upper_contour[i] == 1 && upper_contour[i + 1] == 1) {
				up = 1;
			} else if (up && upper_contour[i] == 1 && upper_contour[i + 1] == dimy) {
				up = 0;
				num_contour_maxima++;
				if (num_contour_maxima == 1) {
					x_contour_maxima = (int*)malloc(sizeof(int));
					y_contour_maxima = (int*)malloc(sizeof(int));
				} else {
					x_contour_maxima = (int*)realloc(x_contour_maxima, sizeof(int)
							*num_contour_maxima);
					y_contour_maxima = (int*)realloc(y_contour_maxima, sizeof(int)
							*num_contour_maxima);
				}
				x_contour_maxima[num_contour_maxima-1] = i;
				y_contour_maxima[num_contour_maxima-1] = h + upper_contour[i];

				if (maxima_x[i] != -1) {
					fprintf(stderr, "mehrere Maxima in Spalte %d\n", i);
				}
				maxima_x[i] = i;
				maxima_y[i] = h + upper_contour[i];
			} else {
				up = 0;
			}
		}
	}
	for (k = 0, y_s = 0.0; k < num_contour_minima ; k++) {
		y_s += y_contour_minima[k];
	}

	if ((num_contour_minima) > 0) {
		y_s = y_s / (double)(num_contour_minima);
	}

	*lower = mx_round(y_s);

	maxima_histogram = mx_histogram_create(0, dimy, 1);
	for (i = 0, norm = 0, val = 0; i < num_contour_maxima; i++) {
		if (*lower - y_contour_maxima[i] > BASELINE_DIST) {
			mx_histogram_update(maxima_histogram, y_contour_maxima[i]);
			val = y_contour_maxima[i];
			norm++;
		}
	}
	if (norm > 1) {
		threshold = mx_histogram_otsu_thresh(&M1, &M2, maxima_histogram);
	} else {
		M2 = val;
	}

	*upper = mx_round(M2);

	for (k = 0, y_s = 0.0, num_min_outlier = 0; k < num_contour_minima; k++) {
		if (y_contour_minima[k] > (*upper + (*lower - *upper)/2))
			y_s += y_contour_minima[k];
		else
			num_min_outlier++;
	}
	if ((num_contour_minima - num_min_outlier) > 0) {
		y_s = y_s / (double)(num_contour_minima - num_min_outlier);
	}

	*lower = mx_round(y_s);

	for (i = 0; i < dimx; i++) {
		if (maxima_x[i] != -1 && (maxima_y[i] > threshold)&&(maxima_y[i] < (*lower
				- (*lower - *upper) / 2))) {
			num_maxima++;
		}
		if (minima_x[i] != -1 && (minima_y[i] < *lower + 2 * BASELINE_DIST) && (minima_y[i]
				> (*upper + (*lower - *upper) / 2))) {
			num_minima++;
		}
	}
	maxx = (mx_real_t*)malloc(sizeof(mx_real_t) * (num_maxima + 2));
	maxy = (mx_real_t*)malloc(sizeof(mx_real_t) * (num_maxima + 2));
	minx = (mx_real_t*)malloc(sizeof(mx_real_t) * (num_minima + 2));
	miny = (mx_real_t*)malloc(sizeof(mx_real_t) * (num_minima + 2));

	/* Erster Punkt: Mittelwert der ersten 5 */
	for (i = 0, k = 0, maxx[0] = 0, maxy[0] = 0; i < dimx; i++) {
		if (maxima_x[i] != -1 && (maxima_y[i] > threshold) && (maxima_y[i] < (*lower
				- (*lower - *upper) / 2))) {
			maxy[0] += maxima_y[i];
			k++;
		}
		if (k >= 5) {
			break;
		}
	}
	maxy[0] /= k;

	/* Letzter Punkt: Mittelwert der letzten 5 */
	for (i = dimx - 1, k = 0, maxx[num_maxima + 1] = dimx - 1, maxy[num_maxima + 1]
			= 0; i >= 0; i--) {
		if (maxima_x[i] != -1 && (maxima_y[i] > threshold) && (maxima_y[i] < (*lower
				- (*lower - *upper) / 2))) {
			maxy[num_maxima + 1] += maxima_y[i];
			k++;
		}
		if (k >= 5) {
			break;
		}
	}
	maxy[num_maxima + 1] /= k;

	/* Erster Punkt: Mittelwert der ersten 5 */
	for (i = 0, k = 0, minx[0] = 0, miny[0] = 0; i < dimx; i++) {
		if (minima_x[i] != -1 && (minima_y[i] < *lower + 2 * BASELINE_DIST) && (minima_y[i]
				> (*upper + (*lower - *upper) / 2))) {
			miny[0] += minima_y[i];
			k++;
		}
		if (k >= 5) {
			break;
		}
	}
	miny[0] /= k;

	/* Letzter Punkt: Mittelwert der letzten 5 */
	for (i = dimx - 1, k = 0, minx[num_minima + 1] = dimx - 1, miny[num_minima + 1]
			= 0; i >= 0; i--) {
		if (minima_x[i] != -1 && (minima_y[i] < *lower + 2 * BASELINE_DIST) && (minima_y[i]
				> (*upper + (*lower - *upper) / 2))) {
			miny[num_minima+1] += minima_y[i];
			k++;
		}
		if (k >= 5) {
			break;
		}
	}
	miny[num_minima + 1] /= k;

	for (i = 0, j = 1; i < dimx; i++) {
		if (maxima_x[i] != -1 && (maxima_y[i] > threshold) && (maxima_y[i] < (*lower
				- (*lower - *upper) / 2))) {
			maxx[j] = maxima_x[i];
			maxy[j++] = maxima_y[i];
		}
	}
	for (i = 0, j = 1; i < dimx; i++) {
		if (minima_x[i] != -1 && (minima_y[i] < *lower + 2 * BASELINE_DIST) && (minima_y[i]
				> (*upper + (*lower - *upper) / 2))) {
			minx[j] = minima_x[i];
			miny[j++] = minima_y[i];
		}
	}

	off_fex_par->upper_base_spline = (mx_real_t*)malloc(sizeof(mx_real_t) * dimx);
	off_fex_par->lower_base_spline = (mx_real_t*)malloc(sizeof(mx_real_t) * dimx);
	smooth_spline(off_fex_par->upper_base_spline, maxx, maxy, num_maxima + 2, 0,
			100000);
	smooth_spline(off_fex_par->lower_base_spline, minx, miny, num_minima + 2, 0,
			100000);

#ifdef DEBUG
	dbg = copy_image(img);
	for (i=0;i<num_maxima+1;i++)
	{
		draw_line(dbg, maxx[i], maxy[i], maxx[i+1], maxy[i+1]);
	}
	for (i=0;i<num_minima+1;i++)
	{
		draw_line(dbg, minx[i], miny[i], minx[i+1], miny[i+1]);
	}
	draw_line(dbg, 0, y_s, bin->width-1, y_s);
	draw_line(dbg, 0, (int)M1, bin->width-1, (int)M1);
	draw_line(dbg, 0, (int)M2, bin->width-1, (int)M2);

	sprintf(index,"%03d.pgm",dbg_num++);
	strcpy(name,"dbg_fex");
	strcat(name,index);
	writeImage(dbg, name);

	free_image(dbg);

	dbg = create_image(img->height+200, img->width, 255);
	for (i=0;i<dbg->width*dbg->height;i++) dbg->pixel[i]=240;
	for (i=0;i<img->width;i++)
	for (j=0;j<img->height;j++)
	dbg->pixel[(j+100)*dbg->width+i] = img->pixel[j*img->width+i];

	for (i=0;i<bin->width-1;i++)
	{
		draw_line(dbg, i, off_fex_par->lower_base_spline[i]+100, i+1, off_fex_par->lower_base_spline[i+1]+100);
		draw_line(dbg, i, off_fex_par->upper_base_spline[i]+100, i+1, off_fex_par->upper_base_spline[i+1]+100);
	}
	sprintf(index,"%03d.pgm",dbg_num++);
	strcpy(name,"dbg_fex_spline");
	strcat(name,index);
	writeImage(dbg, name);

	free_image(dbg);
#endif

	if (maxima_histogram) {
		mx_histogram_destroy(maxima_histogram);
	}

	if (maxima_x) {
		free(maxima_x);
	}

	if (maxima_y) {
		free(maxima_y);
	}

	if (minima_x) {
		free(minima_x);
	}

	if (minima_y) {
		free(minima_y);
	}

	if (maxx) {
		free(maxx);
	}

	if (maxy) {
		free(maxy);
	}

	if (minx) {
		free(minx);
	}

	if (miny) {
		free(miny);
	}

	if (upper_contour) {
		free(upper_contour);
	}

	if (lower_contour) {
		free(lower_contour);
	}

	if (x_contour_maxima) {
		free(x_contour_maxima);
	}

	if (x_contour_minima) {
		free(x_contour_minima);
	}

	if (y_contour_maxima) {
		free(y_contour_maxima);
	}

	if (y_contour_minima) {
		free(y_contour_minima);
	}

	return *lower;
	}

void debug_output(
		im_std(plane,_t) plane,
		char *filename,
		int lb,
		int ub,
		im_std(plane,_t) bin)
	{
	char name1[80];
	char name2[80];
	im_std(plane,_t) dbg, dbg_bw;
	int i, j;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(plane);
	dimy = im_std(plane,_dimy)(plane);

	strcpy(name1, filename);
	strcat(name1, "_dbg.pgm");
	strcpy(name2, filename);
	strcat(name2, "_dbg_bw.pgm");

	fprintf(stderr, "debug: write image %s, %s\n", name1, name2);

	dbg = im_std(plane,_dup)(plane);
	dbg_bw = im_std(plane,_create)(dimx, dimy);

	for (j = 0; j < dimy; j++) {
		for (i = 0; i < dimx; i++) {
			if (bin[j][i] == 0) {
				dbg_bw[j][i] = 0;
			} else {
				dbg_bw[j][i] = IM_GREY_MAX;
			}

			if (j == lb || j == ub) {
				dbg[j][i] = 0;
				dbg[j][i] = 0;

				dbg_bw[j][i] = 0;
				dbg_bw[j][i] = 0;
			}
		}
	}
	/*
	writeImage(dbg, name1);
	writeImage(dbg_bw, name2);
	*/

	im_std(plane,_destroy)(dbg);
	im_std(plane,_destroy)(dbg_bw);
	}

int compute_slope(
		mx_real_t *f,
		int num,
		mx_real_t *y,
		mx_real_t *slope)
	{
	float t_s, t_qs, y_s, y_qs, u, beta = 0.0, zaehler, nenner;
	int i, j, k, m;
	int cnt = 0;

	t_s = 0.0;
	y_s = 0.0;
	t_qs = 0.0;
	y_qs = 0.0;
	u = 0.0;

	for (k = 0; k < num; k++) {
		if (f[k] != -1) {
			t_s += k;
			y_s += f[k];
			t_qs += k * k;
			y_qs += f[k] * f[k];
			u += k * f[k];
			cnt++;
		}
	}

	if (cnt > 0) {
		t_s = t_s / (mx_real_t)(cnt);
		y_s = y_s / (mx_real_t)(cnt);
		t_qs = t_qs / (mx_real_t)(cnt);
		y_qs = y_qs / (mx_real_t)(cnt);
		u = u / (mx_real_t)(cnt);
	}

	zaehler = u - t_s * y_s;
	nenner = t_qs - t_s * t_s;

	if (fabs(nenner) > 0.000001) {
		beta = zaehler / nenner;
	}

	*slope = atand(beta);
	*y = y_s;

	return 0;
	}

int compute_slope2(
		mx_real_t *f,
		int num,
		mx_real_t *y,
		mx_real_t *sin_slope,
		mx_real_t *cos_slope)
	{
	mx_real_t t_s, t_qs, y_s, y_qs, u, beta = 0.0, dy, dx, len;
	int i, j, k, m;
	int cnt = 0;

	t_s = 0.0;
	y_s = 0.0;
	t_qs = 0.0;
	y_qs = 0.0;
	u = 0.0;

	for (k = 0; k < num; k++) {
		if (f[k] != -1) {
			t_s += k;
			y_s += f[k];
			t_qs += k * k;
			y_qs += f[k] * f[k];
			u += k * f[k];
			cnt++;
		}
	}

	if (cnt > 0) {
		t_s = t_s / (mx_real_t)(cnt);
		y_s = y_s / (mx_real_t)(cnt);
		t_qs = t_qs / (mx_real_t)(cnt);
		y_qs = y_qs / (mx_real_t)(cnt);
		u = u / (mx_real_t)(cnt);
	}

	dy = u - t_s * y_s;
	dx = t_qs - t_s * t_s;

	len = sqrt(dx * dx + dy * dy);
	if (len > 0) {
		dx = dx / len;
		dy = dy / len;
	} else {
		dx = 0;
		dy = 0;
	}

	*sin_slope = dy;
	*cos_slope = dx;
	*y = y_s;

	return 0;
	}

pen_fv_t compute_bitmap_fex(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame)
	{
	int i, j, k, top, bottom;
	int cell_size = 25;
	float intensity = 0.0;
	float mean_intensity;
	int num_pixel = 0;
	int dim = 0;
	pen_fv_t fv;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame);
	dimy = im_std(plane,_dimy)(frame);

	dim = off_fex_par->num_cells;
	if (off_fex_par->use_derivatives) {
		dim *= 2;
	}

	fv = pen_fv_create(dim);

	cell_size = dimy / off_fex_par->num_cells;

	top = 0;
	bottom = cell_size;

	for (i = 0; i < off_fex_par->num_cells; i++) {
		num_pixel = dimx * cell_size;
		intensity = 0.0;

		if (bottom > dimy) {
			bottom = dimy;
			num_pixel = dimx * (bottom - top);
		}

		for (j = 0; j < dimx; j++) {
			for (k = top; k < bottom; k++) {
				intensity += frame[k][j];
			}
		}
		mean_intensity = intensity / num_pixel;

		fv->feature[i] = mean_intensity / (mx_real_t)IM_GREY_MAX * off_fex_par->fv_max;

		top = bottom;
		bottom += cell_size;
	}

	return fv;
	}

/* bw transitions */
mx_real_t calc_feature_bw_transitions_mean(
		im_std(plane,_t) frame_bin)
	{
	int i, j;
	int num_bw_transitions;
	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame_bin);
	dimy = im_std(plane,_dimy)(frame_bin);

	for (i = 0, num_bw_transitions = 0; i < dimx; i++) {
		for (j = 0; j < dimy - 1; j++) {
			if (frame_bin[j][i] == 0
					&& frame_bin[j + 1][i] != 0) {
				num_bw_transitions++;
			}
		}
	}

	return num_bw_transitions / dimx;
	}

/* second line to baseline  */
void calc_feature_second_transition_mean(
        mx_real_t* mean_vals,
        im_std(plane,_t) frame_bin)
        {
        int i, j, k;
        int num_bw_transitions;
        int dimx, dimy;
        mx_real_t sum_w = 0.0;
        mx_real_t mean = 0.0;

        dimx = im_std(plane,_dimx)(frame_bin);
        dimy = im_std(plane,_dimy)(frame_bin);

        for (i = 0, num_bw_transitions = 0; i < dimx; i++)
        {
                sum_w = 0.0;
                mean = 0.0;

                for (j = 0, num_bw_transitions = 0; j < dimy - 1; j++) 
                {
                        if (frame_bin[j][i] == 0 && frame_bin[j + 1][i] != 0) 
                        {
                                num_bw_transitions++;
                        }
                }

                if (num_bw_transitions >= 3)
                {

		for (j = 0, num_bw_transitions = 0; j < dimy - 1; j++)
                {
                        if (frame_bin[j][i] == 0 && frame_bin[j + 1][i] != 0)
                        {
                                num_bw_transitions++;
                                if (num_bw_transitions == 2)
                                {
                                        for ( k = j; frame_bin[k][i]==0 && k>0; k--)
                                        {
                                                sum_w++;
                                                mean = mean + k;
                                        }
                                        break;
                                }
                        }
                }
                }

                if ( sum_w > 0.0 )
                {
                        mean_vals[i] = mean / sum_w;
                } else {
                        mean_vals[i] = -1;
                }
        }
        }

mx_real_t calc_width_proportion(im_std(plane,_t) frame_bin)
        {
        int i, j;
        int num_pixel;
		int min_width;
		int max_width;
        int dimx, dimy;
		mx_real_t ratio = 0;

        dimx = im_std(plane,_dimx)(frame_bin);
        dimy = im_std(plane,_dimy)(frame_bin);

		min_width = dimx;
		max_width = 0;
		for (i=0; i < dimy; i++)
		{
			for(j=0,num_pixel = 0; j<dimx; j++)
			{
				if(frame_bin[i][j] == 0)
					num_pixel++;
			}
			if(num_pixel > max_width)
				max_width = num_pixel;
			if(num_pixel < min_width && num_pixel != 0 )
				min_width = num_pixel;
		}
		if(min_width <= max_width)
		ratio = min_width / max_width;
		return ratio;
        }

/* bw transitions 2 */
mx_real_t calc_feature_min_max_bw_transitions_mean(
		im_std(plane,_t) frame_bin_thin)
	{
	int num_bw_transitions;
	int smallest_num_bw_transitions;
	int largest_num_bw_transitions;
	int total_num_bw_transitions;

	int i, j;
	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame_bin_thin);
	dimy = im_std(plane,_dimy)(frame_bin_thin);

	smallest_num_bw_transitions = dimy;
	largest_num_bw_transitions = 0;
	total_num_bw_transitions = 0;
	for (i = 0; i < dimx; i++) {
		for (j = 0, num_bw_transitions = 0; j < dimy - 1; j++) {
			if (frame_bin_thin[j][i] == 0 && frame_bin_thin[j + 1][i] != 0) {
				num_bw_transitions++;
			}
			if (frame_bin_thin[j][i] != 0 && frame_bin_thin[j + 1][i] == 0) {
				num_bw_transitions++;
			}
		}
		total_num_bw_transitions += num_bw_transitions;
		if (num_bw_transitions < smallest_num_bw_transitions) {
			smallest_num_bw_transitions = num_bw_transitions;
		}
		if (num_bw_transitions > largest_num_bw_transitions) {
			largest_num_bw_transitions = num_bw_transitions;
		}
	}

	return (mx_real_t)total_num_bw_transitions / dimx;
	}

/* mean intensity */
mx_real_t calc_feature_mean_intensity(
		im_std(plane,_t) frame,
		im_std(grey,_t) background,
		pen_fextractpar_t off_fex_par)
	{
	int i, j;
	int num_bw_transitions;
	mx_real_t intensity = 0.0;
	mx_real_t mean_intensity;


	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame);
	dimy = im_std(plane,_dimy)(frame);

	for (i = 0, mean_intensity = 0.0; i < dimx; i++) {
		for (intensity = 0.0, j = 0; j < dimy; j++) {
			intensity += frame[j][i];
		}
		mean_intensity += intensity / dimy;
	}

	if (dimx == 0 || background == 0) {
		return 0;
	} else {
		return (mean_intensity / dimx) * off_fex_par->fv_max / background;
	}
	}

/* mean values */
mx_real_t *calc_feature_mean_values(
		mx_real_t* mean_vals,
		im_std(plane,_t) frame_bin)
	{
	int i, j;
	mx_real_t *w;
	mx_real_t mean = 0.0;
	mx_real_t sum_w = 0.0;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame_bin);
	dimy = im_std(plane,_dimy)(frame_bin);

	w = (mx_real_t*)calloc(dimy, sizeof(mx_real_t));

	for (i = 0; i < dimx; i++) {
		sum_w = 0.0;
		for (j = 0; j < dimy; j++) {
			if (frame_bin[j][i] == 0) {
				w[j] = 1.0;
			}	else {
				w[j] = 0.0;
			}
			sum_w += w[j];
		}

		for (mean = 0.0, j = 0; j < dimy; j++) {
			mean += j * w[j];
		}

		if (sum_w > 0.0) {
			mean_vals[i] = mean / sum_w;
		}	else {
			mean_vals[i] = -1.0;
		}
	}
	free (w);
	}

/* minimum values */
void calc_feature_min_values(
		mx_real_t* min_vals,
		int *smallest_min_y,
		im_std(plane,_t) frame_bin)
	{
	int i, j;
	int min_y;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame_bin);
	dimy = im_std(plane,_dimy)(frame_bin);

	if (smallest_min_y) {
		*smallest_min_y = dimy;
	}
	for (i = 0; i < dimx; i++) {
		min_y = dimy;
		for (j = 0; j < dimy; j++) {
			if (frame_bin[j][i] == 0) {
				if (j < min_y) {
					min_y = j;
				}
				if (smallest_min_y) {
					if (j < *smallest_min_y) {
						*smallest_min_y = j;
					}
				}
			}
		}

		if (min_y < dimy) {
			min_vals[i] = min_y;
		}	else {
			min_vals[i] = -1;
		}
	}
	}

/* maximum values */
void calc_feature_max_values(
		mx_real_t* max_vals,
		int *largest_max_y,
		im_std(plane,_t) frame_bin)
	{
	int i, j;
	int max_y;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame_bin);
	dimy = im_std(plane,_dimy)(frame_bin);

	if (largest_max_y) {
		*largest_max_y = 0;
	}
	for (i = 0; i< dimx; i++) {
		max_y = 0;
		for (j = 0; j < dimy; j++) {
			if (frame_bin[j][i] == 0) {
				if (j > max_y) {
					max_y = j;
				}
				if (largest_max_y) {
					if (j > *largest_max_y) {
						*largest_max_y = j;
					}
				}
			}
		}
		if (max_y > 0) {
			max_vals[i] = max_y;
		}	else {
			max_vals[i] = -1.0;
		}
	}
	}

/* min max mena intensity */
mx_real_t calc_feature_min_max_mean_intensity(
		im_std(plane,_t) frame,
		int smallest_min_y,
		int largest_max_y,
		im_std(grey,_t) background,
		pen_fextractpar_t off_fex_par)
	{
	int i, j;
	int num_bw_transitions;
	mx_real_t intensity = 0.0;
	mx_real_t mean_intensity;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame);
	dimy = im_std(plane,_dimy)(frame);

	for (i = 0, mean_intensity = 0.0; i < dimx; i++) {
		for (intensity = 0.0, j = smallest_min_y; j < largest_max_y; j++) {
			intensity += frame[j][i];
		}
		if (largest_max_y > smallest_min_y) {
			mean_intensity += intensity / (largest_max_y - smallest_min_y);
		}	else {
			mean_intensity = 0.0;
		}
	}

	return (mean_intensity / dimx) * off_fex_par->fv_max / background;
	}

/* min max mena intensity 2 */
mx_real_t calc_feature_min_max_mean_intensity2(
		im_std(plane,_t) frame,
		int smallest_min_y,
		int largest_max_y)
	{
	int i, j;
	int num_bw_transitions;
	mx_real_t intensity = 0.0;
	mx_real_t mean_intensity;

	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame);
	dimy = im_std(plane,_dimy)(frame);

	for (i = 0, mean_intensity = 0.0; i < dimx; i++) {
		for (intensity = 0.0, j = smallest_min_y; j < largest_max_y; j++) {
			if (frame[j][i] == 0) {
				intensity ++;
			}
		}
		if (largest_max_y > smallest_min_y) {
			mean_intensity += intensity / (largest_max_y - smallest_min_y);
		}	else {
			mean_intensity = 0.0;
		}
	}

	return mean_intensity;
	}

/* mean intensity */
mx_real_t calc_feature_mean_intensity2(
		im_std(plane,_t) frame)
	{
	int i, j;
	int num_bw_transitions;
	mx_real_t intensity = 0.0;
	mx_real_t mean_intensity;


	int dimx, dimy;

	dimx = im_std(plane,_dimx)(frame);
	dimy = im_std(plane,_dimy)(frame);

	for (i = 0, mean_intensity = 0.0; i < dimx; i++) {
		for (intensity = 0.0, j = 0; j < dimy; j++) {
			if (frame[j][i] == 0) {
				intensity++;
			}
		}
		mean_intensity += intensity / dimy;
	}

	return mean_intensity;
	}

/* NEW NEW NEW NEW */
pen_fv_t compute_geom_fex(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame,
		im_std(plane,_t) frame_bin,
		im_std(plane,_t) frame_thin,
		int left,
		int baseline,
		int lower_baseline,
		int upper_baseline,
		im_std(grey,_t) background,
		int mode)
	{
	int dim;

	int i, j, k, l, f_num = 0;
	mx_real_t intensity = 0.0, mean_intensity;
	mx_real_t f = 0.0;
	mx_real_t mean = 0.0;
	int num_sw_trans = 0, total_num_sw_trans = 0, smallest_num_sw_trans = 0,
			largest_num_sw_trans = 0;
	int num_black_pixels = 0;
	int min_y = 0, smallest_min_y = 0;
	int max_y, largest_max_y = 0;
	mx_real_t *w = NULL;
	mx_real_t *mean_vals = NULL;
	mx_real_t *min_vals = NULL;
	mx_real_t *max_vals = NULL;
	mx_real_t y_val, slope_sin, slope_cos;
	mx_real_t slope_val;
	mx_real_t sum_w = 0.0;
	pen_fv_t fv;
	mx_real_t upper_base;
	mx_real_t lower_base;

	int dimx, dimy;
	int thin_dimx, thin_dimy;

	if (frame) {
		dimx = im_std(plane,_dimx)(frame);
		dimy = im_std(plane,_dimy)(frame);
	} else if (frame_bin) {
		dimx = im_std(plane,_dimx)(frame_bin);
		dimy = im_std(plane,_dimy)(frame_bin);
	}

	if (frame_thin) {
		thin_dimx = im_std(plane,_dimx)(frame_thin);
		thin_dimy = im_std(plane,_dimy)(frame_thin);
	}

	max_y = dimy - 1;

	mean_vals = (mx_real_t*)calloc(dimx, sizeof(mx_real_t));
	min_vals = (mx_real_t*)calloc(dimx, sizeof(mx_real_t));
	max_vals = (mx_real_t*)calloc(dimx, sizeof(mx_real_t));

	/* dependent on the mode we have different dimensions */
	switch(mode) {
		case 1:
			dim = 12;
			break;
		case 2:
			dim = 10;
			break;
		case 3:
			dim = 10;
			break;
		case 4:
			dim = 9;
			break;
		case 5:
			dim = 9;
			break;
		case 6:
			dim = 12;
			break;
	}

	if (off_fex_par->use_derivatives) {
		dim *= 2;
	}

	fv = pen_fv_create(dim);

	/* local geometric features */
	f_num = 0;

	switch(mode) {
		case 1:
			break;
		case 2:
			break;
		case 3:
		case 4:
		case 5:
		case 6:
			if (off_fex_par->spline) {
				upper_base = (off_fex_par->upper_base_spline[left]
						+ off_fex_par->upper_base_spline[left + dimx - 1]) / 2;
				lower_base = (off_fex_par->lower_base_spline[left]
						+ off_fex_par->lower_base_spline[left + dimx - 1]) / 2;
			} else {
				upper_base = upper_baseline;
				lower_base = lower_baseline;
			}
			if (lower_base - upper_base <= 0) {
				upper_base = upper_baseline;
				lower_base = lower_baseline;
			}
			break;
	}

	switch(mode) {
		/*case 1:
			fv->feature[0] = calc_feature_bw_transitions_mean(frame_bin);

			fv->feature[1] = calc_feature_mean_intensity(frame, background, off_fex_par);

			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			fv->feature[2] = (mx_real_t) (baseline - y_val);
			fv->feature[3] = slope_val;

			calc_feature_min_values(min_vals, NULL, frame_bin);
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			fv->feature[4] = (mx_real_t) (baseline - y_val);
			fv->feature[5] = slope_val;

			calc_feature_max_values(max_vals, NULL, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			fv->feature[6] = (mx_real_t)(baseline - y_val);
			fv->feature[7] = slope_val;
			break; */
		case 2:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_bin);

			/* 2 & 3 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			fv->feature[1] = (mx_real_t) (baseline - y_val);
			fv->feature[2] = slope_val;

			/* 4 & 5 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			fv->feature[3] = (mx_real_t) (baseline - y_val);
			fv->feature[4] = slope_val;

			/* 6 & 7 & 8 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			fv->feature[5] = (mx_real_t) (baseline - y_val);
			fv->feature[6] = slope_val;
			fv->feature[7] = (mx_real_t) (largest_max_y - smallest_min_y);

			if (largest_max_y == 0) {
				largest_max_y = dimy;
			}

			/* 9 */
			fv->feature[8] = calc_feature_min_max_mean_intensity(frame, smallest_min_y, largest_max_y, background, off_fex_par);

			/* 10 */
			fv->feature[9] = calc_feature_mean_intensity(frame, background, off_fex_par);
			
			/* 11 12 von Klaus 
			calc_feature_second_transition_mean(mean_vals, frame_bin);
           	compute_slope(mean_vals, dimx, &y_val, &slope_val);
           	fv->feature[10] = (mx_real_t) (baseline - y_val);
            fv->feature[11] = slope_val;
			*/

			/* new feature from Yuanchen 
			fv->feature[9] = calc_width_proportion(frame_bin);
			*/
			
            break;

		case 3:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_bin);

			/* 2 & 3 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[1] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[1] = 0;
			}
			fv->feature[2] = slope_val;

			/* 4 & 5 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			if (smallest_min_y == dimy) {
				smallest_min_y = 0;
			}
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[3] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[3] = 0;
			}
			fv->feature[4] = slope_val;

			/* 6 & 7 & 8 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[5] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[5] = 0;
			}
			fv->feature[6] = slope_val;
			if (lower_base - upper_base != 0) {
				fv->feature[7] = (mx_real_t) (mx_real_t) (largest_max_y - smallest_min_y) / (lower_base - upper_base);
			} else {
				fv->feature[7] = 0;
			}

			if (largest_max_y == 0) {
				largest_max_y = dimy;
			}

			/* 9 */
			fv->feature[8] = calc_feature_min_max_mean_intensity(frame, smallest_min_y, largest_max_y, background, off_fex_par);

			/* 10 */
			fv->feature[9] = calc_feature_mean_intensity(frame, background, off_fex_par);
			break;
		case 4:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_thin);

			/* 2 & 3 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[1] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[1] = 0;
			}
			fv->feature[2] = slope_val;

			/* 4 & 5 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			if (smallest_min_y == dimy) {
				smallest_min_y = 0;
			}
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[3] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[3] = 0;
			}
			fv->feature[4] = slope_val;

			/* 6 & 7 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[5] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[5] = 0;
			}
			fv->feature[6] = slope_val;

			if (largest_max_y == 0) {
				largest_max_y = dimy;
			}

			/* 8 */
			fv->feature[7] = calc_feature_min_max_mean_intensity(frame, smallest_min_y, largest_max_y, background, off_fex_par);

			/* 9 */
			fv->feature[8] = calc_feature_mean_intensity(frame, background, off_fex_par);
			break;
		case 5:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_thin);

			/* 2 & 3 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[1] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[1] = 0;
			}
			fv->feature[2] = slope_val;

			/* 4 & 5 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			if (smallest_min_y == dimy) {
				smallest_min_y = 0;
			}
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[3] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[3] = 0;
			}
			fv->feature[4] = slope_val;

			/* 6 & 7 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			if (lower_base - upper_base != 0) {
				fv->feature[5] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[5] = 0;
			}
			fv->feature[6] = slope_val;

			if (largest_max_y == 0) {
				largest_max_y = dimy;
			}

			/* 8 */
			fv->feature[7] = calc_feature_min_max_mean_intensity2(frame_bin, smallest_min_y, largest_max_y);

			/* 9 */
			fv->feature[8] = calc_feature_mean_intensity2(frame_bin);
			break;
		case 6:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_thin);

			/* 2 & 3 & 4 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope2(mean_vals, dimx, &y_val, &slope_sin, &slope_cos);
			if (lower_base - upper_base != 0) {
				fv->feature[1] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[1] = 0;
			}
			fv->feature[2] = slope_sin;
			fv->feature[3] = slope_cos;

			/* 5 & 6 & 7 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			if (smallest_min_y == dimy) {
				smallest_min_y = 0;
			}
			compute_slope2(min_vals, dimx, &y_val, &slope_sin, &slope_cos);
			if (lower_base - upper_base != 0) {
				fv->feature[4] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[4] = 0;
			}
			fv->feature[5] = slope_sin;
			fv->feature[6] = slope_cos;

			/* 8 & 9 & 10 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope2(max_vals, dimx, &y_val, &slope_sin, &slope_cos);
			if (lower_base - upper_base != 0) {
				fv->feature[7] = (mx_real_t) (lower_base - y_val) / (lower_base - upper_base);
			} else {
				fv->feature[7] = 0;
			}
			fv->feature[8] = slope_sin;
			fv->feature[9] = slope_cos;

			/* 11 */
			fv->feature[10] = calc_feature_min_max_mean_intensity(frame_bin, smallest_min_y, largest_max_y, background, off_fex_par);

			/* 12 */
			fv->feature[11] = calc_feature_mean_intensity(frame_bin, background, off_fex_par);
			break;

		case 1:
			/* 1 */
			fv->feature[0] = calc_feature_min_max_bw_transitions_mean(frame_bin);

			/* 2 & 3 */
			calc_feature_mean_values(mean_vals, frame_bin);
			compute_slope(mean_vals, dimx, &y_val, &slope_val);
			fv->feature[1] = (mx_real_t) (baseline - y_val);
			fv->feature[2] = slope_val;

			/* 4 & 5 */
			calc_feature_min_values(min_vals, &smallest_min_y, frame_bin);
			compute_slope(min_vals, dimx, &y_val, &slope_val);
			fv->feature[3] = (mx_real_t) (baseline - y_val);
			fv->feature[4] = slope_val;

			/* 6 & 7 & 8 */
			calc_feature_max_values(max_vals, &largest_max_y, frame_bin);
			compute_slope(max_vals, dimx, &y_val, &slope_val);
			fv->feature[5] = (mx_real_t) (baseline - y_val);
			fv->feature[6] = slope_val;
			fv->feature[7] = (mx_real_t) (largest_max_y - smallest_min_y);

			if (largest_max_y == 0) {
				largest_max_y = dimy;
			}

			/* 9 */
			fv->feature[8] = calc_feature_min_max_mean_intensity(frame, smallest_min_y, largest_max_y, background, off_fex_par);

			/* 10 */
			fv->feature[9] = calc_feature_mean_intensity(frame, background, off_fex_par);
			
			/* 11 12 von Klaus */
			calc_feature_second_transition_mean(mean_vals, frame_bin);
           	compute_slope(mean_vals, dimx, &y_val, &slope_val);
           	fv->feature[10] = (mx_real_t) (baseline - y_val);
            fv->feature[11] = slope_val;

			/* new feature from Yuanchen 
			fv->feature[9] = calc_width_proportion(frame_bin);
			*/
            break;
	}

	/* ... free stuff ... */
	if (w) {
		free (w);
	}

	free(mean_vals);
	free(min_vals);
	free(max_vals);

	return fv;
	}

pen_fv_t compute_geom_fex_1(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame,
		im_std(plane,_t) frame_bin,
		int baseline,
		im_std(grey,_t) background)
	{
	return compute_geom_fex(
			off_fex_par,
			frame,
			frame_bin,
			NULL,
			0,
			baseline,
			0,
			0,
			background,
			1);
	}

pen_fv_t compute_geom_fex_2(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame,
		im_std(plane,_t) frame_bin,
		int baseline,
		int background)
	{
	return compute_geom_fex(
			off_fex_par,
			frame,
			frame_bin,
			NULL,
			0,
			baseline,
			0,
			0,
			background,
			2);
}

pen_fv_t compute_geom_fex_3(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame,
		im_std(plane,_t) frame_bin,
		int left,
		int lower_baseline,
		int upper_baseline,
		im_std(grey,_t) background)
	{
	return compute_geom_fex(
			off_fex_par,
			frame,
			frame_bin,
			NULL,
			left,
			0,
			lower_baseline,
			upper_baseline,
			background,
			3);
}

pen_fv_t compute_geom_fex_4(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame,
		im_std(plane,_t) frame_bin,
		im_std(plane,_t) frame_thin,
		int left,
		int lower_baseline,
		int upper_baseline,
		int background)
	{
	return compute_geom_fex(
			off_fex_par,
			frame,
			frame_bin,
			frame_thin,
			left,
			0,
			lower_baseline,
			upper_baseline,
			background,
			4);
}

pen_fv_t compute_geom_fex_5(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame_bin,
		im_std(plane,_t) frame_thin,
		int left,
		int lower_baseline,
		int upper_baseline,
		int background)
	{
	return compute_geom_fex(
			off_fex_par,
			NULL,
			frame_bin,
			frame_thin,
			left,
			0,
			lower_baseline,
			upper_baseline,
			background,
			5);
}

/* geometric feature extraction number 6 - already translated code */
pen_fv_t compute_geom_fex_6(
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) frame_bin,
		im_std(plane,_t) frame_thin,
		int left,
		int lower_baseline,
		int upper_baseline,
		im_std(grey,_t) background)
	{
	return compute_geom_fex(
			off_fex_par,
			NULL,
			frame_bin,
			frame_thin,
			left,
			0,
			lower_baseline,
			upper_baseline,
			background,
			6);
}


int move_along_line(
		pen_fvlist_t *result,
		im_std(planelist,_t) *frames,
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) plane,
		im_std(plane,_t) bin,
		im_std(plane,_t) thin,
		im_std(grey,_t) background)
	{

	im_std(plane,_t) frame = NULL;
	im_std(plane,_t) frame_bin = NULL;
	im_std(plane,_t) frame_thin = NULL;
	im_std(plane,_t) frame_norm = NULL;

	pen_fvbuffer_t fvbuffer = NULL;
	pen_fv_t fv = NULL;

	int i, i1, j, j1;
	int left = 0;
	int right = off_fex_par->frame_width;
	int baseline = 0;
	int coresize = 0;
	int lower_baseline = 0;
	int upper_baseline = 0;
	int dim = 0, dim_without_diffs = 0;

	char name[80];
	char index[10];

	int baseline_dest_idx, crop_offset;
	int src_y0, dest_y0, src_y_size;
	int grey_max;
	int jy;
	mx_real_t y, dy;
	mx_real_t scale;
	mx_real_t up_scale, core_scale, under_scale;

	int dimx, dimy;
	int frame_dimx, frame_dimy;
	int frame_norm_dimx, frame_norm_dimy;

	dimx = im_std(plane,_dimx)(plane);
	dimy = im_std(plane,_dimy)(plane);

	fvbuffer = pen_fvbuffer_create(3);

	if (off_fex_par->spline) {
		baseline = compute_baselines_SPLINE(off_fex_par, bin, plane,
				&lower_baseline, &upper_baseline);
	} else {
	  if (off_fex_par->horizontal_density_histogram) {
      baseline = compute_baselines_histogram(bin, &lower_baseline, &upper_baseline);
	  } else {
	    baseline = compute_baselines_contour(bin, &lower_baseline, &upper_baseline);
	  }
	}
	coresize = upper_baseline - lower_baseline;

#ifdef DEBUG
	if (upper_baseline == lower_baseline)
	fprintf(stderr," ACHTUNG: BASELINES IDENTISCH !!!!!\n");
#endif

	pen_msg("baselines: u= %3d(%3.1f), l= %3d, s= %3d(%3.1f), core = %3d",
			upper_baseline, (mx_real_t)upper_baseline / (mx_real_t)(-coresize), lower_baseline,
			dimy, (mx_real_t)dimy / (mx_real_t)(-coresize), -coresize);

#ifdef DEBUG_OUTPUT
	debug_output(image, output_filename, baseline, upper_baseline, bin);
#endif

	while (right < dimx) {
		frame = im_std(plane,_create)(off_fex_par->frame_width, dimy);
		frame_bin = im_std(plane,_create)(off_fex_par->frame_width, dimy);
		frame_thin = im_std(plane,_create)(off_fex_par->frame_width, dimy);

		frame_dimx = im_std(plane,_dimx)(frame);
		frame_dimy = im_std(plane,_dimy)(frame);

		/*
		 * used for normalized frame sizes
		 */
		/* method -1: */
#define NORM_UPPER	100
#define NORM_LOWER	150	/* i.e. NORM_CORE = 50 */
#define NORM_CORE	(NORM_LOWER - NORM_UPPER)
#define NORM_HEIGHT	200
#define NORM_UNDER	(NORM_HEIGHT - NORM_LOWER)
		/* method -2: */
#define NORM_FRM_SIZE	64
#define NORM_FRM_ISIZE	128
#define NORM_FRM_REL_BASELINE	0.75

		/*** settings for former n1.2* frame normalizations ******
		 #define CROP_FRM_SIZE	256
		 #define CROP_FRM_REL_BASELINE	0.67
		 *********************************************/
#define CROP_FRM_SIZE	150
#define CROP_FRM_REL_BASELINE	0.67

		if (off_fex_par->fex_method == -1) {
			frame_norm = im_std(plane,_create)(NORM_HEIGHT,	off_fex_par->frame_width);
		} else if (off_fex_par->fex_method == -2) {
			frame_norm = im_std(plane,_create)(NORM_FRM_ISIZE, off_fex_par->frame_width);
		} else if (off_fex_par->fex_method == -3) {
			frame_norm = im_std(plane,_create)(CROP_FRM_SIZE, off_fex_par->frame_width);
		}

		for (i = left, i1 = 0; i1 < frame_dimx; i++, i1++) {
			for (j = 0; j < frame_dimy; j++) {
				frame[j][i1] = plane[j][i];
				frame_bin[j][i1] = bin[j][i];
				frame_thin[j][i1] = thin[j][i];
			}
		}

		switch (off_fex_par->fex_method) {
			case -3:
				/* do noting, frame has already been created */
				grey_max = 0; /* max grey value per frame 4 filling */

				/* determine position of baseline and cropping region ... */
				baseline_dest_idx = CROP_FRM_SIZE * CROP_FRM_REL_BASELINE;
				crop_offset = lower_baseline - baseline_dest_idx;

				if (crop_offset < 0) {
					src_y0 = 0;
					dest_y0 = -crop_offset;
				} else {
					src_y0 = crop_offset;
					dest_y0 = 0;
				}

				src_y_size = frame_dimy - src_y0;
				if (dest_y0 + src_y_size > CROP_FRM_SIZE) {
					src_y_size = CROP_FRM_SIZE - dest_y0;
				}

				/* determine brightest grey value (background intensity) ... */
				for (j = 0; j < frame_dimy; j++) {
					for (i = 0; i < frame_dimx; i++) {
						if (frame[j][i]> grey_max)
							grey_max = frame[j][i];
					}
				}

				/* fill normalized frame ... */
				/* ... area "above" cropping offset is background ... */
				for (j = 0; j < dest_y0; j++) {
					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j][i] = grey_max;
					}
				}

				/* ... fill in cropped data from original frame ... */
				for (j = 0; j < src_y_size; j++) {
					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j + dest_y0][i] = frame[j + src_y0][i];
					}
				}

				/* ... area "below" cropped frame is background, too */
				for (j = dest_y0 + src_y_size; j < CROP_FRM_SIZE; j++) {
					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j][i] = grey_max;
					}
				}
				break;

			case -2:
				/* do noting, frame has already been created */
				grey_max = 0; /* max grey value per frame 4 filling */

				/* normalize relative position of lower baseline */
				scale = (mx_real_t)lower_baseline
						/(mx_real_t)(NORM_FRM_REL_BASELINE * NORM_FRM_ISIZE);

				for (j = 0, y = 0.0; j < NORM_FRM_ISIZE; j++, y += scale) {
					jy = __trunc(y);
					dy = y - jy;

					if (jy == (frame_dimy - 1)) {
						if (dy < 0.2) {
							dy = 0.0;
						} else {
							jy += 1;
						}
					}

					for (i = 0; i < frame_dimx; i++) {
						if (jy < frame_dimy) {
							frame_norm[j][i] = interpolate_bilin(frame, 0.0,	dy, i, jy);

							if (frame_norm[j][i] > grey_max) {
								grey_max = frame_norm[j][i];
							}
						} else {
							frame_norm[j][i] = grey_max;
						}
					}
				}
				break;

			case -1:
				/* do noting, frame has already been created */

				/* normalize area above upper baseline */
				up_scale = (mx_real_t)upper_baseline / NORM_UPPER;

				for (j = 0, y = 0.0; j < NORM_UPPER; j++, y += up_scale) {
					jy = __trunc(y);
					dy = y - jy;

					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j][i] = interpolate_bilin(frame, 0.0,	dy, i, jy);
					}
				}

				/* normalize core area */
				core_scale = (mx_real_t)(-coresize) / NORM_CORE;

				for (j = NORM_UPPER, y = (mx_real_t)upper_baseline; j < NORM_LOWER; j++, y
						+= core_scale) {
					jy = __trunc(y);
					dy = y - jy;

					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j][i] = interpolate_bilin(frame, 0.0,	dy, i, jy);
					}
				}

				/* normalize "Unterlaengen" area */
				under_scale = (mx_real_t)(frame_dimy - lower_baseline) / NORM_UNDER;

				for (j = NORM_LOWER, y = (mx_real_t)lower_baseline; j < NORM_HEIGHT; j++, y
						+= under_scale) {
					jy = __trunc(y);
					dy = y - jy;

					for (i = 0; i < frame_dimx; i++) {
						frame_norm[j][i] = interpolate_bilin(frame, 0.0,	dy, i, jy);
					}
				}
				break;

			case 0:
				fv = compute_bitmap_fex(off_fex_par, frame);
				break;

			case 1:
				fv = compute_geom_fex_1(off_fex_par, frame, frame_bin, baseline,
						background);
				break;

			case 2:
				fv = compute_geom_fex_2(off_fex_par, frame, frame_bin, baseline,
						background);
				break;

			case 7:
				fv = compute_geom_fex_2(off_fex_par, frame, frame_bin, baseline,
						background);
				break;

			case 3:
				fv = compute_geom_fex_3(off_fex_par, frame, frame_bin, left, baseline,
						upper_baseline, background);
				break;

			case 4:
				fv = compute_geom_fex_4(off_fex_par, frame, frame_bin, frame_thin, left,
						baseline, upper_baseline, background);
				break;

			case 5:
				fv = compute_geom_fex_5(off_fex_par, frame_bin, frame_thin, left,
						baseline, upper_baseline, background);
				break;

			case 6:
				fv = compute_geom_fex_6(off_fex_par, frame_bin, frame_thin, left,
						baseline, upper_baseline, background);
				break;

			default:
				rs_error("illegal feature extraction method %d!", off_fex_par->fex_method);
		}

		/* if real features are to be extracted ...*/
		if (off_fex_par->fex_method >= 0) {
			pen_fvbuffer_addfv(fvbuffer, fv);

			if (off_fex_par->use_derivatives) {
				dim_without_diffs = fv->dim / 2;
				dim = fv->dim;

				if (fvbuffer->current == 2) {
					/* now we can calculate horizontal deviations */
					for (i = 0; i < dim_without_diffs; i++) {
						fvbuffer->fv[1]->feature[dim_without_diffs + i] = (fvbuffer->fv[2]->feature[i]
								- fvbuffer->fv[0]->feature[i]) / 2.0;
					}
					pen_fvlist_addfv(*result, pen_fv_dup(fvbuffer->fv[1]));
				}
			} else {
				pen_fvlist_addfv(*result, pen_fv_dup(fvbuffer->fv[fvbuffer->current]));
			}
		} else {
			//pen_fv_destroy(fv);
			//fv = NULL;
		}

		left = right - off_fex_par->frame_overlap;
		right += (off_fex_par->frame_width - off_fex_par->frame_overlap);

		if (frames) {
			im_std(planelist,_addplane)(*frames, frame);
		} else {
			im_std(plane,_destroy)(frame);
		}

		frame = NULL;
		im_std(plane,_destroy)(frame_bin);
		frame_bin = NULL;
		im_std(plane,_destroy)(frame_thin);
		frame_thin = NULL;

		fv = NULL;
	}

	pen_fvbuffer_destroy(fvbuffer);
	fvbuffer = NULL;

	return PEN_SUCCESS;
	}

int im_std(plane,_fex)(
		pen_fvlist_t *result,
		im_std(planelist,_t) *frames,
		pen_fextractpar_t off_fex_par,
		im_std(plane,_t) plane)
	{
	im_std(plane,_t) binary_plane = NULL;
	im_std(plane,_t) skeleton_plane = NULL;
	im_std(grey,_t) background;

	int dimx, dimy;
	int status;
	int i, j, x , y;
	int sum;

	int lower_bound, upper_bound;

	/* ... check parameters ... */
	if (!plane) {
		return IM_ERROR_ARG;
	}

	/* ... determine plane dimensions ... */
	dimx = im_std(plane,_dimx)(plane);
	dimy = im_std(plane,_dimy)(plane);
	if (dimx <= 0 || dimy <= 0) {
		return IM_ERROR_DATA;
	}

	/* ... binarize ... */
	status = im_std(plane,_binarize_otsu)(&binary_plane, plane, NULL, &background);

	/* ... skeletonize ... */
	status = im_std(plane,_skeletonizeC)(&skeleton_plane, binary_plane);

	/* ... extract features ... */
	move_along_line(
			result,
			frames,
			off_fex_par,
			plane,
			binary_plane,
			skeleton_plane,
			background);

	/* free stuff */
	im_std(plane,_destroy)(binary_plane);
	im_std(plane,_destroy)(skeleton_plane);

	return status;
	}

int im_std(image,_fex)(
		pen_fvlist_t *result,
		im_std(planelist,_t) *frames,
		pen_fextractpar_t off_fex_par,
		im_std(image,_t) image)
	{

	im_std(image,_t) intensity_image = NULL;
	im_std(plane,_t) plane = NULL;

	int dimx, dimy;
	int mask_dimx, mask_dimy;
	int n_planes;
	int status;
	int i, j;

	/* ... check parameters ... */
	if (!image) {
		return IM_ERROR_ARG;
	}

	/* ... determine image dimensions ... */
	dimx = im_std(image,_dimx)(image);
	dimy = im_std(image,_dimy)(image);
	n_planes = im_std(image,_n_planes)(image);
	if (dimx <= 0 || dimy <= 0 || n_planes <= 0) {
		return IM_ERROR_DATA;
	}

	/* ... if more than 1 planes given or input is unsigned ... */
#ifdef IM_GREY_SIGNED
	if (n_planes > 1)
#endif /* IM_GREY_SIGNED */
	{
		/* ... calculate average grey image ... */
		status = im_std(image,_rgb2i)(&intensity_image, image);

		if (status != IM_SUCCESS) {
			return status;
		}

		/* ... create new plane ... */
		plane = intensity_image[0];
	}
#ifdef IM_GREY_SIGNED
	else {
		plane = image[0];
	}
#endif /* IM_GREY_SIGNED */

	if (*result) {
		pen_fvlist_destroy(*result);
		*result = NULL;
	}

	*result = pen_fvlist_create();

	status = im_std(plane,_fex)(
			result,
			frames,
			off_fex_par,
			plane);

	/* free stuff */
	im_std(image,_destroy)(intensity_image);

	return status;
	}
