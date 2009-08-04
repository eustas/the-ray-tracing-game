double doubmax; // approx square root of max double number
double doubmin; // smallest double number
double doubtol; // tolerance of double numbers

// setup constants
void setcns() {
	doubtol = 1.0;
	while (1.0 + doubtol > 1.0) {doubtol *= 0.5;}
	doubtol = sqrt(doubtol);

	doubmin = 0.5;
	for (int j = 1; j <= 100; ++j) {
		doubmin = doubmin * doubmin;
		if ((doubmin * doubmin) <= (doubmin * doubmin * 0.5)) {break;}
	}
	doubmax = 0.7 / sqrt(doubmin);
}

// solve quartic equation using either quadratic, Ferrari's or Neumark's algorithm
int quartic(double a, double b, double c, double d, double *rts) {
	int j, k, nq, nr;
	double roots[4];

	double odd = (a < 0.0) ? -a : a;
	odd += (c < 0.0) ? -c : c;
	double even = (b < 0.0) ? -b : b;
	even += (d < 0.0) ? -d : d;

	if (odd < even * doubtol) {
		nq = qudrtc(b, d, roots, b * b - 4.0 * d);
		j = 0;
		for (k = 0; k < nq; ++k) {
			if (roots[k] > 0.0) {
				rts[j] = sqrt(roots[k]);
				rts[j + 1] = -rts[j];
				++j;
				++j;
			}
		}
		nr = j;
		return nr;
	}
	k = (a < 0.0) ? 1 : 0;
	k += (b < 0.0) ? k + 1 : k;
	k += (c < 0.0) ? k + 1 : k;
	k += (d < 0.0) ? k + 1 : k;
	switch (k) {
		case 0:nr = ferrari(a, b, c, d, rts);break;
		case 1:nr = neumark(a, b, c, d, rts);break;
		case 2:nr = neumark(a, b, c, d, rts);break;
		case 3:nr = ferrari(a, b, c, d, rts);break;
		case 4:nr = ferrari(a, b, c, d, rts);break;
		case 5:nr = neumark(a, b, c, d, rts);break;
		case 6:nr = ferrari(a, b, c, d, rts);break;
		case 7:nr = ferrari(a, b, c, d, rts);break;
		case 8:nr = neumark(a, b, c, d, rts);break;
		case 9:nr = ferrari(a, b, c, d, rts);break;
		case 10:nr = ferrari(a, b, c, d, rts);break;
		case 11:nr = neumark(a, b, c, d, rts);break;
		case 12:nr = ferrari(a, b, c, d, rts);break;
		case 13:nr = ferrari(a, b, c, d, rts);break;
		case 14:nr = ferrari(a, b, c, d, rts);break;
		case 15:nr = ferrari(a, b, c, d, rts);break;
	}
	return nr;
}

// solve the quartic equation; method: Ferrari-Lagrange
int ferrari(double a, double b, double c, double d, double *rts) {
	double v1[4], v2[4], e, f, h, hh;

	double asq = a * a;
	double p = b;
	double q = a * c - 4.0 * d;
	double r = (asq - 4.0 * b) * d + c * c;
	double y = cubic(p, q, r);

	double esq = 0.25 * asq - b - y;
	if (esq < 0.0) {return (0);}

	double fsq = 0.25 * y * y - d;
	if (fsq < 0.0) {return (0);}
	double ef = -(0.25 * a * y + 0.5 * c);
	if (((a > 0.0) && (y > 0.0) && (c > 0.0)) || ((a > 0.0) && (y < 0.0) && (c
			< 0.0)) || ((a < 0.0) && (y > 0.0) && (c < 0.0)) || ((a < 0.0)
			&& (y < 0.0) && (c > 0.0)) || (a == 0.0) || (y == 0.0)
			|| (c == 0.0)) {
		if ((b < 0.0) && (y < 0.0) && (esq > 0.0)) {
			e = sqrt(esq);
			f = ef / e;
		} else if ((d < 0.0) && (fsq > 0.0)) {
			f = sqrt(fsq);
			e = ef / f;
		} else {
			e = sqrt(esq);
			f = sqrt(fsq);
			if (ef < 0.0) {
				f = -f;
			}
		}
	} else {
		e = sqrt(esq);
		f = sqrt(fsq);
		if (ef < 0.0) {
			f = -f;
		}
	}
	double ainv2 = a * 0.5;
	double g = ainv2 - e;
	double gg = ainv2 + e;
	if (((b > 0.0) && (y > 0.0)) || ((b < 0.0) && (y < 0.0))) {
		if ((a > 0.0) && (e != 0.0)) {
			g = (b + y) / gg;
		} else if (e != 0.0) {
			gg = (b + y) / g;
		}
	}
	if ((y == 0.0) && (f == 0.0)) {
		h = 0.0;
		hh = 0.0;
	} else if (((f > 0.0) && (y < 0.0)) || ((f < 0.0) && (y > 0.0))) {
		hh = -0.5 * y + f;
		h = d / hh;
	} else {
		h = -0.5 * y - f;
		hh = d / h;
	}
	int n1 = qudrtc(gg, hh, v1, gg * gg - 4.0 * hh);
	int n2 = qudrtc(g, h, v2, g * g - 4.0 * h);
	int nquar = n1 + n2;
	rts[0] = v1[0];
	rts[1] = v1[1];
	rts[n1 + 0] = v2[0];
	rts[n1 + 1] = v2[1];
	return nquar;
}

// solve the quartic equation; method: Neumark
int neumark(double a, double b, double c, double d, double *rts) {
	double gdisrt, hdisrt, g2, h2;
	double v1[4], v2[4];

	double asq = a * a;

	double p = -b * 2.0;
	double q = b * b + a * c - 4.0 * d;
	double r = (c - a * b) * c + asq * d;
	double y = cubic(p, q, r);

	double bmy = b - y;
	double y4 = y * 4.0;
	double d4 = d * 4.0;
	double bmysq = bmy * bmy;
	double gdis = asq - y4;
	double hdis = bmysq - d4;
	if ((gdis < 0.0) || (hdis < 0.0)) {
		return (0);
	}

	double g1 = a * 0.5;
	double h1 = bmy * 0.5;
	double gerr = asq + y4;
	double herr = hdis;
	if (d > 0.0) {
		herr = bmysq + d4;
	}

	if ((y < 0.0) || (herr * gdis > gerr * hdis)) {
		gdisrt = sqrt(gdis);
		g2 = gdisrt * 0.5;
		if (gdisrt != 0.0) {
			h2 = (a * h1 - c) / gdisrt;
		} else {
			h2 = 0.0;
		}
	} else {
		hdisrt = sqrt(hdis);
		h2 = hdisrt * 0.5;
		if (hdisrt != 0.0) {
			g2 = (a * h1 - c) / hdisrt;
		} else {
			g2 = 0.0;
		}
	}
	double h = h1 - h2;
	double hh = h1 + h2;
	double hmax = hh;
	if (hmax < 0.0) {
		hmax = -hmax;
	}
	if (hmax < h) {
		hmax = h;
	}
	if (hmax < -h) {
		hmax = -h;
	}
	if ((h1 > 0.0) && (h2 > 0.0)) {
		h = d / hh;
	}
	if ((h1 < 0.0) && (h2 < 0.0)) {
		h = d / hh;
	}
	if ((h1 > 0.0) && (h2 < 0.0)) {
		hh = d / h;
	}
	if ((h1 < 0.0) && (h2 > 0.0)) {
		hh = d / h;
	}
	if (h > hmax) {
		h = hmax;
	}
	if (h < -hmax) {
		h = -hmax;
	}
	if (hh > hmax) {
		hh = hmax;
	}
	if (hh < -hmax) {
		hh = -hmax;
	}

	double g = g1 - g2;
	double gg = g1 + g2;
	double gmax = gg;
	if (gmax < 0.0) {
		gmax = -gmax;
	}
	if (gmax < g) {
		gmax = g;
	}
	if (gmax < -g) {
		gmax = -g;
	}
	if ((g1 > 0.0) && (g2 > 0.0)) {
		g = y / gg;
	}
	if ((g1 < 0.0) && (g2 < 0.0)) {
		g = y / gg;
	}
	if ((g1 > 0.0) && (g2 < 0.0)) {
		gg = y / g;
	}
	if ((g1 < 0.0) && (g2 > 0.0)) {
		gg = y / g;
	}
	if (g > gmax) {
		g = gmax;
	}
	if (g < -gmax) {
		g = -gmax;
	}
	if (gg > gmax) {
		gg = gmax;
	}
	if (gg < -gmax) {
		gg = -gmax;
	}

	int n1 = qudrtc(gg, hh, v1, gg * gg - 4.0 * hh);
	int n2 = qudrtc(g, h, v2, g * g - 4.0 * h);
	int nquar = n1 + n2;
	rts[0] = v1[0];
	rts[1] = v1[1];
	rts[n1 + 0] = v2[0];
	rts[n1 + 1] = v2[1];
	return (nquar);
}

// solve the quadratic equation: x**2+b*x+c = 0
int qudrtc(double b, double c, double *rts, double dis) {
	if (dis >= 0.0) {
		double rtdis = sqrt(dis);
		if (b > 0.0) {
			rts[0] = (-b - rtdis) * 0.5;
		} else {
			rts[0] = (-b + rtdis) * 0.5;
		}
		if (rts[0] == 0.0) {
			rts[1] = -b;
		} else {
			rts[1] = c / rts[0];
		}
		return 2;
	}

	rts[0] = 0.0;
	rts[1] = 0.0;
	return 0;
}

// find the lowest real root of the cubic - x**3 + p*x**2 + q*x + r = 0
double cubic(double p, double q, double r) {
	double po3, po3sq, qo3;
	double uo3, u2o3, uo3sq4, uo3cu4;
	double v, vsq, wsq;
	double mcube, n;
	double muo3, s, scube, t, cosk, sinsqk;
	double root;

	double m = 0.0;
	int nrts = 0;
	if ((p > doubmax) || (p < -doubmax)) {
		return -p;
	}
	if ((q > doubmax) || (q < -doubmax)) {
		if (q > 0.0) {
			return -r / q;
		}
		return -sqrt(-q);
	}
	if ((r > doubmax) || (r < -doubmax)) {
		return -curoot(r);
	}
	po3 = p * 0.3333333333333333333333333333333333333333333333;
	po3sq = po3 * po3;
	if (po3sq > doubmax) {
		return -p;
	}
	v = r + po3 * (po3sq + po3sq - q);
	if ((v > doubmax) || (v < -doubmax)) {
		return -p;
	}
	vsq = v * v;
	qo3 = q * 0.3333333333333333333333333333333333333333333333;
	uo3 = qo3 - po3sq;
	u2o3 = uo3 + uo3;
	if ((u2o3 > doubmax) || (u2o3 < -doubmax)) {
		if (p == 0.0) {
			if (q > 0.0) {
				root = -r / q;
			} else {
				root = -sqrt(-q);
			}
		} else {
			root = -q / p;
		}
	}
	uo3sq4 = u2o3 * u2o3;
	if (uo3sq4 > doubmax) {
		if (p == 0.0) {
		
			if (q > 0.0) {
				root = -r / q;
			} else {
				root = -sqrt(fabs(q));
			}
		} else {
			root = -q / p;
		}
	}
	uo3cu4 = uo3sq4 * uo3;
	wsq = uo3cu4 + vsq;
	if (wsq >= 0.0) {
		nrts = 1;
		if (v <= 0.0) {mcube = (-v + sqrt(wsq)) * 0.5;}
		if (v > 0.0) {mcube = (-v - sqrt(wsq)) * 0.5;}
		m = curoot(mcube);
		if (m != 0.0) {n = -uo3 / m;}
		else {n = 0.0;}
		return m + n - po3;
	}
	nrts = 3;
	if (uo3 >= 0.0) {return curoot(v) - po3;}
	muo3 = -uo3;
	s = sqrt(muo3);
	scube = s * muo3;
	t = -v / (scube + scube);
	cosk = cos(acos(t) * 0.3333333333333333333333333333333333333333333333);
	if (po3 < 0.0) {return (s + s) * cosk - po3;}
	sinsqk = 1.0 - cosk * cosk;
	if (sinsqk < 0.0) {sinsqk = 0.0;}
	return s * (-cosk - 1.7320508075688772935274463415059 * sqrt(sinsqk)) - po3;
}

// find cube root of x
double curoot(double x) {
	return (x < 0.0) ? -exp(log(-x)
			* 0.3333333333333333333333333333333333333333333333) : exp(log(x)
			* 0.3333333333333333333333333333333333333333333333);
}