class CLASS_NAME {
	LPDWORD* const my_lines;
	number const my_x0;
	number const my_y0;
	number const my_z0;

int ferrari(const number a, const number b, const number c, const number d, number *rts) const {
	number v1[4], v2[4], e, f, h, hh;

	number asq = a * a;
	number p = b;
	number q = a * c - 4.0 * d;
	number r = (asq - 4.0 * b) * d + c * c;
	number y = cubic(p, q, r);

	number esq = 0.25 * asq - b - y;
	if (esq < 0.0) {return 0;}

	number fsq = 0.25 * y * y - d;
	if (fsq < 0.0) {return 0;}
	number ef = -(0.25 * a * y + 0.5 * c);
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
			if (ef < 0.0) {f = -f;}
		}
	} else {
		e = sqrt(esq);
		f = sqrt(fsq);
		if (ef < 0.0) {
			f = -f;
		}
	}
	number ainv2 = a * 0.5;
	number g = ainv2 - e;
	number gg = ainv2 + e;
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

int qudrtc(const number b, const number c, number *rts, const number dis) const {
	if (dis >= 0.0) {
		number rtdis = sqrt(dis);
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
number cubic(const number p, const number q, const number r) const {
	number po3, po3sq, qo3;
	number uo3, u2o3, uo3sq4, uo3cu4;
	number v, vsq, wsq;
	number mcube, n;
	number muo3, s, scube, t, cosk, sinsqk;
	number root;

	number m = 0.0;
	int nrts = 0;
	if ((p > doubmax) || (p < -doubmax)) {return -p;}
	if ((q > doubmax) || (q < -doubmax)) {
		if (q > 0.0) {return -r / q;}
		return -sqrt(-q);
	}
	if ((r > doubmax) || (r < -doubmax)) {
		return ((r < 0.0) ? exp(log(-r)
			* 0.3333333333333333333333333333333333333333333333) : -exp(log(r)
			* 0.3333333333333333333333333333333333333333333333));
	}
	po3 = p * 0.3333333333333333333333333333333333333333333333;
	po3sq = po3 * po3;
	if (po3sq > doubmax) {return -p;}
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
		m = (mcube < 0.0) ? -exp(log(-mcube)
			* 0.3333333333333333333333333333333333333333333333) : exp(log(mcube)
			* 0.3333333333333333333333333333333333333333333333);
		if (m != 0.0) {n = -uo3 / m;}
		else {n = 0.0;}
		return m + n - po3;
	}
	nrts = 3;
	if (uo3 >= 0.0) {return ((v < 0.0) ? -exp(log(-v)
			* 0.3333333333333333333333333333333333333333333333) : exp(log(v)
			* 0.3333333333333333333333333333333333333333333333)) - po3;}
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


	bool Solve(const number a, const number b,const number c,const number d,const number tStart,const number tEnd, number &result) const {
		result = LAST_T;

		number rts[4];
		int sol = ferrari(a, b, c, d, rts);
		for (int i = 0; i < sol; i++) {
			double s = rts[i];
			if ((s < result) && (s > tStart) && (s < tEnd)) {result = s;}
		}
		return result < LAST_T;
	}

	bool CanSee(const int light, number x0, number y0, number z0, bool &lensed, number &factor, number &lx, number &ly, number &lz) const {
		factor = 1.0;
		lensed = false;
		number vx = lights[light].x;
		number vy = lights[light].y;
		number vz = lights[light].z;
		number wx,wy,wz;
		number sx,sy,sz, coneC,econeA,einvConeA,a,c,cox,coz;
restart:
		number coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
		number invConeA = 1 / coneA;
		bool cAp = coneA > 0.0;
		number dx, dy, dz, w, k, D, t0,t1, t2, t3, t4, tmp,temp, tShi;
		number x,y,z, b,f,e,x1,y1,z1;
		int mx = (int)x0;
		int mz = (int)z0;
		int d;
		int deepness = 0;

		number tNextX, tmx;
		int dmx;
		if (vx > 0.0) {
			tmx = 1.0 / vx;
			tNextX = (mx + 1.0 - x0) / vx;
			dmx = 1;
		} else if (vx < 0.0) {
			tmx = -1.0 / vx;
			tNextX = (mx - x0) / vx;
			dmx = -1;
		} else {
			tNextX = LAST_T;
		}

		number tNextZ, tmz;
		int dmz;
		if (vz > 0.0) {
			tmz = 1.0 / vz;
			tNextZ = (mz + 1.0 - z0) / vz;
			dmz = 1;
		} else if (vz < 0.0) {
			tmz = -1.0 / vz;
			tNextZ = (mz - z0) / vz;
			dmz = -1;
		} else {
			tNextZ = LAST_T;
		}

		number tStart, tEnd;
		bool boundsFloor;
		CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);
		tStart+=0.00000001;
		number rayFirst = tStart;

		int idx = mx + 15 * mz;
		number nx, ny, nz;
		while (true) {
			if ((mx < 0) || (mx >= 15) || (mz < 0) || (mz >= 9)) {break;}
			char type = rt::cells[idx].type;
			char state = rt::cells[idx].state;
			switch (type) {
				case POLAR:
				case MPOLAR:
					tmp = (vx * poX[state]) + (vz * poZ[state]);
					if (tmp < 0.0) {tmp = -tmp;}
					if (tmp > 0.95 * sqrt((vx * vx) + (vz * vz))) {break;}
					return false;


				case RANDOM:
					if (((int)(((vx * vy) + (vy * vz) + (vz * vx) + x0 + y0 + z0) * (65536 + subFrame)) & 0xFF) > 0x80) {return false;}
					break;

				case COLLECTOR:
tShi = ((mx + 0.5 - x0) * vx) + ((0.5 - y0) * vy) + ((mz + 0.5 - z0) * vz);
x1 = x0 + (vx * tShi); y1 = y0 + (vy * tShi); z1 = z0 + (vz * tShi);
					if (state == 0) {

dx = x1 - mx - 0.03;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx) + 0.0216; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll1;}
return false;
  
coll1:
dx = x1 - mx - 0.3;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (0.04 * dx * dx) + 0.0084; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
temp = e + (vx * vx * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll2;}
return false;

coll2:
dx = x1 - mx - 0.6;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx); // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll3;}
return false;

coll3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.9 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto coll4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto coll4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox > (mx + 1.0)) { goto coll4;}
return false;
coll4:
if (vx == 0.0) {break;}
t1 = (mx + 1.0 - x0) / vx + 0.000001;
if (t1 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					} else if (state == 1) {
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.03;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz) + 0.0216; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}  
colu1:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.3;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (0.04 * dz * dz) + 0.0084; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
temp = e + (vz * vz * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}
colu2:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.60;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz); // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}
colu3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.9 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto colu4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto colu4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz > (mz + 1.0)) { goto colu4;}
return false;
colu4:
if (vz == 0.0) {break;}
t1 = (mz + 1.0 - z0) / vz + 0.000001;
if (t1 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {break;}
return false;
					} else if (state == 2) {
dx = x1 - mx - 0.97;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx) + 0.0216; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr1;}
return false;  
colr1:
dx = x1 - mx - 0.7;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (0.04 * dx * dx) + 0.0084; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
temp = e + (vx * vx * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr2;}
return false;
colr2:
dx = x1 - mx - 0.4;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx); // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr3;}
return false;
colr3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.1 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto colr4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto colr4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox < mx) { goto colr4;}
return false;
colr4:
if (vx == 0.0) {break;}
t1 = (mx - x0) / vx + 0.000001;
if (t1 + 0.00002 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {break;}
return false;
					} else if (state == 3) {
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.97;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz) + 0.0216; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}  
cold1:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.7;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (0.04 * dz * dz) + 0.0084; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
temp = e + (vz * vz * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}
cold2:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.40;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz); // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {return false;}
cold3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.1 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto cold4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto cold4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz < mz) { goto cold4;}
return false;
cold4:
if (vz == 0.0) {break;}
t1 = (mz - z0) / vz + 0.000001;
if (t1 + 0.00002 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					}
					break;

				case EMITTER:
					if (state == 0) {
//#include"emi_l_c.cpp"
sx = mx - 0.1 - x0;
sy = 0.5 - y0;
sz = (0.5 + mz) - z0;
k = (sy * vy) + (sz * vz) - (sx * vx * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sz * sz) - (sx * sx * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vz * vz) - (vx * vx * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emil1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
x = (mx + 1.0 - 0.1) - (x0 + (vx * t1));
if (x > E_CONE_H - 0.1) {
	x = (mx + 1.0 - 0.1) - (x0 + (vx * t2));
	if (x > E_CONE_H - 0.1) { goto emil1; }
	t1 = t2;
}
if (x < 0.1) { goto emil1; }
if (t1 < rayFirst) {goto emil1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * x); // tmp == r
if (tmp <= 0.0) {goto emil1;}
return false;
emil1:
if (vx == 0.0) {goto emil2;}
t1 = (mx - x0) / vx;
if (t1 < rayFirst) {goto emil2;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > RAY_R_2) {goto emil2;}
return false;
emil2:
a = (vx * vx * 16.0) + (vy * vy) + (vz * vz);
x = mx + 0.25 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (16.0 * x * vx) + (y * vy) + (z * vz);
c = (16.0 * x * x) + (y * y) + (z * z) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emil3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emil3;}
return false;
emil3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.9 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emil4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emil4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox > (mx + 1.0)) { goto emil4;}
return false;
emil4:
if (vx == 0.0) {break;}
t1 = (mx + 1.0 - x0) / vx;
if (t1 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					} else if (state == 1) {
//#include<emi_u_p.cpp>
sz = mz - 0.1 - z0;
sy = 0.5 - y0;
sx = (0.5 + mx) - x0;
k = (sy * vy) + (sx * vx) - (sz * vz * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sx * sx) - (sz * sz * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vx * vx) - (vz * vz * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emiu1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
z = (mz + 1.0 - 0.1) - (z0 + (vz * t1));
if (z > E_CONE_H - 0.1) {
	z = (mz + 1.0 - 0.1) - (z0 + (vz * t2));
	if (z > E_CONE_H - 0.1) { goto emiu1; }
	t1 = t2;
}
if (z < 0.1) { goto emiu1; }
if (t1 < rayFirst) {goto emiu1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * z); // tmp == r
if (tmp <= 0.0) {goto emiu1;}
return false;
emiu1:
if (vz == 0.0) {goto emiu2;}
t1 = (mz - z0) / vz;
if (t1 < rayFirst) {goto emiu2;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > RAY_R_2) {goto emiu2;}
return false;
emiu2:
a = (vz * vz * 16.0) + (vy * vy) + (vx * vx);
z = mz + 0.25 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (16.0 * z * vz) + (y * vy) + (x * vx);
c = (16.0 * z * z) + (y * y) + (x * x) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emiu3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emiu3;}
return false;
emiu3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.9 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emiu4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emiu4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz > (mz + 1.0)) { goto emiu4;}
return false;
emiu4:
if (vz == 0.0) {break;}
t1 = (mz + 1.0 - z0) / vz;
if (t1 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					} else if (state == 2) {
//#include<emi_r_p.cpp>
sx = mx + 1.1 - x0;
sy = 0.5 - y0;
sz = (0.5 + mz) - z0;
k = (sy * vy) + (sz * vz) - (sx * vx * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sz * sz) - (sx * sx * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vz * vz) - (vx * vx * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emir1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
x = (mx + 0.1) - (x0 + (vx * t1));
if (x < -E_CONE_H + 0.1) {
	x = (mx + 0.1) - (x0 + (vx * t2));
	if (x < -E_CONE_H + 0.1) { goto emir1; }
	t1 = t2;
}
x = -x;
if (x < 0.1) { goto emir1; }
if (t1 < rayFirst) {goto emir1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * x); // tmp == r
if (tmp <= 0.0) {goto emir1;}
return false;
emir1:
if (vx == 0.0) {goto emir2;}
t1 = (mx + 1.0 - x0) / vx;
if (t1 < rayFirst) {goto emir2;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > RAY_R_2) {goto emir2;}
return false;
emir2:
a = (vx * vx * 16.0) + (vy * vy) + (vz * vz);
x = mx + 0.75 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (16.0 * x * vx) + (y * vy) + (z * vz);
c = (16.0 * x * x) + (y * y) + (z * z) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emir3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emir3;}
return false;
emir3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.1 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emir4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emir4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox < mx) { goto emir4;}
return false;
emir4:
if (vx == 0.0) {break;}
t1 = (mx - x0) / vx;
if (t1 + 0.00002 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					} else if (state == 3) {
//#include<emi_d_p.cpp>
sz = mz + 1.1 - z0;
sy = 0.5 - y0;
sx = (0.5 + mx) - x0;
k = (sy * vy) + (sx * vx) - (sz * vz * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sx * sx) - (sz * sz * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vx * vx) - (vz * vz * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emid1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
z = (mz + 0.1) - (z0 + (vz * t1));
if (z < -E_CONE_H + 0.1) {
	z = (mz + 0.1) - (z0 + (vz * t2));
	if (z < -E_CONE_H + 0.1) { goto emid1; }
	t1 = t2;
}
z = -z;
if (z < 0.1) { goto emid1; }
if (t1 < rayFirst) {goto emid1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * z); // tmp == r
if (tmp <= 0.0) {goto emid1;}
return false;
emid1:
if (vz == 0.0) {goto emid2;}
t1 = (mz + 1.0 - z0) / vz;
if (t1 < rayFirst) {goto emid2;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > RAY_R_2) {goto emid2;}
return false;
emid2:
a = (vz * vz * 16.0) + (vy * vy) + (vx * vx);
z = mz + 0.75 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (16.0 * z * vz) + (y * vy) + (x * vx);
c = (16.0 * z * z) + (y * y) + (x * x) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emid3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emid3;}
return false;
emid3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.1 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emid4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emid4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz < mz) { goto emid4;}
return false;
emid4:
if (vz == 0.0) {break;}
t1 = (mz - z0) / vz;
if (t1 + 0.00002 < rayFirst) {break;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {break;}
return false;

					}
					break;


				case PORTAL:
					temp = LAST_T;
					dx = mx + 0.5 - x0;
					dy =      0.5 - y0;
					dz = mz + 0.5 - z0;
					w = (dx * dx) + (dy * dy) + (dz * dz) - (PORTAL_R * PORTAL_R);
					k = (dx * vx) + (dy * vy) + (dz * vz);
					D = (k * k) - w;
					if (D <= 0.0) { goto portalTorusC; }
					t1 = k - sqrt(D);
					if (t1 < tStart) {
						t1 = k + sqrt(D);
						if (t1 > tStart) {return false;}
						goto portalTorusC;
					}
					temp = t1;
portalTorusC:
					dx = x0 - mx - 0.5;
					dy = y0 - 0.1;
					dz = z0 - mz - 0.5;
					a = 2.0 * ((vx * dx) + (vz * dz)); // r: t
					b = a + (2.0 * vy * dy); // l^1 : t
					c = (dx * dx) + (dz * dz); // r: 1
					f = c + (dy * dy) + 0.15; // l^1: 1
					e = (vx * vx) + (vz * vz); // r: t^2
					//t4 = 1.0;
					t3 = 2.0 * b;
					t2 = (2.0 * f) + (b * b) - (0.64 * e);
					t1 = (2.0 * b * f) - (0.64 * a);
					t0 = (f * f) - (0.64 * c);
					if (!Solve(/*t4,*/t3,t2,t1,t0, tStart, tEnd, tmp)) {
						if (temp<LAST_T) {goto portalPortal;}
						break;
					}
					if (temp<tmp) {goto portalPortal;}
					return false;
portalPortal:
					deepness++;
					if (deepness > 5) {return false;}
					t2 = t1 - 0.00001;
					x = x0 + (vx * t2) - mx;
					y = y0 + (vy * t2);
					z = z0 + (vz * t2) - mz;
					x = 1.0 - x;
					y = 1.0 - y;
					z = 1.0 - z;
					d = 1;
					if ((levelData.px[state * 2] != mx) || (levelData.py[state * 2] != mz)) {d = 0;}
					mx = levelData.px[state * 2 + d];
					mz = levelData.py[state * 2 + d];
					x0 = mx + x;
					y0 = y;
					z0 = mz + z;
					CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);

					tStart += 0.000000001;
					tEnd   += 0.000000001;

					if (vx > 0.0) {tNextX = (mx + 1.0 - x0) / vx;} else if (vx < 0.0) {tNextX = (mx - x0) / vx;}
					if (vz > 0.0) {tNextZ = (mz + 1.0 - z0) / vz;} else if (vz < 0.0) {tNextZ = (mz - z0) / vz;}
					idx = mx + 15 * mz;
					break;



				case BLOCK:
				case MBLOCK:
					t1 = tStart;
					if (IntersectBlock(x0, y0, z0, vx, vy, vz, state, mx, mz, nx, ny, nz, t1)) {
						return false;
					}
				break;

				case QUANT:
				case STAR:
					dx = mx + 0.5 - x0;
					dy = 0.5 - y0;
					dz = mz + 0.5 - z0;
					w = (dx * dx) + (dy * dy) + (dz * dz) - (QUANT_R * QUANT_R);
					k = (dx * vx) + (dy * vy) + (dz * vz);
					D = (k * k) - w;
					if (D <= 0.0) { break; }
					t1 = k - sqrt(D);
					if (t1 > 0.0) {
#ifndef GLASS_POINT
						return false;
#else
						if (type == STAR) {return false;}
						goto lens;
#endif
					}
					t1 = k + sqrt(D);
					if (t1 > 0.0) {
						if (type == STAR) {return false;}
						// strange case at all!
#ifndef GLASS_POINT
						return false;
#endif
					}
					break;
lens:
					wx = dx - (vx * FOCUS); wy = dy - (vy * FOCUS); wz = dz - (vz * FOCUS);
					k = (wx * wx) + (wy * wy) + (wz * wz);
					t2 = 1 / sqrt(k);
					lx = wx = wx * t2;
					ly = wy = wy * t2;
					lz = wz = wz * t2;
					w = (vx * wx) + (vy * wy) + (vz * wz);
					if (w == 0.0) {
						/*can't be*/
						return false;
					}
					t1 = ((dx * vx) + (dy * vy) + (dz * vz)) / w;
					if (t1 < 0.0) {
						lx = -lx; ly = -ly; lz = -lz;
					}
					dx = (wx * t1) - dx; dy = (wy * t1) - dy; dz = (wz * t1) - dz;
					w = (dx * dx) + (dy * dy) + (dz * dz);
					if (w > QUANT_R_2) {return false;}
					lensed = true;
					factor = (FOCUS_FACTOR * factor) / k;

					w = (dx * dx) + (dy * dy) + (dz * dz) - (QUANT_R * QUANT_R);
					k = (dx * vx) + (dy * vy) + (dz * vz);
					D = (k * k) - w;
					if (D <= 0.0) {
						/*can't be*/
						return false;
					}
					t1 = sqrt(D) - k + 0.0001;
					x0 = dx + t1 * vx + mx + 0.5;
					y0 = dy + t1 * vy + 0.5;
					z0 = dz + t1 * vz + mz + 0.5;
					goto restart;
					//return true;
				break;

				case MIRROR:
					number sx = .5 + mx - x0;
					number sy = CONE_H  - y0;
					number sz = .5 + mz - z0;
					number k = (sx * vx) + (sz * vz) - (sy * vy * INV_64_CONE_H_2);
					number coneC = (sx * sx) + (sz * sz) - (sy * sy * INV_64_CONE_H_2);
					D = k * k - coneA * coneC;
					if (D <= 0.0) { goto mirrorPlane; }
					D = sqrt(D);
					if (cAp) {
						t1 = (k - D) * invConeA;
						t2 = (k + D) * invConeA;
					} else {
						t2 = (k - D) * invConeA;
						t1 = (k + D) * invConeA;
					}
					y = y0 + vy * t1;
					if (y > CONE_H) {
						y = y0 + vy * t2;
						if (y > CONE_H) { goto mirrorPlane; }
						t1 = t2;
					}
					if (t1 < 0.0) { goto mirrorPlane;}
					if (y < 0.0) { goto mirrorPlane; }
					return false;
mirrorPlane:
					D = vx * nX[state] + vz * nZ[state];
					if (D == 0.0) {break;}
					t1 = ((mx + 0.5 - x0) * nX[state]) + ((mz + 0.5 - z0) * nZ[state]);
					t1 = t1 / D; // t1 == time
					t2 = y0 + (vy * t1) - 0.5; // t2 == y
					t3 = ((x0 + (vx * t1) - mx - .5) * wX[state]) + ((z0 + (vz * t1) - mz - .5) * wZ[state]); // t3 == x
					tmp = ((t2 * t2) + (t3 * t3));
					//if (tmp < MIRROR_R_2) {if (tmp > MIRROR_WR_2) {return false;}/* TODO: mirroring light? */return false;}
					if (tmp < MIRROR_R_2) {return false;}
// mirror torus
					tShi = ((mx + 0.5 - x0) * vx) + ((0.5 - y0) * vy) + ((mz + 0.5 - z0) * vz);
					x1 = x0 + (vx * tShi); y1 = y0 + (vy * tShi); z1 = z0 + (vz * tShi);
					t1 = x1 - mx - 0.5;
					t2 = y1 - 0.5;
					t3 = z1 - mz - 0.5;
					dx = (wX[state] * t1) + (wZ[state] * t3);
					dy = (wZ[state] * t1) - (wX[state] * t3);
					dz = t2;
					t1 = (wX[state] * vx) + (wZ[state] * vz);
					t2 = (wZ[state] * vx) - (wX[state] * vz);
					a = 2.0 * ((t1 * dx) + (vy * dz)); // r: t
					b = a + (2.0 * t2 * dy); // l^1 : t
					c = (dx * dx) + (dz * dz); // r: 1
					f = c + (dy * dy) + 0.0195; // l^1: 1
					e = (t1 * t1) + (vy * vy); // r: t^2
					//t4 = 1.0;
					t3 = 2.0 * b;
					t2 = (2.0 * f) + (b * b) - (0.0784 * e);
					t1 = (2.0 * b * f) - (0.0784 * a);
					t0 = (f * f) - (0.0784 * c);
					if (Solve(/*t4,*/t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {
						return false;
					}
				break;
			}

			if (tNextX < tNextZ) {
				if (tEnd <= tNextX) {return true;}
				if (dmx > 0) { if (mx == 14) { return true; } mx++; idx++;}
				else {if (mx == 0) { return true; } mx--; idx--;}
				tNextX += tmx;
			} else {
				if (tEnd <= tNextZ) { return true; }
				if (dmz > 0) {if (mz == 8) { return true; } mz++; idx+=15;}
				else {if (mz == 0) { return true; } mz--; idx-=15;}
				tNextZ += tmz;
			}
		}
		return true;
	}

	bool TestInside(const number x0, const number y0, const number z0, int mx, int mz) const {
		if (y0 > 1.0) { return false;}
		if ((x0 < 0.0) || (x0 > 15.0) || (z0 < 0.0) || (z0 > 9.0)) {return false;}
		number dx, dy, dz, r;
		int idx = mx + 15 * mz;
		int mdx, mdz;
		char type = rt::cells[idx].type;
		char state = rt::cells[idx].state;
		switch (type) {
			case BLOCK:
			case MBLOCK:
				mdx = ((x0 - mx) < 0.5) ? 1 : 2;
				mdz = ((z0 - mz) < 0.5) ? 1 : 2;
				if (rt::qcells[(mx * 2) + mdx + 32 * (mdz + (2 * mz))] == type) {return true;}
			break;

			case POLAR:
			case MPOLAR:
				return true;

			case QUANT:
			case STAR:
			case PORTAL:
				dx = x0 - mx - 0.5;
				dy = y0 - 0.5;
				dz = z0 - mz - 0.5;
				r = (dx * dx) + (dy * dy) + (dz * dz);
				if (r <= QUANT_R_2) {return true;}
			break;
		}
		return false;
	}

	bool TestCaps(const number x0, const number y0, const number z0, const number vx,const number vy,const number vz, const int mx, const int mz, number &nx, number &ny, number &nz, number tStart) const {
		if (y0 < 1.0) {return false;}
		double tmp = y0 + vy * tStart;
		if (tmp < 0.9999) { return false;}
		int idx = mx + (15 * mz);
		char type = rt::cells[idx].type;
		char state = rt::cells[idx].state;
		int mdx,mdz;
		char qtype;
		switch (type) {
			case BLOCK:
			case MBLOCK:
				mdx = ((x0 + (vx * tStart) - mx) < 0.5) ? 1 : 2;
				mdz = ((z0 + (vz * tStart) - mz) < 0.5) ? 1 : 2;
				qtype = rt::qcells[mdx + (mx * 2) + (32 * (mdz + (mz * 2)))];
				return qtype == type;

			case POLAR:
			case MPOLAR:
				tmp = (vx * poX[state]) + (vz * poZ[state]);
				if (tmp < 0.0) {tmp = -tmp;}
				return (tmp < (0.95 * sqrt((vx*vx)+(vz*vz))));
		}
		return false;
	}

	bool IntersectBlock(const number x0, const number y0, const number z0, const number vx,const number vy,const number vz, char state, const int mx, const int mz, number &nx, number &ny, number &nz,number &tIntersect) const {
		number t = LAST_T;
		number tmp;
		number tStart = tIntersect;
		int flags = bFlags[state];
		number nuz, nux, nuy;
		bool hasr = false;
		ny = 0.0;

		if (vx != 0.0) {
			number ivx = 1.0 / vx;
			number ivx2 = 0.5 * ivx;
			number vzivx2 = vz * ivx2;
			number vyivx2 = vy * ivx2;
			tmp = (mx - x0) * ivx;
			nuz = z0 + tmp * vz - mz;
			nuy = y0 + tmp * vy;

			if ((flags & VL) != 0) {
				if ((nuz >= 0.0) && (nuz <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nuz > 0.5) ? VLD : VLU)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nx = -1.0; nz = 0.0;}}
				}
			}
			tmp += ivx2; nuz += vzivx2; nuy += vyivx2;
			if ((flags & VM) != 0) {
				if ((nuz >= 0.0) && (nuz <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nuz > 0.5) ? VMD : VMU)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nz = 0.0; nx = (vx > 0.0)? -1.0 : 1.0;}}
				}
			}
			tmp += ivx2; nuz += vzivx2; nuy += vyivx2;
			if ((flags & VR) != 0) {
				if ((nuz >= 0.0) && (nuz <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nuz > 0.5) ? VRD : VRU)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nx= 1.0; nz = 0.0;}}
				}
			}
		}

		if (vz != 0.0) {
			number ivz = 1.0 / vz;
			number ivz2 = 0.5 * ivz;
			number vxivz2 = vx * ivz2;
			number vyivz2 = vy * ivz2;
			tmp = (mz - z0) * ivz;
			nux = x0 + tmp * vx - mx;
			nuy = y0 + tmp * vy;

			if ((flags & HU) != 0) {
				if ((nux >= 0.0) && (nux <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nux < 0.5) ? HUL : HUR)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nx=0.0; nz = -1.0;}}
				}
			}
			tmp += ivz2; nux += vxivz2; nuy += vyivz2;
			if ((flags & HM) != 0) {
				if ((nux >= 0.0) && (nux <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nux < 0.5) ? HML : HMR)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nx = 0.0; nz = (vz > 0.0)? -1.0 : 1.0;}}
				}
			}
			tmp += ivz2; nux += vxivz2; nuy += vyivz2;
			if ((flags & HD) != 0) {
				if ((nux >= 0.0) && (nux <= 1.0) && (nuy >= 0.0) && (nuy <= 1.0)) {
					if (((flags & ((nux < 0.5) ? HDL : HDR)) != 0) && (tmp > tStart)) {hasr = true; if (tmp < t) {t = tmp; nx = 0.0; nz = 1.0;}}
				}
			}
		}
		tIntersect = t;
		return hasr;
	}

	bool CalcTimeBox(const number x0, const number y0, const number z0, const number vx, const number vy, const number vz, number &tStart, number &tEnd, bool &boundsFloor) const {
		tStart = 0.0;
		tEnd = LAST_T;
		// XBOX
		if (vx > 0.0) {
			number t = -x0 / vx;
			if (t > tStart) {tStart = t;}
			t = (15.0 - x0) / vx;
			if (t < tEnd) {tEnd = t;}
		} else if (vx < 0.0) {
			number t = (15.0 - x0) / vx;
			if (t > tStart) {tStart = t;}
			t = -x0 / vx;
			if (t < tEnd) {tEnd = t;}
		} else {
			if ((x0 < 0.0) || (x0 > 15.0)) {return false;}
		}

		// YBOX
		if (vy > 0.0) {
			boundsFloor = false;
			number t = -y0 / vy;
			if (t > tStart) {tStart = t;}
			t = (1.0 - y0) / vy;
			if (t < tEnd) {tEnd = t;}
		} else if (vy < 0.0) {
			boundsFloor = true;
			number t = (1.0 - y0) / vy;
			if (t > tStart) {tStart = t;}
			t = -y0 / vy;
			if (t < tEnd) {tEnd = t;}
		} else {
			boundsFloor = false;
			if ((y0 < 0.0) || (y0 > 1.0)) {return false;}
		}

		// ZBOX
		if (vz > 0.0) {
			number t = -z0 / vz;
			if (t > tStart) {tStart = t;}
			t = (9.0 - z0) / vz;
			if (t < tEnd) {tEnd = t;}
		} else if (vz < 0.0) {
			number t = (9.0 - z0) / vz;
			if (t > tStart) {tStart = t;}
			t = -z0 / vz;
			if (t < tEnd) {tEnd = t;}
		} else {
			if ((z0 < 0.0) || (z0 > 9.0)) {return false;}
		}
		
		if (tEnd <= tStart) {
			return false;
		}
		return true;
	}

	void CalcPoint(const int cx, const int cy, number x0, number y0, number z0, number* clr) const {
#ifdef PROTECT
	try {
#endif
		VEC4D* ray0;
#ifndef ANTI_ALIASING
		ray0 = &(modelRays[cx + (256 * cy)]);
#else
		ray0 = &(modelRaysA[cx + (1024 * cy)]);
#endif
		number vx = ray0->x;
		number vy = ray0->y;
		number vz = ray0->z;
		number coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
		number invConeA = 1 / coneA;
		bool cAp = coneA > 0.0;
		number sx,sy,sz,k,coneC, econeA, einvConeA;
		number D, t0, t1, t2, t3, t4, tmp, temp, tShi, x1, y1, z1, vx1, vz1;
		number nx,ny,nz;
		number a,b,c,e,f,g;
		number tResult;
		DWORD result;
		number rx,ry,rz;
		number lx,ly,lz;
		bool lensed = false;
		number tlr;
		number x, y, z;
		number dx,dy,dz,w, cot;
		int d;
		int deepness = 0;

		clr[0] = clr[1] = clr[2] = 0.0;

		number colorFactor = 1.0;
		VEC4D* mat;

		bool precalcActive = true;

		number intense;

		number tStart, tEnd;
		bool boundsFloor;

		number cox, coy, coz, conx, cony, conz;
		int doWhat;
		number rayT1[32], rayT2[32];
		bool hasRay = false;
		int dir;
		int rayIdx = 0;
		int doThing;


		if (!CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor)) { RETURN_BLACK }
// calculate scene ray time
		tStart += 0.000000001;
		tEnd   += 0.000000001;

		number rayLast = tEnd, rayFirst = tStart;

		int mx = (int)(x0 + (tStart * vx));
		int mz = (int)(z0 + (tStart * vz));

		if (TestInside(x0, y0, z0, mx, mz)) {
			clr[0] = clr[1] = clr[2] = (rand() & 0xFF);
			return;
		}

		number tNextX, tmx;
		int dmx;
		if (vx > 0.0) {
			tmx = 1.0 / vx;
			tNextX = (mx + 1.0 - x0) / vx;
			dmx = 1;
		} else if (vx < 0.0) {
			tmx = -1.0 / vx;
			tNextX = (mx - x0) / vx;
			dmx = -1;
		} else {
			tNextX = LAST_T;
		}

		number tNextZ, tmz;
		int dmz;
		if (vz > 0.0) {
			tmz = 1.0 / vz;
			tNextZ = (mz + 1.0 - z0) / vz;
			dmz = 1;
		} else if (vz < 0.0) {
			tmz = -1.0 / vz;
			tNextZ = (mz - z0) / vz;
			dmz = -1;
		} else {
			tNextZ = LAST_T;
		}

		int idx = mx + 15 * mz;
		if (TestCaps(x0, y0, z0, vx, vy, vz, mx, mz, nx, ny, nz, tStart)) {
			nx = 0.0; nz = 0.0; ny = 1.0;
			if ((rt::cells[idx].type != BLOCK) && (rt::cells[idx].type != POLAR)) {mat = &(materials[MAT_DIF_DEF]);}
			else { mat = &(materials[MAT_DIF_BLO]);}
			for (int i = 0; i < 3; i++) {
				intense = DIFFUSE(i) * colorFactor;
				ADD_COLOR
			}
			RETURN_COLOR
		}
		while (true) {
			APPLY_RAY
			RAY_PART* part = rayThings[mx + 15 * mz];
			rayLast = LAST_T;
			hasRay = false;
			rayIdx = 0;
			while (part != NULL) {
				dir = part->dir;
				a = (roX[dir] * vx) + (roZ[dir] * vz);
				b = part->d - (roX[dir] * x0) - (roZ[dir] * z0);
				tmp = 0.5 - y0;
				k = (a * b) + (vy * tmp);
				c = (b * b) + (tmp * tmp) - RAY_R_2;
				a = (a * a) + (vy * vy);
				D = (k * k) - (a * c);
				if (D > 0.0) {
					a = 1 / a;
					D = sqrt(D);
					t1 = a * (k - D);
					t2 = a * (k + D);
					t3 = rayFirst;
					t4 = rayLast;

					tmp = (vx * part->staX) + (vz * part->staZ);
					if (tmp > 0.0) {tmp = (part->sta - (x0 * part->staX) - (z0 * part->staZ)) / tmp; if (t3<tmp){t3=tmp;}}
					else if (tmp < 0.0) {tmp = (part->sta - (x0 * part->staX) - (z0 * part->staZ)) / tmp; if (t4>tmp){t4=tmp;}}

					tmp = (vx * part->endX) + (vz * part->endZ);
					if (tmp > 0.0) {tmp = (part->end - (x0 * part->endX) - (z0 * part->endZ)) / tmp; if (t3<tmp){t3=tmp;}}
					else if (tmp < 0.0) {tmp = (part->end - (x0 * part->endX) - (z0 * part->endZ)) / tmp; if (t4>tmp){t4=tmp;}}

					if (t1 < t3) {t1 = t3;}
					if (t2 > t4) {t2 = t4;}
					if (t2 > t1 + 0.00001) {
						hasRay = true;
						rayT1[rayIdx] = t1;
						rayT2[rayIdx] = t2;
						rayIdx++;
					}
				}
				part = part->next;
			}

			VEC4D* tgt = &precalc[idx];
			char type = rt::cells[idx].type;
			char state = rt::cells[idx].state;
			tResult = LAST_T;

			switch (type) {
				case POLAR:
				case MPOLAR:
					tmp = (vx * poX[state]) + (vz * poZ[state]);
					if (tmp < 0.0) {tmp = -tmp;}
					if (tmp > 0.95 * sqrt((vx*vx)+(vz*vz))) {break;}

					if (type != POLAR) {mat = &(materials[MAT_DIF_DEF]);}
					else {mat = &(materials[MAT_DIF_BLO]);}
					t1 = tStart - 0.0001;
					IntersectBlock(x0, y0, z0, vx, vy, vz, FULL, mx, mz, nx, ny, nz, t1);
rayLast = t1;
					t2 = t1 - 0.0001;
					x = x0 + (vx * t2);
					y = y0 + (vy * t2);
					z = z0 + (vz * t2);
					for (int i = 0; i < 3; i++) {
						if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
							intense = DIFFUSE(i) * colorFactor;
							ADD_COLOR
						}
					}
					RETURN_COLOR

				case RANDOM:
// TODO: rayLast = t1;
					clr[0] += ((int)((vx * vy + z0) * (65536 + subFrame)) & 0xFF) / 255.0;
					clr[1] += ((int)((vy * vz + x0) * (65536 + subFrame)) & 0xFF) / 255.0;
					clr[2] += ((int)((vz * vx + y0) * (65536 + subFrame)) & 0xFF) / 255.0;
					RETURN_COLOR

				case PORTAL:
					temp = LAST_T;
					if (precalcActive) {
						k = (tgt->x * vx) + (tgt->y * vy) + (tgt->z * vz);
						D = (k * k) - tgt->w;
					} else {
						dx = mx + 0.5 - x0;
						dy =      0.5 - y0;
						dz = mz + 0.5 - z0;
						w = (dx * dx) + (dy * dy) + (dz * dz) - (PORTAL_R * PORTAL_R);
						k = (dx * vx) + (dy * vy) + (dz * vz);
						D = (k * k) - w;
					}
					if (D <= 0.0) { goto portalTorusC; }
					t1 = k - sqrt(D);
					if (t1 < tStart) {
						t1 = k + sqrt(D);
						if (t1 > tStart) {
							// is that right?
							RETURN_COLOR
						}
						// may be break?
						goto portalTorusC;
					}
					temp = t1;
portalTorusC:
					tShi = ((mx + 0.5 - x0) * vx) + ((0.1 - y0) * vy) + ((mz + 0.5 - z0) * vz);
					x1 = x0 + (vx * tShi); y1 = y0 + (vy * tShi); z1 = z0 + (vz * tShi);
					dx = x1 - mx - 0.5;
					dy = y1 - 0.1;
					dz = z1 - mz - 0.5;
					a = 2.0 * ((vx * dx) + (vz * dz)); // r: t
					b = a + (2.0 * vy * dy); // l^1 : t
					c = (dx * dx) + (dz * dz); // r: 1
					f = c + (dy * dy) + 0.15; // l^1: 1
					e = (vx * vx) + (vz * vz); // r: t^2
					//t4 = 1.0;
					t3 = 2.0 * b;
					t2 = (2.0 * f) + (b * b) - (0.64 * e);
					t1 = (2.0 * b * f) - (0.64 * a);
					t0 = (f * f) - (0.64 * c);
					if (!Solve(/*t4,*/t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {
						if (temp < LAST_T) {goto portalPortal;}
						break;
					}
					tmp += tShi;
					if (temp < tmp) {goto portalPortal;}
//portalTorus
					t1 = tmp;
					t2 = t1 - 0.000001;
					x = x0 + (vx * t2);
					y = y0 + (vy * t2);
					z = z0 + (vz * t2);
					dx = x - mx - 0.5;
					dy = y - 0.1;
					dz = z - mz - 0.5;
					f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.17);
					nx = dx * f;
					ny = dy * (f + 1.28);
					nz = dz * f;
					e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
					nx = nx * e;
					ny = ny * e;
					nz = nz * e;
					tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
					rx = tmp * nx - vx;
					ry = tmp * ny - vy;
					rz = tmp * nz - vz;
					for (int i = 0; i < 3; i++) {
						if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
							mat = &(materials[MAT_DIF_TOR]);
							intense = DIFFUSE(i) * colorFactor;
							ADD_COLOR
							intense = SPECULAR(i);
							if (intense < 0.0) {
								intense = intense * intense;
								intense = intense * intense * colorFactor;
								mat = &(materials[MAT_SPE_TOR]);
								ADD_COLOR
							}
						}
					}
					RETURN_COLOR
					break;
portalPortal:
					t1 = temp;
					precalcActive = false;
					deepness++;
					if (deepness > 16) {
						clr[0] += (rand() & 0x3F) / 255.0;
						clr[1] += (rand() & 0x3F) / 255.0;
						clr[2] += (rand() & 0x3F) / 255.0;
						RETURN_COLOR
					}
rayLast = t1;
					t2 = t1 - 0.00001;
					x = x0 + (vx * t2) - mx;
					y = y0 + (vy * t2);
					z = z0 + (vz * t2) - mz;
					x = 1.0 - x;
					y = 1.0 - y;
					z = 1.0 - z;
					d = 1;
					if ((levelData.px[state * 2] != mx) || (levelData.py[state * 2] != mz)) {d = 0;}
					mx = levelData.px[state * 2 + d];
					mz = levelData.py[state * 2 + d];
					x0 = mx + x;
					y0 = y;
					z0 = mz + z;
					CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);

					tStart += 0.000000001;
					tEnd   += 0.000000001;

					if (vx > 0.0) {tNextX = (mx + 1.0 - x0) / vx;} else if (vx < 0.0) {tNextX = (mx - x0) / vx;}
					if (vz > 0.0) {tNextZ = (mz + 1.0 - z0) / vz;} else if (vz < 0.0) {tNextZ = (mz - z0) / vz;}
					idx = mx + 15 * mz;
					temp = LAST_T;
					goto portalTorusC;
					break;

				case BLOCK:
				case MBLOCK:
					t1 = tStart - 0.0001;
					if (IntersectBlock(x0, y0, z0, vx, vy, vz, state, mx, mz, nx, ny, nz, t1)) {
rayLast = t1;
						if (type != BLOCK) {mat = &(materials[MAT_DIF_DEF]);}
						else {mat = &(materials[MAT_DIF_BLO]);}
						t2 = t1 - 0.000001;
						x = x0 + (vx * t2);
						y = y0 + (vy * t2);
						z = z0 + (vz * t2);
						for (int i = 0; i < 3; i++) {
							if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
								intense = DIFFUSE(i) * colorFactor;
								ADD_COLOR
							}
						}
						RETURN_COLOR
					}
				break;

				case QUANT:
				case STAR:
					if (precalcActive) {
						k = (tgt->x * vx) + (tgt->y * vy) + (tgt->z * vz);
						D = (k * k) - tgt->w;
					} else {
						dx = mx + 0.5 - x0;
						dy =      0.5 - y0;
						dz = mz + 0.5 - z0;
						w = (dx * dx) + (dy * dy) + (dz * dz) - (QUANT_R * QUANT_R);
						k = (dx * vx) + (dy * vy) + (dz * vz);
						D = (k * k) - w;
					}
					if (D <= 0.0) { break; }
					t1 = k - sqrt(D);
					if (t1 < tStart) {
						t1 = k + sqrt(D);
						if (t1 > tStart) {
rayLast = t1;
							RETURN_COLOR
						}
						break;
					}
					t2 = t1 - 0.00001;
rayLast = t1;
					x = x0 + (vx * t2);
					y = y0 + (vy * t2);
					z = z0 + (vz * t2);
					nx = (x - mx - .5) * INV_QUANT_R;
					ny = (y      - .5) * INV_QUANT_R;
					nz = (z - mz - .5) * INV_QUANT_R;
					tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
					rx = tmp * nx - vx;
					ry = tmp * ny - vy;
					rz = tmp * nz - vz;
					if (type == QUANT) {
#ifdef GLASS_POINT
							APPLY_RAY
							for (int i = 0; i < 3; i++) {
								mat = &(materials[MAT_SPE_PNT]);
								if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
									intense = SPECULAR(i);
									if (intense < 0.0) {
										intense = intense * intense;
										intense = intense * intense * colorFactor;
										ADD_COLOR
									}
								}
							}
							t1 = -((vx * nx) + (vy * ny) + (vz * nz));
							t2 = 1.0 - (GLASS_K * GLASS_K) * (1 - (t1 * t1));
							if (t2 < 0.0) {
								t1 = t1 * 2.0;
								vx = vx + (t1 * nx);
								vy = vy + (t1 * ny);
								vz = vz + (t1 * nz);
								goto fullReflect;
							}
							t2 = sqrt(t2);
							t1 = t1 * GLASS_K;
							if (t1 < 0.0) {t2 = -t2;}
							vx = (GLASS_K * vx) + (t1 - t2) * nx;
							vy = (GLASS_K * vy) + (t1 - t2) * ny;
							vz = (GLASS_K * vz) + (t1 - t2) * nz;

							dx = mx + 0.5 - x;
							dy =      0.5 - y;
							dz = mz + 0.5 - z;

							w = (dx * dx) + (dy * dy) + (dz * dz) - (QUANT_R * QUANT_R);
							k = (dx * vx) + (dy * vy) + (dz * vz);
							D = (k * k) - w;
							if (D < 0.0) {RETURN_COLOR}
							t1 = k + sqrt(D) + 0.00001;
							x0 = x + (vx * t1);
							y0 = y + (vy * t1);
							z0 = z + (vz * t1);

							nx = -(x0 - mx - .5) * INV_QUANT_R;
							ny = -(y0      - .5) * INV_QUANT_R;
							nz = -(z0 - mz - .5) * INV_QUANT_R;

							t1 = -((vx * nx) + (vy * ny) + (vz * nz));
							t2 = 1.0 - (INV_GLASS_K * INV_GLASS_K) * (1 - (t1 * t1));
							if (t2 < 0.0) {RETURN_COLOR}
							t2 = sqrt(t2);
							t1 = t1 * INV_GLASS_K;
							if (t1 < 0.0) {t2 = -t2;}
							vx = (INV_GLASS_K * vx) + (t1 - t2) * nx;
							vy = (INV_GLASS_K * vy) + (t1 - t2) * ny;
							vz = (INV_GLASS_K * vz) + (t1 - t2) * nz;
fullReflect:
							precalcActive = false;
							deepness++; if (deepness > 16) {
								clr[0] += (rand() & 0x3F) / 255.0;
								clr[1] += (rand() & 0x3F) / 255.0;
								clr[2] += (rand() & 0x3F) / 255.0;
								RETURN_COLOR
							}
							CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);

							tStart += 0.000000001;
							tEnd   += 0.000000001;
							rayFirst = tStart;
							rayLast = tEnd;

							if (vx > 0.0) {tmx = 1.0 / vx; tNextX = (mx + 1.0 - x0) / vx; dmx = 1;}
							else if (vx < 0.0) {tmx = -1.0 / vx;tNextX = (mx - x0) / vx; dmx = -1;}
							else {tNextX = LAST_T;}

							if (vz > 0.0) {tmz = 1.0 / vz; tNextZ = (mz + 1.0 - z0) / vz; dmz = 1;}
							else if (vz < 0.0) {tmz = -1.0 / vz; tNextZ = (mz - z0) / vz; dmz = -1;}
							else { tNextZ = LAST_T;}

							coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
							invConeA = 1 / coneA;
							cAp = coneA > 0.0;

							goto noCycle;
#else
							for (int i = 0; i < 3; i++) {
								if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
									mat = &(materials[MAT_DIF_PNT]);
									intense = DIFFUSE(i) * colorFactor;
									ADD_COLOR
									intense = SPECULAR(i);
									if (intense < 0.0) {
										intense = intense * intense;
										intense = intense * intense * colorFactor;
										mat = &(materials[MAT_SPE_PNT]);
										ADD_COLOR
									}
								}
							}
#endif
						RETURN_COLOR
					} else {
						precalcActive = false;
						deepness++; if (deepness > 16) {
							clr[0] += (rand() & 0x3F) / 255.0;
							clr[1] += (rand() & 0x3F) / 255.0;
							clr[2] += (rand() & 0x3F) / 255.0;
							RETURN_COLOR
						}
						x0 = x;
						y0 = y;
						z0 = z;
						vx = -rx;
						vy = -ry;
						vz = -rz;
						CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);

						tStart += 0.000000001;
						tEnd   += 0.000000001;

						if (vx > 0.0) {tmx = 1.0 / vx; tNextX = (mx + 1.0 - x0) / vx; dmx = 1;}
						else if (vx < 0.0) {tmx = -1.0 / vx;tNextX = (mx - x0) / vx; dmx = -1;}
						else {tNextX = LAST_T;}


						if (vz > 0.0) {tmz = 1.0 / vz; tNextZ = (mz + 1.0 - z0) / vz; dmz = 1;}
						else if (vz < 0.0) {tmz = -1.0 / vz; tNextZ = (mz - z0) / vz; dmz = -1;}
						else { tNextZ = LAST_T;}

						coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
						invConeA = 1 / coneA;
						cAp = coneA > 0.0;
					}
				break;

				case COLLECTOR:
					tResult = LAST_T;
tShi = ((mx + 0.5 - x0) * vx) + ((0.5 - y0) * vy) + ((mz + 0.5 - z0) * vz);
x1 = x0 + (vx * tShi); y1 = y0 + (vy * tShi); z1 = z0 + (vz * tShi);

					if (state == 0) {
dx = x1 - mx - 0.03;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx) + 0.0216; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll1;}
tmp = tmp + tShi;

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.03;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.0234);
nx = dx * (f + 0.18);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;
  
coll1:
dx = x1 - mx - 0.3;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (0.04 * dx * dx) + 0.0084; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
temp = e + (vx * vx * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll2;}
tmp = tmp + tShi;
if (tmp > tResult) {goto coll2;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.3;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dz * dz) + (dy * dy) + (0.04 * dx * dx) - 0.0116);
nx = 0.2 * dx * (f + 0.08);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

coll2:
dx = x1 - mx - 0.6;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx); // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto coll3;}
tmp = tmp + tShi;
if (tmp > tResult) {goto coll3;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.6;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.02);
nx = dx * (f + 0.08);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

coll3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.9 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto coll4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto coll4;}
if (t1 > tResult) {goto coll4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox > (mx + 1.0)) { goto coll4;}
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 3.333333333333333333 * 1.2 * (cox - mx - 0.9);
cony = 3.333333333333333333 * (coy - 0.5);
conz = 3.333333333333333333 * (coz - 0.5 - mz);
cot = t1;
tResult = t1;
coll4:
if (vx == 0.0) {goto col5;}
t1 = (mx + 1.0 - x0) / vx + 0.000001;
if (t1 < rayFirst) {goto col5;}
if (t1 > tResult) {goto col5;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {goto col5;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 1.0;
cony = 0.0;
conz = 0.0;
cot = t1;
tResult = t1;

					} else if (state == 1) {

dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.03;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz) + 0.0216; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colu1;}
tmp = tmp + tShi;

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.03;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.0234);
nx = dx * f;
ny = dy * f;
nz = dz * (f + 0.18);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;
  
colu1:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.3;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (0.04 * dz * dz) + 0.0084; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
temp = e + (vz * vz * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto colu2;}
tmp = tmp + tShi;
if (tmp > tResult) {goto colu2;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.3;
f = 4.0 * ((dx * dx) + (dy * dy) + (0.04 * dz * dz) - 0.0116);
nx = dx * f;
ny = dy * f;
nz = 0.2 * dz * (f + 0.08);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

colu2:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.60;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz); // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colu3;}
tmp = tmp + tShi;
if (tmp > tResult) {goto colu3;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.60;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.02);
nx = dx * f;
ny = dy * f;
nz = dz * (f + 0.08);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

colu3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.9 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto colu4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto colu4;}
if (t1 > tResult) {goto colu4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz > (mz + 1.0)) { goto colu4;}
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 3.333333333333333333 * 1.2 * (coz - mz - 0.9);
cony = 3.333333333333333333 * (coy - 0.5);
conx = 3.333333333333333333 * (cox - 0.5 - mx);
cot = t1;
tResult = t1;
colu4:
if (vz == 0.0) {goto col5;}
t1 = (mz + 1.0 - z0) / vz + 0.000001;
if (t1 < rayFirst) {goto col5;}
if (t1 > tResult) {goto col5;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {goto col5;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 1.0;
cony = 0.0;
conx = 0.0;
cot = t1;
tResult = t1;

					} else if (state == 2) {

dx = x1 - mx - 0.97;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx) + 0.0216; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr1;}
tmp = tmp + tShi;

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.97;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.0234);
nx = dx * (f + 0.18);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;
  
colr1:
dx = x1 - mx - 0.7;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (0.04 * dx * dx) + 0.0084; // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
temp = e + (vx * vx * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr2;}
tmp = tmp + tShi;
if (tmp > tResult) {goto colr2;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.7;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dz * dz) + (dy * dy) + (0.04 * dx * dx) - 0.0116);
nx = 0.2 * dx * (f + 0.08);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

colr2:
dx = x1 - mx - 0.4;
dy = y1 - 0.5;
dz = z1 - mz - 0.5;
a = 2.0 * ((vz * dz) + (vy * dy)); // r: t
b = a + (2.0 * vx * dx); // l^1 : t
c = (dz * dz) + (dy * dy); // r: 1
f = c + (dx * dx); // l^1: 1
e = (vz * vz) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto colr3;}
tmp = tmp + tShi;
if (tmp > tResult) {goto colr3;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.4;
dy = coy - 0.5;
dz = coz - mz - 0.5;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.02);
nx = dx * (f + 0.08);
ny = dy * f;
nz = dz * f;
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

colr3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.1 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto colr4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto colr4;}
if (t1 > tResult) {goto colr4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox < mx) { goto colr4;}
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 3.333333333333333333 * 1.2 * (cox - mx - 0.1);
cony = 3.333333333333333333 * (coy - 0.5);
conz = 3.333333333333333333 * (coz - 0.5 - mz);
cot = t1;
tResult = t1;
colr4:
if (vx == 0.0) {goto col5;}
t1 = (mx - x0) / vx + 0.000001;
if (t1 + 0.00002 < rayFirst) {goto col5;}
if (t1 > tResult) {goto col5;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {goto col5;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = -1.0;
cony = 0.0;
conz = 0.0;
cot = t1;
tResult = t1;

					} else if (state == 3) {
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.97;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz) + 0.0216; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.09 * e);
t1 = (2.0 * b * f) - (0.09 * a);
t0 = (f * f) - (0.09 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto cold1;}
tmp = tmp + tShi;

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.97;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.0234);
nx = dx * f;
ny = dy * f;
nz = dz * (f + 0.18);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

cold1:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.7;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (0.04 * 2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (0.04 * dz * dz) + 0.0084; // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
temp = e + (vz * vz * 0.04);
t3 = 2.0 * b * temp;
t2 = (2.0 * temp * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
temp = 1 / (temp * temp);
if (!Solve(temp * t3,temp * t2,temp * t1,temp * t0, tStart - tShi, tEnd - tShi, tmp)) {goto cold2;}
tmp = tmp + tShi;
if (tmp > tResult) {goto cold2;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.7;
f = 4.0 * ((dx * dx) + (dy * dy) + (0.04 * dz * dz) - 0.0116);
nx = dx * f;
ny = dy * f;
nz = 0.2 * dz * (f + 0.08);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

cold2:
dx = x1 - mx - 0.5;
dy = y1 - 0.5;
dz = z1 - mz - 0.40;
a = 2.0 * ((vx * dx) + (vy * dy)); // r: t
b = a + (2.0 * vz * dz); // l^1 : t
c = (dx * dx) + (dy * dy); // r: 1
f = c + (dz * dz); // l^1: 1
e = (vx * vx) + (vy * vy); // r: t^2
t3 = 2.0 * b;
t2 = (2.0 * f) + (b * b) - (0.04 * e);
t1 = (2.0 * b * f) - (0.04 * a);
t0 = (f * f) - (0.04 * c);
if (!Solve(t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {goto cold3;}
tmp = tmp + tShi;
if (tmp > tResult) {goto cold3;}

tResult = tmp;
t2 = tmp - 0.000001;
coy = y0 + vy * t2;
cox = x0 + vx * t2;
coz = z0 + vz * t2;
dx = cox - mx - 0.5;
dy = coy - 0.5;
dz = coz - mz - 0.40;
f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.02);
nx = dx * f;
ny = dy * f;
nz = dz * (f + 0.08);
e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
conx = nx * e;
cony = ny * e;
conz = nz * e;
cot = tmp;

cold3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.1 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto cold4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto cold4;}
if (t1 > tResult) {goto cold4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz < mz) { goto cold4;}
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 3.333333333333333333 * 1.2 * (coz - mz - 0.1);
cony = 3.333333333333333333 * (coy - 0.5);
conx = 3.333333333333333333 * (cox - 0.5 - mx);
cot = t1;
tResult = t1;
cold4:
if (vz == 0.0) {goto col5;}
t1 = (mz - z0) / vz + 0.000001;
if (t1 + 0.00002 < rayFirst) {goto col5;}
if (t1 > tResult) {goto col5;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {goto col5;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = -1.0;
cony = 0.0;
conx = 0.0;
cot = t1;
tResult = t1;

					}
col5:
					if (tResult < LAST_T) {
rayLast = cot;
						x = cox; y = coy; z = coz; nx = conx; ny = cony; nz = conz;
						tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
						rx = tmp * nx - vx;
						ry = tmp * ny - vy;
						rz = tmp * nz - vz;

						mat = &(materials[MAT_DIF_CON]);
						for (int i = 0; i < 3; i++) {
							if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
//								intense = DIFFUSE(i) * colorFactor;
//								ADD_COLOR
								intense = SPECULAR(i);
								if (intense < 0.0) {
//									intense = intense * intense;
//									intense = intense * intense * colorFactor;
									intense = -2.0 * intense * intense * intense * colorFactor;
									ADD_COLOR
								}
							}
						}
						clr[0] += 0.35; clr[1] += 0.25; clr[2] += 0.15;
						RETURN_COLOR
					}
					break;

				case EMITTER:
					tResult = LAST_T;
					doThing = -1;
					if (state == 0) {
//#include<emi_l_p.cpp>
sx = mx - 0.1 - x0;
sy = 0.5 - y0;
sz = (0.5 + mz) - z0;
k = (sy * vy) + (sz * vz) - (sx * vx * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sz * sz) - (sx * sx * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vz * vz) - (vx * vx * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emil1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
x = (mx + 1.0 - 0.1) - (x0 + (vx * t1));
if (x > E_CONE_H - 0.1) {
	x = (mx + 1.0 - 0.1) - (x0 + (vx * t2));
	if (x > E_CONE_H - 0.1) { goto emil1; }
	t1 = t2;
}
if (x < 0.1) { goto emil1; }
if (t1 < rayFirst) {goto emil1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * x); // tmp == r
if (tmp <= 0.0) {goto emil1;}
tmp = 1 / tmp;
cot = t1;
t2 = t1 - 0.00001;
x = x0 + vx * t2;
y = y0 + vy * t2;
z = z0 + vz * t2;
nx = -E_CONE_DY;
ny =  (y - 0.5) * tmp;
nz = (z - mz - 0.5) * tmp;
tResult = t1; doThing = 0; cox = x; coy = y; coz = z; conx= nx; cony = ny; conz = nz;
emil1:
if (vx == 0.0) {goto emil2;}
t1 = (mx - x0) / vx;
if (t1 < rayFirst) {goto emil2;}
if (t1 > tResult) {goto emil2;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > RAY_R_2) {goto emil2;}
cot = t1;
tResult = t1; doThing = 1;
emil2:
a = (vx * vx * 16.0) + (vy * vy) + (vz * vz);
x = mx + 0.25 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (16.0 * x * vx) + (y * vy) + (z * vz);
c = (16.0 * x * x) + (y * y) + (z * z) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emil3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emil3;}
if (t1 > tResult) {goto emil3;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 5.0 * 4.0 * (cox - mx - 0.25);
cony = 5.0 * (coy - 0.5);
conz = 5.0 * (coz - 0.5 - mz);
cot = t1;
tResult = t1; doThing = 2;
emil3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.9 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emil4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emil4;}
if (t1 > tResult) {goto emil4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox > (mx + 1.0)) { goto emil4;}
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 3.333333333333333333 * 1.2 * (cox - mx - 0.9);
cony = 3.333333333333333333 * (coy - 0.5);
conz = 3.333333333333333333 * (coz - 0.5 - mz);
cot = t1;
tResult = t1; doThing = 3;
emil4:
if (vx == 0.0) {goto emi5;}
t1 = (mx + 1.0 - x0) / vx + 0.000001;
if (t1 < rayFirst) {goto emi5;}
if (t1 > tResult) {goto emi5;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {goto emi5;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 1.0;
cony = 0.0;
conz = 0.0;
cot = t1;
tResult = t1; doThing = 4;

					} else if (state == 1) {
//#include<emi_u_p.cpp>
sz = mz - 0.1 - z0;
sy = 0.5 - y0;
sx = (0.5 + mx) - x0;
k = (sy * vy) + (sx * vx) - (sz * vz * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sx * sx) - (sz * sz * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vx * vx) - (vz * vz * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emiu1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
z = (mz + 1.0 - 0.1) - (z0 + (vz * t1));
if (z > E_CONE_H - 0.1) {
	z = (mz + 1.0 - 0.1) - (z0 + (vz * t2));
	if (z > E_CONE_H - 0.1) { goto emiu1; }
	t1 = t2;
}
if (z < 0.1) { goto emiu1; }
if (t1 < rayFirst) {goto emiu1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * z); // tmp == r
if (tmp <= 0.0) {goto emiu1;}
tmp = 1 / tmp;
cot = t1;
t2 = t1 - 0.00001;
z = z0 + vz * t2;
y = y0 + vy * t2;
x = x0 + vx * t2;
nz = -E_CONE_DY;
ny =  (y - 0.5) * tmp;
nx = (x - mx - 0.5) * tmp;
tResult = t1; doThing = 0; cox = x; coy = y; coz = z; conx= nx; cony = ny; conz = nz;
emiu1:
if (vz == 0.0) {goto emiu2;}
t1 = (mz - z0) / vz;
if (t1 < rayFirst) {goto emiu2;}
if (t1 > tResult) {goto emiu2;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > RAY_R_2) {goto emiu2;}
cot = t1;
tResult = t1; doThing = 1;
emiu2:
a = (vz * vz * 16.0) + (vy * vy) + (vx * vx);
z = mz + 0.25 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (16.0 * z * vz) + (y * vy) + (x * vx);
c = (16.0 * z * z) + (y * y) + (x * x) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emiu3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emiu3;}
if (t1 > tResult) {goto emiu3;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 5.0 * 4.0 * (coz - mz - 0.25);
cony = 5.0 * (coy - 0.5);
conx = 5.0 * (cox - 0.5 - mx);
cot = t1;
tResult = t1; doThing = 2;
emiu3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.9 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emiu4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emiu4;}
if (t1 > tResult) {goto emiu4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz > (mz + 1.0)) { goto emiu4;}
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 3.333333333333333333 * 1.2 * (coz - mz - 0.9);
cony = 3.333333333333333333 * (coy - 0.5);
conx = 3.333333333333333333 * (cox - 0.5 - mx);
cot = t1;
tResult = t1; doThing = 3;
emiu4:
if (vz == 0.0) {goto emi5;}
t1 = (mz + 1.0 - z0) / vz + 0.000001;
if (t1 < rayFirst) {goto emi5;}
if (t1 > tResult) {goto emi5;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {goto emi5;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 1.0;
cony = 0.0;
conx = 0.0;
cot = t1;
tResult = t1; doThing = 4;

					} else if (state == 2) {
//#include<emi_r_p.cpp>
sx = mx + 1.1 - x0;
sy = 0.5 - y0;
sz = (0.5 + mz) - z0;
k = (sy * vy) + (sz * vz) - (sx * vx * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sz * sz) - (sx * sx * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vz * vz) - (vx * vx * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emir1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
x = (mx + 0.1) - (x0 + (vx * t1));
if (x < -E_CONE_H + 0.1) {
	x = (mx + 0.1) - (x0 + (vx * t2));
	if (x < -E_CONE_H + 0.1) { goto emir1; }
	t1 = t2;
}
x = -x;
if (x < 0.1) { goto emir1; }
if (t1 < rayFirst) {goto emir1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * x); // tmp == r
if (tmp <= 0.0) {goto emir1;}
tmp = 1 / tmp;
cot = t1;
t2 = t1 - 0.00001;
x = x0 + vx * t2;
y = y0 + vy * t2;
z = z0 + vz * t2;
nx = E_CONE_DY;
ny =  (y - 0.5) * tmp;
nz = (z - mz - 0.5) * tmp;
tResult = t1; doThing = 0; coz = z; coy = y; cox = x; conz= nz; cony = ny; conx = nx;
emir1:
if (vx == 0.0) {goto emir2;}
t1 = (mx + 1.0 - x0) / vx;
if (t1 < rayFirst) {goto emir2;}
if (t1 > tResult) {goto emir2;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > RAY_R_2) {goto emir2;}
cot = t1;
tResult = t1; doThing = 1;
emir2:
a = (vx * vx * 16.0) + (vy * vy) + (vz * vz);
x = mx + 0.75 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (16.0 * x * vx) + (y * vy) + (z * vz);
c = (16.0 * x * x) + (y * y) + (z * z) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emir3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emir3;}
if (t1 > tResult) {goto emir3;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 5.0 * 4.0 * (cox - mx - 0.75);
cony = 5.0 * (coy - 0.5);
conz = 5.0 * (coz - 0.5 - mz);
cot = t1;
tResult = t1; doThing = 2;
emir3:
a = (vx * vx * 1.44) + (vy * vy) + (vz * vz);
x = mx + 0.1 - x0;
y = 0.5 - y0;
z = mz + 0.5 - z0;
k = (1.44 * x * vx) + (y * vy) + (z * vz);
c = (1.44 * x * x) + (y * y) + (z * z) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emir4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emir4;}
if (t1 > tResult) {goto emir4;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
if (cox < mx) { goto emir4;}
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = 3.333333333333333333 * 1.2 * (cox - mx - 0.1);
cony = 3.333333333333333333 * (coy - 0.5);
conz = 3.333333333333333333 * (coz - 0.5 - mz);
cot = t1;
tResult = t1; doThing = 3;
emir4:
if (vx == 0.0) {goto emi5;}
t1 = (mx - x0) / vx + 0.000001;
if (t1 + 0.00002 < rayFirst) {goto emi5;}
if (t1 > tResult) {goto emi5;}
y = y0 + vy * t1 - 0.5;
z = z0 + vz * t1 - mz - 0.5;
t2 = (y * y) + (z * z);
if (t2 > CUT_EMI_R_2) {goto emi5;}
t1 = t1 - 0.0001;
cox = x0 + vx * t1;
coy = y0 + vy * t1;
coz = z0 + vz * t1;
conx = -1.0;
cony = 0.0;
conz = 0.0;
cot = t1;
tResult = t1; doThing = 4;

					} else if (state == 3) {
//#include<emi_d_p.cpp>
sz = mz + 1.1 - z0;
sy = 0.5 - y0;
sx = (0.5 + mx) - x0;
k = (sy * vy) + (sx * vx) - (sz * vz * E_INV_16_CONE_H_2);
coneC = (sy * sy) + (sx * sx) - (sz * sz * E_INV_16_CONE_H_2);
econeA = (vy * vy) + (vx * vx) - (vz * vz * E_INV_16_CONE_H_2);
D = k * k - econeA * coneC;
einvConeA = 1 / econeA;
if (D <= 0.0) { goto emid1; }
D = sqrt(D);
if (econeA > 0.0) {
	t1 = (k - D) * einvConeA;
	t2 = (k + D) * einvConeA;
} else {
	t2 = (k - D) * einvConeA;
	t1 = (k + D) * einvConeA;
}
z = (mz + 0.1) - (z0 + (vz * t1));
if (z < -E_CONE_H + 0.1) {
	z = (mz + 0.1) - (z0 + (vz * t2));
	if (z < -E_CONE_H + 0.1) { goto emid1; }
	t1 = t2;
}
z = -z;
if (z < 0.1) { goto emid1; }
if (t1 < rayFirst) {goto emid1;}
tmp = E_CONE_YR_P - (E_CONE_YR_Y * z); // tmp == r
if (tmp <= 0.0) {goto emid1;}
tmp = 1 / tmp;
cot = t1;
t2 = t1 - 0.00001;
z = z0 + vz * t2;
y = y0 + vy * t2;
x = x0 + vx * t2;
nz = E_CONE_DY;
ny =  (y - 0.5) * tmp;
nx = (x - mx - 0.5) * tmp;
tResult = t1; doThing = 0; cox = x; coy = y; coz = z; conx= nx; cony = ny; conz = nz;
emid1:
if (vz == 0.0) {goto emid2;}
t1 = (mz + 1.0 - z0) / vz;
if (t1 < rayFirst) {goto emid2;}
if (t1 > tResult) {goto emid2;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > RAY_R_2) {goto emid2;}
cot = t1;
tResult = t1; doThing = 1;
emid2:
a = (vz * vz * 16.0) + (vy * vy) + (vx * vx);
z = mz + 0.75 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (16.0 * z * vz) + (y * vy) + (x * vx);
c = (16.0 * z * z) + (y * y) + (x * x) - 0.04;
D = k * k - a * c;
if (D < 0.0) {goto emid3;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emid3;}
if (t1 > tResult) {goto emid3;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 5.0 * 4.0 * (coz - mz - 0.75);
cony = 5.0 * (coy - 0.5);
conx = 5.0 * (cox - 0.5 - mx);
cot = t1;
tResult = t1; doThing = 2;
emid3:
a = (vz * vz * 1.44) + (vy * vy) + (vx * vx);
z = mz + 0.1 - z0;
y = 0.5 - y0;
x = mx + 0.5 - x0;
k = (1.44 * z * vz) + (y * vy) + (x * vx);
c = (1.44 * z * z) + (y * y) + (x * x) - 0.09;
D = k * k - a * c;
if (D < 0.0) {goto emid4;}
D = sqrt(D);
a = 1 / a;
t1 = (k - D) * a;
if (t1 < rayFirst) {t1 = (k + D) * a;}
if (t1 < rayFirst) {goto emid4;}
if (t1 > tResult) {goto emid4;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
if (coz < mz) { goto emid4;}
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = 3.333333333333333333 * 1.2 * (coz - mz - 0.1);
cony = 3.333333333333333333 * (coy - 0.5);
conx = 3.333333333333333333 * (cox - 0.5 - mx);
cot = t1;
tResult = t1; doThing = 3;
emid4:
if (vz == 0.0) {goto emi5;}
t1 = (mz - z0) / vz + 0.000001;
if (t1 + 0.00002 < rayFirst) {goto emi5;}
if (t1 > tResult) {goto emi5;}
y = y0 + vy * t1 - 0.5;
x = x0 + vx * t1 - mx - 0.5;
t2 = (y * y) + (x * x);
if (t2 > CUT_EMI_R_2) {goto emi5;}
t1 = t1 - 0.0001;
coz = z0 + vz * t1;
coy = y0 + vy * t1;
cox = x0 + vx * t1;
conz = -1.0;
cony = 0.0;
conx = 0.0;
cot = t1;
tResult = t1; doThing = 4;

					}
emi5:

					if (doThing == -1) {break;}
					if ((doThing == 0)||(doThing == 2)||(doThing == 3)||(doThing == 4)) {
rayLast = cot;
						x = cox; y = coy; z = coz; nx = conx; ny = cony; nz = conz;
						tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
						rx = tmp * nx - vx;
						ry = tmp * ny - vy;
						rz = tmp * nz - vz;

						mat = &(materials[MAT_DIF_CON]);
						for (int i = 0; i < 3; i++) {
							if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
//								intense = DIFFUSE(i) * colorFactor;
//								ADD_COLOR
								intense = SPECULAR(i);
								if (intense < 0.0) {
//									intense = intense * intense;
//									intense = intense * intense * colorFactor;
									intense = -2.0 * intense * intense * intense * colorFactor;
									ADD_COLOR
								}
							}
						}
//						clr[0] += 0.25; clr[1] += 0.25; clr[2] += 0.25;
						clr[0] += 0.15; clr[1] += 0.25; clr[2] += 0.35;
					} else if (doThing == 1) {
						clr[2] += 1.0;
					}
					RETURN_COLOR
				break;

				case MIRROR:
					doWhat = DO_NOTHING;
					sx = .5 + mx - x0;
					sy = CONE_H - y0;
					sz = .5 + mz - z0;
					k = (sx * vx) + (sz * vz) - (sy * vy * INV_64_CONE_H_2);
					coneC = (sx * sx) + (sz * sz) - (sy * sy * INV_64_CONE_H_2);
					D = k * k - coneA * coneC;
					if (D <= 0.0) { goto mirrorPlane; }
					D = sqrt(D);
					if (cAp) {
						t1 = (k - D) * invConeA;
						t2 = (k + D) * invConeA;
					} else {
						t2 = (k - D) * invConeA;
						t1 = (k + D) * invConeA;
					}
					y = y0 + vy * t1;
					if (y > CONE_H) {
						y = y0 + vy * t2;
						if (y > CONE_H) { goto mirrorPlane; }
						t1 = t2;
					}
					if (y < 0.0) { goto mirrorPlane; }
					if (t1 < rayFirst) {goto mirrorPlane;}
					tmp = CONE_YR_P - (CONE_YR_Y * y); // tmp == r
					if (tmp <= 0.0) {goto mirrorPlane;}
					tmp = 1 / tmp;
					cot = t1;
					t2 = t1 - 0.00001;
					x = x0 + vx * t2;
					nx = (x - mx - 0.5) * tmp;
					y = y0 + vy * t2;
					ny = CONE_DY;
					z = z0 + vz * t2;
					nz = (z - mz - 0.5) * tmp;
					tResult = t1; doWhat = DO_CONE; cox = x; coy = y; coz = z; conx= nx; cony = ny; conz = nz;

mirrorPlane:
					D = vx * nX[state] + vz * nZ[state];
					if (D == 0.0) {goto mirrorTorus;}
					if (precalcActive) {
						t1 = tgt->w / D;
					} else {
						t1 = ((mx + 0.5 - x0) * nX[state]) + ((mz + 0.5 - z0) * nZ[state]);
						t1 = t1 / D;
					}
					if (t1 > tResult) {goto mirrorTorus;}
					if (t1 < rayFirst) {goto mirrorTorus;}
					t2 = y0 + (vy * t1) - 0.5; // t2 == y
					t3 = ((x0 + (vx * t1) - mx - .5) * wX[state]) + ((z0 + (vz * t1) - mz - .5) * wZ[state]); // t3 == x
					tmp = ((t2 * t2) + (t3 * t3));
					if (tmp > MIRROR_WR_2) {goto mirrorTorus;}
					tResult = t1;
					doWhat = DO_MIRROR;

mirrorTorus:
					tShi = ((mx + 0.5 - x0) * vx) + ((0.5 - y0) * vy) + ((mz + 0.5 - z0) * vz);
					x1 = x0 + (vx * tShi); y1 = y0 + (vy * tShi); z1 = z0 + (vz * tShi);
					t1 = x1 - mx - 0.5;
					t2 = y1 - 0.5;
					t3 = z1 - mz - 0.5;
					dx = (wX[state] * t1) + (wZ[state] * t3);
					dy = (wZ[state] * t1) - (wX[state] * t3);
					dz = t2;
					t1 = (wX[state] * vx) + (wZ[state] * vz);
					t2 = (wZ[state] * vx) - (wX[state] * vz);
					a = 2.0 * ((t1 * dx) + (vy * dz)); // r: t
					b = a + (2.0 * t2 * dy); // l^1 : t
					c = (dx * dx) + (dz * dz); // r: 1
					f = c + (dy * dy) + 0.0195; // l^1: 1
					e = (t1 * t1) + (vy * vy); // r: t^2
					//t4 = 1.0;
					t3 = 2.0 * b;
					t2 = (2.0 * f) + (b * b) - (0.0784 * e);
					t1 = (2.0 * b * f) - (0.0784 * a);
					t0 = (f * f) - (0.0784 * c);
					if (!Solve(/*t4,*/t3,t2,t1,t0, tStart - tShi, tEnd - tShi, tmp)) {
						goto mirrorEnd;
					}
					tmp += tShi;
					if (tmp > tResult) {goto mirrorEnd;}
					tResult = tmp;
					doWhat = DO_TORUS;

mirrorEnd:
					if (doWhat == DO_CONE) {
rayLast = cot;
						x = cox; y = coy; z = coz; nx = conx; ny = cony; nz = conz;
						mat = &(materials[MAT_DIF_CON]);
						for (int i = 0; i < 3; i++) {
							if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
								intense = DIFFUSE(i) * colorFactor;
								ADD_COLOR
							}
						}
						clr[0] += 0.25; clr[1] += 0.25; clr[2] += 0.25;
						RETURN_COLOR
						break;
					} else if (doWhat == DO_NOTHING) {
						break;
					} else if (doWhat == DO_TORUS) {
rayLast = tResult;
					t1 = tResult;
// TODO: t2 = t1 - 0.000001;
					t2 = t1 - 0.000001;
					x = x0 + (vx * t2);
					y = y0 + (vy * t2);
					z = z0 + (vz * t2);
					t1 = x - mx - 0.5;
					t2 = y - 0.5;
					t3 = z - mz - 0.5;
					dx = (wX[state] * t1) + (wZ[state] * t3);
					dy = -(wZ[state] * t1) + (wX[state] * t3);
					dz = t2;
					f = 4.0 * ((dx * dx) + (dy * dy) + (dz * dz) - 0.0197);
					nx = dx * f;
					ny = dy * (f + 0.1568);
					nz = dz * f;
					e = 1.0 / sqrt((nx * nx) + (ny * ny) + (nz * nz));
					t1 = nx * e;
					t2 = ny * e;
					t3 = nz * e;
					nx = (wX[state] * t1) - (wZ[state] * t2);
					nz = (wZ[state] * t1) + (wX[state] * t2);
					ny = t3;
					tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
					rx = tmp * nx - vx;
					ry = tmp * ny - vy;
					rz = tmp * nz - vz;
					for (int i = 0; i < 3; i++) {
						if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
							mat = &(materials[MAT_DIF_MIR]);
							intense = DIFFUSE(i) * colorFactor;
							ADD_COLOR
							intense = SPECULAR(i);
							if (intense < 0.0) {
								intense = intense * intense;
								intense = intense * intense * colorFactor;
								mat = &(materials[MAT_SPE_MIR]);
								ADD_COLOR
							}
						}
					}
						clr[0] += 0.1; clr[1] += 0.35; clr[2] += 0.1;						
						RETURN_COLOR
						break;
					}
// DO_MIRROR
					t1 = tResult;
rayLast = t1;
					// TODO: cone first?
					APPLY_RAY
					ny = 0.0;
					nx = nX[state];
					nz = nZ[state];
					tmp = ((nx * vx) + (ny * vy) + (nz * vz));
					if (tmp > 0) { nx = -nx; ny = -ny; nz= -nz; tmp = -tmp;}
					tmp = vx * nx + vy * ny + vz * nz; tmp = tmp + tmp;
					rx = tmp * nx - vx;
					ry = tmp * ny - vy;
					rz = tmp * nz - vz;
					t2 = t1 - 0.000001;
					x = x0 + vx * t1;
					y = y0 + vy * t1;
					z = z0 + vz * t1;
					precalcActive = false;
					deepness++; if (deepness > 16) {
						clr[0] += (rand() & 0x3F) / 255.0;
						clr[1] += (rand() & 0x3F) / 255.0;
						clr[2] += (rand() & 0x3F) / 255.0;
						RETURN_COLOR
					}
					x0 = x;
					y0 = y;
					z0 = z;
					vx = -rx;
					vy = -ry;
					vz = -rz;
					CalcTimeBox(x0, y0, z0, vx, vy, vz, tStart, tEnd, boundsFloor);
					tStart += 0.000000001;
					tEnd   += 0.000000001;
					rayFirst = tStart;
					rayLast = tEnd;
					if (vx > 0.0) {tmx = 1.0 / vx; tNextX = (mx + 1.0 - x0) / vx; dmx = 1;}
					else if (vx < 0.0) {tmx = -1.0 / vx;tNextX = (mx - x0) / vx; dmx = -1;}
					else {tNextX = LAST_T;}
					if (vz > 0.0) {tmz = 1.0 / vz; tNextZ = (mz + 1.0 - z0) / vz; dmz = 1;}
					else if (vz < 0.0) {tmz = -1.0 / vz; tNextZ = (mz - z0) / vz; dmz = -1;}
					else { tNextZ = LAST_T;}
					coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
					invConeA = 1 / coneA;
					cAp = coneA > 0.0;
					goto noCycle;
					//break;
					
				break;
			}
			if (tNextX < tNextZ) {
				//rayFirst = tNextX;
				if (tEnd <= tNextX) {if (boundsFloor) {goto renderFloor;} else { RETURN_COLOR }}
				if (dmx > 0) { if (mx == 14) { RETURN_COLOR } mx++; idx++;}
				else {if (mx == 0) { RETURN_COLOR } mx--; idx--;}
				tNextX += tmx;
			} else {
				//rayFirst = tNextZ;
				if (tEnd <= tNextZ) {if (boundsFloor) {goto renderFloor;}else { RETURN_COLOR }}
				if (dmz > 0) {if (mz == 8) { RETURN_COLOR } mz++; idx+=15;}
				else {if (mz == 0) { RETURN_COLOR } mz--; idx-=15;}
				tNextZ += tmz;
			}
noCycle:
		}
renderFloor:
		x = x0 + (vx * tEnd);
		y = 0.0;
		z = z0 + (vz * tEnd);
		nx = 0.0;
		ny = 1.0;
		nz = 0.0;
		mat = &(materials[MAT_DIF_FLO]);
		for (int i = 0; i < 3; i++) {
			if (CanSee(i, x, y, z, lensed, colorFactor, lx, ly, lz)) {
				intense = DIFFUSE(i) * colorFactor;
				ADD_COLOR
			}
		}
		RETURN_COLOR
#ifdef PROTECT
	} catch (...) {
		clr[0] = rand() & 0xFF; clr[1] = rand() & 0xFF; clr[2] = rand() & 0xFF;
		return;
	}
#endif
	}

	void Perform(const int from, const int to) const {
		LPDWORD* lines = my_lines;
		number x0 = my_x0;
		number y0 = my_y0;
		number z0 = my_z0;
		number clrs0[3];
		number clrs1[3];
		number clrs2[3];
		number clrs3[3];
		for (int yy = from; yy != to; yy++) {
			int y = 127 - (yy / 2);
			if ((yy & 0x1) == 1) {
				y = 255 - y;
			}
			LPDWORD line = lines[y];
#ifndef ANTI_ALIASING
				for (int x = 0; x < 256; x++) {
					if (objectiveWindow[x + (256 * y)]) {
						CalcPoint(x, y, x0, y0, z0, clrs0);
						line[x] = TO_RGB(clrs0);
					}
				}
#else 
				for (int x = 0; x < 256; x++) {
					if (objectiveWindow[x + (256 * y)]) {
						CalcPoint(x * 4    , y, x0, y0, z0, clrs0);
						CalcPoint(x * 4 + 1, y, x0, y0, z0, clrs1);
						CalcPoint(x * 4 + 2, y, x0, y0, z0, clrs2);
						CalcPoint(x * 4 + 3, y, x0, y0, z0, clrs3);
						clrs0[0] = 0.25 * (clrs0[0] + clrs1[0] + clrs2[0] + clrs3[0]);
						clrs0[1] = 0.25 * (clrs0[1] + clrs1[1] + clrs2[1] + clrs3[1]);
						clrs0[2] = 0.25 * (clrs0[2] + clrs1[2] + clrs2[2] + clrs3[2]);
						line[x] = TO_RGB(clrs0);
					}
				}
#endif
		}
	}
public:
	void operator()(const blocked_range<int>& r) const {
		Perform(r.begin(), r.end());
	}
	CLASS_NAME (LPDWORD* lines, number x0, number y0, number z0) :
		my_lines(lines), my_x0(x0), my_y0(y0), my_z0(z0)
	{}
};
