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
t1 = (mz + 1.0 - z0) / vz;
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
