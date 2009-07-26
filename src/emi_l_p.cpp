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
t1 = (mx + 1.0 - x0) / vx;
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
