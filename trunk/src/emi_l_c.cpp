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
