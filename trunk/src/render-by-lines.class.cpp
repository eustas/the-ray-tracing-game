class CLASS_NAME {
	LPDWORD* const my_lines;
	number const my_x0;
	number const my_y0;
	number const my_z0;

	bool Ferrari(number A, number B,number C,number D,number E,number tStart, number tEnd, number &result) const {
		if (A == 0.0) {fwprintf_s(debugFile, L"A\n"); return false; /* todo: cordano */}
		number iA = 1.0 / A;
		number iA2 = iA * iA;
		number iA3 = iA * iA2;
		number iA4 = iA2 * iA2;
		number B2 = B * B;
		number B3 = B * B2;
		number B4 = B2* B2;
		number alpha = (C * iA) - (0.375 * B2 * iA2);
		number beta = (0.125 * B3 * iA3) - (0.5 * B * C * iA2) + (D * iA);
		number gamma = (0.0625 * C * B2 * iA3) - (0.01171875 * B4 * iA4) - (0.25 * B * D * iA2) + (E * iA);
		number D1, D2, tmp;
		result = LAST_T;
		number t;
		if (beta == 0.0) {
			fwprintf_s(debugFile, L"beta\n");
			D1 = (alpha * alpha) - (4.0 * gamma);
			if (D1 < 0.0) {fwprintf_s(debugFile, L"B\n");return false;}
			D1 = sqrt(D1);
			tmp = -B * iA;
			D2 = 0.5 * (-alpha - D1);
			if (D2 >= 0.0) {
				D2 = sqrt(D2);
				t = tmp - D2; if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
				t = tmp + D2; if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
			}
			D2 = 0.5 * (-alpha - D1);
			if (D2 >= 0.0) {
				D2 = sqrt(D2);
				t = tmp - D2; if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
				t = tmp + D2; if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
			}
			return (result < LAST_T);
		}
		number alpha2 = alpha * alpha;
		number P = (-0.083333333333333333333333333333333 * alpha2) - gamma;
		number Q = (0.33333333333333333333333333333333 * alpha * gamma)-(0.0092592592592592592592592592592593 * alpha2 * alpha) - (0.125 * beta * beta);
		D1 = (0.25 * Q * Q) + (0.037037037037037037037037037037037 * P * P * P);
		std::complex <number> R0 = (0.5 * Q);
		std::complex <number> R1 = (D1);
		std::complex <number> R = std::sqrt(R1) - R0;
		std::complex <number> U;
		std::complex <number> y;
		if (R == 0.0) {
			U = std::pow(-Q, 0.33333333333333333333333333333333);
			y = U - (0.83333333333333333333333333333333 * alpha);
		} else {
			U = std::pow(R, 0.33333333333333333333333333333333);
			y = U - (0.83333333333333333333333333333333 * alpha) - ((0.33333333333333333333333333333333 * P) / U);
		}
		number W = alpha + (2.0 * std::real(y));
		if (W <= 0.0) {fwprintf_s(debugFile, L"D %f| %f %f %f %f %f\n",W,A,B,C,D,E);return false;}
		W = sqrt(W);
		number iW = 1 / W;
		number D0 = -(3.0 * alpha) - (2.0 * std::real(y));
		number D4 = 2.0 * beta * iW;
		tmp = -0.25 * B * iA;
		D1 = D0 - D4;
		if (D1 > 0.0) {
			D1 = sqrt(D1);
			t = tmp + (0.5 * (W - D1)); if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
			t = tmp + (0.5 * (W + D1)); if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
		}
		D1 = D0 + D4;
		if (D1 > 0.0) {
			D1 = sqrt(D1);
			t = tmp + (0.5 * (-W - D1)); if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
			t = tmp + (0.5 * (-W + D1)); if ((t > tStart)&&(t < tEnd)&&(t < result)) {result = t;}
		}
		return (result < LAST_T);
	}

	bool CanSee(const int light, number x0, number y0, number z0, bool &lensed, number &factor, number &lx, number &ly, number &lz) const {
		factor = 1.0;
		lensed = false;
		number vx = lights[light].x;
		number vy = lights[light].y;
		number vz = lights[light].z;
		number wx,wy,wz;
		number sx,sy,sz, coneC,econeA,einvConeA,a,c,cox;
restart:
		number coneA = (vx * vx) + (vz * vz) - (vy * vy * INV_64_CONE_H_2);
		number invConeA = 1 / coneA;
		bool cAp = coneA > 0.0;
		number dx, dy, dz, w, k, D, t1, t2, tmp;
		number x,y,z, t3;
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


				case EMITTER:
					if (state == 0) {
#include"emi_l_c.cpp"
					}
					break;


				case PORTAL:
					dx = mx + 0.5 - x0;
					dy =      0.5 - y0;
					dz = mz + 0.5 - z0;
					w = (dx * dx) + (dy * dy) + (dz * dz) - (PORTAL_R * PORTAL_R);
					k = (dx * vx) + (dy * vy) + (dz * vz);
					D = (k * k) - w;
					if (D <= 0.0) { break; }
					t1 = k - sqrt(D);
					if (t1 < tStart) {t1 = k + sqrt(D); if (t1 > tStart) {return false;} break;}
					deepness++;
					if (deepness > 3) {return false;}
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
					if (tmp < MIRROR_R_2) {
						if (tmp > MIRROR_WR_2) {
							return false;
						}
						// TODO: mirroring light?
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
		number D, t0, t1, t2, t3, t4, tmp, temp;
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
		bool doCone;
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
					if (D <= 0.0) { goto portalTorus; }
					t1 = k - sqrt(D);
					if (t1 < tStart) {
						t1 = k + sqrt(D);
						if (t1 > tStart) {
							// is that right?
							RETURN_COLOR
						}
						// may be break?
						goto portalTorus;
					}
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
					break;
portalTorus:
					dx = x0 - mx - 0.5;
					dy = y0;
					dz = z0 - mz - 0.5;
					a = 2.0 * ((vx * dx) + (vz * dz)); // r: t
					b = a + (2.0 * vy * dy); // l^1 : t
					c = (dx * dx) + (dz * dz); // r: 1
					f = c + (dy * dy) + 0.15; // l^1: 1
					e = (vx * vx) + (vz * vz); // r: t^2
					t4 = 1.0;
					t3 = 2.0 * b;
					t2 = (2.0 * f) + (b * b) - (0.64 * e);
					t1 = (2.0 * b * f) - (0.64 * a);
					t0 = (f * f) - (0.64 * c);
					if (Ferrari(t4,t3,t2,t1,t0, tStart, tEnd, tmp)) {
						//fwprintf_s(debugFile, L"t=%f\n", t0 + tmp * ( t1 + tmp * (t2 + tmp * (t3 + tmp * t4))));
						clr[0] += (rand() & 0x3F) / 255.0;
						clr[1] += (rand() & 0x3F) / 255.0;
						clr[2] += (rand() & 0x3F) / 255.0;
						RETURN_COLOR
					}
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

				case EMITTER:
					tResult = LAST_T;
					doThing = -1;
					if (state == 0) {
#include<emi_l_p.cpp>
					} else if (state == 1) {
#include<emi_u_p.cpp>
					} else if (state == 2) {
#include<emi_r_p.cpp>
					} else if (state == 3) {
#include<emi_d_p.cpp>
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
					doCone = false;
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
					tResult = t1; doCone = true; cox = x; coy = y; coz = z; conx= nx; cony = ny; conz = nz;

mirrorPlane:
					D = vx * nX[state] + vz * nZ[state];
					if (D == 0.0) {goto mirrorEnd;}
					if (precalcActive) {
						t1 = tgt->w / D;
					} else {
						t1 = ((mx + 0.5 - x0) * nX[state]) + ((mz + 0.5 - z0) * nZ[state]);
						t1 = t1 / D;
					}
					if (t1 > tResult) {goto mirrorEnd;}
					if (t1 < rayFirst) {goto mirrorEnd;}
					t2 = y0 + (vy * t1) - 0.5; // t2 == y
					t3 = ((x0 + (vx * t1) - mx - .5) * wX[state]) + ((z0 + (vz * t1) - mz - .5) * wZ[state]); // t3 == x
					tmp = ((t2 * t2) + (t3 * t3));
					if (tmp < MIRROR_R_2) {
rayLast = t1;
						if (tmp > MIRROR_WR_2) {
							clr[1] += 1.0;
							RETURN_COLOR
						} else {
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
						}
					}
mirrorEnd:
					if (doCone) {
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
					}
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
