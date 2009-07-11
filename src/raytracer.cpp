#define SCREEN_Z 2.0
#define LAST_T 200000.0

#define QUANT_R 0.2
#define QUANT_R_2 0.04
#define INV_QUANT_R 5.0
#define CONE_H 0.35
#define MIRROR_R_2 0.0225
#define MIRROR_WR_2 0.0169
#define INV_64_CONE_H_2 0.12755102040816326530612244897959
#define CONE_DY 0.24253562503633297351890646211612
#define RAY_R_2 0.0006

#define GLASS_K     1.1764705882352941176470588235294
#define INV_GLASS_K 0.85
#define FOCUS       0.6666666666666666666666666666666
#define FOCUS_FACTOR 0.075

//#define AA_X 0.16666666666666666666666666666666
//#define AA_Y 0.28867513459481288225457439025098
#define AA_X 0.25
#define AA_Y 0.25

#define PORTAL_R 0.3

#define DIFFUSE(light)	( lensed ? \
				((nx * lx) + (ny * ly) + (nz * lz)) : \
				((nx * lights[light].x) + (ny * lights[light].y) + (nz * lights[light].z)) \
			)

#define SPECULAR(light) ((rx * lights[light].x) + (ry * lights[light].y) + (rz * lights[light].z))
#define ADD_COLOR if (intense > 0.0) { clr[0] += lightsClr[i].x * mat->x * intense; clr[1] += lightsClr[i].y * mat->y * intense; clr[2] += lightsClr[i].z * mat->z * intense;}

#define APPLY_RAY \
	if (hasRay) { \
		for (int ridx = 0; ridx< rayIdx; ridx++) { \
			if (rayLast < rayT2[ridx]) {rayT2[ridx] = rayLast;} \
			if (rayFirst > rayT1[ridx]) {rayT1[ridx] = rayFirst;} \
			if (rayT2[ridx] > rayT1[ridx]) { \
				clr[2] += 6.0 * (rayT2[ridx] - rayT1[ridx]); \
			} \
		} \
	} \
	hasRay = false;

#define RETURN_COLOR \
	APPLY_RAY \
	clr[0] *= 255.99; clr[1] *= 255.99; clr[2] *= 255.99; \
	if (clr[0] > 255.99) {clr[0] = 255.99;} \
	if (clr[1] > 255.99) {clr[1] = 255.99;} \
	if (clr[2] > 255.99) {clr[2] = 255.99;} \
	return;

#define RETURN_BLACK \
	clr[0] = 0.0; \
	clr[1] = 0.0; \
	clr[2] = 0.0; \
	return;

#define TO_RGB(colors) (((DWORD)colors[0]) | (((DWORD)colors[1]) << 8) | (((DWORD)colors[2]) << 16))

#define SIN_LA 0.5
#define COS_LA 0.86602540378443864676372317075294

// 1 / (2 * sqrt(17))
#define CONE_YR_P 0.12126781251816648675945323105806

// 1 / (2 * h * sqrt(17))
#define CONE_YR_Y 0.34647946433761853359843780302303

VEC4D lights[3];
VEC4D lightsClr[3];
VEC4D materials[5];

#define MAT_DIF_DEF 0
#define MAT_DIF_BLO 1
#define MAT_DIF_FLO 2
#define MAT_DIF_PNT 3
#define MAT_SPE_PNT 4

#define MAT_DIF_CON 2

int bFlags[128];
#define VLD 0x000001
#define VLU 0x000002
#define VL  0x000003
#define VMD 0x000010
#define VMU 0x000020
#define VM  0x000030
#define VRD 0x000100
#define VRU 0x000200
#define VR  0x000300
#define HDL 0x001000
#define HDR 0x002000
#define HD  0x003000
#define HML 0x010000
#define HMR 0x020000
#define HM  0x030000
#define HUL 0x100000
#define HUR 0x200000
#define HU  0x300000
#define AVH 0x333333

//#define PROTECT
RAY_PART* rayThings[15 * 9];

namespace rt {
	int phase;
	CELL* cells;
	char* qcells;
	RAY_PART* rayParts;
	int rayLen;
	bool glassPoint;
	bool antiAliasing;
};

double nX[16], nZ[16], wX[16], wZ[16], poX[16], poZ[16];
double p[16], roX[16], roZ[16], doX[16], doZ[16];

HANDLE rayTracerThread;
DWORD  rayTracerThreadId;

VEC4D* basicRays;
VEC4D* modelRays;
VEC4D* basicRaysA;
VEC4D* modelRaysA;

VEC4D precalc[9 * 15];

number lastModelAlpha	= 0.0;
number lastModelBeta	= 0.0;

char microInfo[32];
char microFps[32];
char microFrm[32];

#define CLASS_NAME RenderByLinesSimple
#include"render-by-lines.class.cpp"
#undef CLASS_NAME

#define CLASS_NAME RenderByLinesAa
#define ANTI_ALIASING
#include"render-by-lines.class.Aa.cpp"
#undef ANTI_ALIASING
#undef CLASS_NAME

#define CLASS_NAME RenderByLinesGp
#define GLASS_POINT
#include"render-by-lines.class.Gp.cpp"
#undef GLASS_POINT
#undef CLASS_NAME

#define CLASS_NAME RenderByLinesAaGp
#define GLASS_POINT
#define ANTI_ALIASING
#include"render-by-lines.class.AaGp.cpp"
#undef ANTI_ALIASING
#undef GLASS_POINT
#undef CLASS_NAME

double stohasticData[200];

void ResetResults() {
	if (stohasticIdx < 0) {return;}
	stohasticIdx = 0;
}

void AddResult(double sec) {
	if (stohasticIdx < 0) {return;}
	stohasticData[stohasticIdx] = sec;
	stohasticIdx++;

	if (stohasticIdx != 200) {return;}
	stohasticIdx = 0;

	double tmp;
	for (int i = 0; i < 200;i++) {
		for (int j = 0; j < 199; j++) {
			if (stohasticData[j] > stohasticData[j+1]) {
				tmp = stohasticData[j];
				stohasticData[j] = stohasticData[j+1];
				stohasticData[j+1] = tmp;
			}
		}
	}
	tmp = 0.0;
	for (int k = 70; k < 130; k++) {
		tmp += stohasticData[k];
	}
	fwprintf_s(debugFile, L"stohastic mode=%d gp=%d aa=%d td=%f\n", demoMode, rt::glassPoint, rt::antiAliasing, tmp);
}

void DumpToCanvas() {
	for (int yy = 0; yy < 256; yy++) {
		LPDWORD line = microLines[yy];
		int sh2 = 21 + 300 * (300 - 23 - yy - 1);
		for (int xx = 0; xx < 256; xx++) {
			if (objectiveWindow[xx + 256 * yy]) {
				((LPDWORD)mbits)[sh2 + xx] = line[xx];
			}
		}
	}
}

void RecalculateModelRays(number newModelAlpha, number newModelBeta, bool newAntiAliasing) {
	lastModelAlpha = newModelAlpha;
	lastModelBeta = newModelBeta;
	rt::antiAliasing = newAntiAliasing;

	number sina = sin(lastModelAlpha);
	number sinb = sin(lastModelBeta);
	number cosa = cos(lastModelAlpha);
	number cosb = cos(lastModelBeta);
	//VEC4D zRay = {-cosa * cosb, -sinb, -sina * cosb, 1.0f};
	//VEC4D yRay = {-cosa * sinb,  cosb, -sina * sinb, 1.0f};
	//VEC4D xRay = { sina       ,  0.0f, -cosa       , 1.0f};
	VEC4D xRayT = { sina       ,  -cosa * sinb, -cosa * cosb, 0.0f};
	VEC4D yRayT = {        0.0f,  cosb        , -sinb       , 0.0f};
	VEC4D zRayT = {-cosa       , -sinb * sina , -sina * cosb, 0.0f};
	//VEC4D wRayT = {        0.0f,  0.0f        ,         0.0f, 1.0f};
	int idx = 0, idx2 = 0;
	if (!newAntiAliasing) {for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			if (objectiveWindow[idx]) {
				modelRays[idx].x = (xRayT.x * basicRays[idx].x) + (xRayT.y * basicRays[idx].y) + (xRayT.z * basicRays[idx].z);
				modelRays[idx].y = (yRayT.x * basicRays[idx].x) + (yRayT.y * basicRays[idx].y) + (yRayT.z * basicRays[idx].z);
				modelRays[idx].z = (zRayT.x * basicRays[idx].x) + (zRayT.y * basicRays[idx].y) + (zRayT.z * basicRays[idx].z);
				//modelRays[idx].w = 1.0f;
			}
			idx++;
		}
	}} else {for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			if (objectiveWindow[idx]) {
				idx2 = idx * 4;
				modelRaysA[idx2].x = (xRayT.x * basicRaysA[idx2].x) + (xRayT.y * basicRaysA[idx2].y) + (xRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].y = (yRayT.x * basicRaysA[idx2].x) + (yRayT.y * basicRaysA[idx2].y) + (yRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].z = (zRayT.x * basicRaysA[idx2].x) + (zRayT.y * basicRaysA[idx2].y) + (zRayT.z * basicRaysA[idx2].z);
				//modelRays[idx2].w = 1.0f;
				idx2++;

				modelRaysA[idx2].x = (xRayT.x * basicRaysA[idx2].x) + (xRayT.y * basicRaysA[idx2].y) + (xRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].y = (yRayT.x * basicRaysA[idx2].x) + (yRayT.y * basicRaysA[idx2].y) + (yRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].z = (zRayT.x * basicRaysA[idx2].x) + (zRayT.y * basicRaysA[idx2].y) + (zRayT.z * basicRaysA[idx2].z);
				//modelRays[idx2].w = 1.0f;
				idx2++;

				modelRaysA[idx2].x = (xRayT.x * basicRaysA[idx2].x) + (xRayT.y * basicRaysA[idx2].y) + (xRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].y = (yRayT.x * basicRaysA[idx2].x) + (yRayT.y * basicRaysA[idx2].y) + (yRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].z = (zRayT.x * basicRaysA[idx2].x) + (zRayT.y * basicRaysA[idx2].y) + (zRayT.z * basicRaysA[idx2].z);
				//modelRays[idx2].w = 1.0f;
				idx2++;

				modelRaysA[idx2].x = (xRayT.x * basicRaysA[idx2].x) + (xRayT.y * basicRaysA[idx2].y) + (xRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].y = (yRayT.x * basicRaysA[idx2].x) + (yRayT.y * basicRaysA[idx2].y) + (yRayT.z * basicRaysA[idx2].z);
				modelRaysA[idx2].z = (zRayT.x * basicRaysA[idx2].x) + (zRayT.y * basicRaysA[idx2].y) + (zRayT.z * basicRaysA[idx2].z);
				//modelRays[idx2].w = 1.0f;
			}
			idx++;
		}
	}}
	
}

void Precalculate(number x0, number y0, number z0) {
	VEC4D* tgt;
	number dx,dy,dz;
	char nstate;
	for (int y = 0; y < 9; y++) {
		for (int x = 0; x < 15; x++) {
			int idx = x + (y * 15);
			char type = rt::cells[idx].type;
			char state = rt::cells[idx].state;
			tgt = &(precalc[idx]);
			switch (type) {
				case QUANT:
				case STAR:
					dx = x + 0.5 - x0; tgt->x = dx;
					dy = 0.5 - y0; tgt->y = dy;
					dz = y + 0.5 - z0; tgt->z = dz;
					tgt->w = (dx * dx) + (dy * dy) + (dz * dz) - (QUANT_R * QUANT_R);
				break;

				case PORTAL:
					dx = x + 0.5 - x0; tgt->x = dx;
					dy = 0.5 - y0; tgt->y = dy;
					dz = y + 0.5 - z0; tgt->z = dz;
					tgt->w = (dx * dx) + (dy * dy) + (dz * dz) - (PORTAL_R * PORTAL_R);
				break;

				case MIRROR:
				case ROTOR:
					nstate = (state + ((type == ROTOR) ? (rt::phase) : 0)) & 0xF;
					rt::cells[idx].type = MIRROR;
					rt::cells[idx].state = nstate;
					tgt->w = ((x + 0.5 - x0) * nX[nstate]) + ((y + 0.5 - z0) * nZ[nstate]);
				break;

				case ROLAR:
				case MROLAR:
					nstate = (state + rt::phase) & 0xF;
					if (type == ROLAR) {rt::cells[idx].type = POLAR;} else {rt::cells[idx].type = MPOLAR;}
					rt::cells[idx].state = nstate;
				break;
			}
		}
	}

	for (int i = 0; i < 16; i++) {
		p[i] = x0 * nX[i] + z0 * nZ[i];
	}
}

void DoRayTrace() {
	number newModelAlpha = rayAlpha;
	number newModelBeta = rayBeta;
	number newRayDist = rayDist;
	bool newAntiAliasing = antiAliasing;
	bool newGlassPoint = glassPoint;

	if ((newModelAlpha != lastModelAlpha) || (newModelBeta != lastModelBeta) || (newAntiAliasing != rt::antiAliasing)) {
		RecalculateModelRays(newModelAlpha, newModelBeta, newAntiAliasing);
		ResetResults();
	}
	if (newGlassPoint != rt::glassPoint) {
		ResetResults();
	}
	number sina = sin(lastModelAlpha);
	number sinb = sin(lastModelBeta);
	number cosa = cos(lastModelAlpha);
	number cosb = cos(lastModelBeta);

	number x = controlX + .5f;
	number y = controlY + .5f;
	if (microX >= 0.0) {
		x = microX;
		y = microY;
	}

	lightAlpha += lightAlphaDelta;

	for (int i = 0; i < 3; i++) {
		VEC4D* light = &(lights[i]);
		number ang = PI * (1 + 8 * i) / 12.0;
		light->x = cos(ang + lightAlpha) * COS_LA;
		light->y = SIN_LA;
		light->z = sin(ang + lightAlpha) * COS_LA;
	}

	rt::glassPoint = newGlassPoint;

	number x0 = x   + newRayDist * cosa * cosb;
	number y0 = .5f + newRayDist * sinb;
	number z0 = y   + newRayDist * sina * cosb;

	memcpy(rt::cells, levelData.cells, 15 * 9 * sizeof(CELL));
	memcpy(rt::qcells, levelData.qcells, 32 * 20 * sizeof(char));
	rt::phase = phase;
	Precalculate(x0, y0, z0);

	while (castingRay.compare_and_swap(true,false) != false) {Sleep(0);}
	rt::rayLen = rayLen;
	memcpy_s(rt::rayParts, (rt::rayLen) * sizeof(RAY_PART), rayParts, (rt::rayLen) * sizeof(RAY_PART));
	while (castingRay.compare_and_swap(false,true) != true) {Sleep(0);}

	BuildRay();

	tick_count t0 = tick_count::now();
	if (!rt::antiAliasing) {
		if (!rt::glassPoint) {
			parallel_for(blocked_range<int>(0, 256, 1), RenderByLinesSimple(microLines, x0, y0, z0));
		} else {
			parallel_for(blocked_range<int>(0, 256, 1), RenderByLinesGp(microLines, x0, y0, z0));
		}
	} else {
		if (!rt::glassPoint) {
			parallel_for(blocked_range<int>(0, 256, 1), RenderByLinesAa(microLines, x0, y0, z0));
		} else {
			parallel_for(blocked_range<int>(0, 256, 1), RenderByLinesAaGp(microLines, x0, y0, z0));
		}
	}
	tick_count t1 = tick_count::now();
	double sec = (t1-t0).seconds();
	AddResult(sec);
}

void BuildRay() {
	memset(rayThings, 0, 15 * 9 * sizeof(RAY_PART*));
	number mx, mz,mxd,mzd, dx,dz;
	int xd,zd, dir;
	char type, state;
	for (int i = 0; i < rt::rayLen; i++) {
		RAY_PART* part = &(rt::rayParts[i]);
		dir = part->dir;
		mx = (part->sx - 2.0) / 4.0;
		mz = (part->sy - 2.0) / 4.0;
		mxd = mx + (deltaX[dir] * 0.01);
		mzd = mz + (deltaY[dir] * 0.01);
		if ((mxd < 0.0) || (mxd > 15.0) || (mzd < 0.0) || (mzd > 9.0)) {continue;}
		xd = (int)mxd;
		zd = (int)mzd;
		part->next = rayThings[xd + 15 * zd];
		part->d = (roX[dir] * mx) + (roZ[dir] * mz);

		type = rt::cells[xd + 15 * zd].type;
		state = rt::cells[xd + 15 * zd].state;
		dx = doX[dir]; dz = doZ[dir];
		if (type == MIRROR) {
			if (((part->sx & 0x3) == 0) && ((part->sy & 0x3) == 0)) {
				dx = nX[state]; dz = nZ[state];
				if (((deltaX[dir] * dx) + (deltaY[dir] * dz)) < 0.0) {dx = -dx; dz = -dz;}
			}
		}

		part->sta = (dx * mx) + (dz * mz);
		part->staX = dx; part->staZ = dz;


		mx += deltaX[dir] * 0.25;
		mz += deltaY[dir] * 0.25;
		dx = -doX[dir]; dz = -doZ[dir];
		if (type == MIRROR) {
			if ((((part->sx + deltaX[dir]) & 0x3) == 0) && (((part->sy + deltaY[dir]) & 0x3) == 0)) {
				dx = nX[state]; dz = nZ[state];
				if (((deltaX[dir] * dx) + (deltaY[dir] * dz)) > 0.0) {dx = -dx; dz = -dz;}
			}
		}
		part->end = (dx * mx) + (dz * mz);
		part->endX = dx; part->endZ = dz;

		rayThings[xd + 15 * zd] = part;
	}
}

DWORD WINAPI RayTracerThreadFunction(LPVOID) {
	task_scheduler_init tbb_init;

	bool memOnBattery = false;

	tick_count s_tick = tick_count::now();
	int frameCount = 0;

	while (true) {
		bool newOnBattery = onBattery;
		if (!memOnBattery && newOnBattery) {
			DrawUnplugged();
			SetDIBitsToDevice(mdc, 0, 0, 300, 300, 0, 0, 0, 300, mbits, &mbmi, DIB_RGB_COLORS);
		}
		memOnBattery = newOnBattery;
		if (memOnBattery) {
			Sleep(250);
			s_tick = tick_count::now();
			frameCount = 0;
			microFps[0] = 0;
			continue;
		}

		if (!IsWindowVisible(mwnd)) {
			Sleep(250);
			s_tick = tick_count::now();
			frameCount = 0;
			microFps[0] = 0;
			continue;
		}

		tick_count t0 = tick_count::now();
		{
			DoRayTrace();
			frameCount++;
			DumpToCanvas();
		}
		tick_count t1 = tick_count::now();
		double sec = (t1-t0).seconds();
		if (sec <= 0.099999) {
			long mks = (int)((t1-t0).seconds() * 1000000.0);
			long msR = (mks % 1000);
			long msL = (mks / 1000);
			sprintf_s(microInfo, 32, "dt=%d.%03dms", (int)msL, (int)msR);
		} else {
			long ms = (int)((t1-t0).seconds() * 1000.0);
			long sR = (ms % 1000);
			long sL = (ms / 1000);
			sprintf_s(microInfo, 32, "dt=%d.%03ds", (int)sL, (int)sR);
		}

		sec = (t1 - s_tick).seconds();
		if (sec > fpsGap) {
			s_tick = t1;
			sec = frameCount / sec;
			long fps = (sec * 1000);
			long l = fps / 1000;
			long r = fps % 1000;
			sprintf_s(microFps, 32, "FPS=%d.%03d", (int)l, (int)r);
			frameCount = 0;
		}

		DrawMicroText(microInfo, 1, 300 - 12 - 1);
		DrawMicroText(microFps, 1, 1);

		if (stohasticIdx >= 0) {
			sprintf_s(microFrm, 32, "       %03d", stohasticIdx);
			DrawMicroText(microFrm, 300 - 88 - 1, 300 - 12 - 1);
		}

		SetDIBitsToDevice(mdc, 0, 0, 300, 300, 0, 0, 0, 300, mbits, &mbmi, DIB_RGB_COLORS);
	}
	mallocThreadShutdownNotification(NULL);
	return 0;
}

void InitRayTracer() {
	rt::cells = new CELL[15 * 9];
	rt::qcells = new char[32 * 20];
	basicRays = new VEC4D[256 * 256];
	modelRays = new VEC4D[256 * 256];
	basicRaysA = new VEC4D[256 * 1024];
	modelRaysA = new VEC4D[256 * 1024];
	rt::rayParts = new RAY_PART[10240];
	int idx = 0, idx2 = 0;
	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			if (objectiveWindow[idx]) {
				number x0 = (x - 127.5f) / 127.5f;
				number y0 = (127.5f - y) / 127.5f;
				number z0 = SCREEN_Z;
				number l = sqrt( (x0 * x0) + (y0 * y0) + (z0 * z0) );
				basicRays[idx].x = x0 / l;
				basicRays[idx].y = y0 / l;
				basicRays[idx].z = z0 / l;
				basicRays[idx].w = 1.0f;

				x0 = (x + AA_X - 127.5f) / 127.5f;
				y0 = (127.5f - y - AA_Y) / 127.5f;
				z0 = SCREEN_Z;
				l = sqrt( (x0 * x0) + (y0 * y0) + (z0 * z0) );
				basicRaysA[idx2].x = x0 / l;
				basicRaysA[idx2].y = y0 / l;
				basicRaysA[idx2].z = z0 / l;
				basicRaysA[idx2].w = 1.0f;

				x0 = (x - AA_Y - 127.5f) / 127.5f;
				y0 = (127.5f - y - AA_X) / 127.5f;
				z0 = SCREEN_Z;
				l = sqrt( (x0 * x0) + (y0 * y0) + (z0 * z0) );
				basicRaysA[idx2 + 1].x = x0 / l;
				basicRaysA[idx2 + 1].y = y0 / l;
				basicRaysA[idx2 + 1].z = z0 / l;
				basicRaysA[idx2 + 1].w = 1.0f;

				x0 = (x - AA_X - 127.5f) / 127.5f;
				y0 = (127.5f - y + AA_Y) / 127.5f;
				z0 = SCREEN_Z;
				l = sqrt( (x0 * x0) + (y0 * y0) + (z0 * z0) );
				basicRaysA[idx2 + 2].x = x0 / l;
				basicRaysA[idx2 + 2].y = y0 / l;
				basicRaysA[idx2 + 2].z = z0 / l;
				basicRaysA[idx2 + 2].w = 1.0f;

				x0 = (x + AA_Y - 127.5f) / 127.5f;
				y0 = (127.5f - y + AA_X) / 127.5f;
				z0 = SCREEN_Z;
				l = sqrt( (x0 * x0) + (y0 * y0) + (z0 * z0) );
				basicRaysA[idx2 + 3].x = x0 / l;
				basicRaysA[idx2 + 3].y = y0 / l;
				basicRaysA[idx2 + 3].z = z0 / l;
				basicRaysA[idx2 + 3].w = 1.0f;
			}
			idx++;
			idx2 += 4;
		}
	}
	for (int i = 0; i < 16; i++) {
		double nx = - sin((-PI * i) / 16.0);
		double nz = cos((-PI * i) / 16.0);
		double n = sqrt((nx * nx) + (nz * nz));
		nx /= n;
		nz /= n;
		nX[i] = nx;
		nZ[i] = nz;
	}

	for (int i = 0; i < 16; i++) {
		double wx = cos((-PI * i) / 16.0);
		double wz = sin((-PI * i) / 16.0);
		double w = sqrt((wx * wx) + (wz * wz));
		wx /= w;
		wz /= w;
		wX[i] = wx;
		wZ[i] = wz;
	}

	for (int i = 0; i < 16; i++) {
		double pox = deltaX[(16 - i) & 0xF];
		double poz = deltaY[(16 - i) & 0xF];
		double po = sqrt((pox * pox) + (poz * poz));
		pox /= po;
		poz /= po;
		poX[i] = pox;
		poZ[i] = poz;
	}

	for (int i = 0; i < 16; i++) {
		double rox = -deltaY[i];
		double roz = deltaX[i];
		double ro = sqrt((rox * rox) + (roz * roz));
		rox /= ro;
		roz /= ro;
		roX[i] = rox;
		roZ[i] = roz;
	}

	for (int i = 0; i < 16; i++) {
		double dox = deltaX[i];
		double doz = deltaY[i];
		double doo = sqrt((dox * dox) + (doz * doz));
		dox /= doo;
		doz /= doo;
		doX[i] = dox;
		doZ[i] = doz;
	}

//	lightsClr[0].x = 1.0;lightsClr[0].y = 0.0;lightsClr[0].z = 0.0;
//	lightsClr[1].x = 0.0;lightsClr[1].y = 1.0;lightsClr[1].z = 0.0;
//	lightsClr[2].x = 0.0;lightsClr[2].y = 0.0;lightsClr[2].z = 1.0;
	lightsClr[0].x = .666;lightsClr[0].y = .333;lightsClr[0].z = .333;
	lightsClr[1].x = .333;lightsClr[1].y = .666;lightsClr[1].z = .333;
	lightsClr[2].x = .333;lightsClr[2].y = .333;lightsClr[2].z = .666;

	VEC4D* mat;
	mat = &(materials[MAT_DIF_DEF]); mat->x = 0.50; mat->y = 0.75; mat->z = 1.00;
	mat = &(materials[MAT_DIF_BLO]); mat->x = 1.00; mat->y = 0.75; mat->z = 0.50;
	mat = &(materials[MAT_DIF_FLO]); mat->x = 0.75; mat->y = 0.75; mat->z = 0.75;
	mat = &(materials[MAT_DIF_PNT]); mat->x = 0.75; mat->y = 0.75; mat->z = 0.75;
	//mat = &(materials[MAT_SPE_PNT]); mat->x = 0.50; mat->y = 0.50; mat->z = 0.50;
	mat = &(materials[MAT_SPE_PNT]); mat->x = 0.75; mat->y = 0.75; mat->z = 0.75;

	memset(bFlags, 0, 128 * sizeof(int));
	bFlags[FULL      ] = VL  | VR  | HD  | HU ;
	bFlags[LEFT      ] = VL  | VM  | HDL | HUL;
	bFlags[RIGHT     ] = VR  | VM  | HDR | HUR;
	bFlags[UP        ] = VLU | VRU | HU  | HM ;
	bFlags[DOWN      ] = VLD | VRD | HM  | HD ;
	bFlags[LEFT_UP   ] = VLU | VMU | HML | HUL;
	bFlags[LEFT_DOWN ] = VLD | VMD | HML | HDL;
	bFlags[RIGHT_UP  ] = VRU | VMU | HMR | HUR;
	bFlags[RIGHT_DOWN] = VRD | VMD | HMR | HDR;

	bFlags[EMPTY_LEFT_UP   ] = AVH ^ (VLU | HUL);
	bFlags[EMPTY_LEFT_DOWN ] = AVH ^ (VLD | HDL);
	bFlags[EMPTY_RIGHT_UP  ] = AVH ^ (VRU | HUR);
	bFlags[EMPTY_RIGHT_DOWN] = AVH ^ (VRD | HDR);

	RecalculateModelRays(0.0f, 0.0f, false);
	rayTracerThread = CreateThread(NULL, 0, RayTracerThreadFunction, NULL, 0, &rayTracerThreadId);
}
