bool hasAboutMusic;
DWORD palette[256];

void AboutClick() {
	ShowMenu();
}

namespace aboutNS {
#define dq 0.40
#define sq256 15.999999999999

double xi[512], fxi[512];
double yi[384], fyi[384];
double q;

void DoRenderAbout() {
	for (int y = 0; y < 384; y++) {
		double v = (y / 192.0) - 1.0;
		yi[y] = sq256 * v;
		fyi[y] = sq256 * cos(q * v);
	}

	for (int x = 0; x < 512; x++) {
		double v = (x / 256.0) - 1.0;
		xi[x] = sq256 * v;
		fxi[x] = sq256 * sin(q * v);
	}

	int t = 0;
	for (int y = 0; y < 384; y++) {
		int sh1 = 512 * (384 - 1 - y);
		double cyi = yi[y], cfyi = fyi[y];
		for (int x = 0; x < 512; x++) {
			DWORD clr = about[x + sh1];
			if (clr != 0) {
				((LPDWORD)bits)[t] = clr;
			} else {
				((LPDWORD)bits)[t] = palette[(int)((fxi[x] - cyi) * (cfyi - xi[x])) & 0xFF];
			}
			t++;
		}
	}

	q += dq;
}

#undef mq
#undef sq256
} // namespace about

void SetupAbout() {
	int k = 0;
	for (int i = 0; i < 64; i++) {
		DWORD clr = (DWORD)((255.99 * i) / 64.0);
		palette[k] = clr << 16;
		k++;
	}
	for (int i = 0; i < 64; i++) {
		DWORD clr = (DWORD)((255.99 * i) / 64.0);
		palette[k] = (clr << 8) | 0xFF0000;
		k++;
	}
	for (int i = 0; i < 64; i++) {
		DWORD clr = 255 - (DWORD)((255.99 * i) / 64.0);
		palette[k] = (clr << 16) | 0xFF00;
		k++;
	}
	for (int i = 0; i < 64; i++) {
		DWORD clr = 255 - (DWORD)((255.99 * i) / 64.0);
		palette[k] = clr << 8;
		k++;
	}
}

void StartAbout() {
	aboutNS::q = -25.0;
}

void PlayAbout() {
	if (hasAboutMusic) {StartSound(FILE_ABOUT_MUSIC);}
}

void RenderAbout() {
	aboutNS::DoRenderAbout();
}
