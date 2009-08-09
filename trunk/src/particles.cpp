void ClearParticles() {
	for (int i = 0; i < N_PARTICLES; i++) {
		particles[i].fuel = -1;
	}
}

void PrepareBlowParticles() {
	ClearParticles();
	int x0 = ((5 + triesLeft - 1) * 32 + 8) * 2;
	int y0 = (10 * 32 + 8) * 2;
	for (int i = 0; i < N_PARTICLES; i++) {
		int dy = rand() % 31 - 15;
		int dx = rand() % 31 - 15;
		double r = sqrt((dx * dx) + (dy * dy) + 0.0);
		if (r > 15.0) { continue; }
		int clr = (int)(r / 5.0); if (clr < 0) {clr = 0;} if (clr > 2) {clr = 2;}
		int x = x0 + (rand() & 0x1F);
		int y = y0 + (rand() & 0x1F);
		SetupParticle(i, 255, 4, x, y, dx, dy, burnR[clr], burnG[clr], burnB[clr]);
	}
}

int SetupParticle(int place, int fuel, int burnOut, int x, int y, int dx, int dy, int r0, int g0, int b0) {
	place = NextFreeParticle(place);
	particles[place].fuel = fuel;
	particles[place].burnOut = burnOut;
	particles[place].x = x;
	particles[place].y = y;
	particles[place].dx = dx;
	particles[place].dy = dy;
	particles[place].r0 = r0;
	particles[place].g0 = g0;
	particles[place].b0 = b0;
	return place;
}

int NextFreeParticle(int from) {
	while (from < N_PARTICLES) {
		if (particles[from].fuel <= 0) {
			return from;
		}
		from++;
	}
	return N_PARTICLES - 1;
}

void PremultiplyParticles() {
	multParticleHolder = new DWORD[256 * 25];
	for (int i = 0; i < 256; i++) {
		LPDWORD particle = &(multParticleHolder[25 * i]);
		multParticle[i] = particle;
		double k = (i + 0.9999);
		for (int j = 0; j < 25; j++) {
			particle[j] = (DWORD)(k * idealParticle[j] / 36.0);
			//particle[j] = (DWORD)(k * boldParticle[j] / 255.0);
		}
	}
}

void ProcessParticles() {
	for (int i = 0; i < N_PARTICLES; i++) {
		if (particles[i].fuel > 0) {
			particles[i].fuel -= particles[i].burnOut;
			particles[i].x += particles[i].dx;
			particles[i].y += particles[i].dy;
			if ((particles[i].burnOut < 0) && (particles[i].fuel > 255)) {
				particles[i].burnOut = -particles[i].burnOut;
				particles[i].fuel -= particles[i].burnOut;
			}
		}
	}
}
