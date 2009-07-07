bool hasHiScoresMusic;

int diagX, diagY, diagStartX, diagStartY;
bool hasMoreDiagBits;
int nextParticle;

void HiScoreClick() {
	ShowMenu();
}

void PlayHiScores() {
	if (hasHiScoresMusic) {StartSound(FILE_RECORDS_MUSIC);}
}

inline void NextDiagBit() {
	diagX++;
	diagY--;
	if ((diagX >= 512) || (diagY < 24)) {
		if (diagStartX > 0) {
			diagStartX--;
		} else if (diagStartY > 24) {
			diagStartY--;
		} else {
			hasMoreDiagBits = false;
			return;
		}
		diagX = diagStartX;
		diagY = diagStartY;
	}
}

void ProcessHiScores() {
	int pCount = 0;
	while (hasMoreDiagBits && (pCount < 14)) {
		int p = diagX + (512 * diagY);
		DWORD clr = buffer[p];
		if (clr != 0) {
			DWORD b0 = clr & 0xFF;
			DWORD g0 = (clr >> 8) & 0xFF;
			DWORD r0 = (clr >> 16) & 0xFF;
			int dx, dy;
			for (int cnt = 0; cnt < 4; cnt++) {
				while (true) {
					dx = (rand() % 9) - 4;
					dy = (rand() % 9) - 4;
					if ((dx * dx) + (dy * dy) > 16) { continue; }
					break;
				}
				particles[nextParticle].fuel = -1;
				SetupParticle(nextParticle, 255, 3, diagX * 2, diagY * 2, dx, dy, r0, g0, b0);
				nextParticle++;	if (nextParticle == N_PARTICLES) { nextParticle = 0; }
			}

			buffer[p] = (r0 << 8) | (g0 << 16) | b0;
			pCount++;
		}
		NextDiagBit();
	}
}

void RenderHiScores() {
	if (subStage == 0) {
		memcpy_s(buffer, 512 * 384 * 4, hiScores, 512 * 384 * 4);
		diagX = diagStartX = 511;
		diagY = diagStartY = 383;
		hasMoreDiagBits = true;
		ClearParticles();
		nextParticle = 0;
	}
	if (subStage == 980) {
		ShowMenu();
		return;
	}

	subStage++;
	if (subStage > 40) {
		ProcessHiScores();
	}
	for (int y = 0; y < 384; y++) {
		int sh1 = 512 * y;
		int sh2 = 512 * (384 - 1 - y);
		for (int x = 0; x < 512; x++) { bits[x + sh1] = buffer[x + sh2]; }
	}
	DrawParticles();
	ProcessParticles();
}
