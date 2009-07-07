void DrawMicroText(char* text, int x, int y) {
	int i = 0;
	int l = strlen(text);
	for (int i = 0; i < 11; i++) {
		int c = ' ';
		if (i < l) {
			c = text[i] & 0xFF;
		}
		if (c < 128) {
				c = microTrans[c];
		} else {
			c = 19;
		}
		for (int yy = 0; yy < 12; yy++) {
			for (int xx = 0; xx < 8; xx++) {
				mbits[(i * 8) + x + xx + 300 * (300 - 1 - y - yy)] = microFont[xx + (c * 8) + (160 * yy)];
			}
		}
	}
}

void DrawObjective() {
	for (int yy = 0; yy < 300; yy++) {
		int sh1 = yy * 300;
		int sh2 = 300 * (300 - yy - 1);
		for (int xx = 0; xx < 300; xx++) {
			((LPDWORD)mbits)[sh2 + xx] = objective[sh1 + xx];
		}
	}
}

void DrawUnplugged() {
	for (int yy = 0; yy < 256; yy++) {
		int sh1 = yy * 256;
		int sh2 = 21 + 300 * (300 - 23 - yy - 1);
		for (int xx = 0; xx < 256; xx++) {
			if (objectiveWindow[xx + 256 * yy]) {
				((LPDWORD)mbits)[sh2 + xx] = unplugged[sh1 + xx];
			}
		}
	}
}

void DrawLevelNN() {
	int n1 = levelNN / 10;
	int n2 = levelNN % 10;
	DrawBigNum(13 * 8, 10 * 32, n1);
	DrawBigNum(13 * 8 + 24, 10 * 32, n2);
}

void DrawScore() {
	int scr = score;
	for (int i = 5; i>=0; i--) {
		int n = scr % 10;
		scr /= 10;
		DrawSmallNum(512 - (16 * 6) + 16 * i, 320 + 17, n);
	}
}

void DrawMessage(LPDWORD image) {
	for (int yy = 0; yy < 48; yy++) {
		int sh1 = yy * 128;
		int sh2 = 192 + 512 * (384 - yy - 168 - 1);
		for (int xx = 0; xx < 128; xx++) {
			((LPDWORD)bits)[sh2 + xx] = image[sh1 + xx];
		}
	}
}

void DrawBigNum(int x, int y, int n) {
	int sh0 = n * 24;
	for (int yy = 0; yy < 24; yy++) {
		int sh1 = sh0 + (yy * 240);
		int sh2 = x + 512 * (384 - yy - y - 1);
		for (int xx = 0; xx < 24; xx++) {
			((LPDWORD)bits)[sh2 + xx] = bigNum[sh1 + xx];
		}
	}
}

void DrawSmallNum(int x, int y, int n) {
	int sh0 = n * 16;
	for (int yy = 0; yy < 15; yy++) {
		int sh1 = sh0 + (yy * 160);
		int sh2 = x + 512 * (384 - yy - y - 1);
		for (int xx = 1; xx < 16; xx++) {
			((LPDWORD)bits)[sh2 + xx] = smallNum[sh1 + xx];
		}
	}
}

void DrawFadeOut() {
	for (int t = 0; t < 384 * 512; t++) {
		DWORD v = ((LPDWORD)bits)[t];
		bits[t] = ((v & 0xF0F0F0) >> 4) + ((v & 0xF8F8F8) >> 3) + ((v & 0xFCFCFC) >> 2) + ((v & 0xFEFEFE) >> 1);
	}
}

void DrawDiagnose() {
	int sh0 = 256 * 16 * diagnose;
	for (int y = 0; y < 16; y++) {
		int sh1 = sh0 + y * 256;
		int sh2 = (15 - y) * 512;
		for (int x = 0; x < 256; x++) {
			DWORD v = status[sh1 + x];
			if (v != 0) {
				DWORD vF0F = ((((v & 0xFF00FF) * diagnoseState)) >> 8) & 0xFF00FF;
				DWORD v0F0 = ((((v & 0x00FF00) * diagnoseState)) >> 8) & 0x00FF00;
				bits[sh2 + x] = vF0F | v0F0;
			} else {
				bits[sh2 + x] = 0;
			}
		}
	}
}

void DrawPanel() {
	for (int y = 0; y < 64; y++) {
		int sh1 = (512 * (64 - 1 - y));
		int sh2 = y * 512;
		for (int x = 0; x < 512; x++) {
			((LPDWORD)bits)[x + sh1] = panel[x + sh2];
		}
	}
	for (int i = 0; i < triesLeft; i++) {
		DrawTile((5 + i) * 32, 10 * 32, 4, 3, (LPDWORD)bits, tiles[0]);
	}

	int emax = energy / ENERGY_SCALE;
	emax = (emax < 128) ? emax : 128;
	emax = (emax > 2) ? emax : 2;
	int esh2 = 48 * 512 + 17 * 16;
	int esh1 = 4 * 512 + 17 * 16;
	for (int x = 0; x < emax; x++) {
		for (int y = 0; y < 9; y++) {
			((LPDWORD)bits)[x + (512 * y) + esh1] = panel[x + esh2];
		}
	}

	int tmax = temperature / TEMPERATURE_SCALE;
	tmax = (tmax < 128) ? tmax : 128;
	tmax = (tmax > 2) ? tmax : 2;
	int tsh2 = 16 * 512 + 17 * 16;
	int tsh1 = 36 * 512 + 17 * 16;
	for (int x = 0; x < tmax; x++) {
		for (int y = 0; y < 9; y++) {
			((LPDWORD)bits)[x + (512 * y) + tsh1] = panel[x + tsh2];
		}
	}

	DrawLevelNN();
	DrawScore();

	int mVis = IsWindowVisible(mwnd) ? 5 : 4;
	DrawTile(512 - 32, 384 - 32, mVis, 2, (LPDWORD)bits, tiles[0]);

	if (diagnoseState > 0) {
		DrawDiagnose();
	}
}

void DrawParticle4x4rd(int x, int y, DWORD *canvas, DWORD b, DWORD g, DWORD r, DWORD* particle) {
	g <<= 8;
	b <<= 16;
	for (int xx = 0; xx < 4; xx++) {
		if (((x + xx) < 0) || ((x + xx) >= 512)) {
			continue;
		}
		for (int yy = 0; yy < 4; yy++) {
			if (((y + yy) < 0) || ((y + yy) >= 384)) {
				continue;
			}
			DWORD before = canvas[x + xx + (512 * (384 - y - 1 - yy))];
			DWORD p = particle[xx + (4 * yy)];
			DWORD pp = 255 - p;
			DWORD rr = (before & 0xFF) * pp + (p * r);
			DWORD gg = (before & 0xFF00) * pp + (p * g);
			DWORD bb = (before & 0xFF0000) * pp + (p * b);
			DWORD after = (rr & 0xFF00) | (gg & 0xFF0000) | (bb & 0xFF000000);
			canvas[x + xx + (512 * (384 - y - 1 - yy))] = (after >> 8);
		}
	}
}

void DrawParticles() {
	for (int i = 0; i < N_PARTICLES; i++) {
		if (particles[i].fuel > 0) {
			DrawParticle5x5(particles[i].x / 2, particles[i].y / 2, bits,
				particles[i].r0, particles[i].g0, particles[i].b0, multParticle[particles[i].fuel]);
		}
	}
}

void DrawParticle5x5(int x, int y, DWORD *canvas, DWORD r, DWORD g, DWORD b, DWORD* particle) {
	g <<= 8;
	r <<= 16;
	x -= 2;
	y -= 2;
	for (int xx = 0; xx < 5; xx++) {
		if (((x + xx) < 0) || ((x + xx) >= 512)) {
			continue;
		}
		for (int yy = 0; yy < 5; yy++) {
			if (((y + yy) < 0) || ((y + yy) >= 384)) {
				continue;
			}
			DWORD before = canvas[x + xx + (512 * (384 - y - 1 - yy))];
			DWORD p = particle[xx + (5 * yy)];
			DWORD pp = 256 - p;
			DWORD rr = (before & 0xFF0000) * pp + (p * r);
			DWORD gg = (before & 0xFF00) * pp + (p * g);
			DWORD bb = (before & 0xFF) * pp + (p * b);
			DWORD after = (rr & 0xFF000000) | (gg & 0xFF0000) | (bb & 0xFF00);
			canvas[x + xx + (512 * (384 - y - 1 - yy))] = (after >> 8);
		}
	}
}




void DrawMicro() {
	if (microX < 0.0) {return;}
	if (((subFrame >> 3) & 0x1) == 0) {return;}
	int x = 16 + ((int)(32 * microX));
	int y = 16 + ((int)(32 * microY));
	DrawParticle5x5(x, y, bits, 255,255,255, (DWORD *)boldParticle);
}

void DrawControl() {
	int dx = 32 + 32 * controlX;
	int dy = 32 + 32 * controlY;

	double shift = phase / 16.0;

	double a0 = 2.0 * 3.1415 * shift;
	int x0 = (int)(15.0 * cos(a0));
	int y0 = (int)(15.0 * sin(a0));
	DrawParticle5x5(dx + x0, dy + y0, bits, 255,128,128, (DWORD *)boldParticle);

	double a1 = 2.0 * 3.1415 * (0.33333 + shift);
	int x1 = (int)(15.0 * cos(a1));
	int y1 = (int)(15.0 * sin(a1));
	DrawParticle5x5(dx + x1, dy + y1, bits, 128,255,128, (DWORD *)boldParticle);

	double a2 = 2.0 * 3.1415 * (0.66666 + shift);
	int x2 = (int)(15.0 * cos(a2));
	int y2 = (int)(15.0 * sin(a2));
	DrawParticle5x5(dx + x2, dy + y2, bits, 128,128,255, (DWORD *)boldParticle);

	double a3 = 2.0 * 3.1415 * (0.16666 + shift);
	int x3 = (int)(15.0 * cos(a3));
	int y3 = (int)(15.0 * sin(a3));
	DrawParticle5x5(dx + x3, dy + y3, bits, 255,255,0, (DWORD *)boldParticle);

	double a4 = 2.0 * 3.1415 * (0.5 + shift);
	int x4 = (int)(15.0 * cos(a4));
	int y4 = (int)(15.0 * sin(a4));
	DrawParticle5x5(dx + x4, dy + y4, bits, 0,255,255, (DWORD *)boldParticle);

	double a5 = 2.0 * 3.1415 * (0.83333 + shift);
	int x5 = (int)(15.0 * cos(a5));
	int y5 = (int)(15.0 * sin(a5));
	DrawParticle5x5(dx + x5, dy + y5, bits, 255,0,255, (DWORD *)boldParticle);
}

void DrawRay() {
	for (int i = 0; i < rayLen; i++) {
		RAY_PART part = rayParts[i];
		int dir = part.dir;
		int x = (part.sx * 16) + (2 * rayshX[dir]);
		int y = (part.sy * 16) + (2 * rayshY[dir]);
		for (int j = 0; j < 16; j++) {
			if (rayDensity[j & 0x3]) {
				DrawParticle4x4rd(x / 2, y / 2, bits, 255, 0, 0, (DWORD *)lightParticle);
			}
			x += deltaX[dir];
			y += deltaY[dir];
		}
	}
}

void DrawCells() {
	CELL* cells = levelData.cells;
	for (int x = 0; x < 15; x++) {
		for (int y = 0; y < 9; y++) {
			int idx = x + (15 * y);
			char type = cells[idx].type;
			if (type == EMPTY) {
				continue;
			}
			int state = cells[idx].state;
			int px = 16 + (x * 32);
			int py = 16 + (y * 32);
			switch (type) {
				case MIRROR:
					DrawTile(px, py, state & 0x7, 1, bits, tiles[4 * (state >> 3)]);
				break;

				case POLAR:
					state = state & 0x7;
					DrawTile(px, py, state & 0x3, 2, bits, tiles[state & 0x4]);
				break;

				case MPOLAR:
					state = state & 0x7;
					DrawTile(px, py, state & 0x3, 3, bits, tiles[state & 0x4]);
				break;

				case QUANT:
					if (state == QUANT) {
						DrawQTile(px + 8, py + 8, 4, 1, bits);
					}
				break;

				case STAR:
					DrawQTile(px + 8, py + 8, 5, 1, bits);
				break;

				case ROTOR:
					state = (state + phase) & 0xF;
					DrawTile(px, py, state & 0x7, 1, bits, tiles[4 * (state >> 3)]);
				break;

				case ROLAR:
					state = (state + phase) & 0x7;
					DrawTile(px, py, state & 0x3, 2, bits, tiles[state & 0x4]);
				break;

				case MROLAR:
					state = (state + phase) & 0x7;
					DrawTile(px, py, state & 0x3, 3, bits, tiles[state & 0x4]);
				break;

				case PORTAL:
					DrawTile(px, py, 3 + state, 0, bits, tiles[0]);
				break;

				case RANDOM:
					DrawTile(px, py, 7, 0, bits, tiles[0]);
				break;

				case EMITTER:
					DrawTile(px, py, 0, 0, bits, tiles[state]);
				break;

				case COLLECTOR:
					DrawTile(px, py, 1, 0, bits, tiles[state]);
				break;
			}
		}
	}
}

void DrawCells2() {
	CELL* cells = levelData.cells;
	for (int x = 0; x < 15; x++) {
		for (int y = 0; y < 9; y++) {
			int idx = x + (15 * y);
			char type = cells[idx].type;
			if (type == EMPTY) {
				continue;
			}
			char state = cells[idx].state;
			int px = 16 + (x * 32);
			int py = 16 + (y * 32);
			switch (type) {
				case RANDOM:
					DrawTile(px, py, 7, 0, bits, tiles[0]);
				break;
			}
		}
	}

	char* qcells = levelData.qcells;
	int idx = 0;
	for (int y = 0; y < 20; y++) {
		for (int x = 0; x < 32; x++) {
			char type = qcells[idx];
			idx++;
			if (type == 0) {
				continue;
			}
			DrawQTile(x * 16, y * 16, (type == BLOCK) ? 5 : 4, 0, bits);
		}
	}
}

void DrawQTile(int dx, int dy, int sx, int sy, LPDWORD canvas) {
	for (int y = 0; y < 16; y++) {
		int sh1 = dx + (512 * (384 - dy - 1 - y));
		int sh2 = (sx * 16) + ((y + (sy * 16)) << 8);
		for (int x = 0; x < 16; x++) {
			canvas[x + sh1] = tiles[0][x + sh2];
		}
	}
}

void DrawTile(int dx, int dy, int sx, int sy, LPDWORD canvas, DWORD* bank) {
	for (int y = 0; y < 32; y++) {
		int sh1 = dx + (512 * (384 - dy - 1 - y));
		int sh2 = (sx * 32) + ((y + (sy * 32)) << 8);
		for (int x = 0; x < 32; x++) {
			canvas[x + sh1] = bank[x + sh2];
		}
	}
}

