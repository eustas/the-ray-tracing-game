
int alarm, phaseF, nextQuant = 0;

long rDownX, rDownY;
bool rDown = false;

int rWin;

void KillGremlin() {
	for (int i = 0; i < levelData.maxG; i++) {
		if (gremlins[i].alive) {
			int x = (gremlins[i].x - 32) / 64;
			int y = (gremlins[i].y - 32) / 64;
			if ((controlX == x) && (controlY == y)) {
				gremlins[i].alive = false;
				score += 60;
				diagnose = DIAGNOSE_BUG_FIXED; diagnoseState = 256;
			}
		}
	}
}


void GameRDown(int win) {
	POINT pnt;

	rWin = win;
	pnt.x = (rWin == 0) ? 256 : 150;
	pnt.y = (rWin == 0) ? 192 : 150;
	ClientToScreen((rWin == 0) ? wnd : mwnd, &pnt);
	SetCursorPos(pnt.x, pnt.y);

	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	rDownX = pnt.x;
	rDownY = pnt.y;
	rDown = true;
}

void GameRUp() {
	GameMMove();
	rDown = false;
}

void GameMMove() {
	if (!rDown) {return;}
	POINT pnt;
	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	long x = pnt.x;
	long y = pnt.y;
	long dx = x - rDownX;
	long dy = y - rDownY;
	//rDownX = x;
	//rDownY = y;

	double tmp = rayBeta + (dy * 0.006);
	if (tmp > PI_2) {tmp = PI_2;}
	if (tmp < 0.0) {tmp = 0.0;}
	rayBeta = tmp;

	tmp = rayAlpha + (dx * 0.01);
	while (tmp > PIx2) {tmp -= PIx2;}
	while (tmp < 0.0) {tmp += PIx2;}
	rayAlpha = tmp;

	pnt.x = (rWin == 0) ? 256 : 150;
	pnt.y = (rWin == 0) ? 192 : 150;
	ClientToScreen((rWin == 0) ? wnd : mwnd, &pnt);
	SetCursorPos(pnt.x, pnt.y);
	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	x = pnt.x;
	y = pnt.y;
	rDownX = x;
	rDownY = y;

}

void PlayQuant() {
	Media(quantsR[nextQuant]);
	Media(quantsP[nextQuant]);
	nextQuant = (nextQuant + 1) & 0x7;
}

void ShiftMirrors() {
	rayLen = 0;
	CELL* cells = levelData.cells;
	for (int i = 0; i < 9 * 15; i++) {
		if (cells[i].type == MIRROR) {
			int r = (rand() % 3) - 1;
			cells[i].state = (char)((cells[i].state + 16 + r) & 0xF);
		}
	}
}

void NextLevel() {
	int last = wcscmp(LAST_LEVEL, levelData.nextLevel);
	if (last == 0) {
		subFrame = 0;
		subStage = GAME_END;
	} else {
		levelScore = score;
		StartLevel(levelData.nextLevel);
	}
}

void LooseTry() {
	if (triesLeft == 1) {
		subFrame = 0;
		subStage = GAME_OVER;
	} else {
		triesLeft--;
		InitLevel();
	}
}

void InitLevel() {
	LoadLevel(levelFile, &levelData);
	subStage = GAME_RECHARGE;
	subFrame = temperature = energy = phase = phaseF = 0;
	rayDensity[0] = rayDensity[1] = rayDensity[2] = rayDensity[3] = true;
	controlX = levelData.controlX; controlY = levelData.controlY;
	score = levelScore;
	diagnose = 7; diagnoseState = 0;
	ClearParticles();
	StartSound(FILE_RECHARGE);
	alarm = ALARM_OK;
	for (int i=0;i<5;i++) {
		gremlins[i].alive = false;
	}
}

void StartGame() {
	srand(GetTickCount());
	rayDist = 2.5;
	rayAlpha = PI_2;
	rayBeta = PI_2 / 4;

	DockMicroscope();
	//ShowWindow(mwnd, SW_SHOWDEFAULT); UpdateWindow(mwnd);
	levelNN = levelScore = 0;
	StartLevel(firstLevel);
}

void OfferGameExit() {
	int oldSubStage = subStage;
	subStage = GAME_EXIT;
	int but = MessageBox(NULL,l10nStr[L10N_EXIT_TEXT], l10nStr[L10N_EXIT_CAPTION], MB_YESNO | MB_ICONEXCLAMATION);
	if (but == IDYES) {ShowMenu();}
	else {subStage = oldSubStage;}
}

void StartLevel(WCHAR *fileName) {
	microX = -1.0;
	microY = -1.0;
	stage = GAME;
	levelNN++;
	wcscpy_s(levelFile, 256, fileName);
	triesLeft = 3;
	InitLevel();
}

void GameClick() {
	POINT pnt;
	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	long x0 = pnt.x, y0 = pnt.y;

	if ((y0 >= 16) && (x0 >= 16 ) && (x0 < (512 - 16)) && (y0 < (384 - 64 - 16))) {
		microX = (x0 - 16) / 32.0;
		microY = (y0 - 16) / 32.0;
		return;
	}

	if (y0 < 384 - 32) {return;}

	if (x0 >= 512 - 32) {
		if (IsWindowVisible(mwnd)) { ShowWindow(mwnd, SW_HIDE); }
		else { ShowWindow(mwnd, SW_SHOWDEFAULT); }
		return;
	}

	if (x0 >= 512 - 64) {
		antiAliasing = !antiAliasing;
		return;
	}

	if (x0 >= 512 - 96) {
		glassPoint = !glassPoint;
		return;
	}

}

void SetAlarm(int newAlarm) {
	if (newAlarm == alarm) {return;}
	alarm = newAlarm;
	if (alarm == ALARM_OK) { StopSound(); }
	else if (alarm == ALARM_WARN) { StartSound(FILE_ALARM); }
	else { StartSound(FILE_FATAL); }
}

void HitCollector(int dir) {
	if (levelData.cDir == dir) {
		subFrame = 0;
		subStage = GAME_WIN;
	}
}

void ProcessQuants() {
	for (int i = 0; i < 15 * 9; i++) {
		if (levelData.cells[i].type == QUANT) {
			if (levelData.cells[i].state != QUANT) {
				levelData.cells[i].state--;
				if (levelData.cells[i].state == 0) {levelData.cells[i].type = EMPTY;}
			}
		}
	}
}

void KillQuant(int x, int y) {
	if (levelData.cells[x + 15 * y].state == QUANT) {
		PlayQuant();
		int pcount = 0, place = 0;
		while (pcount != 64) {
			int x1 = rand() % 17 - 8, y1 = rand() % 17 - 8;
			if ((x1 * x1) + (y1 * y1) > 256) { continue; }
			int x2 = rand() % 9 - 4, y2 = rand() % 9 - 4;
			if ((x2 * x2) + (y2 * y2) > 16) { continue; }
			int color = rand() % 3;
			place = SetupParticle(place, 255, 10, 2 * (x1 + 32 + 32 * x),2 * (y1 + 32 + 32 * y),x2,y2, blowR[color], blowG[color], blowB[color]);
			pcount++;
		}
		levelData.quantCount--;
		score += 120;
		levelData.cells[x + 15 * y].state = 14;
		if (levelData.quantCount == 0) {
			score += 250;
			int ox = levelData.ox;
			int oy = levelData.oy;
			levelData.cells[ox + 15 * oy].type = EMPTY;
			int p = 33 + (ox * 2) + (2 * 32 * oy);
			levelData.qcells[p] = levelData.qcells[p +  1] = levelData.qcells[p + 32] = levelData.qcells[p + 33] = 0;

			pcount = 0; place = 0;
			while (pcount != 128) {
				int x2 = rand() % 7 - 3, y2 = rand() % 7 - 3;
				if ((x2 * x2) + (y2 * y2) > 9) { continue; }
				int x1 = rand() & 0x1F, y1 = rand() & 0x1F;
				int color = rand() % 3;
				place = SetupParticle(place, 255, 4, 2 * (x1 + 16 + 32 * ox), 2 * (y1 + 16 + 32 * oy),x2,y2, blowR[color], blowG[color], blowB[color]);
				pcount++;
			}
		}
	}
}

void CastRay() {
	int idx = 0;
	int ndir = 0;
	int nstate = 0;

	CELL* cells = levelData.cells;
	char* qcells = levelData.qcells;

	int dir = levelData.eDir;
	int x = (levelData.ex * 4) + 4 + deltaX[dir];
	int y = (levelData.ey * 4) + 4 + deltaY[dir];

	rayParts[idx].sx = x;
	rayParts[idx].sy = y;
	rayParts[idx].dir = dir;
	idx++;

	while (true) {
		x += deltaX[dir];
		y += deltaY[dir];

		if ((x <= 1) || (y <= 1) || (x >= 63) || (y >= 39)) {
			rayParts[idx].sx = x;
			rayParts[idx].sy = y;
			rayParts[idx].dir = dir;
			idx++;
			break;
		}

		int xb = 4 * x - deltaX[dir];
		int yb = 4 * y - deltaY[dir];
		if (xb < 0) {xb = -1;} else {xb = xb / 8;}
		if (yb < 0) {yb = -1;} else {yb = yb / 8;}
		int xa = 4 * x + deltaX[dir];
		int ya = 4 * y + deltaY[dir];
		if (xa < 0) {xa = -1;} else {xa = xa / 8;}
		if (ya < 0) {ya = -1;} else {ya = ya / 8;}

		if (((x & 0x3) == 0) && ((y & 0x3) == 0)) {
			int cx = x - 2;
			int cy = y - 2;
			if (cx < 0) {cx = -1;} else {cx = cx / 4;}
			if (cy < 0) {cy = -1;} else {cy = cy / 4;}
			char type = cells[cx + 15 * cy].type;
			int state = cells[cx + 15 * cy].state;
			if (type == EMITTER) {break;}
			if (type == COLLECTOR) {HitCollector(dir); break;}
			if (type == QUANT) {KillQuant(cx, cy); break;}
			if (type == STAR) { overloadPlus = DUST_BURN; break;}
			if ((type == MIRROR) || (type == ROTOR)) {
				nstate = (state + ((type == ROTOR) ? (phase) : 0)) & 0xF;
				ndir = (48 - dir - nstate) & 0xF;
				if (((dir + 8) & 0xF) == ndir) { overloadPlus = FULL_REFLECT; break;}
				dir = ndir;
				rayParts[idx].sx = x;
				rayParts[idx].sy = y;
				rayParts[idx].dir = dir;
				idx++;
				continue;
			} else if (type == PORTAL) {
				int d = 1;
				if ((levelData.px[state * 2] != cx) || (levelData.py[state * 2] != cy)) {d = 0;}
				x = (levelData.px[state * 2 + d] * 4) + 4;
				y = (levelData.py[state * 2 + d] * 4) + 4;
				rayParts[idx].sx = x;
				rayParts[idx].sy = y;
				rayParts[idx].dir = dir;
				idx++;
				continue;
			}
		}

		if (((x & 0x3) == 2) || ((y & 0x3) == 2)) {
			int xc = 4 * x - 8 - deltaX[dir];
			int yc = 4 * y - 8 - deltaY[dir];
			if (xc < 0) {xc = -1;} else {xc = xc / 16;}
			if (yc < 0) {yc = -1;} else {yc = yc / 16;}
			int xd = 4 * x - 8 + deltaX[dir];
			int yd = 4 * y - 8 + deltaY[dir];
			if (xd < 0) {xd = -1;} else {xd = xd / 16;}
			if (yd < 0) {yd = -1;} else {yd = yd / 16;}
			char type = EMPTY;
			char state = 0;
			if ((xd >= 0) && (xd < 15) && (yd >= 0) && (yd < 9)) {
				type = cells[xd + 15 * yd].type;
				state = cells[xd + 15 * yd].state;
			}
			if (type == EMPTY) {}
			else if (type == RANDOM) {
				int dir0;
				if (yc > yd) {dir0 = 0;}
				else if (yc < yd) {dir0 = 8;}
				else if (xc < xd) {dir0 = 4;}
				else {dir0 = 12;}
				x = xd * 4 + 4;
				y = yd * 4 + 4;
				DWORD r1 = 31 * (phaseF + x);
				DWORD r2 = (7 * phaseF + y + 3);
				DWORD r3 = (11 * phaseF + x * y + 7);
				DWORD rnd = ((r1 * r2 * r3) / 71) % 9;
				dir = (dir0 + rnd) & 0xF;
				rayParts[idx].sx = x;
				rayParts[idx].sy = y;
				rayParts[idx].dir = dir;
				idx++;
				continue;
			} else if ((type == POLAR) || (type == ROLAR)) {
				nstate = (8 - state + ((type == ROLAR) ? phase : 0)) & 0x7;
				if ((dir & 0x7) != nstate) {break;}
			} else if ((type == MPOLAR) || (type == MROLAR)) {
				nstate = (8 - state + ((type == MROLAR) ? phase : 0)) & 0x7;
				if ((dir & 0x7) != nstate) {
					if (xc == xd) {ndir = (48 - dir) & 0xF;}
					else if (yc == yd) {ndir = (40 - dir) & 0xF;}
					else {
						if ((xa - xb) == (ya - yb)) {ndir = (44 - dir) & 0xF;}
						else {ndir = (36 - dir) & 0xF;}
					}
					if (((dir + 8) & 0xF) == ndir) {overloadPlus = PARTIAL_REFLECT; break;}
					dir = ndir;
					rayParts[idx].sx = x;
					rayParts[idx].sy = y;
					rayParts[idx].dir = dir;
					idx++;
					continue;
				}
			}
		}

		char qtype = qcells[xa + 32 * ya];
		if (qtype == BLOCK) {
			if (dir == 0) {
				if (qcells[xa + (32 * ya) - 32] != 0) {break;}
			}
			else {break;}
		}
		else if (qtype == MBLOCK) {
			ndir = (dir + 8) & 0xF;
			if (xa == xb) {
				ndir = (48 - dir) & 0xF;
				if (dir == 4) {
					if (qcells[xa - 1 + 32 * ya] != MBLOCK) {ndir = dir;}
				}
			} else if (ya == yb) {ndir = (40 - dir) & 0xF;}
			else {
				if (qcells[xb + 32 * ya] == MBLOCK) {
					if (qcells[xa + 32 * yb] == MBLOCK) {
						if ((dir & 0x1) == 0) {
							overloadPlus = PARTIAL_REFLECT; break;
						}
						ndir = edgeDir[dir];
					} else {
						ndir = (48 - dir) & 0xF;
					}
				} else if (qcells[xa + 32 * yb] == MBLOCK) {ndir = (40 - dir) & 0xF;}
				else {
					if ((xa - xb) == (ya - yb)) {ndir = (44 - dir) & 0xF;}
					else {ndir = (36 - dir) & 0xF;}
				}
			}
			if (((dir + 8) & 0xF) == ndir) {overloadPlus = PARTIAL_REFLECT; break;}
			dir = ndir;
			rayParts[idx].sx = x;
			rayParts[idx].sy = y;
			rayParts[idx].dir = dir;
			idx++;
			continue;
		}

		rayParts[idx].sx = x;
		rayParts[idx].sy = y;
		rayParts[idx].dir = dir;
		idx++;
	}

	rayLen = idx;
}

void RenderGame() {
	if (subStage == GAME_EXIT) { return; }

	subFrame++;
	diagnoseState = (diagnoseState > 4) ? (diagnoseState - 4) : 0;

	if (subStage == GAME_RECHARGE) {
		energy += 40;
		if (energy > 128 * ENERGY_SCALE) {
			energy = 128 * ENERGY_SCALE;
			subStage = GAME_PLAY;
			StopSound();
		}
		if ((subFrame & 0x1) == 0) { ShiftMirrors(); }
		memset((LPBYTE)bits, 0, 512 * 384 * 4);

		diagnose = DIAGNOSE_RECHARGE; diagnoseState = 256;

		DrawCells();
		DrawCells2();
		DrawPanel();
		return;
	} else if (subStage == GAME_BLOW) {
		ProcessParticles();
		DrawParticles();
		if (subFrame == 120) { LooseTry(); }
		return;
	} else if (subStage == GAME_DISCHARGE) {
		if ((subFrame & 0x1) == 0) { DrawFadeOut(); }
		if (subFrame == 120) { LooseTry(); }
		return;
	} if (subStage == GAME_WIN) {
		diagnose = DIAGNOSE_CONNECT; diagnoseState = 256;
		if (energy > 0) {
			energy -= (int)((128 * ENERGY_SCALE) / 50);
			score += 10;
		} else {
			if (subFrame >= 60) {
				if ((subFrame % 15) == 0) {
					if (triesLeft > 0) {
						triesLeft--;
						score += 500;
					}
				}
			}
		}
		DrawPanel();
		if (subFrame == 120) { NextLevel(); }
		return;
	} if (subStage == GAME_END) {
		DrawMessage(gameWin);
		if (subFrame == 160) { ShowMenu(); }
		return;
	} if (subStage == GAME_OVER) {
		DrawMessage(gameOver);
		if (subFrame == 80) { ShowMenu(); }
		return;
	}

	GameMMove();

	ProcessQuants();
	ProcessParticles();

	if ((subFrame & 0x1) == 0) {
		phase = (phase + 15) & 0xF;
		ProcessGameInput();
	}
	if (((subFrame % 2) < 1) && ((subFrame % 5) < 3) && ((subFrame % 7) < 5)) {
		phaseF++;
	}
	if (levelData.maxG > 0) {
		if ((rand() % 500) <= levelData.gBorn) {
			int g = rand() % levelData.maxG;
			if (!gremlins[g].alive) {
				gremlins[g].alive = true;
				gremlins[g].x = 32 + rand() % 960;
				gremlins[g].y = 32 + rand() % 576;
			}
		}
	}
	int dx, dy, live, burn, idx;
	int parti = NextFreeParticle(0);
	for (int i = 0; i < levelData.maxG; i++) {
		if (gremlins[i].alive) {
			for (int j = 0; j < 3; j++) {
				while (true) {
					dx = rand() % 13 - 6;
					dy = rand() % 13 - 6;
					if (((dx * dx) + (dy * dy)) <= 36) {break;}
				}
				live = (rand() % 15) + 6;
				int color = (rand() % 3);
				parti = SetupParticle(parti, 1, -(255 / live), gremlins[i].x - live * dx, gremlins[i].y - live * dy, dx, dy, blowR[color], blowG[color], blowB[color]);
			}
			dx = gremlins[i].x + (rand() % (levelData.gSpeed * 2 + 1)) - levelData.gSpeed;
			if (dx < 32) {dx = 32;} if (dx >= 960 + 32) {dx = 960 + 31;}
			dy = gremlins[i].y + (rand() % (levelData.gSpeed * 2 + 1)) - levelData.gSpeed;
			if (dy < 32) {dy = 32;} if (dy >= 576 + 32) {dy = 576 + 31;}
			gremlins[i].x = dx;
			gremlins[i].y = dy;
			dx = (dx - 32) / 64;
			dy = (dy - 32) / 64;
			idx = dx + (15 * dy);
			if (levelData.cells[idx].type == MIRROR) {
				if ((rand() % 60) <= levelData.gEvil) {
					levelData.cells[idx].state = (levelData.cells[idx].state + 15 + (rand() % 3)) & 0xF;
				}
			}
		}
	}

	overloadPlus = -COOL_DOWN;

	while (castingRay.compare_and_swap(true,false) != false) {
		Sleep(0);
	}
	CastRay();
	while (castingRay.compare_and_swap(false,true) != true) {
		Sleep(0);
	}

	temperature += overloadPlus;
	if (temperature < 0) {
		temperature = 0;
	}

	if (overloadPlus > 0) {
		SetAlarm((overloadPlus > PARTIAL_REFLECT) ? ALARM_ERROR : ALARM_WARN);
		diagnose = (temperature > 96 * TEMPERATURE_SCALE) ? DIAGNOSE_MELTDOWN : DIAGNOSE_WARMUP;
		diagnoseState = 256;
	} else {
		SetAlarm(ALARM_OK);
	}

	if (energy > 0) {
		energy--;
	} else {
		subFrame = 0;
		SetAlarm(ALARM_OK);
		StartSound(FILE_DISCHARGE);
		subStage = GAME_DISCHARGE;
	}

	if (energy < 16 * ENERGY_SCALE) {
		bool hid = ((subFrame % 7) < 4) && ((subFrame % 11) < 5);
		rayDensity[0] = (((subFrame + 4) & 0x3) == 0) && hid;
		rayDensity[1] = (((subFrame + 3) & 0x3) == 0) && hid;
		rayDensity[2] = (((subFrame + 2) & 0x3) == 0) && hid;
		rayDensity[3] = (((subFrame + 1) & 0x3) == 0) && hid;
		diagnose = DIAGNOSE_LOW_ENERGY; diagnoseState = 256;
	}

	if (temperature > 128 * TEMPERATURE_SCALE) {
		temperature = 128 * TEMPERATURE_SCALE;
		subFrame = 0;
		PrepareBlowParticles();
		subStage = GAME_BLOW;
		SetAlarm(ALARM_OK);
		StartSound(FILE_BOOM);
		diagnose = DIAGNOSE_BLOW; diagnoseState = 256;
	}

	memset((LPBYTE)bits, 0, 512 * 384 * 4);
	DrawCells();
	DrawRay();
	DrawCells2();
	DrawPanel();
	DrawParticles();
	DrawControl();
	DrawMicro();
}
