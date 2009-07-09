bool hasMenuMusic;
int buttonZ[5];

void PlayMenu() {
	if (hasMenuMusic) {StartMusic(FILE_MENU_MUSIC);}
}

void MenuClick() {
	int button = CalcActiveMenuButton();
	if (button == BTN_GAME) {
		StartGame();
	} else if (button == BTN_MANUAL) {
		PlayManual();
		stage = MANUAL;
		subStage = 0;
	} else if (button == BTN_HISCORES) {
		stage = HISCORE;
		subStage = 0;
		PlayHiScores();
	} else if (button == BTN_ABOUT) {
		stage = ABOUT;
		subStage = 0;
		StartAbout();
		PlayAbout();
	} else if (button == BTN_EXIT) {
		PostQuitMessage(0);
	}
}

int CalcActiveMenuButton() {
	POINT pnt;
	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	long x0 = pnt.x;
	long y0 = pnt.y;

	int Y = (y0 - 12) / 72;
	if ( (Y < 0) || (Y >= 5) ) { return -1; }
	if ((Y == 0) && (!hasManual)) { return -1;}
	if ((Y == 2) && (!hasHiScores)) { return -1;}

	long y1 = y0 - 48 - (Y * 72);
	long x1 = x0 - 256;
	long r = (x1 * x1 * 24 * 24) + (y1 * y1 * 128 * 128);
	if (r <= 24 * 24 * 128 * 128) {
		return Y;
	}
	return -1;
}

void RenderMenu() {
	for (int i = 0; i < 5; i++) {
		int z = buttonZ[i] - 6;
		if (z < 0) { z = 0; }
		buttonZ[i] = z;
	}
	int active = CalcActiveMenuButton();
	if (active >= 0) {
		int z = buttonZ[active] + 16;
		if (z > 128) { z = 128; }
		buttonZ[active] = z;
	}

	for (int y = 0; y < 384; y++) {
		int sh1 = y * 512;
		int sh2 = (383 - y) * 512;
		int Y = (y - 12) / 72;
		if ((Y < 0) || (Y >= 5)) { Y = 0; }
		DWORD z = buttonZ[Y] + 32;
		DWORD zz = 256 - z;
		for (int x = 0; x < 512; x++) {
			DWORD v1 = menu[sh1 + x];
			DWORD v3 = menuHighlight[sh1 + x];
			if (v3 != 0x0000FF) {
				DWORD vF0F = ((((v1 & 0xFF00FF) * zz) + ((v3 & 0xFF00FF) * z)) >> 8) & 0xFF00FF;
				DWORD v0F0 = ((((v1 & 0xFF00) * zz) + ((v3 & 0xFF00) * z)) >> 8) & 0xFF00;
				DWORD v2 = menuText[sh1 + x];
				if (v2 != 0) {
					DWORD v = vF0F | v0F0;
					bits[sh2 + x] = (((v & 0xFF0000) > (v2 & 0xFF0000)) ? (v & 0xFF0000) : (v2 & 0xFF0000)) |
						(((v & 0xFF00) > (v2 & 0xFF00)) ? (v & 0xFF00) : (v2 & 0xFF00)) |
						(((v & 0xFF) > (v2 & 0xFF)) ? (v & 0xFF) : (v2 & 0xFF));
				} else {
					bits[sh2 + x] = vF0F | v0F0;
				}
			} else {
				bits[sh2 + x] = v1;
			}
		}
	}
}

void ShowMenu() {
	ShowWindow(mwnd, SW_HIDE);
	UpdateWindow(mwnd);
	stage = MENU;
	buttonZ[0] = buttonZ[1] = buttonZ[2] = buttonZ[3] = buttonZ[4] = 0;
	PlayMenu();
}
