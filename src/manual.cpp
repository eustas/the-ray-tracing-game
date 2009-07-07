bool hasManualMusic;

int forwBtn = 0, backBtn = 0;

void PlayManual() {
	if (hasManualMusic) {StartMusic(FILE_MANUAL_MUSIC);}
}

void ManualClick() {
	int btn = CalcManualButton();
	subStage += btn;
	if ((subStage < 0) || (subStage > 4)) {
		ShowMenu();
	}
}

int CalcManualButton() {
	POINT pnt;
	GetCursorPos(&pnt);
	ScreenToClient(wnd, &pnt);
	long x0 = pnt.x;
	long y0 = pnt.y;
	int result;

	if ((y0 >= 24) && (y0 + 24 < 384)) {
		return 0;
	}
	if (y0 < 24) {
		x0 -= 16;
		y0 -= 12;
		result = -1;
	} else {
		x0 -= (512 - 16);
		y0 -= (384 - 12);
		result = 1;
	}

	long r = (x0 * x0 * 12 * 12) + (y0 * y0 * 16 * 16);
	if (r <= 12 * 12 * 16 * 16) {
		return result;
	}
	return 0;
}

void RenderManual() {
	int btn = CalcManualButton();
	forwBtn = (forwBtn - 6 > 0) ? forwBtn - 6 : 0;
	forwBtn += (btn == BTN_FORW) ? 16 : 0;
	forwBtn = (forwBtn < 128) ? forwBtn : 128;

	backBtn = (backBtn - 6 > 0) ? backBtn - 6 : 0;
	backBtn += (btn == BTN_BACK) ? 16 : 0;
	backBtn = (backBtn < 128) ? backBtn : 128;

	LPDWORD man = manual[subStage];

	for (int y = 0; y < 384; y++) {
		int sh1 = y * 512;
		int sh2 = (383 - y) * 512;
		if ((y > 25) && (y < 359)) {
			for (int x = 0; x < 512; x++) {
				bits[sh2 + x] = man[sh1 + x];
			}
		} else {
			DWORD z;
			if (y < 24) {
				z = backBtn;
			} else {
				z = forwBtn;
			}
			z += 96;
			DWORD zz = 256 - z;
			for (int x = 0; x < 512; x++) {
				DWORD v1 = man[sh1 + x];
				DWORD v3 = manualButtons[sh1 + x];
				if (v3 != 0) {
					DWORD vF0F = ((((v1 & 0xFF00FF) * zz) + ((v3 & 0xFF00FF) * z)) >> 8) & 0xFF00FF;
					DWORD v0F0 = ((((v1 & 0xFF00) * zz) + ((v3 & 0xFF00) * z)) >> 8) & 0xFF00;
					bits[sh2 + x] = vF0F | v0F0;
				} else {
					bits[sh2 + x] = v1;
				}
			}
		}
	}
}
