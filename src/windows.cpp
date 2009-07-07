bool isLaptop = false;
DWORD lastOnBatteryQuestTick = 0;

void DoWheel(int zDelta) {
	double z = zDelta / WHEEL_DELTA;
	double tmp = rayDist * exp(-0.12 * z);
	if (tmp <  0.5) {tmp =  0.5;}
	if (tmp > 15.0) {tmp = 15.0;}
	rayDist = tmp;
}

bool IsLaptop() {
	SYSTEM_POWER_STATUS sysPwrStat;
	if (!GetSystemPowerStatus(&sysPwrStat)) { return false; }
	return (sysPwrStat.BatteryFlag != 128);
}

bool OnBattery() {
	if (!isLaptop) {
		return false;
	}
	SYSTEM_POWER_STATUS sysPwrStat;
	if (!GetSystemPowerStatus(&sysPwrStat)) { return false; }
	return (sysPwrStat.ACLineStatus == 0);
}

void DockMicroscope() {
	RECT rcMicroscope = {0, 0, 300, 300};
	AdjustWindowRectEx(&rcMicroscope, WS_CAPTION | WS_SYSMENU, false, WS_EX_TOOLWINDOW);
	int w = rcMicroscope.right - rcMicroscope.left;
	int h = rcMicroscope.bottom - rcMicroscope.top;

	RECT rcMain;
	GetWindowRect(wnd, &rcMain);
	int mh = rcMain.bottom - rcMain.top;
	MoveWindow(mwnd, rcMain.right, rcMain.top + ((mh - h) / 2), w, h, TRUE);
}

LRESULT WINAPI MicroscopeMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
	int zDelta;
	LRESULT result;

	switch(msg) {
		case WM_CREATE:
			mwnd = hWnd;
			mdc = GetDC(hWnd);
			InitRayTracer();
			return 0;

		case WM_CLOSE:
			ShowWindow(hWnd, SW_HIDE);
			return 0;

		case WM_PAINT:
			ValidateRect(hWnd, NULL);
			return 0;

		case WM_NCHITTEST:
			result = DefWindowProc(hWnd, msg, wParam, lParam);
			if (result == HTCAPTION) {
				result = HTCLIENT;
			}
			return result;
		break;

		case WM_KEYDOWN:
			if (wParam == VK_F2) {antiAliasing = !antiAliasing;}
			break;
		case WM_HELP:
			glassPoint = !glassPoint;
			break;
		case WM_MOUSEWHEEL:
			zDelta = GET_WHEEL_DELTA_WPARAM(wParam);
			DoWheel(zDelta);
			return 0;
		case WM_RBUTTONDOWN:
			if (stage == GAME) {GameRDown(1);}
			break;
		case WM_RBUTTONUP:
			if (stage == GAME) {GameRUp();}
			break;

	}

	return DefWindowProc(hWnd, msg, wParam, lParam);
}


LRESULT WINAPI MainMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
	DWORD tick;
	int zDelta;

	switch(msg) {
		case WM_CREATE:
			wnd = hWnd;
			dc = GetDC(hWnd);
			SetTimer(hWnd, TIMER_ID, 25, NULL);
			SetCursor(cursor);
			ShowMenu();
			return 0;

		case MM_MCINOTIFY:
			if (wParam == MCI_NOTIFY_SUCCESSFUL) {
				if (stage == MENU) {PlayMenu();}
			}
			return 0;

		case WM_TIMER:
			Render();
			tick = GetTickCount();
			if ((tick < lastOnBatteryQuestTick) || (tick >= (lastOnBatteryQuestTick + 1000))) {
				onBattery = OnBattery();
				lastOnBatteryQuestTick = tick;
			}
			break;

		case WM_WINDOWPOSCHANGED:
			DockMicroscope();
			break;

		case WM_LBUTTONDOWN:
			if (stage == MENU) {
				MenuClick();
			} else if (stage == MANUAL) {
				ManualClick();
			} else if (stage == GAME) {
				GameClick();
			} else if (stage == ABOUT) {
				AboutClick();
			} else if (stage == HISCORE) {
				HiScoreClick();
			}
			break;

		case WM_KEYDOWN:
			if (wParam == VK_F1) {glassPoint = !glassPoint;}
			if (wParam == VK_F2) {antiAliasing = !antiAliasing;}
			break;
		case WM_MOUSEWHEEL:
			zDelta = GET_WHEEL_DELTA_WPARAM(wParam);
			DoWheel(zDelta);
			return 0;
		case WM_RBUTTONDOWN:
			if (stage == GAME) {GameRDown(0);}
			break;
		case WM_RBUTTONUP:
			if (stage == GAME) {GameRUp();}
			break;

		case WM_DESTROY:
			PostQuitMessage(0);
			return 0;

		case WM_PAINT:
			ValidateRect(hWnd, NULL);
			return 0;

	}

	return DefWindowProc(hWnd, msg, wParam, lParam);
}

void ProcessGameInput() {
#define Key(code) (GetAsyncKeyState(code) < 0)
	number tmp;
	if (Key(VK_ESCAPE)) {OfferGameExit();}
	if ((Key(0x53) || Key(VK_DOWN)) && (controlY < 8)) { controlY++; }
	if ((Key(0x57) || Key(VK_UP)) && (controlY > 0)) { controlY--; }
	if ((Key(0x41) || Key(VK_LEFT)) && (controlX > 0)) { controlX--; }
	if ((Key(0x44) || Key(VK_RIGHT)) && (controlX < 14)) { controlX++; }
	if ((Key(0x51) || Key(VK_SHIFT)) && (levelData.cells[controlX + 15 * controlY].type == MIRROR)) {
		//if ((subFrame & 0x3) == 0) {
			levelData.cells[controlX + 15 * controlY].state = (char)((levelData.cells[controlX + 15 * controlY].state + 1) & 0xF);
		//}
	}
	if ((Key(0x45) || Key(VK_CONTROL)) && (levelData.cells[controlX + 15 * controlY].type == MIRROR)) {
		//if ((subFrame & 0x3) == 0) {
			levelData.cells[controlX + 15 * controlY].state = (char)((levelData.cells[controlX + 15 * controlY].state + 15) & 0xF);
		//}
	}
	if (Key(0x47) || Key(VK_DELETE)) { tmp = rayAlpha - 0.09f; if (tmp < 0.0f) {tmp += PIx2;} rayAlpha = tmp;}
	if (Key(0x4A) || Key(VK_NEXT)) { tmp = rayAlpha + 0.09f; if (tmp > PIx2) {tmp -= PIx2;} rayAlpha = tmp;}
	if (Key(0x59) || Key(VK_HOME)) { tmp = rayBeta + 0.06f; if (tmp > PI_2) {tmp = PI_2;} rayBeta = tmp;}
	if (Key(0x48) || Key(VK_END)) { tmp = rayBeta - 0.06f; if (tmp < 0.0f) {tmp = 0.0f;} rayBeta = tmp;}
	if (Key(0x54) || Key(VK_INSERT)) { tmp = rayDist / 1.04f; if (tmp < 0.5f) {tmp = 0.5f;} rayDist = tmp;}
	if (Key(0x55) || Key(VK_PRIOR)) { tmp = rayDist * 1.04f; if (tmp > 15.0f) {tmp = 15.0f;} rayDist = tmp;}
#undef Key
}
