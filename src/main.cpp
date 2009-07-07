#include "main.h"

#include "constants.cpp"
#include "variables.cpp"

#include "raytracer.cpp"

#include "drawers.cpp"
#include "loaders.cpp"
#include "particles.cpp"
#include "windows.cpp"
#include "game.cpp"
#include "manual.cpp"
#include "menu.cpp"
#include "hiscores.cpp"
#include "about.cpp"

void Render() {
	switch (stage) {
		case GAME:	RenderGame();		break;
		case MENU:	RenderMenu();		break;
		case MANUAL:	RenderManual();		break;
		case HISCORE:	RenderHiScores();	break;
		case ABOUT:	RenderAbout();		break;
	}

	SetDIBitsToDevice(dc, 0, 0, 512, 384, 0, 0, 0, 384, (LPBYTE)bits, &bmi, DIB_RGB_COLORS);
}

INT WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, INT) {
	_wfopen_s(&debugFile, FILE_DEBUG, L"w");
	for (int i = 0; i < 8; i++) {
		Media(quantsL[i]);
	}
	castingRay = false;

	isLaptop = IsLaptop();

	hasManualMusic = GetFileAttributes(FILE_MANUAL_MUSIC) != INVALID_FILE_ATTRIBUTES;
	hasMenuMusic = GetFileAttributes(FILE_MENU_MUSIC) != INVALID_FILE_ATTRIBUTES;
	hasHiScoresMusic = GetFileAttributes(FILE_RECORDS_MUSIC) != INVALID_FILE_ATTRIBUTES;
	hasAboutMusic = GetFileAttributes(FILE_ABOUT_MUSIC) != INVALID_FILE_ATTRIBUTES;

	hasManual = GetFileAttributes(FILE_MANUAL_01) != INVALID_FILE_ATTRIBUTES;
	hasHiScores = GetFileAttributes(FILE_HISCORES) != INVALID_FILE_ATTRIBUTES;

	cache_aligned_allocator<DWORD> dwordAllocator;
	for (int line = 0; line < 256; line++) {microLines[line] = dwordAllocator.allocate(256);}

	HINSTANCE inst = GetModuleHandle(NULL);
	HANDLE icon = LoadImage(inst, MAKEINTRESOURCE(MAIN_ICON), IMAGE_ICON, 0, 0, LR_DEFAULTSIZE);
	cursor = LoadCursor(NULL, IDC_ARROW);
	HBRUSH brush = CreateSolidBrush(0);

	bits = new DWORD[512 * 384];
	buffer = new DWORD[512 * 384];
	mbits = new DWORD[300 * 300];
	objectiveWindow = new bool[256 * 256];
	l10nRaw = new WCHAR[2048];
	particles = new PARTICLE[N_PARTICLES];
	memset(particles, 0, sizeof(PARTICLE) * N_PARTICLES);
	for (int i = 0; i< 128; i++) {microTrans[i] = 19;}
	microTrans['0'] = 0; microTrans['1'] = 1; microTrans['2'] = 2; microTrans['3'] = 3; microTrans['4'] = 4;
	microTrans['5'] = 5; microTrans['6'] = 6; microTrans['7'] = 7; microTrans['8'] = 8; microTrans['9'] = 9;
	microTrans['F'] =10; microTrans['P'] =11; microTrans['S'] =12; microTrans['m'] =13; microTrans['.'] =14;
	microTrans['s'] =15; microTrans['t'] =16; microTrans['d'] =17; microTrans['='] =18; microTrans[' '] =19;
	PremultiplyParticles();

	levelData.cells = new CELL[9 * 15];
	levelData.qcells = new char[20 * 32];
	rayParts = new RAY_PART[10240];

	RECT rcMain = {0, 0, 512, 384};
	RECT rcMicroscope = {0, 0, 300, 300};
	AdjustWindowRect(&rcMain, WS_CAPTION | WS_SYSMENU, false);
	AdjustWindowRectEx(&rcMicroscope, WS_CAPTION | WS_SYSMENU, false, WS_EX_TOOLWINDOW);

	bmi.bmiHeader.biSize = sizeof(bmi.bmiHeader);
	bmi.bmiHeader.biBitCount = 32;
	bmi.bmiHeader.biCompression = BI_RGB;
	bmi.bmiHeader.biPlanes = 1;
	bmi.bmiHeader.biHeight = 384;
	bmi.bmiHeader.biWidth = 512;

	mbmi.bmiHeader.biSize = sizeof(bmi.bmiHeader);
	mbmi.bmiHeader.biBitCount = 32;
	mbmi.bmiHeader.biCompression = BI_RGB;
	mbmi.bmiHeader.biPlanes = 1;
	mbmi.bmiHeader.biHeight = 300;
	mbmi.bmiHeader.biWidth = 300;

	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			int r = ((x * 2) - 255) * ((x * 2) - 255) + ((y * 2) - 255) * ((y * 2) - 255);
			objectiveWindow[x + 256 * y] = (r <= (255 * 255));
		}
	}

	LoadAll();
	SetupAbout();
	DrawObjective();

	WNDCLASSEX wcMain = {
		sizeof(WNDCLASSEX), CS_CLASSDC, MainMsgProc, 0L, 0L,
		inst, (HICON)icon, cursor, brush, NULL,
		MAIN_CLS, NULL
	};
	RegisterClassEx(&wcMain);

	WNDCLASSEX wcMicroscope = {
		sizeof(WNDCLASSEX), CS_CLASSDC, MicroscopeMsgProc, 0L, 0L,
		inst, (HICON)icon, cursor, brush, NULL,
		MICROSCOPE_CLS, NULL
	};
	RegisterClassEx(&wcMicroscope);

	HWND hWnd = CreateWindow(
		MAIN_CLS, l10nStr[L10N_MAIN_WINDOW],
		WS_CAPTION | WS_SYSMENU, 32 , 32, rcMain.right - rcMain.left, rcMain.bottom - rcMain.top,
		NULL, NULL, inst, NULL
	);

	CreateWindowEx(WS_EX_TOOLWINDOW,
		MICROSCOPE_CLS, l10nStr[L10N_MICROSCOPE],
		WS_CAPTION | WS_SYSMENU, 560 , 72, rcMicroscope.right - rcMicroscope.left, rcMicroscope.bottom - rcMicroscope.top,
		NULL, NULL, inst, NULL
	);

	ShowWindow(hWnd, SW_SHOWDEFAULT);
	UpdateWindow(hWnd);

	bool bGotMsg;
	MSG msg;
	msg.message = WM_NULL;

	while (true) {
		GetMessage(&msg, NULL, 0U, 0U);
		if (WM_QUIT == msg.message) {
			break;
		}
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}


	StopSound();

	UnregisterClass(MAIN_CLS, inst);
	UnregisterClass(MICROSCOPE_CLS, inst);

	mallocProcessShutdownNotification();

	fclose(debugFile);

	return 0;
}
