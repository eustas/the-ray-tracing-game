#ifndef MAIN_H
#define MAIN_H

#include "files.h"
#include<complex>
#include <windows.h>
#include <stdio.h>
#include <math.h>
#define PI   3.1415926f
#define PIx2 6.2831853f
#define PI_2 1.5707963f

#include <png.h>

#include <tbb\task_scheduler_init.h>
#include <tbb\cache_aligned_allocator.h>
#include <tbb\parallel_for.h>
#include <tbb\blocked_range.h>
#include <tbb\tick_count.h>
#include <tbb\atomic.h>

#include <tbbmalloc\Customize.h>
using namespace tbb;

#define StartSound(fileName) PlaySound(fileName, NULL, SND_FILENAME | SND_ASYNC)
#define StartMusic(fileName) PlaySound(fileName, NULL, SND_FILENAME | SND_ASYNC | SND_LOOP)
#define StopSound() PlaySound(NULL, NULL, SND_ASYNC)
#define Media(command) mciSendString(command, NULL, 0, NULL)

//
// TYPES
//

//#define number float
#define number double

typedef struct _GREMIN {
	bool alive;
	int x;
	int y;
} GREMLIN;

typedef struct _VEC4D {
	number x;
	number y;
	number z;
	number w;
} VEC4D;

typedef struct _PARTICLE {
	int fuel;
	int burnOut;
	DWORD r0, g0, b0;
	int x, y;
	int dx, dy;
} PARTICLE;

typedef struct _CELL {
	char type;
	char state;
} CELL;

typedef struct _RAY_PART {
	int sx;
	int sy;
	int dir;
	number d;
	number sta;
	number end;
	number staX;
	number staZ;
	number endX;
	number endZ;
	_RAY_PART* next;
} RAY_PART;

typedef struct _LEVEL {
	CELL* cells;
	char* qcells;
	int ex;
	int ey;
	int eDir;
	int cx;
	int cy;
	int cDir;
	int px[8];
	int py[8];
	int ox;
	int oy;
	int quantCount;
	int controlX;
	int controlY;
	int maxG;
	int gBorn;
	int gSpeed;
	int gEvil;
	WCHAR nextLevel[256];
} LEVEL;

// RAYTRACER
void PerformParallelRaytrace();
DWORD WINAPI RayTracerThreadFunction(LPVOID lParam);
void InitRayTracer();
void DumpToCanvas();
void DoRayTrace();
void RecalculateModelRays(number newModelAlpha, number newModelBeta, bool newAntiAliasing);
void Precalculate(number x0, number y0, number z0);
void BuildRay();

// LOADERS
void LoadAll();
void LoadTiles();
//LPDWORD LoadBmp(LPCWSTR fileName, int w, int h);
LPDWORD LoadPng(LPCWSTR fileName, int w, int h);
void LoadLevel(WCHAR* fileName, LEVEL* level);
void LoadL10n();

// DRAWERS
void DrawMicroText(char * text, int x, int y);
void DrawParticle4x4rd(int x, int y, DWORD *canvas, DWORD b, DWORD g, DWORD r, DWORD* particle);
void DrawParticle5x5(int x, int y, DWORD *canvas, DWORD b, DWORD g, DWORD r, DWORD* particle);
void DrawControl();
void DrawRay();
void DrawCells();
void DrawCells2();
void DrawTile(int dx, int dy, int sx, int sy, LPDWORD canvas, DWORD* bank);
void DrawQTile(int dx, int dy, int sx, int sy, LPDWORD canvas);
void DrawParticles();
void DrawPanel();
void DrawFadeOut();
void DrawBigNum(int x, int y, int n);
void DrawLevelNN();
void DrawMessage(LPDWORD image);
void DrawSmallNum(int x, int y, int n);
void DrawScore();
void DrawUnplugged();
void DrawObjective();
void DrawDiagnose();
void DrawMicro();

// WINDOWS
LRESULT WINAPI MicroscopeMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
LRESULT WINAPI MainMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
bool IsLaptop();
bool OnBattery();
void ProcessGameInput();
void DockMicroscope();
void DoWheel(int zDelta);

// PARTICLES
int  NextFreeParticle(int from);
void ProcessParticles();
void PremultiplyParticles();
void PrepareBlowParticles();
int  SetupParticle(int place, int fuel, int burnOut, int x, int y, int dx, int dy, int r0, int g0, int b0);
void ClearParticles();

// GAME
void StartLevel(WCHAR* fileName);
void InitLevel();
void LooseTry();
void NextLevel();
void ShiftMirrors();
void StartGame();
void GameClick();
void CastRay();
void HitCollector(int dir);
void KillQuant(int x, int y);
void RenderGame();
void ProcessQuants();
void SetAlarm(int newAlarm);
void PlayQuant();
void OfferGameExit();
void GameRDown(int win);
void GameRUp();
void GameMMove();
void KillGremlin();

// MAIN
void Render();

// MANUAL
void PlayManual();
void ManualClick();
void RenderManual();
int  CalcManualButton();

// HISCORES
void PlayHiScores();
void RenderHiScores();
void ProcessHiScores();
void NextDiagBit();
void HiScoreClick();

// MENU
void PlayMenu();
void RenderMenu();
int  CalcActiveMenuButton();
void MenuClick();
void ShowMenu();

// ABOUT
void PlayAbout();
void SetupAbout();
void StartAbout();
void RenderAbout();
void AboutClick();

// LEVEL ITEMS
#define EMPTY  '.'
#define MIRROR 'M'
#define QUANT  '@'
#define BLOCK  'B'
#define MBLOCK 'D'
#define ROTOR  'R'
#define RANDOM 'A'
#define PORTAL 'T'
#define POLAR  'P'
#define MPOLAR 'F'
#define ROLAR  'L'
#define MROLAR 'W'
#define STAR   '*'
#define EMITTER   'E'
#define COLLECTOR 'C'

// LEVEL ITEM STATES
#define LEFT  'L'
#define RIGHT 'R'
#define UP    'U'
#define DOWN  'D'
#define FULL  'F'
#define EMPTY_LEFT_UP    '1'
#define EMPTY_RIGHT_UP   '2'
#define EMPTY_RIGHT_DOWN '3'
#define EMPTY_LEFT_DOWN  '4'
#define LEFT_UP    '5'
#define RIGHT_UP   '6'
#define RIGHT_DOWN '7'
#define LEFT_DOWN  '8'

// STAGES
#define GAME	0
#define MANUAL	1
#define MENU	2
#define HISCORE	3
#define ABOUT	4

#define GAME_RECHARGE	0
#define GAME_PLAY	1
#define GAME_BLOW	2
#define GAME_OVER	3
#define GAME_WIN	4
#define GAME_DISCHARGE	5
#define GAME_END	6
#define GAME_EXIT	7

// DIAGNOSES
#define DIAGNOSE_RECHARGE	0
#define DIAGNOSE_WARMUP		1
#define DIAGNOSE_MELTDOWN	2
#define DIAGNOSE_CONNECT	3
#define DIAGNOSE_LOW_ENERGY	4
#define DIAGNOSE_BLOW		5
#define DIAGNOSE_BUG_FIXED	6

#define ALARM_OK 0
#define ALARM_WARN 1
#define ALARM_ERROR 2

// CONSTS
#define N_PARTICLES 4096

#define ENERGY_SCALE 32
#define TEMPERATURE_SCALE 100
#define COOL_DOWN 60
#define FULL_REFLECT 120
#define PARTIAL_REFLECT 60
#define DUST_BURN 120
#define OVER_MAX 1000

#define MAIN_CLS L"TheRayGame"
#define MICROSCOPE_CLS L"TheRayGameMicroscope"

#define BTN_MANUAL	0
#define BTN_GAME	1
#define BTN_HISCORES	2
#define BTN_ABOUT	3
#define BTN_EXIT	4

#define BTN_FORW	1
#define BTN_BACK	-1
// IDS
#define MAIN_ICON	101
#define TIMER_ID	102

// L10N
#define L10N_MAIN_WINDOW	0
#define L10N_MICROSCOPE		1
#define L10N_FIRST_LEVEL	2
#define L10N_EXIT_CAPTION	3
#define L10N_EXIT_TEXT		4

#define LAST_LEVEL L"LAST"

#endif
