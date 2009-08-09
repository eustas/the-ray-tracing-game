int stage, subStage, subFrame, phase;

bool hasManual;
bool hasHiScores;

bool glassPoint = false;
bool antiAliasing = false;

number rayDist = 2.5;
number rayAlpha = PI_2;
number rayBeta = PI_2 / 4;
number microX = -1.0;
number microY = -1.0;
int demoMode = 0;
number fpsGap = 1.0;

number doubmax; // approx square root of max number number
number doubmin; // smallest number number
//number doubtol; // tolerance of number numbers

number lightAlphaDelta = 0.005;
double lightAlpha = 0.0;

int stohasticIdx = -1;

atomic<bool> castingRay;

WCHAR debug[256];
FILE* debugFile;

int microTrans[128];

HWND wnd, mwnd;
HDC dc, mdc;
HCURSOR cursor;
BITMAPINFO bmi, mbmi;

LPDWORD bits, mbits, buffer;
LPDWORD tiles[5], unplugged, panel, gameOver, gameWin, bigNum, smallNum, status, hiScores;
LPDWORD menu, menuHighlight, menuText, manualButtons, manual[6], objective, microFont, about;

LPDWORD microLines[256];

LPDWORD multParticle[256];
LPDWORD multParticleHolder;

bool* objectiveWindow;

int triesLeft, temperature, energy, score, levelScore, levelNN;
int controlX, controlY, overloadPlus;

int diagnose = 0;
int diagnoseState = 0;

WCHAR* l10nRaw;
LPWSTR l10nStr[20];
WCHAR levelFile[256], firstLevel[256];

LEVEL levelData;

RAY_PART* rayParts;
int rayLen;
bool rayDensity[4];

PARTICLE* particles;
bool onBattery = false;

GREMLIN gremlins[5];