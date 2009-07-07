//                       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
const int deltaX[16] = {-2,-2,-2,-1, 0, 1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2};
const int deltaY[16] = { 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 2, 2, 2, 2, 2, 1};
const int rayshX[16] = {-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-2,-2,-2};
const int rayshY[16] = {-2,-2,-2,-2,-2,-2,-2,-1,-2,-2,-2,-2,-2,-2,-2,-2};

const DWORD boldParticle[25] = {
	  0, 81,109, 81,  0,
	 81,169,199,169, 81,
	109,199,255,199,109,
	 81,169,199,169, 81,
	  0, 81,109, 81,  0
};
// 1
// 1 1
// 1 2 1
// 1 3 3 1
// 1 3 6 3 1
// 
//    01 03 06 03 01
// 01 01 03 06 03 01
// 03 03 09 18 09 03
// 06 06 18 36 18 06
// 03 03 09 18 09 03
// 01 01 03 06 03 01

const double idealParticle[25] = {
	 1,  3,  6,  3,  1,
	 3,  9, 18,  9,  3,
	 6, 18, 36, 18,  6,
	 3,  9, 18,  9,  3,
	 1,  3,  6,  3,  1
};

const int blowR[3] = { 255,   0, 128};
const int blowG[3] = { 192, 255, 192};
const int blowB[3] = { 128,   0, 255};

const int burnR[3] = { 255,  0, 255};
const int burnG[3] = { 128,  0,   0};
const int burnB[3] = {   0,  0,   0};

const DWORD lightParticle[16] = {
	  0, 64, 64,  0,
	 64,255,255, 64,
	 64,255,255, 64,
	  0, 64, 64,  0
};

const WCHAR* quantsP[8] = {
	L"play Q0",
	L"play Q1",
	L"play Q2",
	L"play Q3",
	L"play Q4",
	L"play Q5",
	L"play Q6",
	L"play Q7"
};
const WCHAR* quantsR[8] = {
	L"seek Q0 to start",
	L"seek Q1 to start",
	L"seek Q2 to start",
	L"seek Q3 to start",
	L"seek Q4 to start",
	L"seek Q5 to start",
	L"seek Q6 to start",
	L"seek Q7 to start"
};
const WCHAR* quantsL[8] = {
	L"open waveaudio!quant.wav alias Q0 wait",
	L"open waveaudio!quant.wav alias Q1 wait",
	L"open waveaudio!quant.wav alias Q2 wait",
	L"open waveaudio!quant.wav alias Q3 wait",
	L"open waveaudio!quant.wav alias Q4 wait",
	L"open waveaudio!quant.wav alias Q5 wait",
	L"open waveaudio!quant.wav alias Q6 wait",
	L"open waveaudio!quant.wav alias Q7 wait"
};
