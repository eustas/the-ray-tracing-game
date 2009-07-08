void LoadAll() {
	LoadL10n();
	LoadTiles();

	menu		= LoadBmp(FILE_MENU		, 512, 384);
	menuText	= LoadBmp(FILE_MENU_TEXT	, 512, 384);
	menuHighlight	= LoadBmp(FILE_MENU_HIGHLIGHT	, 512, 384);
	manualButtons	= LoadBmp(FILE_MANUAL_BUTTONS	, 512, 384);
	//manual[0]	= LoadBmp(FILE_MANUAL_01	, 512, 384);
	manual[0]	= LoadPng(FILE_MANUAL_01PNG	, 512, 384);
	manual[1]	= LoadBmp(FILE_MANUAL_02	, 512, 384);
	manual[2]	= LoadBmp(FILE_MANUAL_03	, 512, 384);
	manual[3]	= LoadBmp(FILE_MANUAL_04	, 512, 384);
	manual[4]	= LoadBmp(FILE_MANUAL_05	, 512, 384);
	hiScores	= LoadBmp(FILE_HISCORES		, 512, 384);
	about		= LoadBmp(FILE_ABOUT		, 512, 384);

	panel = LoadBmp(FILE_PANEL, 512, 64);
	unplugged = LoadBmp(FILE_UNPLUGGED, 256, 256);
	bigNum = LoadBmp(FILE_BIGNUM, 240, 24);
	smallNum = LoadBmp(FILE_SMALLNUM, 160, 15);
	gameOver = LoadBmp(FILE_GAME_OVER, 128, 48);
	gameWin = LoadBmp(FILE_GAME_WIN, 128, 48);
	objective = LoadBmp(FILE_OBJECTIVE, 300, 300);
	status = LoadBmp(FILE_STATUS, 256, 128);
	microFont = LoadBmp(FILE_MICRO_FONT, 160, 12);

	wcscpy_s(firstLevel, 256, l10nStr[L10N_FIRST_LEVEL]);
}

LPDWORD LoadBmp(LPCWSTR fileName, int w, int h) {
	LPDWORD result = new DWORD[w * h];

	BYTE* tmp = new BYTE[w * h * 3 + 54];
	memset(tmp, 0, w * h * 3 + 54);

	DWORD read = 0;
	HANDLE file = CreateFile(fileName, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	ReadFile(file, (LPVOID)tmp, w * h * 3 + 54, &read, NULL);
	CloseHandle(file);

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			int shi = 54 + (3 * (x + (y * w)));
			DWORD v = (tmp[shi]) + (tmp[shi + 1] << 8) + (tmp[shi + 2] << 16);
			result[x + ((h - 1 - y) * w)] = v;
		}
	}

	delete tmp;

	return result;
}

void LoadL10n() {
	FILE* in;
	_wfopen_s(&in, FILE_L10N_STR, L"r, ccs=UTF-8");
	WCHAR* head = l10nRaw;

	for (int i = 0; i <= L10N_EXIT_TEXT; i++) {
		l10nStr[i] = head;
		fgetws(head, 512, in);
		int l = wcslen(head);
		while (((head[l - 1] == 10) || (head[l - 1] == 13)) && (l > 0)) {
			head[l - 1] = 0;
			l--;
		}
		head[l] = 0;
		l++;
		head = &(head[l]);
	}

	fclose(in);
}

void LoadTiles() {
	tiles[0] = LoadBmp(FILE_TILES, 256, 128);

	tiles[1] = new DWORD[256 * 128];
	tiles[2] = new DWORD[256 * 128];
	tiles[3] = new DWORD[256 * 128];
	tiles[4] = new DWORD[256 * 128];

	for (int sx = 0; sx < 8; sx++) {
		for (int sy = 0; sy < 4; sy++) {
			for (int y = 0; y < 32; y++) {
				for (int x = 0; x < 32; x++) {
					tiles[4][(sx << 5) + 31 - y + (((sy << 5)      + x) << 8)] = tiles[0][(sx << 5) + x + (((sy << 5) + y) << 8)];
					tiles[3][(sx << 5) + 31 - y + (((sy << 5) + 31 - x) << 8)] = tiles[0][(sx << 5) + x + (((sy << 5) + y) << 8)];
					tiles[2][(sx << 5) + 31 - x + (((sy << 5) + 31 - y) << 8)] = tiles[0][(sx << 5) + x + (((sy << 5) + y) << 8)];
					tiles[1][(sx << 5)      + y + (((sy << 5)      + x) << 8)] = tiles[0][(sx << 5) + x + (((sy << 5) + y) << 8)];
				}
			}
		}
	}
}

void LoadLevel(WCHAR* fileName, LEVEL* level) {
	FILE* in;
	CELL* cells = level->cells;
	char buf[40];
	_wfopen_s(&in, fileName, L"r");
	for (int y = 0; y < 9; y++) {
		fgets(buf, 40, in);
		for (int x = 0; x < 15; x++) {
			cells[x + (15 * y)].type = buf[2 * x];
			cells[x + (15 * y)].state = buf[2 * x + 1];
		}
	}
	fgets(buf, 40, in);
	level->ox = buf[0] - 'a';
	level->oy = buf[1] - 'a';

	fgets(buf, 40, in);
	level->controlX = buf[0] - 'a';
	level->controlY = buf[1] - 'a';

	fgetws(level->nextLevel, 256, in);

	fclose(in);

	char* qcells = level->qcells;
	memset(qcells, 0, 20 * 32);
	int quantCount = 0;

	int pd[4]; pd[0] = 0; pd[1] = 0; pd[2] = 0; pd[3] = 0;
	int d;

	for (int y = 0; y < 9; y++) {
		for (int x = 0; x < 15; x++) {
			int idx = x + (15 * y);
			char type = cells[idx].type;
			char state = cells[idx].state;
			switch (type) {
				case PORTAL:
					state = (char)(state - 'a');
					cells[idx].state = state;
					d = pd[state];
					level->px[state * 2 + d] = x;
					level->py[state * 2 + d] = y;
					pd[state]++;
				break;

				case EMITTER:
					level->ex = x;
					level->ey = y;
					if (state == LEFT) {
						cells[idx].state = 0; level->eDir = 0;
					} else if (state == RIGHT) {
						cells[idx].state = 2; level->eDir = 8;
					} else if (state == UP) {
						cells[idx].state = 1; level->eDir = 4;
					} else {
						cells[idx].state = 3; level->eDir = 12;
					}
				break;

				case QUANT:
					quantCount++;
				break;

				case COLLECTOR:
					level->cx = x;
					level->cy = y;
					if (state == LEFT) {
						cells[idx].state = 0; level->cDir =  8;
					} else if (state == RIGHT) {
						cells[idx].state = 2; level->cDir =  0;
					} else if (state == UP) {
						cells[idx].state = 1; level->cDir = 12;
					} else {
						cells[idx].state = 3; level->cDir =  4;
					}
				break;

				case MIRROR:
				case POLAR:
				case MPOLAR:
				case ROTOR:
				case ROLAR:
				case MROLAR:
					cells[idx].state = (char)(state - 'a');
				break;

				case BLOCK:
				case MBLOCK:
					//cells[idx].type = 0;
					switch (state) {
						case LEFT:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
						break;

						case FULL:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case RIGHT:
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case UP:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
						break;

						case DOWN:
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case LEFT_UP:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
						break;

						case RIGHT_UP:
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
						break;

						case LEFT_DOWN:
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
						break;

						case RIGHT_DOWN:
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case EMPTY_LEFT_DOWN:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case EMPTY_LEFT_UP:
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;

						case EMPTY_RIGHT_DOWN:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 1] = type;
						break;

						case EMPTY_RIGHT_UP:
							qcells[33 + (2 * y * 32) + (2 * x)] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 32] = type;
							qcells[33 + (2 * y * 32) + (2 * x) + 33] = type;
						break;
					}
				break;
			}
		}
	}
	level->quantCount = quantCount;
}

LPDWORD LoadPng(LPCWSTR fileName, int w, int h) {
  LPDWORD result = new DWORD[w * h];

  png_struct    *png_ptr = NULL;
  png_info	*info_ptr = NULL;
  png_byte      buf[8];
  png_byte      *png_pixels = NULL;
  png_byte      **row_pointers = NULL;
  png_byte      *pix_ptr = NULL;
  png_uint_32   row_bytes;

  png_uint_32   width;
  png_uint_32   height;
  int           bit_depth;
  int           channels;
  int           color_type;
  int           alpha_present;
  int           row, col;
  int           ret;
  int           i;
  long          dep_16;

  FILE* png_file;

 _wfopen_s(&png_file, fileName, L"rb");

  /* read and check signature in PNG file */
  ret = fread (buf, 1, 8, png_file);
  if (ret != 8) {
	  fclose(png_file);
  fwprintf_s(debugFile,L"e1 %d %d %d %s", ret, ferror(png_file), feof(png_file),fileName);
  return result;}

  ret = png_check_sig (buf, 8);
  if (!ret) {fclose(png_file);
  fwprintf_s(debugFile,L"e2");
  return result;}

  /* create png and info structures */
  png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {fclose(png_file);
  fwprintf_s(debugFile,L"e3");
  return result; /* out of memory */}

  info_ptr = png_create_info_struct (png_ptr);
  if (!info_ptr) {
    png_destroy_read_struct (&png_ptr, NULL, NULL);
    fclose(png_file);
  fwprintf_s(debugFile,L"e4");
    return result; /* out of memory */
  }

  if (setjmp (png_jmpbuf(png_ptr))) {
    png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
    fclose(png_file);
  fwprintf_s(debugFile,L"e5");
    return result;
  }

  /* set up the input control for C streams */
  png_init_io (png_ptr, png_file);
  png_set_sig_bytes (png_ptr, 8);  /* we already read the 8 signature bytes */

  /* read the file information */
  png_read_info (png_ptr, info_ptr);

  /* get size and bit-depth of the PNG-image */
  png_get_IHDR (png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, NULL, NULL, NULL);

  /* set-up the transformations */

  /* transform paletted images into full-color rgb */
  if (color_type == PNG_COLOR_TYPE_PALETTE) {png_set_expand (png_ptr);}
  /* expand images to bit-depth 8 (only applicable for grayscale images) */
  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {png_set_expand (png_ptr);}
  /* transform transparency maps into full alpha-channel */
  if (png_get_valid (png_ptr, info_ptr, PNG_INFO_tRNS)) {png_set_expand (png_ptr);}

  /* all transformations have been registered; now update info_ptr data,
   * get rowbytes and channels, and allocate image memory */

  png_read_update_info (png_ptr, info_ptr);

  /* get the new color-type and bit-depth (after expansion/stripping) */
  png_get_IHDR (png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, NULL, NULL, NULL);

  /* calculate new number of channels and store alpha-presence */
  if (color_type == PNG_COLOR_TYPE_GRAY) {channels = 1;}
  else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {channels = 2;}
  else if (color_type == PNG_COLOR_TYPE_RGB) {channels = 3;}
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {channels = 4;}
  else {channels = 0; /* should never happen */}
  alpha_present = (channels - 1) % 2;

  /* check if alpha is expected to be present in file */

  /* row_bytes is the width x number of channels x (bit-depth / 8) */
  row_bytes = png_get_rowbytes (png_ptr, info_ptr);

  if ((png_pixels = (png_byte *) malloc (row_bytes * height * sizeof (png_byte))) == NULL) {
    png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
    fclose(png_file);
  fwprintf_s(debugFile,L"e6");
    return result;
  }

  if ((row_pointers = (png_byte **) malloc (height * sizeof (png_bytep))) == NULL) {
    png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
    free (png_pixels);
    png_pixels = NULL;
    fclose(png_file);
  fwprintf_s(debugFile,L"e7");
    return result;
  }

  /* set the individual row_pointers to point at the correct offsets */
  for (i = 0; i < (height); i++) {row_pointers[i] = png_pixels + i * row_bytes;}

  /* now we can go ahead and just read the whole image */
  png_read_image (png_ptr, row_pointers);

  /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
  png_read_end (png_ptr, info_ptr);

  /* clean up after the read, and free any memory allocated - REQUIRED */
  png_destroy_read_struct (&png_ptr, &info_ptr, (png_infopp) NULL);
  fclose(png_file);

  /* write data to PNM file */
  pix_ptr = png_pixels;

  fwprintf_s(debugFile,L"%d %d %d", height, width, channels);
  for (row = 0; row < height; row++) {
    for (col = 0; col < width; col++) {
      for (i = 0; i < (channels - alpha_present); i++) {((LPBYTE)result)[(((row * width) + col) * 4) + 2 - i] = (int) *pix_ptr++;}
    } /* end for col */
  } /* end for row */

  if (row_pointers != (unsigned char**) NULL) {free (row_pointers);}
  if (png_pixels != (unsigned char*) NULL) {free (png_pixels);}
  return result;
}
