#if defined(WIN32) || defined(_WIN32)

#define COLOR_BLACK		0x0000
#define COLOR_BLUE		0x0001
#define COLOR_GREEN		0x0002
#define COLOR_CYAN		0x0003
#define COLOR_RED		0x0004
#define COLOR_MAGENTA	0x0005
#define COLOR_YELLOW	0x0006
#define COLOR_WHITE		0x0007

#define COLOR_NULL		0x0000
#define COLOR_BRIGHT	0x0008
#define COLOR_REVERSE	0x4000
#define COLOR_UNDERLINE	0x8000
#define COLOR_RESET		0x0007

#else

#define COLOR_BLACK		0
#define COLOR_RED		1
#define COLOR_GREEN		2
#define COLOR_YELLOW	3
#define COLOR_BLUE		4
#define COLOR_MAGENTA	5
#define COLOR_CYAN		6
#define COLOR_WHITE		7

#define COLOR_NULL		0
#define COLOR_BRIGHT	1
#define COLOR_REVERSE	7
#define COLOR_UNDERLINE	3
#define COLOR_RESET		0

#endif

void ConsoleSetTextColor( short int textAttribute , short int foregroundTextColor , short int backgroundTextColor );
void ConsoleMoveCursorPosition( short int horizontalShift , short int verticalShift );
void ConsoleResetTextColor();
void ConsoleClearScreen();
void ConsoleClearLine();
