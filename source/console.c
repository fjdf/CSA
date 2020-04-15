#include "console.h"

#if defined(WIN32) || defined(_WIN32)

#include <Windows.h>

void ConsoleSetTextColor( short int textAttribute, short int foregroundTextColor, short int backgroundTextColor){
	SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE) , (WORD)( (textAttribute) | (foregroundTextColor) | (backgroundTextColor << 4) ) );
}

void ConsoleResetTextColor(){
	SetConsoleTextAttribute( GetStdHandle(STD_OUTPUT_HANDLE) , (WORD) COLOR_RESET );
}

void ConsoleMoveCursorPosition( short int horizontalShift, short int verticalShift){
	HANDLE hConsoleOutput;
	COORD dwCursorPosition;
	CONSOLE_SCREEN_BUFFER_INFO lpConsoleScreenBufferInfo;
	hConsoleOutput = GetStdHandle( STD_OUTPUT_HANDLE );
	GetConsoleScreenBufferInfo( hConsoleOutput , &lpConsoleScreenBufferInfo );
	dwCursorPosition = lpConsoleScreenBufferInfo.dwCursorPosition;
	dwCursorPosition.X += (SHORT)horizontalShift;
	dwCursorPosition.Y += (SHORT)verticalShift;
	SetConsoleCursorPosition( hConsoleOutput , dwCursorPosition );
}

void ConsoleClearScreen(){
	HANDLE hConsoleOutput;
	CONSOLE_SCREEN_BUFFER_INFO lpConsoleScreenBufferInfo;
	COORD dwSize;
	DWORD nLength;
	COORD dwWriteCoord;
	DWORD lpNumberOfCharsWritten;
	hConsoleOutput = GetStdHandle( STD_OUTPUT_HANDLE );
	GetConsoleScreenBufferInfo( hConsoleOutput , &lpConsoleScreenBufferInfo );
	dwSize = lpConsoleScreenBufferInfo.dwSize;
	nLength = (DWORD)( dwSize.X * dwSize.Y );
	dwWriteCoord.X = (SHORT) 0;
	dwWriteCoord.Y = (SHORT) 0;
	FillConsoleOutputCharacter( hConsoleOutput, (TCHAR) ' ' , nLength , dwWriteCoord , &lpNumberOfCharsWritten );
	FillConsoleOutputAttribute( hConsoleOutput, (WORD) COLOR_RESET , nLength , dwWriteCoord , &lpNumberOfCharsWritten );
	SetConsoleCursorPosition( hConsoleOutput , dwWriteCoord );
}

void ConsoleClearLine(){
	HANDLE hConsoleOutput;
	CONSOLE_SCREEN_BUFFER_INFO lpConsoleScreenBufferInfo;
	COORD dwSize;
	DWORD nLength;
	COORD dwWriteCoord;
	DWORD lpNumberOfCharsWritten;
	COORD dwCursorPosition;
	hConsoleOutput = GetStdHandle( STD_OUTPUT_HANDLE );
	GetConsoleScreenBufferInfo( hConsoleOutput , &lpConsoleScreenBufferInfo );
	dwSize = lpConsoleScreenBufferInfo.dwSize;
	nLength = (DWORD)( dwSize.X ); // number of columns in a row
	dwCursorPosition = lpConsoleScreenBufferInfo.dwCursorPosition;
	dwWriteCoord.X = (SHORT) 0; // first column
	dwWriteCoord.Y = dwCursorPosition.Y; // same row
	FillConsoleOutputCharacter( hConsoleOutput, (TCHAR) ' ' , nLength , dwWriteCoord , &lpNumberOfCharsWritten );
	FillConsoleOutputAttribute( hConsoleOutput, (WORD) COLOR_RESET , nLength , dwWriteCoord , &lpNumberOfCharsWritten );
	SetConsoleCursorPosition( hConsoleOutput , dwWriteCoord );
}

#else

#include <stdio.h>

void ConsoleSetTextColor( short int textAttribute, short int foregroundTextColor, short int backgroundTextColor){
	printf( "%c[%hd;%hd;%hdm" , 0x1B , textAttribute , ( 30 + foregroundTextColor ) , ( 40 + backgroundTextColor ) );
}

void ConsoleMoveCursorPosition( short int horizontalShift, short int verticalShift){
	if( verticalShift > 0 ) printf( "%c[%hdA" , 0x1B , verticalShift ); // up
	else if( verticalShift < 0 ) printf( "%c[%hdB" , 0x1B , (-verticalShift) ); // down
	if( horizontalShift > 0 ) printf( "%c[%hdC" , 0x1B , horizontalShift ); // right
	else if( horizontalShift < 0 ) printf( "%c[%hdD" , 0x1B , (-horizontalShift) ); // left
}

void ConsoleResetTextColor(){
	printf( "%c[0m" , 0x1B );
}

void ConsoleClearScreen(){
	printf( "%c[2J" , 0x1B );
}

void ConsoleClearLine(){
	printf( "%c[2K" , 0x1B );
}

#endif
