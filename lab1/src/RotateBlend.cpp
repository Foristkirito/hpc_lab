// RotateStretchBlend.cpp : 
//

#include "stdafx.h"

#define WINDOW_CLASS "ROTATE_BLEND"


LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);

// ------------------------------------------------------------------------
int APIENTRY WinMain(HINSTANCE hInstance,
                     HINSTANCE hPrevInstance,
                     LPSTR     lpCmdLine,
                     int       nCmdShow)
{
	WNDCLASSEX wcex;
	wcex.cbSize			= sizeof(WNDCLASSEX); 
	wcex.style			= 0;
	wcex.lpfnWndProc	= (WNDPROC)WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= NULL; //LoadIcon(hInstance, (LPCTSTR)IDI_ROTATEBLEND);
	wcex.hCursor		= LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground	= NULL; //(HBRUSH)(COLOR_WINDOW+1);
	wcex.lpszMenuName	= NULL;
	wcex.lpszClassName	= WINDOW_CLASS;
	wcex.hIconSm		= NULL; // LoadIcon(wcex.hInstance, (LPCTSTR)IDI_SMALL);
	RegisterClassEx(&wcex);

	HWND hWnd = CreateWindow(WINDOW_CLASS, 
						"Rotate and blend sample program", 
						WS_CAPTION | WS_SYSMENU,
						CW_USEDEFAULT, 0, 
						CW_USEDEFAULT, 0, 
						NULL, NULL, hInstance, NULL);

	if (hWnd == NULL)
		return 0;

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	if (AppInit(hWnd) != 0)
		return 0;

	MSG msg;
	while (1)
	{         
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
                break;

            DispatchMessage(&msg);
        }
		else
			AppIdle(hWnd);
    }

	AppTerm();

	return msg.wParam;
}


// ------------------------------------------------------------------------
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message) 
	{
		case WM_DESTROY:
			PostQuitMessage(0);
			break;
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
   }
   return 0;
}