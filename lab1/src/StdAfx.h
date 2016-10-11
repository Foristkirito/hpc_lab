// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_)
#define AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

//#include <dvec.h>
#include <windows.h>
#include <stdlib.h>
#include <malloc.h>
#include <ddraw.h>
#include <math.h>
#include <mmsystem.h>

// in app.cpp
int AppInit (HWND hWnd);
void AppIdle (HWND hWnd);
void AppTerm (void);

// in bitmap.cpp
DWORD MyLoadBitmap32(void **ppBitmap, BITMAPINFOHEADER *pBIH, char filename[]);
DWORD MyLoadBitmap32Pow2Pitch(void **ppBitmap, DWORD * Pitch, BITMAPINFOHEADER *pBIH, char filename[]);


// in dd.cpp
HRESULT dd_Init(HWND hWnd, DWORD width, DWORD height);
void dd_Term(void);
HRESULT dd_BltBackToPrimary(HWND hWnd, RECT region);
HRESULT dd_BltBackToPrimaryFull(HWND hWnd);
HRESULT dd_LockSurface(BYTE **pBuf, DWORD *pdwWidth, DWORD *pdwHeight, long *plPitch);
void dd_UnlockSurface(BYTE *pBuf);


//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__A9DB83DB_A9FD_11D0_BFD1_444553540000__INCLUDED_)
