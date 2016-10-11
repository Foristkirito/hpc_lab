// app.cpp : Contains generic Windows functions and calls AppInit, AppIdle, and AppTerm
//

#include "stdafx.h"

// background image
#define BKGBMP_FILENAME "background.bmp"
DWORD *pBkgBmp = NULL;
DWORD BkgBmpWidth, BkgBmpHeight, BkgBmpPitch;

// foreground image
#define FGBMP_FILENAME "foreground.bmp"
DWORD *pFgBmp = NULL;
DWORD FgBmpWidth, FgBmpHeight, FgBmpPitch;
POINT FgPos = { 0, 0};
POINT FgDir = {+2, 0};
DWORD *pRotatedImage = NULL;	// same size as foreground
float Angle = 0.0f, AngleDir = -0.04f;

//---------------------------------------------------------------------------
void UpdatePositionAndRotation (void)
{
	Angle += AngleDir;

	FgPos.x += FgDir.x;

	if (FgPos.x < 0)
	{
		FgPos.x = 0;
		FgDir.x = -FgDir.x;
		AngleDir *= -1;
	} else if (FgPos.x + FgBmpWidth >= (int)BkgBmpWidth)
	{
		FgPos.x = BkgBmpWidth - FgBmpWidth;
		FgDir.x = -FgDir.x;
		AngleDir *= -1;
	}
}


// --------------------------------------------------------------------------
void Rotate (DWORD *pDst, long DPitch, DWORD *pSrc, long SPitch, long width, long  height, float angle)
{
	DPitch /= sizeof(DWORD);
	SPitch /= sizeof(DWORD);

	for (int y=-height/2; y<height/2; y++)
	{
		for (int x=-width/2; x<width/2; x++)
		{
			float fSrcX = (float)(width /2.0 + x*cos(angle) - y*sin(angle));
			float fSrcY = (float)(height/2.0 + x*sin(angle) + y*cos(angle));

			long dwSrcX = (long)fSrcX;
			long dwSrcY = (long)fSrcY;
			
			if (dwSrcX > 0 && dwSrcY > 0 && dwSrcX < width-1 && dwSrcY < height-1)
			{
				DWORD *pTopLeft, *pTopRight, *pBottomLeft, *pBottomRight;
				pTopLeft		= pSrc + dwSrcY * SPitch + dwSrcX;
				pTopRight		= pTopLeft + 1;
				pBottomLeft		= pTopLeft + SPitch;
				pBottomRight	= pBottomLeft + 1;

				float fx = fSrcX - dwSrcX;
				float fy = fSrcY - dwSrcY;

				// alpha interpolation
				DWORD TopLeftAlpha     = (*pTopLeft)    >> 24;
				DWORD TopRightAlpha    = (*pTopRight)   >> 24;
				DWORD BottomLeftAlpha  = (*pBottomLeft) >> 24;
				DWORD BottomRightAlpha = (*pBottomRight)>> 24;
				float alphaFP = TopLeftAlpha	* (1.0f-fx) * (1.0f-fy) +
								TopRightAlpha	* fx * (1.0f-fy) +
								BottomLeftAlpha	* (1.0f-fx) * fy +
								BottomRightAlpha* fx * fy;
				DWORD  alpha_value = (DWORD)(alphaFP + 0.5);

				// red interpolation
				DWORD TopLeftRed     = ((*pTopLeft)    >> 16) & 0xff;
				DWORD TopRightRed    = ((*pTopRight)   >> 16) & 0xff;
				DWORD BottomLeftRed  = ((*pBottomLeft) >> 16) & 0xff;
				DWORD BottomRightRed = ((*pBottomRight)>> 16) & 0xff;
				float redFP =	TopLeftRed		* (1.0f-fx) * (1.0f-fy) +
								TopRightRed		* fx * (1.0f-fy) +
								BottomLeftRed	* (1.0f-fx) * fy +
								BottomRightRed	* fx * fy;
				DWORD  red_value = (DWORD)(redFP + 0.5);
				
				// green interpolation
				DWORD TopLeftGreen     = ((*pTopLeft)    >> 8) & 0xff;
				DWORD TopRightGreen    = ((*pTopRight)   >> 8) & 0xff;
				DWORD BottomLeftGreen  = ((*pBottomLeft) >> 8) & 0xff;
				DWORD BottomRightGreen = ((*pBottomRight)>> 8) & 0xff;
				float greenFP = TopLeftGreen	  * (1.0f-fx) * (1.0f-fy) +  
								 TopRightGreen	  * fx * (1.0f-fy) +         
								 BottomLeftGreen  * (1.0f-fx) * fy +         
								 BottomRightGreen * fx * fy;                 
				DWORD  green_value = (DWORD)(greenFP + 0.5);

				// blue interpolation
				DWORD TopLeftBlue     = *pTopLeft & 0xff;
				DWORD TopRightBlue    = *pTopRight & 0xff;
				DWORD BottomLeftBlue  = *pBottomLeft & 0xff;
				DWORD BottomRightBlue = *pBottomRight & 0xff;
				float blueFP = TopLeftBlue		* (1.0f-fx) * (1.0f-fy) +  
							   TopRightBlue		* fx * (1.0f-fy) +         
							   BottomLeftBlue	* (1.0f-fx) * fy +         
							   BottomRightBlue	* fx * fy;                 
				DWORD  blue_value = (DWORD)(blueFP + 0.5);

				pDst[x+width/2] = (alpha_value << 24) + (red_value << 16) + 
								  (green_value << 8)  + (blue_value);
			} else {
				pDst[x+width/2] = 0;
			}
		}
		pDst += DPitch;
	}
}
// --------------------------------------------------------------------------
void AlphaBlend32to32(DWORD *pDst, long dPitch,
				DWORD *pSrc, long sPitch,
				DWORD width, DWORD height)
{
	for (DWORD row=0; row < height; row++)
	{
		for (DWORD col=0; col < width; col++)
		{
			DWORD sColor = pSrc[col];
			DWORD sAlpha = sColor >> 24;
			if (sAlpha != 0)
			{
				if (sAlpha == 255)
					pDst[col] = sColor;
				else
				{
					DWORD dColor = pDst[col];
					DWORD dAlpha = 255 - sAlpha;

					DWORD sBlue = sColor & 0xff;
					DWORD dBlue = dColor & 0xff;
					DWORD Blue  = (sBlue * sAlpha + dBlue * dAlpha) / 256;

					DWORD sGreen = (sColor >> 8) & 0xff;
					DWORD dGreen = (dColor >> 8) & 0xff;
					DWORD Green  = (sGreen * sAlpha + dGreen * dAlpha) / 256;

					DWORD sRed = (sColor >> 16) & 0xff;
					DWORD dRed = (sColor >> 16) & 0xff;
					DWORD Red  = (sRed * sAlpha + dRed * dAlpha) / 256;

					pDst[col] = (Red << 16) + (Green << 8) + Blue;
				}
			}
		}
		pSrc += sPitch / 4;
		pDst += dPitch / 4;
	}
}

// --------------------------------------------------------------------------
void MemCopyRect(BYTE *pDst, long dPitch, BYTE *pSrc, long sPitch, 
		DWORD width, DWORD height)
{
	for (DWORD row=0; row < height; row ++)
	{
		CopyMemory(pDst, pSrc, width);
		pDst += dPitch;
		pSrc += sPitch;
	}
}

// ------------------------------------------------------------------------
int AppInit (HWND hWnd)
{
	int RetVal;

	// load the bitmap for the background
	BITMAPINFOHEADER BIH;
	RetVal = MyLoadBitmap32((void**)&pBkgBmp, &BIH, BKGBMP_FILENAME);
	if (RetVal != 0)
		return -1;
	BkgBmpWidth	 = BIH.biWidth;
	BkgBmpHeight = BIH.biHeight;
	BkgBmpPitch	 = BIH.biWidth*4;

	// adjust the window to match the size of the bitmap
	RECT ClientRect;
	SetRect(&ClientRect, 0, 0, BIH.biWidth, BIH.biHeight);
	AdjustWindowRect(&ClientRect, WS_CAPTION | WS_SYSMENU, FALSE);
	SetWindowPos(hWnd, 0, 0, 0, ClientRect.right - ClientRect.left, ClientRect.bottom - ClientRect.top, SWP_NOMOVE | SWP_NOZORDER);

	// load the bitmap for the foreground
	RetVal = MyLoadBitmap32((void **)&pFgBmp, &BIH, FGBMP_FILENAME);
	if (RetVal != 0)
	{
		free (pBkgBmp);
		return -1;
	}
	FgBmpWidth  = BIH.biWidth;
	FgBmpHeight = BIH.biHeight;
	FgBmpPitch  = BIH.biWidth*4;

	// create the bitmap for the rotated foreground image
	pRotatedImage = (DWORD *)malloc(FgBmpWidth * FgBmpHeight * sizeof(DWORD));
	if (pRotatedImage == NULL)
	{
		free (pBkgBmp);
		return -1;
	}

	// place the ball vertically centered 
	FgPos.y = BkgBmpHeight / 2 - FgBmpHeight / 2;

	RetVal = dd_Init(hWnd, BkgBmpWidth, BkgBmpHeight);


	return RetVal;
}


// ------------------------------------------------------------------------
void AppIdle (HWND hWnd)
{
	BYTE *pBackBuf;
	DWORD BackBufWidth, BackBufHeight;
	long  BackBufPitch;

	DWORD StartTime, ElapTime;
	_asm {
		RDTSC
		mov StartTime, eax
	}

	UpdatePositionAndRotation();   // updates FgPos.x and Angle

	// rotates pFgBmp by Angle radians and places result in pRotatedImage
	Rotate(pRotatedImage, FgBmpPitch, pFgBmp, FgBmpPitch, FgBmpWidth, FgBmpHeight, Angle);

	// retrieve a write pointer to the off-screen surface from DirectDraw
	HRESULT dd_RetVal = dd_LockSurface(&pBackBuf, &BackBufWidth, &BackBufHeight, &BackBufPitch);
	if(dd_RetVal != DD_OK) 
		return;
	// erase the old image by copying the pBkgBmp to the off-screen surface
	MemCopyRect(pBackBuf, BackBufPitch, (BYTE*)pBkgBmp, BkgBmpPitch, BkgBmpWidth*4, BkgBmpHeight);

	// calculate the offset of the upper left pixel of the foreground image on the off-screen surface
	DWORD *pBallPosOnSurface = (DWORD*)(pBackBuf + FgPos.x*4 + FgPos.y*BackBufPitch);

	// blend the rotated foreground image onto the background image
	// note: the blend process reads from the destination
	AlphaBlend32to32((DWORD *)pBallPosOnSurface, BackBufPitch, pRotatedImage, FgBmpPitch, FgBmpWidth, FgBmpHeight);

	dd_UnlockSurface(pBackBuf);

	// copy the off-screen surface to the visible client window
	dd_BltBackToPrimaryFull(hWnd);

	_asm {
		RDTSC
		sub eax, StartTime
		mov ElapTime, eax
	}
	static DWORD MinElapsedTime = 0xffffffff;
	if (MinElapsedTime > ElapTime)
	{
		MinElapsedTime = ElapTime;
		char WinTitle[80];
		wsprintf(WinTitle, "%d ", MinElapsedTime);
		SetWindowText(hWnd, WinTitle);
	}
}


// ------------------------------------------------------------------------
void AppTerm (void)
{
	if (pBkgBmp != NULL) free (pBkgBmp);
	if (pFgBmp != NULL) free (pFgBmp);
	if (pRotatedImage != NULL) free (pRotatedImage);
	
	dd_Term();
}
