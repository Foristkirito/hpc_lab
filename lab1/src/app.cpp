// app.cpp : Contains generic Windows functions and calls AppInit, AppIdle, and AppTerm
//

#include "stdafx.h"
#include <omp.h>
#include <intrin.h>

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
void Rotate(DWORD *pDst, long DPitch, DWORD *pSrc, long SPitch, long width, long  height, float angle)
{
	DPitch /= sizeof(DWORD);
	SPitch /= sizeof(DWORD);
	float sin_angle = sin(angle);
	float cos_angle = cos(angle);

	for (int y = -height / 2; y<height / 2; y++)
	{
		float tmp_y_sin = width / 2.0 - y * sin_angle;
		float tmp_y_cos = height / 2.0 + y * cos_angle;
		#pragma omp parallel for
		for (int x = -width / 2; x<width / 2; x++)
		{
			float fSrcX = (float)(x*cos_angle + tmp_y_sin);
			float fSrcY = (float)(x*sin_angle + tmp_y_cos);

			long dwSrcX = (long)fSrcX;
			long dwSrcY = (long)fSrcY;

			if (dwSrcX > 0 && dwSrcY > 0 && dwSrcX < width - 1 && dwSrcY < height - 1)
			{
				DWORD *pTopLeft, *pTopRight, *pBottomLeft, *pBottomRight;
				pTopLeft = pSrc + dwSrcY * SPitch + dwSrcX;
				pTopRight = pTopLeft + 1;
				pBottomLeft = pTopLeft + SPitch;
				pBottomRight = pTopLeft + SPitch + 1;

				float fx = fSrcX - dwSrcX;
				float fy = fSrcY - dwSrcY;
				float fxy = fx * fy;

				DWORD topLeftValue = (*pTopLeft);
				DWORD topRightValue = (*pTopRight);
				DWORD bottomLeftValue = (*pBottomLeft);
				DWORD bottomRightValue = (*pBottomRight);

				DWORD TopLeftAlpha = topLeftValue >> 24;
				DWORD TopLeftRed = (topLeftValue >> 16) & 0xff;
				DWORD TopLeftGreen = (topLeftValue >> 8) & 0xff;
				DWORD TopLeftBlue = topLeftValue & 0xff;

				//top right

				DWORD TopRightAlpha = topRightValue >> 24;
				DWORD TopRightRed = (topRightValue >> 16) & 0xff;
				DWORD TopRightGreen = (topRightValue >> 8) & 0xff;
				DWORD TopRightBlue = topRightValue & 0xff;


				//bottom left

				DWORD BottomLeftAlpha = bottomLeftValue >> 24;
				DWORD BottomLeftRed = (bottomLeftValue >> 16) & 0xff;
				DWORD BottomLeftGreen = (bottomLeftValue >> 8) & 0xff;
				DWORD BottomLeftBlue = bottomLeftValue & 0xff;


				//bottom right

				DWORD BottomRightAlpha = bottomRightValue >> 24;
				DWORD BottomRightRed = (bottomRightValue >> 16) & 0xff;
				DWORD BottomRightGreen = (bottomRightValue >> 8) & 0xff;
				DWORD BottomRightBlue = bottomRightValue & 0xff;


				float tmp_1 = 1.0f - fx - fy + fxy;
				float tmp_2 = fx - fxy;
				float tmp_3 = fy - fxy;

				// 使用向量化

				__m128 tmp;
				__m128 alphaV;
				__m128 redV;
				__m128 greenV;
				__m128 blueV;
				tmp = _mm_set_ps(tmp_1, tmp_2, tmp_3, fxy);
				alphaV = _mm_set_ps(TopLeftAlpha, TopRightAlpha, BottomLeftAlpha, BottomRightAlpha);
				redV = _mm_set_ps(TopLeftRed, TopRightRed, BottomLeftRed, BottomRightRed);
				greenV = _mm_set_ps(TopLeftGreen, TopRightGreen, BottomLeftGreen, BottomRightGreen);
				blueV = _mm_set_ps(TopLeftBlue, TopRightBlue, BottomLeftBlue, BottomRightBlue);
				__m128 lastSum = _mm_set_ps(0.5f, 0.5f, 0.5f, 0.5f);
				__m128 tmpAlpha = _mm_mul_ps(alphaV, tmp);
				tmpAlpha = _mm_hadd_ps(tmpAlpha, tmpAlpha);
				tmpAlpha = _mm_hadd_ps(tmpAlpha, tmpAlpha);
				tmpAlpha = _mm_add_ps(tmpAlpha, lastSum);
				int alpha_value = (int)(tmpAlpha.m128_f32[0]);
				__m128 tmpRed = _mm_mul_ps(redV, tmp);
				tmpRed = _mm_hadd_ps(tmpRed, tmpRed);
				tmpRed = _mm_hadd_ps(tmpRed, tmpRed);
				tmpRed = _mm_add_ps(tmpRed, lastSum);
				int red_value = (int)(tmpRed.m128_f32[0]);

				__m128 tmpGreen = _mm_mul_ps(greenV, tmp);
				tmpGreen = _mm_hadd_ps(tmpGreen, tmpGreen);
				tmpGreen = _mm_hadd_ps(tmpGreen, tmpGreen);
				tmpGreen = _mm_add_ps(tmpGreen, lastSum);
				int green_value = (int)(tmpGreen.m128_f32[0]);

				__m128 tmpBlue = _mm_mul_ps(blueV, tmp);
				tmpBlue = _mm_hadd_ps(tmpBlue, tmpBlue);
				tmpBlue = _mm_hadd_ps(tmpBlue, tmpBlue);
				tmpBlue = _mm_add_ps(tmpBlue, lastSum);
				int blue_value = (int)(tmpBlue.m128_f32[0]);

				pDst[x + width / 2] = (alpha_value << 24) | (red_value << 16) |
					(green_value << 8) | (blue_value);
			}
			else {
				pDst[x + width / 2] = 0;
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
					DWORD Blue  = (sBlue * sAlpha + dBlue * dAlpha) >> 8;

					DWORD sGreen = (sColor >> 8) & 0xff;
					DWORD dGreen = (dColor >> 8) & 0xff;
					DWORD Green  = (sGreen * sAlpha + dGreen * dAlpha) >> 8;

					DWORD sRed = (sColor >> 16) & 0xff;
					DWORD dRed = (dColor >> 16) & 0xff;
					DWORD Red  = (sRed * sAlpha + dRed * dAlpha) >> 8;

					pDst[col] = (Red << 16) + (Green << 8) + Blue;
				}
			}
		}
		pSrc += sPitch >> 2;
		pDst += dPitch >> 2;
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
