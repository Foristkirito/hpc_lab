// bitmap.cpp :
//

#include "stdafx.h"


// ------------------------------------------------------------------------
static void _cdecl ErrorMsg (DWORD LastError, char * szFormat, ...)
{
	char MsgBuf[1024];

	int chars = wvsprintf(MsgBuf, szFormat, (char *)(&szFormat+1));
	strcat(MsgBuf, "\n");
	if (LastError != 0)
		FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL,
			LastError,
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
			MsgBuf+chars+1,
			1024-chars-1,
			NULL);

	MessageBox( NULL, MsgBuf, "Fatal Error", MB_OK|MB_ICONEXCLAMATION|MB_TASKMODAL);
}


// ------------------------------------------------------------------------
DWORD MyLoadBitmap32(void **ppBitmap, BITMAPINFOHEADER *pBIH, char filename[])
{
	HANDLE hFile;
	BITMAPFILEHEADER BFH;
	BITMAPINFOHEADER BIH;
	DWORD dwBytesRead;
	int i, BitmapSize, NewBitmapSize, BytesPerRow;
	BYTE *pData = NULL;
	DWORD *pData32 = NULL;
	DWORD RetVal = 0xffffffff;

	*ppBitmap = NULL;
	ZeroMemory(pBIH, sizeof (BITMAPINFOHEADER));

	hFile = CreateFile (filename, GENERIC_READ, FILE_SHARE_READ,
						NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE)
	{
		RetVal = GetLastError();
		ErrorMsg(RetVal, "Failed opening file '%s'.", filename);
		goto end;
	}

	if (ReadFile(hFile, &BFH, sizeof(BITMAPFILEHEADER), &dwBytesRead, NULL) == 0)
	{
		RetVal = GetLastError();
		ErrorMsg(RetVal, "Failed reading bitmap file header in file '%s'.", filename);
		goto end;
	}
	if (dwBytesRead != sizeof(BITMAPFILEHEADER))
	{
		ErrorMsg(0, "Failed reading bitmap file header in file '%s', only read %d bytes, asked for %d.", filename, dwBytesRead, sizeof(BITMAPFILEHEADER));
		goto end;
	}

	if (BFH.bfType != 0x4d42) // "BM"
	{
		ErrorMsg(0, "'%s' is not a valid bitmap.", filename);
		goto end;
	}

	if (ReadFile(hFile, &BIH, sizeof(BITMAPINFOHEADER), &dwBytesRead, NULL) == 0)
	{
		RetVal = GetLastError();
		ErrorMsg(RetVal, "Failed reading bitmap info header in file '%s'.", filename);
		goto end;
	}
	if (dwBytesRead != sizeof(BITMAPINFOHEADER))
	{
		ErrorMsg(0, "Failed reading bitmap info header in file '%s', only read %d bytes, asked for %d.", filename, dwBytesRead, sizeof(BITMAPINFOHEADER));
		goto end;
	}

	if (BIH.biSize < sizeof(BITMAPINFOHEADER))
	{
		ErrorMsg(0, "Unexpected Bitmap info header 'biSize' field.\nExpected %d bytes, got %d.", sizeof(BITMAPINFOHEADER), BIH.biSize);
		goto end;
	}

	if (BIH.biPlanes != 1)
	{
		ErrorMsg(0, "Bitmap info header claims %d planes, can only read 1 plane.", BIH.biPlanes);
		goto end;
	}

	if (BIH.biBitCount != 24 && BIH.biBitCount != 32)
	{
		ErrorMsg(0, "Can only read 24 and 32 bits per pixel images. Bitmap info header claims a %d bit count.", BIH.biBitCount);
		goto end;
	}

	if (BIH.biCompression != 0)
	{
		ErrorMsg(0, "Bitmap is compressed.");
		goto end;
	}

	if (BFH.bfOffBits != sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER))
	{
		ErrorMsg(0, "Unusual offset of bits. Expecting %d, bitmap file header claims %d.", sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER), BFH.bfOffBits);
		goto end;
	}


	BitmapSize = BIH.biWidth * BIH.biHeight * BIH.biBitCount / 8;
	pData = (BYTE *)malloc(BitmapSize);
	if (pData == NULL)
	{
		ErrorMsg(0, "Failed to malloc %d bytes for the data.", BitmapSize);
		goto end;
	}

	BytesPerRow = BIH.biWidth * BIH.biBitCount / 8;
   	pData += BIH.biHeight * BytesPerRow;

   	// read the data right side up
   	for (i = 0; i < BIH.biHeight; i ++)
   	{
		pData -= BytesPerRow;
		if (ReadFile(hFile, pData, BytesPerRow, &dwBytesRead, NULL) == 0)
		{
			RetVal = GetLastError();
			ErrorMsg(RetVal, "Failed to read row %d of %d of the bitmap data.", i, BIH.biHeight);
			goto end;
		}
	}

	// check if color conversion is required
	if (BIH.biBitCount == 32)
	{
		*ppBitmap = pData;
	}
	else if (BIH.biBitCount == 24)
	{
		NewBitmapSize = BIH.biWidth * BIH.biHeight * 4;
		pData32 = (DWORD *)malloc(NewBitmapSize);
		if (pData32 == NULL)
		{
			ErrorMsg(0, "Failed to malloc %d bytes for the 32 bit per pixel data.", NewBitmapSize);
			goto end;
		}
		for (i=0; i<BIH.biWidth * BIH.biHeight; i++)
		{
			pData32[i] = (DWORD)((pData[i*3]) | (pData[i*3+1] << 8) | (pData[i*3+2]) << 16);
		}
		*ppBitmap = pData32;
		free (pData);
		BIH.biBitCount = 32;
		BIH.biSizeImage = NewBitmapSize;
	}


	CopyMemory(pBIH, &BIH, sizeof(BITMAPINFOHEADER));
	RetVal = 0;


end:
	if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
	if (RetVal != 0)
	{
		if (pData != NULL) free (pData);
	}

	return RetVal;
}
