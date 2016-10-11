/**************************************************************************
 Optimizing C for the Pentium(r) II Processor
 XDC'97

  Routinues to deal with DirectDraw

Exported functions in this module:
HRESULT dd_Init(HWND hWnd, DWORD width, DWORD height);
void dd_Term(void);
HRESULT dd_BltBackToPrimary(HWND hWnd);
HRESULT dd_LockSurface(BYTE **pBuf, DWORD *pdwWidth, DWORD *pdwHeight, long *plPitch);
void dd_UnlockSurface(BYTE *pBuf);
**************************************************************************/
#include "stdafx.h"

static	LPDIRECTDRAW		lpDD			= NULL;
static	LPDIRECTDRAWCLIPPER	lpDDClipper		= NULL;
static	RGNDATA				*pClipList		= NULL;
static	DWORD				MaxClipListSize;
static	HINSTANCE DDLibInst;

static LPDIRECTDRAWSURFACE	lpDDSPrimary = NULL;
static LPDIRECTDRAWSURFACE	lpDDSBack	 = NULL;
static DDSURFACEDESC PrimarySurfaceDesc;

static DWORD SurfWidth  = 0;
static DWORD SurfHeight = 0;

static void dd_Error(char * pMsg, HRESULT rval);

// --------------------------------------------------------------------------
HRESULT dd_Init(HWND hWnd, DWORD width, DWORD height)
{
	HRESULT        ddrval;
	DDSURFACEDESC  ddsd;

	ddrval = DirectDrawCreate( NULL, &lpDD, NULL );
	if( ddrval != DD_OK )
	{
		dd_Error("DirectDrawCreate", ddrval);
		dd_Term();
		return ddrval;
	}

	ddrval = lpDD->SetCooperativeLevel( hWnd, DDSCL_NORMAL );
	if( ddrval != DD_OK )
	{
		dd_Error("SetCooperativeLevel", ddrval);
		dd_Term();
		return ddrval;
	}

	//  Create the primary surface
	ZeroMemory( &ddsd, sizeof( ddsd ) );
	ddsd.dwSize         = sizeof( ddsd );
	ddsd.dwFlags        = DDSD_CAPS;
	ddsd.ddsCaps.dwCaps = DDSCAPS_PRIMARYSURFACE;

	ddrval = lpDD->CreateSurface( &ddsd, &lpDDSPrimary, NULL );
	if( ddrval != DD_OK )
	{
		dd_Error("CreateSurface Primary", ddrval);
		dd_Term();
		return ddrval;
	}

	//  Create the offscreen surface
	ZeroMemory( &ddsd, sizeof( ddsd ) );
	ddsd.dwSize  = sizeof( ddsd );
	ddsd.dwFlags = DDSD_CAPS | DDSD_HEIGHT | DDSD_WIDTH;

	ddsd.ddsCaps.dwCaps = DDSCAPS_OFFSCREENPLAIN | DDSCAPS_SYSTEMMEMORY;

	ddsd.dwHeight = height;
	ddsd.dwWidth  = width;

	ddrval = lpDD->CreateSurface( &ddsd, &lpDDSBack, NULL );
	if( ddrval != DD_OK )
	{
		dd_Error("CreateSurface Offscreen", ddrval);
		dd_Term();
		return ddrval;
	}

	//  Create a clipper
	ddrval = lpDD->CreateClipper( 0, &lpDDClipper, NULL);
	if( ddrval != DD_OK )
	{
		dd_Error("IDirectDraw_CreateClipper", ddrval);
		dd_Term();
		return ddrval;
	}

	ddrval = lpDDClipper->SetHWnd( 0, hWnd);
	if( ddrval != DD_OK )
	{
		dd_Error("IDirectDrawClipper_SetHWnd", ddrval);
		dd_Term();
		return ddrval;
	}

	ddrval = lpDDSPrimary->SetClipper( lpDDClipper );
	if( ddrval != DD_OK )
	{
		dd_Error("SetClipper", ddrval);
		dd_Term();
		return ddrval;
	}

	ZeroMemory(&PrimarySurfaceDesc, sizeof(PrimarySurfaceDesc));
	PrimarySurfaceDesc.dwSize = sizeof(PrimarySurfaceDesc);
	ddrval = lpDDSPrimary->GetSurfaceDesc(&PrimarySurfaceDesc);
	if( ddrval != DD_OK )
	{
		dd_Error("GetSurfaceDesc", ddrval);
		dd_Term();
		return ddrval;
	}
	if (PrimarySurfaceDesc.ddpfPixelFormat.dwRGBBitCount != 32)
	{
		MessageBox(NULL, 
			"Application requires a 32-bit per pixel surface.\n"
			"Change your display properties.\n", 
			"Fatal Error" , MB_OK|MB_ICONEXCLAMATION|MB_TASKMODAL);
		dd_Term();
		return -1;
	}

	SurfWidth = width; SurfHeight = height;
	return DD_OK;
}


// --------------------------------------------------------------------------
void dd_Term(void)
{
    if( lpDD != NULL )
    {
        if( lpDDSPrimary != NULL )
        {
            lpDDSPrimary->Release();
            lpDDSPrimary = NULL;
        }
        if( lpDDSBack != NULL )
        {
            lpDDSBack->Release();
            lpDDSBack = NULL;
        }
        if( lpDDClipper != NULL )
        {
            lpDDClipper->Release();
            lpDDClipper = NULL;
        }

        lpDD->Release();
        lpDD = NULL;
    }
	SurfWidth = 0; SurfHeight = 0;
}


// --------------------------------------------------------------------------
HRESULT dd_BltBackToPrimary(HWND hWnd, RECT rcSrc)
{
	static BOOL PaintOK = TRUE;
	HRESULT     ddrval;
	BOOL        bChanged;

	RECT rcDst = rcSrc;

	if( !lpDDSBack )
		return DDERR_GENERIC;

	if(PaintOK == TRUE)
	{
		PaintOK = FALSE;

		ddrval = lpDDClipper->IsClipListChanged(&bChanged);
		if (bChanged)
		{
			ddrval = lpDDClipper->GetClipList(NULL, pClipList, &MaxClipListSize);
			if(ddrval!=DD_OK)
			{
				dd_Error("IDirectDrawSurface_Blt", ddrval);
				return ddrval;
			}
		}

//		GetClientRect(hWnd, &rcDst);
		ClientToScreen(hWnd, (POINT*)&rcDst);
		ClientToScreen(hWnd, ((POINT*)&rcDst)+1);

//		lpDD->WaitForVerticalBlank(DDWAITVB_BLOCKBEGIN, NULL);
		ddrval = lpDDSPrimary->Blt(&rcDst, lpDDSBack, &rcSrc, DDBLT_WAIT, 0);
		if(ddrval!=DD_OK)
		{
			dd_Error("IDirectDrawSurface_Blt", ddrval);
			return ddrval;
		}
	}
	PaintOK = TRUE;

	return DD_OK;
}

// --------------------------------------------------------------------------
HRESULT dd_BltBackToPrimaryFull(HWND hWnd)
{
	static BOOL PaintOK = TRUE;
	HRESULT     ddrval;
	BOOL        bChanged;

	RECT rcSrc = {0,0, SurfWidth, SurfHeight};
	RECT rcDst = rcSrc;

	if( !lpDDSBack )
		return DDERR_GENERIC;

	if(PaintOK == TRUE)
	{
		PaintOK = FALSE;

		ddrval = lpDDClipper->IsClipListChanged(&bChanged);
		if (bChanged)
		{
			ddrval = lpDDClipper->GetClipList(NULL, pClipList, &MaxClipListSize);
			if(ddrval!=DD_OK)
			{
				dd_Error("IDirectDrawSurface_Blt", ddrval);
				return ddrval;
			}
		}

//		GetClientRect(hWnd, &rcDst);
		ClientToScreen(hWnd, (POINT*)&rcDst);
		ClientToScreen(hWnd, ((POINT*)&rcDst)+1);

//		lpDD->WaitForVerticalBlank(DDWAITVB_BLOCKBEGIN, NULL);
		ddrval = lpDDSPrimary->Blt(&rcDst, lpDDSBack, &rcSrc, DDBLT_WAIT, 0);
		if(ddrval!=DD_OK)
		{
			dd_Error("IDirectDrawSurface_Blt", ddrval);
			return ddrval;
		}
	}
	PaintOK = TRUE;

	return DD_OK;
}

// --------------------------------------------------------------------------
HRESULT dd_LockSurface(BYTE **pBuf, DWORD *pdwWidth, DWORD *pdwHeight, long *plPitch)
{
	DDSURFACEDESC        ddsd;
	HRESULT              ddrval;

	if( !lpDDSBack )
		return DDERR_GENERIC;

	ZeroMemory(&ddsd, sizeof(ddsd));
	ddsd.dwSize   = sizeof(ddsd);
	ddsd.dwFlags  = DDSD_HEIGHT | DDSD_WIDTH;

	//  Lock the back buffer
	ddrval = lpDDSBack->Lock( NULL, &ddsd, DDLOCK_SURFACEMEMORYPTR | DDLOCK_WAIT, NULL);
	if( ddrval != DD_OK )
	{
		dd_Error("IDirectDrawSurface_Lock", ddrval);
		return ddrval;
	}

	*pBuf		= (BYTE*)ddsd.lpSurface;
	*pdwWidth	= ddsd.dwWidth;
	*pdwHeight	= ddsd.dwHeight;
	*plPitch	= ddsd.lPitch;
	return DD_OK;
}

// --------------------------------------------------------------------------
void dd_UnlockSurface(BYTE *pBuf)
{
   lpDDSBack->Unlock(pBuf);
}

// --------------------------------------------------------------------------
static void dd_Error(char * pMsg, HRESULT rval)
{
	char ErrorString[50];
	switch(rval)
	{
	case DD_OK :
		wsprintf (ErrorString, "DD_OK");
		break;
	case  DDERR_ALREADYINITIALIZED :
		wsprintf (ErrorString, "DDERR_ALREADYINITIALIZED ");
		break;
	case  DDERR_BLTFASTCANTCLIP :
		wsprintf (ErrorString, "DDERR_BLTFASTCANTCLIP ");
		break;
	case  DDERR_CANNOTATTACHSURFACE :
		wsprintf (ErrorString, "DDERR_CANNOTATTACHSURFACE ");
		break;
	case  DDERR_CANNOTDETACHSURFACE :
		wsprintf (ErrorString, "DDERR_CANNOTDETACHSURFACE ");
		break;
	case  DDERR_CANTCREATEDC :
		wsprintf (ErrorString, "DDERR_CANTCREATEDC ");
		break;
	case  DDERR_CANTDUPLICATE :
		wsprintf (ErrorString, "DDERR_CANTDUPLICATE ");
		break;
	case  DDERR_CANTPAGELOCK :
		wsprintf (ErrorString, "DDERR_CANTPAGELOCK ");
		break;
	case  DDERR_CANTPAGEUNLOCK :
		wsprintf (ErrorString, "DDERR_CANTPAGEUNLOCK ");
		break;
	case  DDERR_CLIPPERISUSINGHWND :
		wsprintf (ErrorString, "DDERR_CLIPPERISUSINGHWND ");
		break;
	case  DDERR_COLORKEYNOTSET :
		wsprintf (ErrorString, "DDERR_COLORKEYNOTSET ");
		break;
	case  DDERR_CURRENTLYNOTAVAIL :
		wsprintf (ErrorString, "DDERR_CURRENTLYNOTAVAIL ");
		break;
	case  DDERR_DCALREADYCREATED :
		wsprintf (ErrorString, "DDERR_DCALREADYCREATED ");
		break;
	case  DDERR_DIRECTDRAWALREADYCREATED :
		wsprintf (ErrorString, "DDERR_DIRECTDRAWALREADYCREATED ");
		break;
	case  DDERR_EXCEPTION :
		wsprintf (ErrorString, "DDERR_EXCEPTION ");
		break;
	case  DDERR_EXCLUSIVEMODEALREADYSET :
		wsprintf (ErrorString, "DDERR_EXCLUSIVEMODEALREADYSET ");
		break;
	case  DDERR_GENERIC :
		wsprintf (ErrorString, "DDERR_GENERIC ");
		break;
	case  DDERR_HEIGHTALIGN :
		wsprintf (ErrorString, "DDERR_HEIGHTALIGN ");
		break;
	case  DDERR_HWNDALREADYSET :
		wsprintf (ErrorString, "DDERR_HWNDALREADYSET ");
		break;
	case  DDERR_HWNDSUBCLASSED :
		wsprintf (ErrorString, "DDERR_HWNDSUBCLASSED ");
		break;
	case  DDERR_IMPLICITLYCREATED :
		wsprintf (ErrorString, "DDERR_IMPLICITLYCREATED ");
		break;
	case  DDERR_INCOMPATIBLEPRIMARY :
		wsprintf (ErrorString, "DDERR_INCOMPATIBLEPRIMARY ");
		break;
	case  DDERR_INVALIDCAPS :
		wsprintf (ErrorString, "DDERR_INVALIDCAPS ");
		break;
	case  DDERR_INVALIDCLIPLIST :
		wsprintf (ErrorString, "DDERR_INVALIDCLIPLIST ");
		break;
	case  DDERR_INVALIDDIRECTDRAWGUID :
		wsprintf (ErrorString, "DDERR_INVALIDDIRECTDRAWGUID ");
		break;
	case  DDERR_INVALIDMODE :
		wsprintf (ErrorString, "DDERR_INVALIDMODE ");
		break;
	case  DDERR_INVALIDOBJECT :
		wsprintf (ErrorString, "DDERR_INVALIDOBJECT ");
		break;
	case  DDERR_INVALIDPARAMS :
		wsprintf (ErrorString, "DDERR_INVALIDPARAMS ");
		break;
	case  DDERR_INVALIDPIXELFORMAT :
		wsprintf (ErrorString, "DDERR_INVALIDPIXELFORMAT ");
		break;
	case  DDERR_INVALIDPOSITION :
		wsprintf (ErrorString, "DDERR_INVALIDPOSITION ");
		break;
	case  DDERR_INVALIDRECT :
		wsprintf (ErrorString, "DDERR_INVALIDRECT ");
		break;
	case  DDERR_INVALIDSURFACETYPE :
		wsprintf (ErrorString, "DDERR_INVALIDSURFACETYPE ");
		break;
	case  DDERR_LOCKEDSURFACES :
		wsprintf (ErrorString, "DDERR_LOCKEDSURFACES ");
		break;
	case  DDERR_NO3D :
		wsprintf (ErrorString, "DDERR_NO3D ");
		break;
	case  DDERR_NOALPHAHW :
		wsprintf (ErrorString, "DDERR_NOALPHAHW ");
		break;
	/*
   case  DDERR_NOANTITEARHW :
		wsprintf (ErrorString, "DDERR_NOANTITEARHW ");
		break;
      */
	case  DDERR_NOBLTHW :
		wsprintf (ErrorString, "DDERR_NOBLTHW ");
		break;
   /*
	case  DDERR_NOBLTQUEUEHW :
		wsprintf (ErrorString, "DDERR_NOBLTQUEUEHW ");
		break;
   */
	case  DDERR_NOCLIPLIST :
		wsprintf (ErrorString, "DDERR_NOCLIPLIST ");
		break;
	case  DDERR_NOCLIPPERATTACHED :
		wsprintf (ErrorString, "DDERR_NOCLIPPERATTACHED ");
		break;
	case  DDERR_NOCOLORCONVHW :
		wsprintf (ErrorString, "DDERR_NOCOLORCONVHW ");
		break;
	case  DDERR_NOCOLORKEY :
		wsprintf (ErrorString, "DDERR_NOCOLORKEY ");
		break;
	case  DDERR_NOCOLORKEYHW :
		wsprintf (ErrorString, "DDERR_NOCOLORKEYHW ");
		break;
	case  DDERR_NOCOOPERATIVELEVELSET :
		wsprintf (ErrorString, "DDERR_NOCOOPERATIVELEVELSET ");
		break;
	case  DDERR_NODC :
		wsprintf (ErrorString, "DDERR_NODC ");
		break;
	case  DDERR_NODDROPSHW :
		wsprintf (ErrorString, "DDERR_NODDROPSHW ");
		break;
	case  DDERR_NODIRECTDRAWHW :
		wsprintf (ErrorString, "DDERR_NODIRECTDRAWHW ");
		break;
	case  DDERR_NOEMULATION :
		wsprintf (ErrorString, "DDERR_NOEMULATION ");
		break;
	case  DDERR_NOEXCLUSIVEMODE :
		wsprintf (ErrorString, "DDERR_NOEXCLUSIVEMODE ");
		break;
	case  DDERR_NOFLIPHW :
		wsprintf (ErrorString, "DDERR_NOFLIPHW ");
		break;
	case  DDERR_NOGDI :
		wsprintf (ErrorString, "DDERR_NOGDI ");
		break;
	case  DDERR_NOHWND :
		wsprintf (ErrorString, "DDERR_NOHWND ");
		break;
	case  DDERR_NOMIPMAPHW :
		wsprintf (ErrorString, "DDERR_NOMIPMAPHW ");
		break;
	case  DDERR_NOMIRRORHW :
		wsprintf (ErrorString, "DDERR_NOMIRRORHW ");
		break;
	case  DDERR_NOOVERLAYDEST :
		wsprintf (ErrorString, "DDERR_NOOVERLAYDEST ");
		break;
	case  DDERR_NOOVERLAYHW :
		wsprintf (ErrorString, "DDERR_NOOVERLAYHW ");
		break;
	case  DDERR_NOPALETTEATTACHED :
		wsprintf (ErrorString, "DDERR_NOPALETTEATTACHED ");
		break;
	case  DDERR_NOPALETTEHW :
		wsprintf (ErrorString, "DDERR_NOPALETTEHW ");
		break;
	case  DDERR_NORASTEROPHW :
		wsprintf (ErrorString, "DDERR_NORASTEROPHW ");
		break;
	case  DDERR_NOROTATIONHW :
		wsprintf (ErrorString, "DDERR_NOROTATIONHW ");
		break;
	case  DDERR_NOSTRETCHHW :
		wsprintf (ErrorString, "DDERR_NOSTRETCHHW ");
		break;
	case  DDERR_NOT4BITCOLOR :
		wsprintf (ErrorString, "DDERR_NOT4BITCOLOR ");
		break;
	case  DDERR_NOT4BITCOLORINDEX :
		wsprintf (ErrorString, "DDERR_NOT4BITCOLORINDEX ");
		break;
	case  DDERR_NOT8BITCOLOR :
		wsprintf (ErrorString, "DDERR_NOT8BITCOLOR ");
		break;
	case  DDERR_NOTAOVERLAYSURFACE :
		wsprintf (ErrorString, "DDERR_NOTAOVERLAYSURFACE ");
		break;
	case  DDERR_NOTEXTUREHW :
		wsprintf (ErrorString, "DDERR_NOTEXTUREHW ");
		break;
	case  DDERR_NOTFLIPPABLE :
		wsprintf (ErrorString, "DDERR_NOTFLIPPABLE ");
		break;
	case  DDERR_NOTFOUND :
		wsprintf (ErrorString, "DDERR_NOTFOUND ");
		break;
	case  DDERR_NOTLOCKED :
		wsprintf (ErrorString, "DDERR_NOTLOCKED ");
		break;
	case  DDERR_NOTPAGELOCKED :
		wsprintf (ErrorString, "DDERR_NOTPAGELOCKED ");
		break;
	case  DDERR_NOTPALETTIZED :
		wsprintf (ErrorString, "DDERR_NOTPALETTIZED ");
		break;
	case  DDERR_NOVSYNCHW :
		wsprintf (ErrorString, "DDERR_NOVSYNCHW ");
		break;
	case  DDERR_NOZBUFFERHW :
		wsprintf (ErrorString, "DDERR_NOZBUFFERHW ");
		break;
	case  DDERR_NOZOVERLAYHW :
		wsprintf (ErrorString, "DDERR_NOZOVERLAYHW ");
		break;
	case  DDERR_OUTOFCAPS :
		wsprintf (ErrorString, "DDERR_OUTOFCAPS ");
		break;
	case  DDERR_OUTOFMEMORY :
		wsprintf (ErrorString, "DDERR_OUTOFMEMORY ");
		break;
	case  DDERR_OUTOFVIDEOMEMORY :
		wsprintf (ErrorString, "DDERR_OUTOFVIDEOMEMORY ");
		break;
	case  DDERR_OVERLAYCANTCLIP :
		wsprintf (ErrorString, "DDERR_OVERLAYCANTCLIP ");
		break;
	case  DDERR_OVERLAYCOLORKEYONLYONEACTIVE :
		wsprintf (ErrorString, "DDERR_OVERLAYCOLORKEYONLYONEACTIVE ");
		break;
	case  DDERR_OVERLAYNOTVISIBLE :
		wsprintf (ErrorString, "DDERR_OVERLAYNOTVISIBLE ");
		break;
	case  DDERR_PALETTEBUSY :
		wsprintf (ErrorString, "DDERR_PALETTEBUSY ");
		break;
	case  DDERR_PRIMARYSURFACEALREADYEXISTS :
		wsprintf (ErrorString, "DDERR_PRIMARYSURFACEALREADYEXISTS ");
		break;
	case  DDERR_REGIONTOOSMALL :
		wsprintf (ErrorString, "DDERR_REGIONTOOSMALL ");
		break;
	case  DDERR_SURFACEALREADYATTACHED :
		wsprintf (ErrorString, "DDERR_SURFACEALREADYATTACHED ");
		break;
	case  DDERR_SURFACEALREADYDEPENDENT :
		wsprintf (ErrorString, "DDERR_SURFACEALREADYDEPENDENT ");
		break;
	case  DDERR_SURFACEBUSY :
		wsprintf (ErrorString, "DDERR_SURFACEBUSY ");
		break;
	case  DDERR_CANTLOCKSURFACE :
		wsprintf (ErrorString, "DDERR_CANTLOCKSURFACE ");
		break;
	case  DDERR_SURFACEISOBSCURED :
		wsprintf (ErrorString, "DDERR_SURFACEISOBSCURED ");
		break;
	case  DDERR_SURFACELOST :
		wsprintf (ErrorString, "DDERR_SURFACELOST ");
		break;
	case  DDERR_SURFACENOTATTACHED :
		wsprintf (ErrorString, "DDERR_SURFACENOTATTACHED ");
		break;
	case  DDERR_TOOBIGHEIGHT :
		wsprintf (ErrorString, "DDERR_TOOBIGHEIGHT ");
		break;
	case  DDERR_TOOBIGSIZE :
		wsprintf (ErrorString, "DDERR_TOOBIGSIZE ");
		break;
	case  DDERR_TOOBIGWIDTH :
		wsprintf (ErrorString, "DDERR_TOOBIGWIDTH ");
		break;
	case  DDERR_UNSUPPORTED :
		wsprintf (ErrorString, "DDERR_UNSUPPORTED ");
		break;
	case  DDERR_UNSUPPORTEDFORMAT :
		wsprintf (ErrorString, "DDERR_UNSUPPORTEDFORMAT ");
		break;
	case  DDERR_UNSUPPORTEDMASK :
		wsprintf (ErrorString, "DDERR_UNSUPPORTEDMASK ");
		break;
	case  DDERR_UNSUPPORTEDMODE :
		wsprintf (ErrorString, "DDERR_UNSUPPORTEDMODE ");
		break;
	case  DDERR_VERTICALBLANKINPROGRESS :
		wsprintf (ErrorString, "DDERR_VERTICALBLANKINPROGRESS ");
		break;
	case  DDERR_WASSTILLDRAWING :
		wsprintf (ErrorString, "DDERR_WASSTILLDRAWING ");
		break;
	case  DDERR_WRONGMODE :
		wsprintf (ErrorString, "DDERR_WRONGMODE ");
		break;
	case  DDERR_XALIGN :
		wsprintf (ErrorString, "DDERR_XALIGN ");
		break;
	default:
		wsprintf (ErrorString, "???");
		break;
	}
	MessageBox(NULL, ErrorString, pMsg, MB_OK|MB_ICONEXCLAMATION|MB_TASKMODAL);
}
