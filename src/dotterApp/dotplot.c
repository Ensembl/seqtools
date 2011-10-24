/*  File: dotplot.c
 *  Author: Gemma Barson, 2010-09-08
 *  Copyright (c) 2010 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: The dotplot widget is the main part of the Dotter window. It
 *              calculates and draws the dot-matrix plot.
 *----------------------------------------------------------------------------
 */

#include <dotterApp/dotter_.h>
#include <dotterApp/seqtoolsExonView.h>
#include <gtk/gtk.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#define PIXELS_PER_MARK_X                           100   /* number of pixels between each major tick mark on the x scale */
#define PIXELS_PER_MARK_Y                           50    /* number of pixels between each major tick mark on the y scale */
#define CROSSHAIR_TEXT_PADDING                      5     /* padding between the crosshair and the coord display text */ 
#define ANNOTATION_LABEL_PADDING		    5	  /* padding around annotation labels, if shown */
#define ANNOTATION_LABEL_LEN			    8	  /* number of chars to allow to show for annotation labels */


int atob_0[]	/* NEW (starting at 0) ASCII-to-binary translation table */
= {
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,23,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,
NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,

NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR 
};

//static int atob[]	/* OLD (starting at 1) ASCII-to-binary translation table  (Inherited from blast) */
//= {
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,24,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA, /* A, B, C, D, ..., P */
//15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA, /* Q, R, S, T, ... */
//NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA, /* a, b, c, d, ..., p */
//15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA, /* q, r, s, t, ... */
//
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
//NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
//};

char aa_btoa[]	/* binary-to-ASCII translation table */
= "-ARNDCQEGHILKMFPSTWYVBZX*" ;


#define NN 5

int ntob[] = {
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

int ntob_compl[] = {
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

double aafq[20]	/* Amino acid residue frequencies used by S. Altschul */
= {.081, .057, .045, .054, .015, .039, .061, .068, .022, .057,
.093, .056, .025, .040, .049, .068, .058, .013, .032, .067 } ;


/* Local function declarations */
static void                       drawDotterScale(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawImage(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawHsps(GtkWidget *dotplot, GdkDrawable *drawable);
static void			  drawBreaklines(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawRubberBand(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       calculateDotplotBorders(GtkWidget *dotplot, DotplotProperties *properties);
static GdkColormap*               insertGreyRamp (DotplotProperties *properties);
static void                       transformGreyRampImage(GdkImage *image, unsigned char *pixmap, DotplotProperties *properties);
static void                       initPixmap(unsigned char **pixmap, const int width, const int height);
static const char*                getShortMspName(const MSP const *msp);
static void                       calculateImageHsps(int strength, int sx, int sy, int ex, int ey, DotplotProperties *properties);
static void                       getMspScreenCoords(const MSP const *msp, DotplotProperties *properties, int *sx, int *ex, int *sy, int *ey);
static void                       setCoordsFromPos(GtkWidget *dotplot, const int x, const int y);
static void                       getCoordsFromPos(GtkWidget *dotplot, const int x, const int y, int *refCoord, int *matchCoord);
static void			  getPosFromSelectedCoords(GtkWidget *dotplot, int *x, int *y);
static void                       getPosFromCoords(DotplotProperties *properties, int qCoord, int sCoord, int *x, int *y);
static void                       setPoint(GdkPoint *point, const int x, const int y, GdkRectangle *rect);
static gdouble                    getScaleFactor(DotplotProperties *properties, const gboolean horizontal);
static void                       initWindow(const char *winsizeIn, DotplotProperties *properties);
static void                       calculateImage(DotplotProperties *properties);
static void                       drawDotplot(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       dotplotDrawCrosshair(GtkWidget *dotplot, GdkDrawable *drawable);

#ifdef ALPHA
static void                       reversebytes(void *ptr, int n);
#endif

/***********************************************************
 *                          Properties                     *
 ***********************************************************/

DotplotProperties* dotplotGetProperties(GtkWidget *widget)
{
  return widget ? (DotplotProperties*)(g_object_get_data(G_OBJECT(widget), "DotplotProperties")) : NULL;
}

static void onDestroyDotplot(GtkWidget *widget)
{
  DotplotProperties *properties = dotplotGetProperties(widget);
  
  if (properties)
    {
      if (properties->pixelmap)
	{
	  g_free(properties->pixelmap);
	  properties->pixelmap = NULL;
	}

    if (properties->hspPixmap)
      {
	g_free(properties->hspPixmap);
	properties->hspPixmap = NULL;
      }
    
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "DotplotProperties", NULL);
    }
}


/* Create the properties object for the dotplot. Note that we still create the properties
 * even if widget is null, because we could be running in batch mode. */
static DotplotProperties* dotplotCreateProperties(GtkWidget *widget,
                                                  DotterWindowContext *dwc,
                                                  const gboolean hspsOn,
                                                  const gboolean breaklinesOn,
                                                  const char *exportFileName)
{
  DotplotProperties *properties = g_malloc(sizeof *properties);
  
  properties->dotterWinCtx = dwc;
  properties->hozExons1 = NULL;
  properties->hozExons2 = NULL;
  properties->vertExons1 = NULL;
  properties->vertExons2 = NULL;
  
  properties->plotRect.x = 0;
  properties->plotRect.y = 0;
  properties->plotRect.width = 0;
  properties->plotRect.height = 0;
  
  if (widget) /* only create colormap if we have a widget (i.e. we're not in batch/non-graphical mode) */
    {
      properties->colorMap = insertGreyRamp(properties);
      gtk_widget_set_default_colormap(properties->colorMap);
    }
  
  properties->image = NULL;
  properties->lineLen = UNSET_INT;

  properties->pixelmap = NULL;
  properties->hspPixmap = NULL;

  properties->crosshairOn = TRUE;
  properties->crosshairCoordsOn = TRUE;
  properties->crosshairFullscreen = TRUE;

  properties->pixelmapOn = !hspsOn;
  properties->hspMode = hspsOn ? DOTTER_HSPS_LINE : DOTTER_HSPS_OFF;
  
  properties->gridlinesOn = FALSE;
  properties->breaklinesOn = breaklinesOn;
  properties->hozLabelsOn = TRUE;
  properties->vertLabelsOn = TRUE;

  properties->dragStart.x = UNSET_INT;
  properties->dragStart.y = UNSET_INT;
  properties->dragEnd.x = UNSET_INT;
  properties->dragEnd.y = UNSET_INT;
      
  properties->exportFileName = exportFileName;

  if (widget)
    {
      g_object_set_data(G_OBJECT(widget), "DotplotProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyDotplot), NULL); 
    }
  
  return properties;
}


DotterHspMode dotplotGetHspMode(GtkWidget *dotplot)
{
  DEBUG_ENTER("dotplotGetHspMode");

  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  DEBUG_EXIT("dotplotGetHspMode returning ");
  return properties->hspMode;
}

int dotplotGetSlidingWinSize(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  return properties->slidingWinSize;
}

gboolean dotplotSetSlidingWinSize(GtkWidget *dotplot, const int newValue, GError **error)
{
  gboolean changed = FALSE;
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  if (newValue <= 0)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_INVALID_WIN_SIZE, "Sliding window size must be greater than 0.\n");
    }
  else if (newValue != properties->slidingWinSize)
    {
      properties->slidingWinSize = newValue;
      changed = TRUE;
    }
  
  return changed;
}

void dotplotSetBreaklinesOn(GtkWidget *dotplot, const gboolean breaklinesOn)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->breaklinesOn = breaklinesOn;
  calculateDotplotBorders(dotplot, properties);
}

void dotplotSetHozLabelsOn(GtkWidget *dotplot, const gboolean labelsOn)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->hozLabelsOn = labelsOn;
  calculateDotplotBorders(dotplot, properties);
}

void dotplotSetVertLabelsOn(GtkWidget *dotplot, const gboolean labelsOn)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->vertLabelsOn = labelsOn;
  calculateDotplotBorders(dotplot, properties);
}

int dotplotGetImageWidth(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  return properties->imageWidth;
}

int dotplotGetImageHeight(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  return properties->imageHeight;
}

double dotplotGetExpResScore(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  return properties->expResScore;
}

int dotplotGetPixelFac(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  return properties->pixelFac;
}

void dotplotToggleBumpExons(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  exonViewToggleBumped(properties->hozExons1);
  exonViewToggleBumped(properties->hozExons2);
  exonViewToggleBumped(properties->vertExons1);
  exonViewToggleBumped(properties->vertExons2);
}

/* Returns true if the pixelmap image should be displayed */
static gboolean showImage(DotplotProperties *properties)
{
  /* Show the greyscale if hsp mode is greyscale or the dot-plot pixelmap is enabled */
  return (properties->pixelmapOn || properties->hspMode == DOTTER_HSPS_GREYSCALE);
}


/* This returns true if the dot-plot should be shown. */
static gboolean showDotplot(DotplotProperties *properties)
{
  /* The HSP greyscale mode overrides the dotplot, so don't show it if that is enabled. */
  return (properties->pixelmapOn && !(properties->hspMode == DOTTER_HSPS_GREYSCALE && properties->hspPixmap));
}


/* Initialise the given pixmap */
static void initPixmap(unsigned char **pixmap, const int width, const int height)
{
  DEBUG_ENTER("initPixmap");

  if (*pixmap)
    g_free(*pixmap);
  
  const int pixelmapLen = width  * height;
  *pixmap = (unsigned char *)g_malloc(sizeof(unsigned char*) * pixelmapLen);
  
  int i = 0;
  for (i=0; i < pixelmapLen; i++)
    {
      (*pixmap)[i] = 0;
    }
  
  DEBUG_EXIT("initPixmap returning ");
}


/* Utility to reset all values in a pixmap to 0 */
static void resetPixmapBackground(unsigned char *pixmap, DotplotProperties *properties)
{
  const int pixelmapLen = properties->image->width * properties->image->height;
  
  int i = 0;
  for (i = 0; i < pixelmapLen; i++) 
    {
      pixmap[i] = 0;
    }
}


/* Set the HSP viewing mode */
void setHspMode(GtkWidget *dotplot, DotterHspMode hspMode)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  properties->hspMode = hspMode;
  
  if (properties->hspMode)
    {
      if (!properties->hspPixmap)
        {
          /* The HSP pixelmap doesn't exist yet so create it */
          initPixmap(&properties->hspPixmap, properties->image->width, properties->image->height);
        }
      
      /* For greyscale mode, loop through all match MSPs and set the pixel strength */
      if (properties->hspMode == DOTTER_HSPS_GREYSCALE) 
        {
          resetPixmapBackground(properties->hspPixmap, properties);
          MSP *msp = dc->mspList;
          
          for ( ; msp; msp = msp->next)
            {
              const char *mspName = getShortMspName(msp);
              if (!mspName || strcmp(mspName, dc->matchSeqName))
                {
                  /* Not an MSP from our match sequence */
                  continue;
                }
              
              int sx, ex, sy, ey;
              getMspScreenCoords(msp, properties, &sx, &ex, &sy, &ey);
              
              /* Draw as greyscale pixmap */
              const int strength = (int)msp->score;
              calculateImageHsps(strength, sx, sy, ex, ey, properties);
            }
          
          /* Overwrite the image with the HSP pixmap */
          transformGreyRampImage(properties->image, properties->hspPixmap, properties);
        }
    }

  if (showDotplot(properties))
    {
      /* Make sure the dot-plot pixmap is set (we may have overwritten it if the previous mode
       * was HSP_GREYSCALE). */
      
      if (!properties->pixelmap)
        {
          /* The dot-plot pixelmap doesn't exist yet so create it */
          initPixmap(&properties->pixelmap, properties->image->width, properties->image->height);
          calculateImage(properties);
        }
      
      transformGreyRampImage(properties->image, properties->pixelmap, properties);
    }
  
  widgetClearCachedDrawable(dotplot, NULL);
  gtk_widget_queue_draw(dotplot);
}


/* Toggle grid lines on/off */
void dotplotToggleGrid(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->gridlinesOn = !properties->gridlinesOn;

  widgetClearCachedDrawable(dotplot, NULL);
  gtk_widget_queue_draw(dotplot);
}


/* Toggle the pixelmap on/off */
void dotplotTogglePixelmap(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->pixelmapOn = !properties->pixelmapOn;
  
  if (!properties->pixelmap)
    {
      /* The dot-plot pixelmap doesn't exist yet so create it */
      initPixmap(&properties->pixelmap, properties->image->width, properties->image->height);
      calculateImage(properties);
    }
  
  if (properties->hspMode != DOTTER_HSPS_GREYSCALE)
    {
      /* Make sure the image has the dot-plot pixels, not the HSP pixels, which might previously
       * have overwritten it */
      transformGreyRampImage(properties->image, properties->pixelmap, properties);
    }
  
  widgetClearCachedDrawable(dotplot, NULL);
  gtk_widget_queue_draw(dotplot);
}


/* Toggle functions for the crosshair properties */
void toggleCrosshairOn(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->crosshairOn = !properties->crosshairOn;
  
  exonViewSetShowCrosshair(properties->hozExons1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->hozExons2, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->vertExons1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->vertExons2, properties->crosshairOn && properties->crosshairFullscreen);
  
  refreshDotplot(dotplot);
}

void toggleCrosshairCoordsOn(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->crosshairCoordsOn = !properties->crosshairCoordsOn;
  gtk_widget_queue_draw(dotplot);
}

void toggleCrosshairFullscreen(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  properties->crosshairFullscreen = !properties->crosshairFullscreen;
  
  exonViewSetShowCrosshair(properties->hozExons1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->hozExons2, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->vertExons1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->vertExons2, properties->crosshairOn && properties->crosshairFullscreen);

  refreshDotplot(dotplot);
}

/***********************************************************
 *			       Events                      *
 ***********************************************************/

/* Expose handler for dot-plot window */
static gboolean onExposeDotplot(GtkWidget *dotplot, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(dotplot)->bin_window;
  DotplotProperties *properties = dotplotGetProperties(dotplot);

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(dotplot);

      if (!bitmap)
        {
          if (properties->exportFileName)
            {
              /* Set the background colour to something sensible for printing */
              GdkColor *defaultBgColor = getGdkColor(DOTCOLOR_BACKGROUND, properties->dotterWinCtx->dotterCtx->defaultColors, FALSE, TRUE);
              setWidgetBackgroundColor(gtk_widget_get_toplevel(dotplot), defaultBgColor);
            }

          bitmap = createBlankPixmap(dotplot);
          drawDotplot(dotplot, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
          
          /* Draw anything else that needs to be refreshed on each expose */
          dotplotDrawCrosshair(dotplot, window);
          drawRubberBand(dotplot, window);
        }
    }

  if (properties->exportFileName)
    {
      /* We're in batch mode and have been asked to export the plot to pdf; print it, 
      * and then quit */
      GError *error = NULL;
      exportPlot(dotplot, GTK_WINDOW(gtk_widget_get_toplevel(dotplot)), properties->exportFileName, &error);
      gtk_main_quit();
    }

  
  return TRUE;
}

/* mouse-press handler for the dot plot window */
static gboolean onButtonPressDotplot(GtkWidget *dotplot, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      setCoordsFromPos(dotplot, event->x, event->y);
      handled = FALSE; /* let the parent update any other widgets that need refreshing */
    }
  else if (event->type == GDK_BUTTON_PRESS && event->button == 2) /* middle click */
    {
      /* Remember start point for drag */
      DotplotProperties *properties = dotplotGetProperties(dotplot);
      setPoint(&properties->dragStart, event->x, event->y, &properties->plotRect);
      handled = TRUE;
    }
  
  return handled;
}


/* mouse-release handler for the dot plot window */
static gboolean onButtonReleaseDotplot(GtkWidget *dotplot, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_RELEASE && event->button == 2) /* middle click */
    {
      /* Call dotter on the selected region */
      DotplotProperties *properties = dotplotGetProperties(dotplot);
      DotterContext *dc = properties->dotterWinCtx->dotterCtx;
      
      if (properties->dragStart.x != UNSET_INT && properties->dragEnd.x != UNSET_INT)
        {
          /* Get the sequence coords for the start and end of the drag rectangle. */
	  const gdouble zoomFactor = 0; 
	
	  int qStart, qEnd, sStart, sEnd;
	  getCoordsFromPos(dotplot, properties->dragStart.x, properties->dragStart.y, &qStart, &sStart);
	  getCoordsFromPos(dotplot, properties->dragEnd.x, properties->dragEnd.y, &qEnd, &sEnd);

	  IntRange qRange, sRange;
	  intrangeSetValues(&qRange, qStart, qEnd);
	  intrangeSetValues(&sRange, sStart, sEnd);
	
          /* Ignore small mouse moves as they are likely to be accidental or cancelled clicks */
          if (qRange.max - qRange.min > 10 && sRange.max - sRange.min > 10)
            {
	      g_debug("Calling dotter internally with the range: q=%d %d, s=%d %d\n", qRange.min, qRange.max, sRange.min, sRange.max);
              callDotterInternal(dc, &qRange, &sRange, zoomFactor, properties->breaklinesOn) ;
            }
        }
      
      /* Clear drag */
      setPoint(&properties->dragStart, UNSET_INT, UNSET_INT, NULL);
      setPoint(&properties->dragEnd, UNSET_INT, UNSET_INT, NULL);
      
      gtk_widget_queue_draw(dotplot);
      handled = TRUE;
    }
  
  return handled;
}


/* Mouse move handler */
static gboolean onMouseMoveDotplot(GtkWidget *dotplot, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->state & GDK_BUTTON1_MASK)  /* left-drag */
    {
      setCoordsFromPos(dotplot, event->x, event->y);
      handled = FALSE; /* let the parent update any other widgets that need refreshing */
    }
  else if (event->state & GDK_BUTTON2_MASK)  /* middle-drag */
    {
      DotplotProperties *properties = dotplotGetProperties(dotplot);
      setPoint(&properties->dragEnd, event->x, event->y, &properties->plotRect);
      gtk_widget_queue_draw(dotplot);
      handled = TRUE;
    }
  
  return handled;
}


/***********************************************************
 *                        Initialisation                   *
 ***********************************************************/

/* Calculate the width/height of the image */
static int getImageDimension(DotplotProperties *properties, const gboolean horizontal)
{
  DEBUG_ENTER("getImageDimension(horizontal=%d)", horizontal);

  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  const IntRange const *seqRange = horizontal ? &dwc->refSeqRange : &dwc->matchSeqRange;
  const int seqLen = getRangeLength(seqRange);
  DEBUG_OUT("Sequence length = %d\n", seqLen);
  
  int imageLen = (int)ceil((double)seqLen / getScaleFactor(properties, horizontal));
  DEBUG_OUT("Image length = %d\n", imageLen);
  
  if (imageLen % 4)
    {
      imageLen += 4 - (imageLen % 4);
    }

  DEBUG_OUT("Image length mod4 = %d\n", imageLen);

  if (seqLen / dc->numFrames > imageLen * dwc->zoomFactor)
    {
      g_critical("seqLen/resfac > imageLen * zoom (%d > %d (%d*%f))", seqLen / dc->numFrames, (int)(imageLen * dwc->zoomFactor), imageLen, dwc->zoomFactor);
    }
  
  DEBUG_EXIT("getImageDimension returning ");
  return imageLen;
}  


/* Set the coords for the initial crosshair position. Use the given coords if in range, or
 * use the center of each seq range otherwise */
static void initCrosshairCoords(const int qcenter, const int scenter, DotterWindowContext *dwc)
{
  if (valueWithinRange(qcenter, &dwc->refSeqRange))
    {
      dwc->refCoord = qcenter;
    }
  else
    {
      dwc->refCoord = getRangeCentre(&dwc->refSeqRange);
    }
  
  if (valueWithinRange(qcenter, &dwc->matchSeqRange))
    {
      dwc->matchCoord = qcenter;
    }
  else
    {
      dwc->matchCoord = getRangeCentre(&dwc->matchSeqRange);
    }
}


/* Create the actual drawing area for the dotplot */
static GtkWidget* createDotplotDrawingArea(DotterWindowContext *dwc, 
                                           const char *loadFileName,
                                           const char *saveFileName,
                                           const char *exportFileName,
                                           const gboolean hspsOn,
                                           const gboolean breaklinesOn,
                                           const char *initWinsize,
                                           const int pixelFacIn,
                                           const int zoomFacIn,
                                           const int qcenter,
                                           const int scenter)
{
  DEBUG_ENTER("createDotplotDrawingArea");

  /* If we're in batch mode and we're exporting the matrix only (i.e. we have a save file
   * but we don't have an export file), then we don't need to create the plot */
  gboolean batch = (saveFileName || exportFileName);
  gboolean showPlot = (exportFileName || !batch);

  GtkWidget *dotplot = (showPlot ? gtk_layout_new(NULL, NULL) : NULL);
  
  DotplotProperties *properties = dotplotCreateProperties(dotplot, dwc, hspsOn, breaklinesOn, exportFileName);
  
  if (loadFileName)
    {
      /* Load existing dotplot file (initialises and populates the pixmap) */
      GError *error = NULL;
      loadPlot(dotplot, loadFileName, &error);
      
      prefixError(error, "Error loading dot plot. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  else 
    {
      initWindow(initWinsize, properties);
      
      /* Set pixelFac so that expResScore is at 1/5 of the range. 
       * This positions expResScore at 51.2 */
      properties->pixelFac = pixelFacIn ? pixelFacIn : 0.2 * NUM_COLORS / properties->expResScore;
      
      /* Calculate the image size */
      properties->imageWidth = getImageDimension(properties, TRUE);
      properties->imageHeight = getImageDimension(properties, FALSE);
      properties->lineLen = properties->imageWidth;
      DEBUG_OUT("Set image w=%d, h=%d, line len=%d\n", properties->imageWidth, properties->imageHeight, properties->lineLen);
      
      unsigned char **pixmap = NULL; /* which pixelmap we're displaying at the start */
      
      if (properties->hspMode == DOTTER_HSPS_GREYSCALE)
        {
          pixmap = &properties->hspPixmap;
          initPixmap(pixmap, properties->imageWidth, properties->imageHeight);
        }
      else if (properties->pixelmapOn)
        {
          pixmap = &properties->pixelmap;
          initPixmap(pixmap, properties->imageWidth, properties->imageHeight);
          calculateImage(properties);
        }

      if (showPlot)
        {
          /* Create the image and push the pixelmap to it */
          properties->image = gdk_image_new(GDK_IMAGE_NORMAL, gdk_visual_get_system(), properties->imageWidth, properties->imageHeight);
          
          if (pixmap)
            transformGreyRampImage(properties->image, *pixmap, properties);
        }
    }
  
  if (saveFileName)
    {
      /* Batch mode: save the dot matrix;  */
      GError *error = NULL;
      savePlot(dotplot, properties, saveFileName, &error);
      
      prefixError(error, "Error saving dot plot. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  if (showPlot)
    {
      /* Set the crosshair coords (only applicable if not in batch mode) */
      if (!batch)
        initCrosshairCoords(qcenter, scenter, properties->dotterWinCtx);

      /* Calculate the size of the dot-plot */
      calculateDotplotBorders(dotplot, properties);
    }
  
  DEBUG_EXIT("createDotplotDrawingArea returning ");
  return dotplot;
}


/* Put all the dotplot and exon widgets in a table. The table gets put in a
 * scrolled window and the scrolled window is returned. */
static GtkWidget* createDotplotTable(GtkWidget *dotplotCont, 
                                     GtkWidget *hozExons1,
                                     GtkWidget *hozExons2,
                                     GtkWidget *vertExons1,
                                     GtkWidget *vertExons2)
{
  const int xpad = 0;
  const int ypad = 0;
  const int numRows = 4;
  const int numCols = 4;
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  
  gtk_table_attach(table, dotplotCont, 1, 2, 1, 2, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, vertExons1, 2, 3, 1, 2, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, vertExons2, 3, 4, 1, 2, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, hozExons1, 1, 2, 2, 3, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, hozExons2, 1, 2, 3, 4, GTK_FILL, GTK_FILL, xpad, ypad);
  
  /* Put the table in a scrolled window with a horizontal scrollbar */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  
  return GTK_WIDGET(scrollWin);
}


/* Return the reverse strand if the given strand is forward and vice versa */
static BlxStrand getOppositeStrand(BlxStrand strand)
{
  return (strand == BLXSTRAND_FORWARD ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
}


/* Create the exons view widgets for dotter */
static void createDotterExonViews(GtkWidget *dotplot, 
                                  DotterWindowContext *dwc, 
                                  GtkWidget **hozExons1,
                                  GtkWidget **hozExons2,
                                  GtkWidget **vertExons1,
                                  GtkWidget **vertExons2)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterContext *dc = dwc->dotterCtx;
  const gboolean showCrosshair = properties->crosshairOn && properties->crosshairFullscreen;

  /* The strand that was passed to dotter is the main strand, so show that first */
  *hozExons1 = createDotterExonView(dotplot, 
                                    NULL,
                                    TRUE,
                                    dc->refSeqStrand, 
                                    dwc, 
                                    properties->imageWidth,
                                    properties->imageHeight,
                                    &dwc->refSeqRange, 
                                    showCrosshair,
                                    &properties->hozExons1);
  
  *hozExons2 = createDotterExonView(dotplot, 
                                    NULL,
                                    TRUE,
                                    getOppositeStrand(dc->refSeqStrand), 
                                    dwc, 
                                    properties->imageWidth, 
                                    properties->imageHeight,
                                    &dwc->refSeqRange, 
                                    showCrosshair, 
                                    &properties->hozExons2);
  
  *vertExons1 = createDotterExonView(dotplot, 
                                     NULL, 
                                     FALSE,
                                     dc->matchSeqStrand, 
                                     dwc, 
                                     properties->imageWidth,
                                     properties->imageHeight,
                                     &dwc->matchSeqRange,
                                     showCrosshair,
                                     &properties->vertExons1);

  *vertExons2 = createDotterExonView(dotplot, 
                                     NULL, 
                                     FALSE,
                                     getOppositeStrand(dc->matchSeqStrand),
                                     dwc, 
                                     properties->imageWidth, 
                                     properties->imageHeight,
                                     &dwc->matchSeqRange, 
                                     showCrosshair,
                                     &properties->vertExons2);
}


/* Create the drawing area for the dot plot. The 'dotplot' output arg is populated with the
 * actual drawing area; the return arg is the overall container widget for the dot plot that will
 * be packed into the dotter window UNLESS we're running in batch mode, in which case the return 
 * will be no widgets displayed to the user. */
GtkWidget* createDotplot(DotterWindowContext *dwc, 
                         const char *loadFileName,
                         const char *saveFileName,
                         const char *exportFileName,
                         const gboolean hspsOn,
                         const gboolean breaklinesOn,
                         const char *initWinsize,
                         const int pixelFacIn,
                         const int zoomFacIn,
                         const int qcenter,
                         const int scenter,
                         GtkWidget **dotplot)
{
  DEBUG_ENTER("createDotplot");

  /* Create the actual drawing area for the dot plot */
  *dotplot = createDotplotDrawingArea(dwc, loadFileName, saveFileName, exportFileName, hspsOn, breaklinesOn, initWinsize, pixelFacIn, zoomFacIn, qcenter, scenter);

  /* Create and realise the container widgets etc. only if not in batch mode, or if
   * in batch mode and exporting to PDF (which requires the widget to draw itself). */
  const gboolean batch = (exportFileName || saveFileName);
  const gboolean createWindow = (!batch || exportFileName);

  if (!createWindow)
    {
      return NULL;
    }

  /* Create the 4 exon views: one for each strand, for both sequences */
  GtkWidget *hozExons1 = NULL;
  GtkWidget *hozExons2 = NULL;
  GtkWidget *vertExons1 = NULL;
  GtkWidget *vertExons2 = NULL;
  createDotterExonViews(*dotplot, dwc, &hozExons1, &hozExons2, &vertExons1, &vertExons2);

  /* Put everything in a table */
  GtkWidget *parent = createDotplotTable(*dotplot, hozExons1, hozExons2, vertExons1, vertExons2);
  
  gtk_widget_add_events(*dotplot, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(*dotplot, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(*dotplot, GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(*dotplot), "expose-event", G_CALLBACK(onExposeDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "button-press-event", G_CALLBACK(onButtonPressDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "button-release-event", G_CALLBACK(onButtonReleaseDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "motion-notify-event",   G_CALLBACK(onMouseMoveDotplot), NULL);

  DEBUG_EXIT("createDotplot returning ");
  return parent;
}


/* Delete image and pixmaps. Needed if we have to re-create them at a different size. */
static void clearPixmaps(DotplotProperties *properties)
{
  if (properties->image)
    {
      gdk_image_unref(properties->image);
      properties->image = NULL;
    }
  
  if (properties->pixelmap)
    {
      g_free(properties->pixelmap);
      properties->pixelmap = NULL;
    }
  
  if (properties->hspPixmap)
    {
      g_free(properties->hspPixmap);
      properties->hspPixmap = NULL;
    }
}


/* Return the alphabet size for the given sequence type (4 for DNA, 20 for protein) */
static int getAlphabetSize(const BlxSeqType seqType)
{
  return (seqType == BLXSEQ_DNA ? 4 : 20);
}



static void initWindow(const char *winsizeIn, DotplotProperties *properties)
{
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  double exp1, exp2, exp3, lambda;
  int win1, win2, win3;
  
  const int alphabetSize = getAlphabetSize(dc->displaySeqType);
  
  /* Call winsizeFromlambdak even if we don't want to set the window
   size in order to get the other parameters (properties->expResScore) */
  
  if (dc->blastMode == BLXMODE_BLASTX) 
    {
      win1 = winsizeFromlambdak(dc->matrix, atob_0, alphabetSize, dc->peptideSeqs[0], dc->matchSeq, &exp1, &lambda);
      win2 = winsizeFromlambdak(dc->matrix, atob_0, alphabetSize, dc->peptideSeqs[1], dc->matchSeq, &exp2, &lambda);
      win3 = winsizeFromlambdak(dc->matrix, atob_0, alphabetSize, dc->peptideSeqs[2], dc->matchSeq, &exp3, &lambda);
      properties->expResScore = (exp1 + exp2 + exp3)/3.0;
      properties->slidingWinSize = (win1 + win2 + win3)/3.0;
    }
  else if (dc->blastMode == BLXMODE_BLASTN)
    {
      properties->slidingWinSize = winsizeFromlambdak(dc->matrix, ntob, alphabetSize, dc->refSeq, dc->matchSeq, &properties->expResScore, &lambda);
    }
  else
    {
      properties->slidingWinSize = winsizeFromlambdak(dc->matrix, atob_0, alphabetSize, dc->refSeq, dc->matchSeq, &properties->expResScore, &lambda);
    }
  
  if (!winsizeIn || toupper(*winsizeIn) == 'K') 
    {
      if (properties->slidingWinSize < 3) 
        {
          g_critical("Karlin/Altschul estimate of window size = %d ignored. Using 10 instead.\n", properties->slidingWinSize);
          properties->slidingWinSize = 10;
        }
      
      if (properties->slidingWinSize > 50) 
        {
          g_critical("Karlin/Altschul estimate of window size = %d ignored. Using 50 instead.\n", properties->slidingWinSize);
          properties->slidingWinSize = 50;
        }
    }
  else
    {
      if (!atoi(winsizeIn))
        {
          g_error("Bad window size specification: %s\n", winsizeIn);
        }
      
      properties->slidingWinSize = atoi(winsizeIn);
    }
}


/* Get a base in the horizontal (reference) sequence. The given index is zero-based within our
 * current display range. */
static char getHozSeqBase(DotterWindowContext *dwc, const int idx, const int frame, const int offset)
{
  DotterContext *dc = dwc->dotterCtx;

  char result;

  if (dc->blastMode == BLXMODE_BLASTX)
    {
      result = dc->peptideSeqs[frame][idx + offset];
    }
  else
    {
      /* Reverse the sequence if the scale is reversed */
      const int coord = dc->hozScaleRev ? dwc->refSeqRange.max - idx : dwc->refSeqRange.min + idx;
      
      /* Complement the sequence if it's the reverse strand */
      const gboolean complement = (dc->refSeqStrand == BLXSTRAND_REVERSE && dc->refSeqType == BLXSEQ_DNA && dc->hozScaleRev);

      result = getSequenceIndex(dc->refSeq, coord, complement, &dc->refSeqFullRange, dc->refSeqType);
    }
  
  return result;
}


/* Get a base in the vertical (match) sequence. The given index is zero-based within our
 * current display range. */
static char getVertSeqBase(DotterWindowContext *dwc, const int idx)
{
  DotterContext *dc = dwc->dotterCtx;

  /* Reverse the sequence if the scale is reversed. */
  const int coord = dc->vertScaleRev ? dwc->matchSeqRange.max - idx : dwc->matchSeqRange.min + idx;

  const gboolean complement = (dc->matchSeqStrand == BLXSTRAND_REVERSE && dc->refSeqType == BLXSEQ_DNA && dc->vertScaleRev);
  
  return getSequenceIndex(dc->matchSeq, coord, complement, &dc->matchSeqFullRange, dc->matchSeqType);
}


/* Populate the array of binary values for the match sequence (used by calculateImage) */
static void populateMatchSeqBinaryVals(DotterWindowContext *dwc, const int slen, const int translationTable[], int *sIndex)
{
  int sIdx = 0;
  
  for ( ; sIdx < slen; ++sIdx) 
    {
      const char sBase = getVertSeqBase(dwc, sIdx);
      const int asciiVal = (int)sBase;
      const int aminoAcidId = translationTable[asciiVal];
      sIndex[sIdx] = aminoAcidId;
    }
}


/* Create the score vector that is used by calculateImage. Allocates memory but does not populate it
 * except for the not-a-residue entry) */
static void createScoreVec(DotterWindowContext *dwc, const int vecLen, const int qlen, BlxHandle *handle, gint32 ***scoreVecPtr)
{
  DotterContext *dc = dwc->dotterCtx;
  
  *scoreVecPtr = handleAlloc(handle, vecLen * sizeof(int*));
  gint32 **scoreVec = *scoreVecPtr;

  /* Allocate memory for each row in the score vector. Each row is the length of the ref sequence */
  int rowIdx = 0;
  for ( ; rowIdx < vecLen; ++rowIdx)
    {
      scoreVec[rowIdx] = handleAlloc(handle, qlen * sizeof(gint32));
    }

  if (dc->blastMode != BLXMODE_BLASTN)
    {
      /* Populate non-protein symbols in scorevector */
      int qIdx = 0;
      for (qIdx = 0; qIdx < qlen; ++qIdx) 
        {
          scoreVec[vecLen - 1][qIdx] = dc->matrix[vecLen - 2][vecLen - 2];
        }
    }
}


static void populateScoreVec(DotterWindowContext *dwc, const int vecLen, const int qlen, const int frame, const int offset, const int translationTable[], int **scoreVec)
{
  DotterContext *dc = dwc->dotterCtx;
  
  /* Loop through each row */
  int rowIdx = 0;
  for (rowIdx = 0; rowIdx < vecLen - 1; ++rowIdx)
    {
      /* Loop through each index on the row */
      int qIdx = 0;
      for (qIdx = 0; qIdx < qlen; ++qIdx)
        {
          const char qBase = getHozSeqBase(dwc, qIdx, frame, offset);
          const int asciiVal = (int)qBase;
          const int aminoAcidId = translationTable[asciiVal];
          const gint32 score = dc->matrix[rowIdx][aminoAcidId]; /* score of this base in the q seq wrt the current amino acid ID 'i' */
          
          scoreVec[rowIdx][qIdx] = score;
        }
    }
}


/* Return the translation table required for a sequence of the given type and strand */
static int* getTranslationTable(const BlxSeqType seqType, const BlxStrand strand)
{
  int *result = NULL;
  
  if (seqType == BLXSEQ_PEPTIDE)
    {
      result = atob_0;
    }
  else
    {
      result = (strand == BLXSTRAND_REVERSE ? ntob_compl : ntob);
    }
     
  return result;
}


static int* getRowToDelete(const int sIdx, const int slen, const int *sIndex, int **scoreVec, int *zero, const int slidingWinSize, const BlxStrand strand)
{
  int *delrow = NULL;
  
  if (strand == BLXSTRAND_REVERSE)
    {
      if (sIdx < slen - slidingWinSize) 
        {
          delrow = scoreVec[sIndex[sIdx + slidingWinSize]];
        }
      else
        {
          delrow = zero;
        }
    }
  else
    {
      if (sIdx >= slidingWinSize) 
        {
          delrow = scoreVec[sIndex[sIdx - slidingWinSize]];
        }
      else
        {
          delrow = zero;
        }
    } 
  
  return delrow;
}


/* This does the work for calculateImage, for a particular strand and reading frame of the reference sequence */
static void doCalculateImage(const BlxStrand qStrand, 
                             const int incrementVal, 
                             const int sStart, 
                             const int frame,
                             DotterWindowContext *dwc,
                             DotplotProperties *properties,
                             const int pepQSeqLen,
                             const int pepQSeqOffset,
                             const int slen,
                             const int qlen,
                             const int vecLen,
                             const int win2,
                             int **scoreVec,
                             int *sIndex,
                             int *sum1,
                             int *sum2,
                             int *zero)
{
  DotterContext *dc = dwc->dotterCtx;
  const int pixelmapLen = properties->imageWidth * properties->imageHeight;

  register int qIdx, sIdx, qmax, dotpos, dotposq, dotposs;
  
  register int *newsum;	/* The current row of scores being calculated */
  register int *oldsum;	/* Remembers the previous row of calculated scores */
  register int *delrow;	/* Pointer to the row in scoreVec to subtract */
  register int *addrow;	/* Pointer to the row in scoreVec to add */
  
  /* Reset the sum vectors */
  int idx = 0;
  for ( ; idx < pepQSeqLen; ++idx) 
    {
      sum1[idx] = 0;
      sum2[idx] = 0;
    }
  
  /* Get the range of valid calculations (excluding the initial sliding window size, where we don't have enough 
   * info to calculate the average properly - exclude the winsize at the start if fwd or the end if reverse) */
  IntRange validRange;
  validRange.min = (qStrand == BLXSTRAND_REVERSE ? 0 : properties->slidingWinSize);
  validRange.max = (qStrand == BLXSTRAND_REVERSE ? slen - properties->slidingWinSize : slen);
  
  /* Re-populate the score vector for this reading frame */
  populateScoreVec(dwc, vecLen, pepQSeqLen, frame, pepQSeqOffset, getTranslationTable(dc->displaySeqType, qStrand), scoreVec);
  
  /* Loop through each base in the match sequence */
  for (sIdx = sStart ; sIdx >= 0 && sIdx < slen; sIdx += incrementVal)
    {   
      /* Set oldsum to the previous row. (newsum will be overwritten, but we re-use the
       * same two vectors (sum1 and sum2) here to save having to keep allocating memory) */
      oldsum = (sIdx & 1) ? sum2 : sum1;
      newsum = (sIdx & 1) ? sum1 : sum2;
      
      delrow = getRowToDelete(sIdx, slen, sIndex, scoreVec, zero, properties->slidingWinSize, qStrand);
      
      /* We add the pre-calculated value from the score vector for the current amino acid */
      addrow = scoreVec[sIndex[sIdx]];
      *newsum = *addrow;
      ++addrow;
      
      qmax = min(properties->slidingWinSize, pepQSeqLen);
      
      for (qIdx = 1; qIdx < qmax ; ++qIdx)
        {
          ++newsum;
          *newsum = *oldsum + *addrow;
          ++oldsum;
          ++addrow;
        }
      
      qmax = (dc->blastMode != BLXMODE_BLASTX && dwc->selfComp ? sIdx + 1 : pepQSeqLen);
      qmax = min(qmax, pepQSeqLen);
      
      for ( ; qIdx < qmax ; ++qIdx) 
        {
          ++newsum;
          *newsum = *oldsum + *addrow - *delrow;
          ++oldsum;
          ++addrow;
          ++delrow;
          
          if (*newsum > 0 && valueWithinRange(sIdx, &validRange)) 
            {
              dotposq = (qIdx - win2)/dwc->zoomFactor;
              dotposs = (sIdx - (incrementVal * win2))/dwc->zoomFactor;
              
              /* Only fill half the submatrix */
              const int qPosLocal = qIdx - win2 - (dotposq * dwc->zoomFactor);  /* query position in local submatrix (of one pixel) */
              int sPosLocal = sIdx - (incrementVal * win2) - (dotposs * dwc->zoomFactor);  /* subject position in local submatrix (of one pixel) */
              
              if (qStrand == BLXSTRAND_REVERSE)
                {
                  /* Set the origin (0,0) to the bottom left corner of submatrix
                   Ugly but correct. Zoom = pixels/submatrix */
                  sPosLocal = dwc->zoomFactor - 1 - sPosLocal;
                }
              
              if (sPosLocal >= qPosLocal)
                {
                  dotpos = properties->imageWidth*dotposs + dotposq;
                  
                  if (dotpos < 0 || dotpos >= pixelmapLen) 
                    {
                      g_critical ( "Pixel %d out of bounds. Pixelmap len=%d, mode =%d, ref sequqnece strand=%s\n", dotpos, pixelmapLen, dc->blastMode, (qStrand == BLXSTRAND_REVERSE ? "reverse" : "forward"));
                    }
                  else
                    {
                      /* Keep the max dot value of all diagonals in this pixel */
                      const int val = *newsum * properties->pixelFac / properties->slidingWinSize;
                      unsigned char dotValue = (val > 255 ? 255 : (unsigned char)val);
                      unsigned char *curDot = &properties->pixelmap[dotpos];
                      
                      if (dotValue > *curDot)
                        {
                          *curDot = dotValue;
                        }
                    }
                }
            }
        } 
    }
}


/* Print some debug info for the calculateImage function to stdout */
static void printCalculateImageStats(DotterWindowContext *dwc, const int qlen, const int slen)
{
  DotterContext *dc = dwc->dotterCtx;
  
  double speed = 17.2;  /* Speed in Mdots/seconds. SGI MIPS R10000 (clobber) */
  /* speed = 5.7;  DEC Alpha AXP 3000/700 */
  /* speed = 3.7;  SGI R4400: */
  
  double numDots = qlen/1e6*slen; /* total number of dots (millions) */
  
  if (dwc->selfComp) 
    numDots /= 2;
  
  if (dc->blastMode == BLXMODE_BLASTN && !(dc->watsonOnly || dc->crickOnly)) 
    numDots *= 2;
  
  if (dc->blastMode == BLXMODE_BLASTX) 
    numDots *= 3;
  
  int min = (int)(numDots/speed/60);
  int sec = (int)(numDots/speed) - min*60;
  
  g_message("%d vs. %d residues => %.2f million dots. ", qlen, slen, numDots);
  
  if (min+sec >= 2) 
    {
      g_message("(Takes ");
      
      if (min)
        g_message("%d:%.2d minutes", min, sec);
      else 
        g_message("%d seconds", sec);
      
      g_message(" on an SGI MIPS R10000)");
    }
  
  g_message("\n");
  fflush(stdout);
}


/* This function calculates the data that will be put into the dotplot image. It puts the max
 * diagonal for each pixel into *pixelmap */
static void calculateImage(DotplotProperties *properties)
{
  DEBUG_ENTER("calculateImage");

  g_assert(properties->slidingWinSize > 0);
  
  register int qIdx, sIdx;     /* Loop variables */
  register int dotpos;
  
  BlxHandle handle = handleCreate();
  
  /* Extract some often-used data */
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  const int qlen = getRangeLength(&dwc->refSeqRange);
  const int slen = getRangeLength(&dwc->matchSeqRange);
  const int win2 = properties->slidingWinSize/2;

  /* Print some statistics about what we're about to do */
  printCalculateImageStats(dwc, qlen, slen);
  
  /* Find the offset of the current display range within the full range of the bit of reference sequence we have */
  const int qOffset = dc->refSeqStrand == BLXSTRAND_REVERSE 
    ? dc->refSeqFullRange.max - dwc->refSeqRange.max
    : dwc->refSeqRange.min - dc->refSeqFullRange.min;
  
  /* Convert from nucleotides to peptides, if applicable */
  const int resFactor = (dc->blastMode == BLXMODE_BLASTX ? dc->numFrames : 1);
  const int pepQSeqLen = qlen / resFactor;
  const int pepQSeqOffset = qOffset / resFactor;
  const int vecLen = (dc->displaySeqType == BLXSEQ_DNA ? 6 : 25);

  /* Initialize lookup tables for faster execution. scoreVec is an array of precalculated 
   * scores for qseq residues. sIndex contains the match sequence forward strand bases as 
   * binary values (i.e. amino-acid IDs 0 -> 23) */
  gint32 **scoreVec = NULL;
  createScoreVec(dwc, vecLen, pepQSeqLen, &handle, &scoreVec);

  gint32 *sIndex = handleAlloc(&handle, slen * sizeof(gint32));
  populateMatchSeqBinaryVals(dwc, slen, getTranslationTable(dc->matchSeqType, BLXSTRAND_FORWARD), sIndex);
  
  /* Allocate some vectors for use in averaging the values for whole rows at a time. Initialise the
   * 'zero' array now but leave the sum arrays because these will be reset in doCalculateWindow. */
  gint32 *zero = handleAlloc(&handle, pepQSeqLen * sizeof(gint32));
  gint32 *sum1 = handleAlloc(&handle, pepQSeqLen * sizeof(gint32));
  gint32 *sum2 = handleAlloc(&handle, pepQSeqLen * sizeof(gint32));

  int idx = 0;
  for (idx = 0; idx < pepQSeqLen; ++idx) 
    {
      zero[idx] = 0;
    }

  if (dc->blastMode == BLXMODE_BLASTX)
    {
      /* Protein -> Nucleotide matches. Calculate the result for each reqding frame of the
       * reference sequence, and use the overall max values. */
      int frame = 0;
      for ( ; frame < dc->numFrames; ++frame)
        {
          doCalculateImage(BLXSTRAND_FORWARD, 1, 0, frame,
                           dwc, properties, pepQSeqLen, pepQSeqOffset, slen, qlen, vecLen, win2,
                           scoreVec, sIndex, sum1, sum2, zero);

        }
    }
  else if (dc->blastMode == BLXMODE_BLASTP) 
    {
      /* Protein -> Protein matches. Straightforward comparison of the two sequences. */
      doCalculateImage(BLXSTRAND_FORWARD, 1, 0, 0,
                       dwc, properties, pepQSeqLen, pepQSeqOffset, slen, qlen, vecLen, win2,
                       scoreVec, sIndex, sum1, sum2, zero);
    }
  else if (dc->blastMode == BLXMODE_BLASTN)
    {
      /* Nucleotide -> Nucleotide matches. Calculate the result for each strand of the reference
       * sequence, and use the overall max values. */
      
      if (!dc->crickOnly)
        {
          doCalculateImage(BLXSTRAND_FORWARD, 1, 0, 0,
                          dwc, properties, pepQSeqLen, pepQSeqOffset, slen, qlen, vecLen, win2,
                          scoreVec, sIndex, sum1, sum2, zero);
        }

      if (!dc->watsonOnly)
        {
          doCalculateImage(BLXSTRAND_REVERSE, -1, slen - 1, 0,
                           dwc, properties, pepQSeqLen, pepQSeqOffset, slen, qlen, vecLen, win2,
                           scoreVec, sIndex, sum1, sum2, zero);
        }
    }

  if (dwc->selfComp && dc->displayMirror) 
    {
      /* Copy mirror image */
      
      int dotposCopy;
      
      for (sIdx = 0; sIdx < properties->imageHeight; ++sIdx) 
        { 
          for (qIdx = 0; qIdx < sIdx ; ++qIdx) 
            {
              dotpos = properties->imageWidth * sIdx + qIdx;
              dotposCopy = properties->imageWidth * qIdx + sIdx;

              const int pixelmapLen = properties->imageWidth * properties->imageHeight;

              if (dotpos < 0 || dotpos >= pixelmapLen)
                g_critical ( "Source pixel out of bounds (%d) in mirrorCopy: %d\n",
                            pixelmapLen-1, dotpos);
              if (dotposCopy < 0 || dotposCopy >= pixelmapLen)
                g_critical ( "Destination pixel out of bounds (%d) in mirrorCopy: %d\n",
                            pixelmapLen-1, dotposCopy);
              properties->pixelmap[dotposCopy] = properties->pixelmap[dotpos];
            }
        }
    }
  
  handleDestroy(&handle);
  
  DEBUG_EXIT("calculateImage returning ");
}


void loadPlot(GtkWidget *dotplot, const char *loadFileName, GError **error)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* Open the file */
  FILE *loadFile = fopen (loadFileName, "rb");
  
  if (!loadFile)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_OPENING_FILE, "Could not open file '%s'", loadFileName);
      return;
    }

  /* 'dotstart' is the sum of how many bytes there are before the dot-plot data starts. We remember this
   * value so that we can look ahead to check the length of the file with fseek but still know where
   * to jump back to in the file to start reading the dot-plot data. */
  int dotstart = 0;

  /* Read in the format, zoom, and image size */
  unsigned char format;
  gboolean ok = fread(&format, 1, sizeof(unsigned char), loadFile) == sizeof(unsigned char);
  dotstart += sizeof(unsigned char);

  if (format != 1 && format != 2 && format != 3) 
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Unknown dotter file format version: %d", format);
      return;
    }
  
  /* Read in the zoom factor (an int in formats 1 and 2, or gdouble in later formats) */
  if (format == 1 || format == 2)
    {
      ok &= fread(&dwc->zoomFactor, 1, sizeof(gint32), loadFile) == sizeof(gint32);
      dotstart += sizeof(gint32);
    }
  else
    {
      ok &= fread(&dwc->zoomFactor, 1, sizeof(gdouble), loadFile) == sizeof(gdouble);
      dotstart += sizeof(gdouble);
    }
  
  /* Read in the image size */
  ok &= fread(&properties->imageWidth, 1, sizeof(gint32), loadFile) == sizeof(gint32);
  dotstart += sizeof(gint32);

  ok &= fread(&properties->imageHeight, 1, sizeof(gint32), loadFile) == sizeof(gint32);
  dotstart += sizeof(gint32);
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Error reading file '%s'", loadFileName);
      return;
    }
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(gint32));
  reversebytes(&properties->imageHeight, sizeof(gint32));
#endif
  
  properties->lineLen = properties->imageWidth;
  
  if (format == 1) 
    {
      /* Don't actually know these variables for sure - guess the most common */
      properties->pixelFac = 50;
      properties->slidingWinSize = 25;
    }
  else 
    {
      ok &= fread(&properties->pixelFac, 1, sizeof(gint32), loadFile) == sizeof(gint32);
      dotstart += sizeof(gint32);

      ok &= fread(&properties->slidingWinSize, 1, sizeof(gint32), loadFile) == sizeof(gint32);
      dotstart += sizeof(gint32);

      /* Read in the matrix name (read the length of the name first) */
      gint32 matrixNameLen;
      ok &= fread(&matrixNameLen, 1, sizeof(gint32), loadFile) == sizeof(gint32);
      dotstart += sizeof(gint32);

#ifdef ALPHA
      if (ok)
        {
          reversebytes(&properties->pixelFac, sizeof(gint32));
          reversebytes(&properties->slidingWinSize, sizeof(gint32));
          reversebytes(&matrixNameLen, sizeof(gint32));
        }
#endif

      ok &= matrixNameLen <= MAX_MATRIX_NAME_LENGTH;
      char matrixName[MAX_MATRIX_NAME_LENGTH + 1];

      ok &= fread(&matrixName, sizeof(char), matrixNameLen, loadFile) == matrixNameLen * sizeof(char);
      dotstart += sizeof(char) * matrixNameLen;

      if (ok)
        {
          /* Terminate the matrix name and copy it into the context. */
          matrixName[matrixNameLen] = 0;
          g_free(dc->matrixName);
          dc->matrixName = g_strdup(matrixName);
      
          /* Read in the matrix data */
          int i = 0;
          int j = 0;
          
          for (i = 0; i < 24; i++)
            {
              for (j = 0; j < 24; j++) 
                {
                  gint32 matrixVal;
                  ok &= fread(&matrixVal, 1, sizeof(gint32), loadFile) == sizeof(gint32); 
                  dotstart += sizeof(gint32);
#ifdef ALPHA
                  reversebytes(&matrixVal, sizeof(gint32));
#endif
                  dc->matrix[i][j] = matrixVal;
                }
            }
        }
    }
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Error reading file '%s'", loadFileName);
      return;
    }
  
  fseek(loadFile, 0, SEEK_END);
  int n = ftell(loadFile);

  const int pixelmapLen = properties->imageHeight*properties->imageWidth;

  if (n - dotstart != pixelmapLen)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Wrong number of pixels in %s: %d. Expected %d * %d = %d\n", 
            loadFileName, n - dotstart, properties->imageWidth, properties->imageHeight, pixelmapLen);
      return;
    }
  
  /* Allocate memory for the pixmap */
  initPixmap(&properties->pixelmap, properties->imageWidth, properties->imageHeight);
  
  fseek(loadFile, dotstart, SEEK_SET);
  
  int i = 0;
  for (i = 0; i < pixelmapLen; i++)
    {
      unsigned char pixelVal;
      ok &= fread(&pixelVal, 1, sizeof(unsigned char), loadFile) == sizeof(unsigned char); 
#ifdef ALPHA
      reversebytes(&pixelVal, sizeof(unsigned char));
#endif
      properties->pixelmap[i] = pixelVal;
    }
  
  
//  if ((n = fread(properties->pixelmap, sizeof(unsigned char), pixelmapLen, loadFile)) != pixelmapLen)
//    {
//      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Read wrong number of pixels from %s: %d. Expected %d * %-d = %d\n", 
//              loadFileName, n, properties->imageWidth, properties->imageHeight, pixelmapLen);
//      return;
//    }
  
  fclose(loadFile);

  g_message("Dotplot file '%s' was loaded successfully.\n", loadFileName);

  if (format == 1) 
    {
      g_message("Old dotplot file format '1' was used, so the windowsize and pixel factor are estimates.\n");
    }

  fflush(stdout);

  /* Create the image */
  if (properties->image)
    gdk_image_unref(properties->image);

  properties->image = gdk_image_new(GDK_IMAGE_NORMAL, gdk_visual_get_system(), properties->imageWidth, properties->imageHeight);
}


/* savePlot writes a 1 byte unsigned-char at the start of the file to indicate the file format. 
 * It then writes various parameters (as determined by that format) and then the dot-matrix itself.
 *
 * The format char is one of the following:
 *   '1'  The original format, with the following fields: format (uchar), zoom (int), image width (int), image height (int) 
 *   '2'  In addition to the format '1' fields: pixelFac (int), win size (int), matrix name len (int), matrix name (MNlen chars) 
 *   '3'  as format 2 but changed zoom from int to gdouble 
 *
 * Note that formats 1 and 2 used to assume that the 'int' date type was always 4 bytes. We now make
 * sure this is the case by using the gint32 data type, which is guaranteed to be 32 bits (4 bytes)
 * on any system.  The standard for chars/uchars (1 byte) and doubles/gdoubles (8 bytes) should be 
 * consistent for all systems so we don't need to worry about those.
 */
void savePlot(GtkWidget *dotplot, DotplotProperties *propertiesIn, const char *saveFileName, GError **error)
{
  /* We may be passed the properties as null if we are given the dotplot. */
  DotplotProperties *properties = propertiesIn ? propertiesIn : dotplotGetProperties(dotplot);
  
  g_assert(properties);
  
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;

  /* This is the latest file format number. Increment this if you change anything that will break
   * compatibility with the current format. */
  unsigned char format = 3;

  gboolean batch = (saveFileName ? TRUE : FALSE); /* whether this is part of a batch process */

  
  /* Open the file. Use the given file name (i.e. if we're in batch mode) or ask the user to 
   * select a file. */
  static const char *fileName = NULL;

  if (batch)
    fileName = saveFileName;
  else if (dotplot)
    fileName = getSaveFileName(dotplot, fileName, NULL, ".dotter", "Save dot-plot in dotter format");
  
  FILE *saveFile = NULL;
  
  if (fileName) 
    {
      saveFile = fopen (fileName, "wb");
      
      if (!saveFile)
        {
          g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_OPENING_FILE, "Failed to open file '%s'.\n", fileName);
          return;
        }
    }
  else
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_OPENING_FILE, "No file name specified.\n");
      return;
    }
  
  /* Get the length of the matrix name */
  const gint32 MNlen = strlen(dc->matrixName);
  gint32 MNlenSave = MNlen;
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(gint32));
  reversebytes(&properties->imageHeight, sizeof(gint32));
  reversebytes(&pixelFac, sizeof(gint32));
  reversebytes(&slidingWinSize, sizeof(gint32));
  reversebytes(&MNlenSave, sizeof(gint32));
#endif
  
  gboolean ok = fwrite(&format, 1, sizeof(unsigned char), saveFile) == sizeof(unsigned char);
  ok &= fwrite(&dwc->zoomFactor, 1, sizeof(gdouble), saveFile) == sizeof(gdouble);
  ok &= fwrite(&properties->imageWidth, 1, sizeof(gint32), saveFile) == sizeof(gint32);
  ok &= fwrite(&properties->imageHeight, 1, sizeof(gint32), saveFile) == sizeof(gint32);
  ok &= fwrite(&properties->pixelFac,  1, sizeof(gint32), saveFile) == sizeof(gint32); /* New feature of format 2  */
  ok &= fwrite(&properties->slidingWinSize, 1, sizeof(gint32), saveFile) == sizeof(gint32); /* New feature of format 2  */
  ok &= fwrite(&MNlenSave, 1, sizeof(gint32), saveFile) == sizeof(gint32); /* New feature of format 2  */
  ok &= fwrite(dc->matrixName, sizeof(char), MNlen, saveFile) == sizeof(char) * MNlen; /* New feature of format 2  */
  
  
  /* Loop through the matrix and write the values to the file */
  int i = 0;
  int j = 0;
  gint32 mtx = 0;
  
  for (i = 0; i < 24; i++)
    {
      for (j = 0; j < 24; j++)
        {
          mtx = dc->matrix[i][j];
#ifdef ALPHA
          reversebytes(&mtx, sizeof(gint32));
#endif
          ok &= fwrite(&mtx, 1, sizeof(gint32), saveFile) == sizeof(gint32); /* New feature of format 2  */
        }
    }
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(gint32));
  reversebytes(&properties->imageHeight, sizeof(gint32));
  reversebytes(&pixelFac, sizeof(gint32));
  reversebytes(&slidingWinSize, sizeof(gint32));
#endif
  
  
  /* Loop through the dotplot values and write them to file */
  const int imgSize = properties->imageWidth * properties->imageHeight;
  
  for (i = 0; i < imgSize; i++)
    {
      unsigned char pixelVal = properties->pixelmap[i];
#ifdef ALPHA
      reversebytes(&pixelVal, sizeof(unsigned char));
#endif
      ok &= fwrite(&pixelVal, 1, sizeof(unsigned char), saveFile) == sizeof(unsigned char); 
    }
  
  
//  ok &= fwrite(properties->pixelmap, sizeof(unsigned char), imgSize, saveFile) == imgSize;
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_SAVING_FILE, "Error writing data to file '%s'.\n", fileName);
    }
  
  fclose(saveFile);
  saveFile = 0;
}


/* Export the plot to a graphical format (currently only PDF is supported by GTK).
 * the drawable can be NULL. */
void exportPlot(GtkWidget *dotplot, 
                GtkWindow *window,
                const char *exportFileName,
                GError **error)
{
  gboolean batch = (exportFileName ? TRUE : FALSE); /* whether this is part of a batch process */
  
  /* Open the file. Use the given file name (i.e. if we're in batch mode) or ask the user to 
   * select a file. */
  static const char *fileName = NULL;

  if (batch)
    fileName = exportFileName;
  else if (dotplot)
    fileName = getSaveFileName(dotplot, fileName, NULL, ".pdf", "Export dot-plot");
  
  GtkPrintSettings *printSettings = gtk_print_settings_new();
  GtkPageSetup *pageSetup = gtk_page_setup_new();
  
  gtk_page_setup_set_orientation(pageSetup, GTK_PAGE_ORIENTATION_LANDSCAPE);
  gtk_print_settings_set_orientation(printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
  gtk_print_settings_set_quality(printSettings, GTK_PRINT_QUALITY_HIGH);
  gtk_print_settings_set_resolution(printSettings, DEFAULT_PRINT_RESOLUTION);

  blxPrintWidget(dotplot, NULL, window, &printSettings, &pageSetup, fileName, TRUE, PRINT_FIT_BOTH);
}


/* Clear previous image and pixmaps and re-create them */
static void recalculateDotplotBorders(GtkWidget *dotplot, DotplotProperties *properties)
{
  DEBUG_ENTER("recalculateDotplotBorders");
  
  /* Clear any previous image and pixmaps */
  clearPixmaps(properties);

  /* Recalculate the image size */
  properties->imageWidth = getImageDimension(properties, TRUE);
  properties->imageHeight = getImageDimension(properties, FALSE);
  properties->lineLen = properties->imageWidth;
  
  DEBUG_OUT("Set image w=%d, h=%d, line len=%d\n", properties->imageWidth, properties->imageHeight, properties->lineLen);

  /* Create the image */
  if (properties->image)
    gdk_image_unref(properties->image);
    
  properties->image = gdk_image_new(GDK_IMAGE_NORMAL, gdk_visual_get_system(), properties->imageWidth, properties->imageHeight);

  /* Create the new pixmap */
  if (properties->hspMode == DOTTER_HSPS_GREYSCALE)
    {
      initPixmap(&properties->hspPixmap, properties->imageWidth, properties->imageHeight);
      transformGreyRampImage(properties->image, properties->hspPixmap, properties);
    }
  else if (properties->pixelmapOn)
    {
      initPixmap(&properties->pixelmap, properties->imageWidth, properties->imageHeight);
      calculateImage(properties);
      transformGreyRampImage(properties->image, properties->pixelmap, properties);
    }
  
  calculateDotplotBorders(dotplot, properties);
  
  DEBUG_EXIT("recalculateDotplotBorders");
}


/* Utility to get the total height required for the dotplot */
int getDotplotHeight(DotplotProperties *properties)
{
  int result = properties->plotRect.y + 
               properties->plotRect.height + 		  
               DEFAULT_Y_PADDING + 
               SCALE_LINE_WIDTH;
  
  if (properties->breaklinesOn && properties->hozLabelsOn)
    {
      /* Add space for breakline labels */
      result += properties->dotterWinCtx->dotterCtx->charHeight + (2 * ANNOTATION_LABEL_PADDING);
    }
  
  return result;
}

/* Utility to get the total width of the dotplot */
int getDotplotWidth(DotplotProperties *properties)
{
  int result = properties->plotRect.x +
               properties->plotRect.width + 
               DEFAULT_X_PADDING + 
               SCALE_LINE_WIDTH;
  
  if (properties->breaklinesOn && properties->vertLabelsOn)
    {
      /* Add space for breakline labels */
      result += (properties->dotterWinCtx->dotterCtx->charWidth * ANNOTATION_LABEL_LEN) + 
                (2 * ANNOTATION_LABEL_PADDING);
    }
  
  return result;
}


static void calculateDotplotBorders(GtkWidget *dotplot, DotplotProperties *properties)
{
  DEBUG_ENTER("calculateDotplotBorders");

  DotterContext *dc = properties->dotterWinCtx->dotterCtx;

  properties->plotRect.x = DEFAULT_X_PADDING + dc->scaleWidth + dc->charHeight;
  properties->plotRect.y = DEFAULT_Y_PADDING + dc->scaleHeight + dc->charHeight;
  properties->plotRect.width = properties->imageWidth;
  properties->plotRect.height = properties->imageHeight;
  
  const int totalWidth = getDotplotWidth(properties);
  
  const int totalHeight = getDotplotHeight(properties);
  
  if (dotplot)
    {
      gtk_layout_set_size(GTK_LAYOUT(dotplot), totalWidth, totalHeight);
      gtk_widget_set_size_request(dotplot, totalWidth, totalHeight);

      widgetClearCachedDrawable(dotplot, NULL);
      gtk_widget_queue_draw(dotplot);
    }
  
  DEBUG_EXIT("calculateDotplotBorders returning ");
}


/***********************************************************
 *                        Drawing                          *
 ***********************************************************/

/* Refresh cached drawables */
void refreshDotplot(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  gtk_widget_queue_draw(dotplot);
  gtk_widget_queue_draw(properties->hozExons1);
  gtk_widget_queue_draw(properties->hozExons2);
  gtk_widget_queue_draw(properties->vertExons1);
  gtk_widget_queue_draw(properties->vertExons2);
}


/* Clear all cached drawables and redraw everything */
void redrawDotplot(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  widgetClearCachedDrawable(properties->hozExons1, NULL);
  widgetClearCachedDrawable(properties->hozExons2, NULL);
  widgetClearCachedDrawable(properties->vertExons1, NULL);
  widgetClearCachedDrawable(properties->vertExons2, NULL);

  callFuncOnAllChildWidgets(dotplot, widgetClearCachedDrawable);

  gtk_widget_queue_draw(dotplot);
}


/* Recalculate and redraw everything */
void recalcDotplot(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  recalculateDotplotBorders(dotplot, properties);
  
  calculateDotterExonViewBorders(properties->hozExons1, properties->imageWidth, properties->imageHeight);
  calculateDotterExonViewBorders(properties->hozExons2, properties->imageWidth, properties->imageHeight);
  calculateDotterExonViewBorders(properties->vertExons1, properties->imageWidth, properties->imageHeight);
  calculateDotterExonViewBorders(properties->vertExons2, properties->imageWidth, properties->imageHeight);
  
  redrawDotplot(dotplot);
}


/* Draws the main dotplot components (does not include transient components
 * like the crosshair or rubberband). */
static void drawDotplot(GtkWidget *dotplot, GdkDrawable *drawable)
{
  drawImage(dotplot, drawable);
  drawDotterScale(dotplot, drawable);
  drawBreaklines(dotplot, drawable);
  drawHsps(dotplot, drawable);
}


/* This prepares the dotplot for printing. It draws the crosshair onto the
 * cached pixmap, because it is this pixmap that will be used in the print.
 * (Normally the crosshair does not get drawn to the pixmap - it just gets 
 * drawn to the screen by the expose function). */
void dotplotPrepareForPrinting(GtkWidget *dotplot)
{
  GdkDrawable *drawable = widgetGetDrawable(dotplot);
  
  if (!drawable)
    {
      drawable = createBlankPixmap(dotplot);
      drawDotplot(dotplot, drawable);
    }
  
  if (drawable)
    {
      dotplotDrawCrosshair(dotplot, drawable);
    }
  
  /* Do the same for each exon view */
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  exonViewPrepareForPrinting(properties->hozExons1);
  exonViewPrepareForPrinting(properties->hozExons2);
  exonViewPrepareForPrinting(properties->vertExons1);
  exonViewPrepareForPrinting(properties->vertExons2);
}


static void findScaleUnit (gdouble cutoff, int *u, int *sub)
{
  gdouble unit = *u ;
  gdouble subunit = *u ;
  
  if (cutoff < 0)
    cutoff = -cutoff ;
  
  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5000001 ;	/* safe rounding */
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }
  subunit /= 10 ;
  if (subunit > *sub)
    *sub = roundNearest(subunit) ;
  *u = roundNearest(unit) ;
}


/* Utility used by drawScaleMarkers to draw grid lines */
static void drawGridline(GdkDrawable *drawable, DotplotProperties *properties, const int xStart, const int yStart, const gboolean horizontal)
{
  if (properties->gridlinesOn)
    {
      GdkGC *gc = gdk_gc_new(drawable);

      DotterWindowContext *dwc = properties->dotterWinCtx;
      DotterContext *dc = properties->dotterWinCtx->dotterCtx;
      
      GdkColor *gridColor = getGdkColor(DOTCOLOR_GRID, dc->defaultColors, FALSE, dwc->usePrintColors);
      gdk_gc_set_foreground(gc, gridColor);
      
      int xEnd = horizontal ? xStart : xStart + properties->plotRect.width;
      int yEnd = horizontal ? yStart + properties->plotRect.height : yStart;
    
      gdk_draw_line(drawable, gc, xStart, yStart, xEnd, yEnd);
      
      g_object_unref(gc);
    }  
}


/* Utility used by drawScaleMarkers to draw labels on tick marks */
static void drawTickmarkLabel(GtkWidget *dotplot, DotterContext *dc, GdkDrawable *drawable, GdkGC *gc, const int coordIn, int x, int y, const gboolean horizontal)
{
  int coord = getDisplayCoord(coordIn, dc, horizontal);
  char *displayText = convertIntToString(coord);
  PangoLayout *layout = gtk_widget_create_pango_layout(dotplot, displayText);
  g_free(displayText);
  
  /* We'll centre the text at half the text width for the horizontal axis, or half the 
   * height for the vertical. */
  int width = UNSET_INT, height = UNSET_INT;
  pango_layout_get_pixel_size(layout, &width, &height);
  
  x -= horizontal ? width / 2 : width;
  y -= horizontal ? height : height / 2;;
  
  gdk_draw_layout(drawable, gc, x, y, layout);
  g_object_unref(layout);  
}


static void drawScaleMarkers(GtkWidget *dotplot, 
                             GdkDrawable *drawable, 
                             GdkGC *gc, 
                             const IntRange const *displayRange,
                             DotplotProperties *properties,
                             const gboolean horizontal)
{
  /* We scale down by the zoom factor, and also divide by the number of frames so that one dot
   * represents one peptide */
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  const gdouble scaleFactor = getScaleFactor(properties, horizontal);
  const gboolean reversedScale = ((horizontal && dc->hozScaleRev) || (!horizontal && dc->vertScaleRev));
  const int direction = reversedScale ? -1 : 1;
  
  int basesPerMark = 1.0;
  int basesPerSubmark = 1.0;
  gdouble cutoff = horizontal ? PIXELS_PER_MARK_X : PIXELS_PER_MARK_Y;
  cutoff *= scaleFactor;
  
  findScaleUnit(cutoff, &basesPerMark, &basesPerSubmark);
  
  const int variableBorder = horizontal ? properties->plotRect.x : properties->plotRect.y;
  const int staticBorder = horizontal 
    ? properties->plotRect.y - SCALE_LINE_WIDTH
    : properties->plotRect.x - SCALE_LINE_WIDTH;
  
  /* Round up the start coord to find the coord and position of the first submark */
  int startCoord = UNSET_INT;
  int endCoord = UNSET_INT;
  int firstMarkCoord = UNSET_INT;

  if (reversedScale)
    {
      /* Horizontal scale is reversed. 	Round the max coord down to the nearest basesPerSubmark */
      startCoord = displayRange->max;
      endCoord = displayRange->min;
      firstMarkCoord = (int)((double)startCoord / (double)basesPerSubmark) * basesPerSubmark;
    }
  else
    {
      /* Round the min coord up to the nearest basesPerSubmark */
      startCoord = displayRange->min;
      endCoord = displayRange->max;
      firstMarkCoord = ceil((double)startCoord / (double)basesPerSubmark) * basesPerSubmark;
    }

  const int firstMarkPos = variableBorder + abs(firstMarkCoord - startCoord) / scaleFactor;
  const int numSubmarks = (int)((gdouble)(abs(endCoord - firstMarkCoord) + 1) / basesPerSubmark) + 1;

  int i = 0;

  for ( ; i < numSubmarks; ++i)
    {
      const int basesFromFirstMark = i * basesPerSubmark;
      const int coord = firstMarkCoord + (direction * basesFromFirstMark);
      const gboolean isMajorTick = (coord % basesPerMark == 0);

      const int currentPos = firstMarkPos + (basesFromFirstMark / scaleFactor);
      const int tickHeight = isMajorTick ? DEFAULT_MAJOR_TICK_HEIGHT : DEFAULT_MINOR_TICK_HEIGHT;
      
      /* Draw the marker line */
      int x1 = horizontal ? currentPos : staticBorder - tickHeight;
      int x2 = horizontal ? currentPos : staticBorder;
      int y1 = horizontal ? staticBorder - tickHeight : currentPos;
      int y2 = horizontal ? staticBorder : currentPos;

      gdk_draw_line(drawable, gc, x1, y1, x2, y2);
      
      /* Draw a grid line, if applicable */
      drawGridline(drawable, properties, x2, y2, horizontal);
          
      /* Draw a lable on major tick marks */
      if (isMajorTick)
        {
	  drawTickmarkLabel(dotplot, dc, drawable, gc, coord, x1, y1, horizontal);
        }
    }
}

static void drawLabel(GtkWidget *dotplot, GdkDrawable *drawable, GdkGC *gc)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;

  PangoLayout *layout = gtk_widget_create_pango_layout(dotplot, dc->refSeqName);

  int textWidth = UNSET_INT, textHeight = UNSET_INT;
  pango_layout_get_pixel_size(layout, &textWidth, &textHeight);

  int x = properties->plotRect.x + (properties->plotRect.width / 2) - (textWidth / 2);
  int y = properties->plotRect.y - (dc->scaleHeight + dc->charHeight);

  gdk_draw_layout(drawable, gc, x, y, layout);

  g_object_unref(layout);

  /* Vertical label */
  layout = gtk_widget_create_pango_layout(dotplot, dc->matchSeqName);

  PangoMatrix pangoMtx = PANGO_MATRIX_INIT;
  pango_matrix_rotate(&pangoMtx, 90.0);

  PangoContext *pangoCtx = pango_layout_get_context(layout);
  pango_context_set_base_gravity(pangoCtx, PANGO_GRAVITY_EAST);
  pango_context_set_matrix(pangoCtx, &pangoMtx);
  
  pango_layout_get_pixel_size(layout, &textWidth, &textHeight);
  
  x = properties->plotRect.x - (dc->scaleWidth + dc->charHeight);
  y = properties->plotRect.y + (properties->plotRect.height / 2) - (textHeight / 2);
  
  gdk_draw_layout(drawable, gc, x, y, layout);
  
  g_object_unref(layout);

  /* revert context */
  pango_context_set_base_gravity(pangoCtx, PANGO_GRAVITY_SOUTH);
  PangoMatrix pangoMtxOrig = PANGO_MATRIX_INIT;
  pango_context_set_matrix(pangoCtx, &pangoMtxOrig);
}


/* Draw the outer rectangle of the dot plot, with a scale along the left and top */
static void drawDotterScale(GtkWidget *dotplot, GdkDrawable *drawable)
{
  GdkGC *gc = gdk_gc_new(drawable);

  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  /* Draw the outline rectangle of the dot plot */
  gdk_gc_set_line_attributes(gc, SCALE_LINE_WIDTH, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  
  gdk_draw_rectangle(drawable, gc, FALSE, 
                     properties->plotRect.x - SCALE_LINE_WIDTH, properties->plotRect.y - SCALE_LINE_WIDTH, 
                     properties->plotRect.width + SCALE_LINE_WIDTH, properties->plotRect.height + SCALE_LINE_WIDTH);

  drawLabel(dotplot, drawable, gc);
  
  drawScaleMarkers(dotplot, drawable, gc, &properties->dotterWinCtx->refSeqRange, properties, TRUE);
  drawScaleMarkers(dotplot, drawable, gc, &properties->dotterWinCtx->matchSeqRange, properties, FALSE);
  
  g_object_unref(gc);
}


/* Draw an individual breakline */
static void drawBreakline(const MSP const *msp, GtkWidget *dotplot, DotplotProperties *properties, GdkDrawable *drawable, GdkGC *gc)
{
  g_assert(msp && msp->type == BLXMSP_FS_SEG);
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;

  /* The q range min and max should be the same coord */
  if (msp->qRange.min != msp->qRange.max)
    g_warning("Breakline coords should be the same but min=%d and max=%d\n", msp->qRange.min, msp->qRange.max);

  /* See if this msp is the vertical or horizontal sequence. It could be both for a self comparison. */
  gboolean horizontal = (msp->qname && strcmp(msp->qname, dc->refSeqName) == 0);
  gboolean vertical = (msp->qname && strcmp(msp->qname, dc->matchSeqName) == 0);
  
  if (horizontal)
    {
      /* Find the x position of this coord and draw a vertical line here */
      int sy = properties->plotRect.y;
      int ey = properties->plotRect.y + properties->plotRect.height + ANNOTATION_LABEL_PADDING;
    
      int sx, ex; /* start/end coordinates [0..n] */
      getMspScreenCoords(msp, properties, &sx, &ex, NULL, NULL);
      gdk_draw_line(drawable, gc, sx, sy, ex, ey);
    
      /* Draw a label at the bottom (if labels enabled) */
      if (properties->hozLabelsOn && msp->desc)
	{
	  PangoLayout *layout = gtk_widget_create_pango_layout(dotplot, msp->desc);
	  gdk_draw_layout(drawable, gc, ex, ey, layout);
	  g_object_unref(layout);
	}
    }
  
  if (vertical)
    {
      /* Find the y position of this coord and draw a horizontal line here */
      int sx = properties->plotRect.x;
      int ex = properties->plotRect.x + properties->plotRect.width + ANNOTATION_LABEL_PADDING;

      int sy, ey; /* start/end coordinates [0..n]. Note that breakline coords are always
		   * stored in the q coords even for the vertical sequence. */
      getMspScreenCoords(msp, properties, &sy, &ey, NULL, NULL);
    
      /* hack the coords because they've been based at the x coord. Really we should store them
       * properly in the s range but for historic reasons they're in the q range. */
      sy += properties->plotRect.y - properties->plotRect.x;
      ey += properties->plotRect.y - properties->plotRect.x;
    
      gdk_draw_line(drawable, gc, sx, sy, ex, ey);
    
      /* Draw a label at the RHS (if labels enabled) */
      if (properties->vertLabelsOn && msp->desc)
	{
	  PangoLayout *layout = gtk_widget_create_pango_layout(dotplot, msp->desc);
	  gdk_draw_layout(drawable, gc, ex, ey, layout);
	  g_object_unref(layout);
	}
    }
}


/* If a sequence file contains multiple sequences, this draws breaklines between the segments.
 * Only does anything if the 'draw breaklines' option is enabled. */
static void drawBreaklines(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = dwc->dotterCtx;
  
  if (properties->breaklinesOn)
    {
      GdkGC *gc = gdk_gc_new(drawable);

      gdk_gc_set_line_attributes(gc, SCALE_LINE_WIDTH, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
      GdkColor *color = getGdkColor(DOTCOLOR_BREAKLINE, dc->defaultColors, FALSE, dwc->usePrintColors);
      gdk_gc_set_foreground(gc, color);

      /* Loop through all MSPs and draw any that are segment ends */
      const MSP const *msp = dc->mspList;
    
      for ( ; msp; msp = msp->next)
	{
	  if (msp->type == BLXMSP_FS_SEG)
	    {
	      drawBreakline(msp, dotplot, properties, drawable, gc);
	    }
	}
      
      g_object_unref(gc);
    }  
}


/* Draw the crosshair that indicates the position of the currently-selected coords */
static void dotplotDrawCrosshair(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  if (properties->crosshairOn)
    {
      GdkGC *gc = gdk_gc_new(drawable);

      DotterWindowContext *dwc = properties->dotterWinCtx;
      DotterContext *dc = dwc->dotterCtx;

      /* Set the line color for the crosshair */
      GdkColor *lineColor = getGdkColor(DOTCOLOR_CROSSHAIR, dc->defaultColors, FALSE, dwc->usePrintColors);
      gdk_gc_set_foreground(gc, lineColor);
      
      int x = UNSET_INT, y = UNSET_INT;
      getPosFromSelectedCoords(dotplot, &x, &y);

      /* Draw the horizontal line (y position is at the match sequence coord position). x coords
       * depend on whether it's across the whole widget or just the dot-plot rectangle */
      int x1 = properties->crosshairFullscreen ? dc->scaleWidth: properties->plotRect.x;
      int width = properties->crosshairFullscreen ? dotplot->allocation.width : properties->plotRect.width;
      gdk_draw_line(drawable, gc, x1, y, x1 + width, y);

      /* Draw the vertical line (x position is at the ref sequence coord position - inverted if the display is reversed) */
      int y1 = properties->crosshairFullscreen ? dc->scaleHeight : properties->plotRect.y;
      int height = properties->crosshairFullscreen ? dotplot->allocation.height : properties->plotRect.height;
      gdk_draw_line(drawable, gc, x, y1, x, y1 + height);
      
      if (properties->crosshairCoordsOn)
        {
          /* Draw the coord text */
          char *displayText = blxprintf("%d, %d", 
					getDisplayCoord(dwc->refCoord, dc, TRUE), 
					getDisplayCoord(dwc->matchCoord, dc, FALSE));
          
          if (displayText && strlen(displayText) > 0)
            {
              PangoLayout *layout = gtk_widget_create_pango_layout(dotplot, displayText);
              
              int textWidth = UNSET_INT, textHeight = UNSET_INT;
              pango_layout_get_pixel_size(layout, &textWidth, &textHeight);

              /* We draw the label below-right by default, or above/left if it won't fit */
              if (x + CROSSHAIR_TEXT_PADDING + textWidth > properties->plotRect.x + properties->plotRect.width)
                {
                  x -= CROSSHAIR_TEXT_PADDING + textWidth;
                }
              else
                {
                  x += CROSSHAIR_TEXT_PADDING;
                }
              
              if (y + CROSSHAIR_TEXT_PADDING + textHeight > properties->plotRect.y + properties->plotRect.height)
                {
                  y -= CROSSHAIR_TEXT_PADDING + textHeight;
                }
              else
                {
                  y += CROSSHAIR_TEXT_PADDING;
                }

              gdk_draw_layout(drawable, gc, x, y, layout);
              g_object_unref(layout);
            }

          g_free(displayText);
        }
      
      g_object_unref(gc);
    }
}


/* Draw the current image */
static void drawImage(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DEBUG_ENTER("drawImage");
  
  DotplotProperties *properties = dotplotGetProperties(dotplot);
 
  if (showImage(properties))
    {  
      GdkGC *gc = gdk_gc_new(drawable);
      
      gdk_draw_image(drawable, gc, properties->image,
                     0, 0, properties->plotRect.x, properties->plotRect.y,
                     properties->image->width, properties->image->height); 

      /* We cannot calculate the image in the sliding-window size border,
       * so shade that in to indicate it is not part of the plot */
      DotterWindowContext *dwc = properties->dotterWinCtx;
      DotterContext *dc = dwc->dotterCtx;

      GdkColor *color = getGdkColor(DOTCOLOR_BORDER, dc->defaultColors, FALSE, dwc->usePrintColors);
      const double alpha = 1;
      const int hBorder = properties->slidingWinSize / (2 * getScaleFactor(properties, TRUE));
      const int vBorder = properties->slidingWinSize / (2 * getScaleFactor(properties, FALSE));

      cairo_operator_t op = CAIRO_OPERATOR_OVER;
      GdkRectangle *r = &properties->plotRect;

      drawRect(drawable, color, r->x, r->y, vBorder, r->height, alpha, op);
      drawRect(drawable, color, r->x + r->width - vBorder, r->y, vBorder, r->height, alpha, op);
      drawRect(drawable, color, r->x, r->y, r->width, hBorder, alpha, op);
      drawRect(drawable, color, r->x, r->y + r->height - hBorder, r->width, hBorder, alpha, op);      

      g_object_unref(gc);
    }
  
  DEBUG_EXIT("drawImage returning ");
}


/* Utility to cut off anything before the ':' in an MSP name. Returns the full name if 
 * there is no colon. Returns a pointer into the original name, so the result should not 
 * be free'd */
static const char* getShortMspName(const MSP const *msp)
{
  const char *result = NULL;
  const char *sName = mspGetSName(msp);
  
  if (sName)
    {
      char *mspName = strchr(sName, ':');  
      result = mspName ? ++mspName : sName;
    }
  
  return result;
}


static void calculateImageHsps(int strength, int sx, int sy, int ex, int ey, DotplotProperties *properties)
{
  /* Get zero-based coords from the edge of the drawing rectangle */
  sx -= properties->plotRect.x;
  ex -= properties->plotRect.x;
  sy -= properties->plotRect.y;
  ey -= properties->plotRect.y;
  
  const int xInc = (sx < ex) ? 1 : -1;
  const int xLen = (ex - sx) * xInc + 1;
  const int yInc = (sy < ey) ? 1 : -1;
  const int yLen = (ey - sy) * yInc + 1;
  
  const int width = properties->image->width;
  const int height = properties->image->height;
  const int pixelmapLen = width * height;
  
  unsigned char dotValue;
  
  if (strength < NUM_COLORS)
    dotValue = (unsigned char)strength;
  else 
    dotValue = NUM_COLORS - 1;
  
  int i = 0;
  for (i = 0; i < xLen && i < yLen; ++i)
    {
      const int x = sx + i * xInc;
      const int y = sy + i * yInc;
      
      if (x >= 0 && x < width && y >= 0 && y < height) 
	{
          const int dotpos = width * (y) + x;

          /* Use the new value if greater than the existing value. Ignore pixels out of range. */
          if (dotpos >= 0 && dotpos < pixelmapLen && dotValue > properties->hspPixmap[dotpos]) 
            {
              properties->hspPixmap[dotpos] = dotValue;
            }
	}
    }
}


/* Utility to get the screen coords for the start and end coords of the given MSP */
static void getMspScreenCoords(const MSP const *msp, DotplotProperties *properties, int *sx, int *ex, int *sy, int *ey)
{
  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));

  const int qStart = msp->qRange.min;
  const int qEnd = msp->qRange.max;
  const int sStart = sameDirection ? msp->sRange.min : msp->sRange.max;
  const int sEnd = sameDirection ? msp->sRange.max : msp->sRange.min;

  getPosFromCoords(properties, qStart, sStart, sx, sy);
  getPosFromCoords(properties, qEnd, sEnd, ex, ey);
}
  

static void drawHsps(GtkWidget *dotplot, GdkDrawable *drawable)
{
  GdkGC *gc = gdk_gc_new(drawable);

  DotplotProperties *properties = dotplotGetProperties(dotplot);

  if (properties->hspMode != DOTTER_HSPS_LINE && properties->hspMode != DOTTER_HSPS_FUNC)
    {
      /* Greyscale hsps are drawn as a pixmap via the drawImage function, so do nothing here */
      return;
    }
  
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* we'll clip the hsp lines to the dotplot drawing area */
  gdk_gc_set_clip_origin(gc, 0, 0);
  gdk_gc_set_clip_rectangle(gc, &properties->plotRect);
  
  /* Loop through all MSPs */
  MSP *msp = dc->mspList;
  int sx, ex, sy, ey; /* start/end coordinates [0..n] */
  
  for (msp = dc->mspList; msp;  msp = msp->next)
    {    
      const char *mspName = getShortMspName(msp);
      if (!mspName || strcmp(mspName, dc->matchSeqName) != 0)
        {
          /* Not an MSP from our sequence */
          continue;
        }
        
      getMspScreenCoords(msp, properties, &sx, &ex, &sy, &ey);
      
      if (properties->hspMode == DOTTER_HSPS_LINE) 
        {
          /* Draw in red */
          GdkColor hspLineColor;
          gdk_color_parse("red", &hspLineColor);
          gdk_colormap_alloc_color(gdk_colormap_get_system(), &hspLineColor, TRUE, TRUE);
          gdk_gc_set_foreground(gc, &hspLineColor);

          gdk_draw_line(drawable, gc, sx, sy, ex, ey);
        }
      else if (properties->hspMode == DOTTER_HSPS_FUNC)
        {
          /* Get color based on score: */
          GdkColor scoreColor;
          
          if (msp->score < 75.0)
            {
              gdk_color_parse("dark red", &scoreColor);
            }
          if (msp->score < 100.0) 
            {
              gdk_color_parse("magenta", &scoreColor);
            }
          else
            {
              gdk_color_parse("red", &scoreColor);
            }
          
          gdk_colormap_alloc_color(gdk_colormap_get_system(), &scoreColor, TRUE, TRUE);
          gdk_gc_set_foreground(gc, &scoreColor);
          
          gdk_draw_line(drawable, gc, sx, sy, ex, ey);
        }
    }
  
  g_object_unref(gc);
}


/* If middle-dragging, draw a rubber-band box around the currently-selected region */
static void drawRubberBand(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  if (properties->dragStart.x != UNSET_INT && properties->dragEnd.x != UNSET_INT)
    {
      GdkGC *gc = gdk_gc_new(drawable);
      
      gdk_gc_set_function(gc, GDK_INVERT);
      
      gdk_draw_line(drawable, gc, properties->dragStart.x, properties->dragStart.y, properties->dragEnd.x, properties->dragStart.y);
      gdk_draw_line(drawable, gc, properties->dragStart.x, properties->dragStart.y, properties->dragStart.x, properties->dragEnd.y);
      gdk_draw_line(drawable, gc, properties->dragStart.x, properties->dragEnd.y, properties->dragEnd.x, properties->dragEnd.y);
      gdk_draw_line(drawable, gc, properties->dragEnd.x, properties->dragStart.y, properties->dragEnd.x, properties->dragEnd.y);
      
      g_object_unref(gc);
    }
}


/***********************************************************
 *                       Internal routines                 *
 ***********************************************************/

static GdkColormap *insertGreyRamp (DotplotProperties *properties)
{
  gboolean success[NUM_COLORS];

  GdkVisual *visual = gdk_visual_get_system();
  GdkColormap *cmap = gdk_colormap_get_system();
  
  if (visual->type == GDK_VISUAL_PSEUDO_COLOR) 
    {
      g_critical("Pseudo-color display not supported.\n");
    }
  
  /* true color display */
  int i = 0;
  
  for (i = 0; i < NUM_COLORS; i++)
    {
      properties->greyRamp[i].red = i<<8;
      properties->greyRamp[i].green =  i<<8;
      properties->greyRamp[i].blue = i<<8;
    }
  
  gdk_colormap_alloc_colors(cmap, properties->greyRamp, NUM_COLORS, FALSE, TRUE, success);
  
  return cmap;
}


static void transformGreyRampImage(GdkImage *image, unsigned char *pixmap, DotplotProperties *properties)
{
  DEBUG_ENTER("transformGreyRampImage");

  if (!pixmap)
    {
      g_warning("Pixelmap is NULL; image will not be drawn.\n");
      return;
    }
    
  /* Note1 : here we stick to client byte-order, and rely on Xlib to swap
   if the server is different */
  
  /* Note2: The exact function of this routine depends on the visual. If it's
   a pseudo color visual, we're just getting the pixel values for the 
   greyramp colormap entries. These values don't change, but the colors
   they represent do. If it's a true-colour visual, the map entries 
   are changed to give the correct values. It so happens that the 
   transformation here is the same in both cases, this is not true 
   for changing the grey-ramp. */
  
  int row, col;
#if G_BYTE_ORDER == G_BIG_ENDIAN
  gboolean byterev = (image->byte_order == GDK_LSB_FIRST);
#else
  gboolean byterev = (image->byte_order == GDK_MSB_FIRST);
#endif
  
  /* Switch on the number of bytes per pixel */
  switch (image->bpp)
  { 
    case 1:
      for (row = 0; row < image->height; row++)
	{ 
	  guint8 *ptr = ((guint8 *)image->mem) + row * image->bpl;
	  guint8 *sptr = ((guint8 *)pixmap) + row * properties->lineLen; 
	  for (col = 0 ; col < image->width; col++)
	    *ptr++ = (guint8) properties->greyMap[*sptr++];
	}
      break;
    case 2:
      for (row = 0; row < image->height; row++)
	{ 
	  guint16 *ptr = (guint16 *)(((guint8 *)(image->mem))+row*image->bpl);
	  guint8 *sptr = ((guint8 *)pixmap) +row*properties->lineLen;
	  if (byterev)
	    for (col = 0 ; col < image->width; col++)
	      { 
		register guint32 pixel = properties->greyMap[*sptr++];
		*ptr++ = ((pixel & 0xff00) >> 8) | ((pixel & 0xff) << 8);
	      }
	  else
	    for (col = 0 ; col < image->width; col++)
	      *ptr++ = (guint16) properties->greyMap[*sptr++];
	}
      break;
    case 3:
      for (row = 0; row < image->height; row++)
	{ 
	  guint8 *ptr = ((guint8 *)image->mem) + row*image->bpl;
	  guint8 *sptr = ((guint8 *)pixmap) +row*properties->lineLen;
	  for (col = 0 ; col < image->width; col++)
	    {
	      register guint32 pixel = properties->greyMap[*sptr++]; 
	      *ptr++ = (guint8)pixel;
	      *ptr++ = (guint8)(pixel>>8);
	      *ptr++ = (guint8)(pixel>>16); 
	    }
	}
      break;
    case 4:
      for (row = 0; row < image->height; row++)
	{ 
	  guint32 *ptr = (guint32 *)(((guint8 *)image->mem) + row*image->bpl);
	  guint8 *sptr = ((guint8 *)pixmap) +row*properties->lineLen;
	  if (byterev)
	    for (col = 0 ; col < image->width; col++)
	      { 
		register guint32 pixel = properties->greyMap[*sptr++];
		*ptr++ = 
                ((pixel & 0xff000000) >> 24) |
                ((pixel & 0xff0000) >> 8) |
                ((pixel & 0xff00) << 8) |
                ((pixel & 0xff) << 24);
	      }
	  else
	    for (col = 0 ; col < image->width; col++)
	      *ptr++ = (guint32) properties->greyMap[*sptr++];
	}
      break;
  }

  DEBUG_EXIT("transformGreyRampImage returning ");
}


/* Update the dot plot crosshair after a change to the selected coords */
void dotplotUpdateCrosshair(GtkWidget *dotplot, DotterWindowContext *dwc)
{
  gtk_widget_queue_draw(dotplot);
}


/* Update the dot-plot image following a change in the greyramp */
void dotplotUpdateGreymap(GtkWidget *dotplot, gpointer data)
{
  DEBUG_ENTER("dotplotUpdateGreymap");

  DotplotProperties *properties = dotplotGetProperties(dotplot);
  unsigned char *ramp = (unsigned char*)data;
  
  /* First calculate the new mapping from weight to pixels */
  int i = 0;
  for (i = 0; i < NUM_COLORS; i++)
    {
      GdkColor *color = &properties->greyRamp[ramp[i]];
      properties->greyMap[i] = color->pixel;
    }
  
  /* Now recreate the image from the pixmap. Use the HSP pixmap if enabled, otherwise
   * use the standard pixelmap. */ 
  if (properties->hspMode == DOTTER_HSPS_GREYSCALE)
    {
      transformGreyRampImage(properties->image, properties->hspPixmap, properties);
    }
  else if (properties->pixelmap)
    {
      transformGreyRampImage(properties->image, properties->pixelmap, properties);
    }
  
  widgetClearCachedDrawable(dotplot, NULL);
  gtk_widget_queue_draw(dotplot);
  
  DEBUG_EXIT("dotplotUpdateGreymap returning ");
}


/* Get the reference/match sequence coords at the given x/y position */
static void getCoordsFromPos(GtkWidget *dotplot, const int x, const int y, 
			     int *refCoord, int *matchCoord)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* For ref seq, convert to number of nucleotides */
  gdouble scaleFactor = getScaleFactor(properties, TRUE);
  const int numBasesHoz = (x - properties->plotRect.x) * scaleFactor;
  
  if (dc->hozScaleRev)
    {
      *refCoord = dwc->refSeqRange.max - numBasesHoz;
    }
  else
    {
      *refCoord = dwc->refSeqRange.min + numBasesHoz;
    }

  /* Round to nearest whole pixel and limit to valid range */
  *refCoord = roundToValue(*refCoord, getResFactor(dc, TRUE));
  boundsLimitValue(refCoord, &dwc->refSeqRange);

  scaleFactor = getScaleFactor(properties, FALSE);
  const int numBasesVert = (y - properties->plotRect.y) * scaleFactor;  
  
  if (dc->vertScaleRev)
    {
      *matchCoord = dwc->matchSeqRange.max - numBasesVert;
    }
  else
    {
      *matchCoord = dwc->matchSeqRange.min + numBasesVert;
    }
  
  /* Round to nearest whole pixel and limit to valid range */
  *matchCoord = roundToValue(*matchCoord, getResFactor(dc, FALSE));
  boundsLimitValue(matchCoord, &dwc->matchSeqRange);
}


/* Set the currently-selected coords from the given x/y position */
static void setCoordsFromPos(GtkWidget *dotplot, const int x, const int y)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  getCoordsFromPos(dotplot, x, y, &dwc->refCoord, &dwc->matchCoord);
  
  refreshDotplot(dotplot);
}  


/* Get the x,y position of the given coords */
static void getPosFromCoords(DotplotProperties *properties, int qCoord, int sCoord, int *x, int *y)
{
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  const gdouble hScaleFactor = getScaleFactor(properties, TRUE);
  const gdouble vScaleFactor = getScaleFactor(properties, FALSE);

  if (x)
    {
      if (dc->hozScaleRev)
	*x = properties->plotRect.x + (dwc->refSeqRange.max - qCoord) / hScaleFactor;
      else
	*x = properties->plotRect.x + (qCoord - dwc->refSeqRange.min) / hScaleFactor;
    }
  
  if (y)
    {
      if (dc->vertScaleRev)
	*y = properties->plotRect.y + (dwc->matchSeqRange.max - sCoord) / vScaleFactor;
      else
	*y = properties->plotRect.y + (sCoord - dwc->matchSeqRange.min) / vScaleFactor;
    }
}


/* Get the x,y position of the currently selected coords */
static void getPosFromSelectedCoords(GtkWidget *dotplot, int *x, int *y)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  getPosFromCoords(properties, dwc->refCoord, dwc->matchCoord, x, y);
}


/* If the given point is outside the given rect, move it to the edge of the rect */
static void boundsLimitPoint(GdkPoint *point, GdkRectangle *rect)
{
  if (point->x < rect->x)
    point->x = rect->x;
  
  if (point->x > rect->x + rect->width)
    point->x = rect->x + rect->width;
  
  if (point->y < rect->y)
    point->y = rect->y;
  
  if (point->y > rect->y + rect->height)
    point->y = rect->y + rect->height;
  
}


/* Utility to set the values of a GdkPoint and optionally bounds limit it to the given rectangle */
static void setPoint(GdkPoint *point, const int x, const int y, GdkRectangle *rect)
{
  point->x = x;
  point->y = y;
  
  if (rect)
    boundsLimitPoint(point, rect);
}


/* Utility to get the scale factor for the horizontal/vertical sequence. Includes zoom and conversion
 * between nucleotide/display coords if applicable. */
static gdouble getScaleFactor(DotplotProperties *properties, const gboolean horizontal)
{
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  gdouble result = properties->dotterWinCtx->zoomFactor * (gdouble)getResFactor(dc, horizontal);
  return result;
}


/* REVERSEBYTES changes order of n bytes at location ptr
 */
#ifdef ALPHA
static void reversebytes(void *ptr, int n)
{ 
  static char copy[256], *cp;
  int  i;
  
  cp = ptr;
  memcpy(copy, ptr, n);  /* Note: strcpy doesn't work - stops at \0 */
  
  for(i=0; i<n; i++) *cp++ = copy[n-i-1];
}
#endif


