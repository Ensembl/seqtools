/*
 *  dotplot.c
 *  blixem
 *
 *  Created by Gemma Barson on 08/09/2010.
 *  Copyright 2010 Sanger Institute. All rights reserved.
 *
 */

#include <SeqTools/dotter_.h>
#include <SeqTools/seqtoolsExonView.h>
#include <gtk/gtk.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <wh/aceio.h>

#define PIXELS_PER_MARK_X                           100   /* number of pixels between each major tick mark on the x scale */
#define PIXELS_PER_MARK_Y                           50    /* number of pixels between each major tick mark on the y scale */
#define CROSSHAIR_TEXT_PADDING                      5     /* padding between the crosshair and the coord display text */ 
#define NUM_COLORS                                  256   /* the number of colors in our greyscale */


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

/* binary-to-ASCII translation table */
char bton[] = "ACGTN*";

double aafq[20]	/* Amino acid residue frequencies used by S. Altschul */
= {.081, .057, .045, .054, .015, .039, .061, .068, .022, .057,
.093, .056, .025, .040, .049, .068, .058, .013, .032, .067 } ;




typedef struct _DotplotProperties
{
  DotterWindowContext *dotterWinCtx;  /* pointer to the dotter context for the window that this tool belongs to */
  GtkWidget *exonView1;		      /* forward strand exon view */
  GtkWidget *exonView2;		      /* reverse strand exon view */
  GdkRectangle plotRect;	      /* the space where the dot plot will be */
  
  gulong greyMap[NUM_COLORS];     /* maps weight -> pixel value. fixed mapping in pseudo colour displays
                                     variable mapping in truecolor displays */
  GdkColor greyRamp[NUM_COLORS];  /* 256 grey colors, black->white, only used in true color displays */
  GdkColormap *colorMap;          /* the greyramp colormap */
  
  int imageWidth;
  int imageHeight;
  GdkImage *image;                /* the greyramp image */
  int lineLen;                    /* line length of the image */
  
  
  double expResScore;
  int pixelFac;
  int slidingWinSize;
  
  /* Dynamic properties: */
  unsigned char *pixelmap;        /* source data for drawing the dot-plot */
  unsigned char *hspPixmap;        /* source data for drawing the HSP dot-plot */

  gboolean crosshairOn;           /* whether to show the crosshair that marks the position of the currently-selected coord */
  gboolean crosshairCoordsOn;     /* whether to display the crosshair label */
  gboolean crosshairFullscreen;   /* whether to show the crosshair over the whole widget or just within the dot-plot rectangle */

  gboolean pixelmapOn;            /* whether to show the dot-plot pixelmap or not */
  DotterHspMode hspMode;          /* how (and whether) to show high-scoring pairs from Blast */
  
  gboolean gridlinesOn;           /* whether to show grid lines */
  
  GdkPoint dragStart;             /* start point for mid-click drag */
  GdkPoint dragEnd;		  /* end point for mid-click drag */
} DotplotProperties;


/* Local function declarations */
static void                       drawDotterScale(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawImage(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawCrosshair(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawHsps(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       drawRubberBand(GtkWidget *dotplot, GdkDrawable *drawable);
static void                       calculateDotplotBorders(GtkWidget *dotplot, DotplotProperties *properties);
static GdkColormap*               insertGreyRamp (DotplotProperties *properties);
static void                       transformGreyRampImage(GdkImage *image, unsigned char *pixmap, DotplotProperties *properties);
static void                       initPixmap(unsigned char **pixmap, const int width, const int height);
static const char*                getShortMspName(const MSP const *msp);
static void                       setHspPixmapStrength(int strength, int sx, int sy, int ex, int ey, DotplotProperties *properties);
static void                       getMspScreenCoords(const MSP const *msp, DotplotProperties *properties, int *sx, int *ex, int *sy, int *ey);
static void                       setCoordsFromPos(GtkWidget *dotplot, const int x, const int y);
static void                       getCoordsFromPos(GtkWidget *dotplot, const int x, const int y, int *refCoord, int *matchCoord);
static void			  getPosFromCoords(GtkWidget *dotplot, int *x, int *y);
static void                       setPoint(GdkPoint *point, const int x, const int y, GdkRectangle *rect);
static gdouble                    getScaleFactor(DotplotProperties *properties, const gboolean horizontal);
static void                       initWindow(const char *winsizeIn, DotplotProperties *properties);
static void                       calculateImage(DotplotProperties *properties);


/***********************************************************
 *                          Properties                     *
 ***********************************************************/

static DotplotProperties* dotplotGetProperties(GtkWidget *widget)
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

static void dotplotCreateProperties(GtkWidget *widget,
				    DotterWindowContext *dwc,
                                    const gboolean hspsOn)
{
  if (widget)
    {
      DotplotProperties *properties = g_malloc(sizeof *properties);
      
      properties->dotterWinCtx = dwc;
      properties->exonView1 = NULL;
      properties->exonView2 = NULL;
      
      properties->plotRect.x = 0;
      properties->plotRect.y = 0;
      properties->plotRect.width = 0;
      properties->plotRect.height = 0;
      
      properties->colorMap = insertGreyRamp(properties);
      gtk_widget_set_default_colormap(properties->colorMap);

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

      properties->dragStart.x = UNSET_INT;
      properties->dragStart.y = UNSET_INT;
      properties->dragEnd.x = UNSET_INT;
      properties->dragEnd.y = UNSET_INT;
      
      g_object_set_data(G_OBJECT(widget), "DotplotProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyDotplot), NULL); 

    }
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

void dotplotSetSlidingWinSize(GtkWidget *dotplot, const int newValue, GError **error)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  if (newValue <= 0)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_INVALID_WIN_SIZE, "Sliding window size must be greater than 0.\n");
    }
  else
    {
      properties->slidingWinSize = newValue;
      widgetClearCachedDrawable(dotplot, NULL);
      gtk_widget_queue_draw(dotplot);
    }
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
              if (strcmp(mspName, dc->matchSeqName))
                {
                  /* Not an MSP from our match sequence */
                  continue;
                }
              
              int sx, ex, sy, ey;
              getMspScreenCoords(msp, properties, &sx, &ex, &sy, &ey);
              
              /* Draw as greyscale pixmap */
              const int strength = (int)msp->score;
              setHspPixmapStrength(strength, sx, sy, ex, ey, properties);
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
  
  exonViewSetShowCrosshair(properties->exonView1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->exonView2, properties->crosshairOn && properties->crosshairFullscreen);
  
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
  
  exonViewSetShowCrosshair(properties->exonView1, properties->crosshairOn && properties->crosshairFullscreen);
  exonViewSetShowCrosshair(properties->exonView2, properties->crosshairOn && properties->crosshairFullscreen);

  refreshDotplot(dotplot);
}

/***********************************************************
 *			       Events                      *
 ***********************************************************/

/* Expose handler for dot-plot window */
static gboolean onExposeDotplot(GtkWidget *dotplot, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(dotplot)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(dotplot);

      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(dotplot);
          drawImage(dotplot, bitmap);
          drawDotterScale(dotplot, bitmap);
          drawHsps(dotplot, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          
          /* Draw anything else that needs to be refreshed on each expose */
          drawCrosshair(dotplot, window);
          drawRubberBand(dotplot, window);
        }
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
              callDotterInternal(dc, &qRange, &sRange, zoomFactor) ;
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
                                           const gboolean hspsOn,
                                           const char *initWinsize,
                                           const int pixelFacIn,
                                           const int zoomFacIn,
                                           const int qcenter,
                                           const int scenter)
{
  DEBUG_ENTER("createDotplotDrawingArea");

  GtkWidget *dotplot = gtk_layout_new(NULL, NULL);
  
  dotplotCreateProperties(dotplot, dwc, hspsOn);
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
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
    
      if (properties->hspMode == DOTTER_HSPS_GREYSCALE)
        {
          /* Initialise the HSPs pixmap */
          initPixmap(&properties->hspPixmap, properties->imageWidth, properties->imageHeight);
        }
      else if (properties->pixelmapOn)
        {
          /* Initialise and populate the pixmap */
          initPixmap(&properties->pixelmap, properties->imageWidth, properties->imageHeight);
          calculateImage(properties);
        }
    }
  
  /* Create the image */
  properties->image = gdk_image_new(GDK_IMAGE_NORMAL, gdk_visual_get_system(), properties->imageWidth, properties->imageHeight);
  
  if (saveFileName)
    {
      GError *error = NULL;
      savePlot(dotplot, saveFileName, &error);
      
      prefixError(error, "Error saving dot plot. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  else
    {
      initCrosshairCoords(qcenter, scenter, properties->dotterWinCtx);
    }      
  
  
  calculateDotplotBorders(dotplot, properties);
  
  DEBUG_EXIT("createDotplotDrawingArea returning ");
  return dotplot;
}


/* Create the drawing area for the dot plot. The 'dotplot' output arg is populated with the
 * actual drawing area; the return arg is the overall container widget for the dot plot that will
 * be packed into the dotter window UNLESS we're running in batch mode, in which case the return 
 * will be no widgets displayed to the user. */
GtkWidget* createDotplot(DotterWindowContext *dwc, 
                         const char *loadFileName,
                         const char *saveFileName,
                         const gboolean hspsOn,
                         const char *initWinsize,
                         const int pixelFacIn,
                         const int zoomFacIn,
                         const int qcenter,
                         const int scenter,
                         GtkWidget **dotplot)
{
  DEBUG_ENTER("createDotplot");

  /* Create the actual drawing area for the dot plot */
  *dotplot = createDotplotDrawingArea(dwc, loadFileName, saveFileName, hspsOn, initWinsize, pixelFacIn, zoomFacIn, qcenter, scenter);

  if (saveFileName)
    {
      /* Don't create the container widgets */
      return NULL;
    }

  /* Put the dotplot in a scrolled window. It only needs the vertical scrollbar because the 
   * container including the exon views etc as well will have a horizontal scrollbar for everything. */
  GtkWidget *dotplotScrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(dotplotScrollWin), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(dotplotScrollWin), *dotplot);
  
  /* Create labels for the axes */
  DotterContext *dc = dwc->dotterCtx;
  GtkWidget *hozLabel = gtk_label_new(dc->refSeqName);
  GtkWidget *vertLabel = gtk_label_new(dc->matchSeqName);
  gtk_label_set_angle(GTK_LABEL(vertLabel), 90);

  /* Create the exon views (one for each strand, with the passed ref seq strand first) */
  DotplotProperties *properties = dotplotGetProperties(*dotplot);
  const gboolean showCrosshair = properties->crosshairOn && properties->crosshairFullscreen;
  
  GtkWidget *exon1 = createDotterExonView(*dotplot, NULL, dc->refSeqStrand, dc, dwc, properties->imageWidth, &dwc->refSeqRange, showCrosshair, &properties->exonView1);
  GtkWidget *exon2 = createDotterExonView(*dotplot, NULL, dc->refSeqStrand == BLXSTRAND_FORWARD ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD, dc, dwc, properties->imageWidth, &dwc->refSeqRange, showCrosshair, &properties->exonView2);
  
  /* Put everything in a table */
  const int xpad = 0;
  const int ypad = 0;
  const int numRows = 4;
  const int numCols = 2;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  
  gtk_table_attach(table, hozLabel, 1, 2, 0, 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, vertLabel, 0, 1, 1, 2, GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, dotplotScrollWin, 1, 2, 1, 2, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, exon1, 1, 2, 2, 3, GTK_FILL, GTK_FILL, xpad, ypad);
  gtk_table_attach(table, exon2, 1, 2, 3, 4, GTK_FILL, GTK_FILL, xpad, ypad);

//  /* Put exon views in a vbox */
//  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
//  gtk_box_pack_start(GTK_BOX(vbox), properties->exonView1, FALSE, FALSE, 0);
//  gtk_box_pack_start(GTK_BOX(vbox), properties->exonView2, FALSE, FALSE, 0);
//  
//  /* put dotplot vs exon view in a paned window */
//  GtkWidget *panedWin = gtk_vpaned_new();
//  gtk_paned_add1(GTK_PANED(panedWin), dotplotScrollWin);
//  gtk_paned_add2(GTK_PANED(panedWin), vbox);
//
//  /* put paned window in the table */
//  gtk_table_attach(table, panedWin, 1, 2, 1, 2, GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, xpad, ypad);

  /* Put the table in a scrolled window with a horizontal scrollbar */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  
  gtk_widget_add_events(*dotplot, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(*dotplot, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(*dotplot, GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(*dotplot), "expose-event", G_CALLBACK(onExposeDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "button-press-event", G_CALLBACK(onButtonPressDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "button-release-event", G_CALLBACK(onButtonReleaseDotplot), NULL);
  g_signal_connect(G_OBJECT(*dotplot), "motion-notify-event",   G_CALLBACK(onMouseMoveDotplot), NULL);

  DEBUG_EXIT("createDotplot returning ");
  return GTK_WIDGET(scrollWin);
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
          g_error("Bad window size specification: %s", winsizeIn);
        }
      
      properties->slidingWinSize = atoi(winsizeIn);
    }
}


/* Utility to allocate memory of the given size with g_malloc. Returns the allocated memory and
 * adds a pointer to it to the given list, so that all memory allocated to this list can be destroyed
 * with destroyListData */
static gpointer allocMemoryToList(GSList **dataList, size_t numBytes)
{
  gpointer result = g_malloc(numBytes);
  
  /* prepend so that we delete data in reverse order */
  *dataList = g_slist_prepend(*dataList, result); 

  return result;
}

/* Utilitiy to call g_free on all elements in the given list. */
static void destroyListData(GSList *dataList)
{
  DEBUG_ENTER("destroyListData");

  GSList *listItem = dataList;
  for ( ; listItem; listItem = listItem->next)
    {
      g_free(listItem->data);
      listItem->data = NULL;
    }
  
  g_slist_free(dataList);
  
  DEBUG_EXIT("destroyListData returning ");
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
      const gboolean complement = (dc->refSeqStrand == BLXSTRAND_REVERSE && dc->refSeqType == BLXSEQ_DNA);

      result = getRefSeqBase(dc->refSeq, coord, complement, &dc->refSeqFullRange, dc->refSeqType);
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

  /* I'm not sure why but calculateImage seems to rely on the match sequence always being uncomplemented. */
  const gboolean complement = FALSE;
  
  return getRefSeqBase(dc->matchSeq, coord, complement, &dc->matchSeqFullRange, dc->matchSeqType);
}


/* Get the array of binary values for the match sequence, usd by calculateImage */
static void getMatchSeqBinaryVals(DotterWindowContext *dwc, const int slen, const int translationTable[], int *sIndex)
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
static void createScoreVec(DotterWindowContext *dwc, const int vecLen, const int qlen, GSList **dataList, int ***scoreVecPtr)
{
  DotterContext *dc = dwc->dotterCtx;
  
  *scoreVecPtr = allocMemoryToList(dataList, vecLen * sizeof(int*));
  int **scoreVec = *scoreVecPtr;

  /* Allocate memory for each row in the score vector. Each row is the length of the ref sequence */
  int rowIdx = 0;
  for ( ; rowIdx < vecLen; ++rowIdx)
    {
      scoreVec[rowIdx] = allocMemoryToList(dataList, qlen * sizeof(int));
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
          const int score = dc->matrix[rowIdx][aminoAcidId]; /* score of this base in the q seq wrt the current amino acid ID 'i' */
          
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


/* This function calculates the data that will be put into the dotplot image. It puts the max
 * diagonal for each pixel into *pixelmap */
static void calculateImage(DotplotProperties *properties)
{
  DEBUG_ENTER("calculateImage");

  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  g_assert(properties->slidingWinSize > 0);
  
  register int 
  qIdx, sIdx, qmax,     /* Loop variables */
  *newsum,	/* New sum pointer */
  *oldsum,	/* Old sum pointer */
  *delrow,	/* Pointer to scoreVec, row to remove */
  *addrow,	/* Pointer to scoreVec, row to add */
  dotpos, dotposq, dotposs;
  
  int 
  i, min, sec, frame, 
  win2 = properties->slidingWinSize/2;
  
  GSList *dataList = NULL;
  
  double speed = 17.2;  /* Speed in Mdots/seconds. SGI MIPS R10000 (clobber) */
  /* speed = 5.7;  DEC Alpha AXP 3000/700 */
  /* speed = 3.7;  SGI R4400: */
  
  const int qlen = getRangeLength(&dwc->refSeqRange);
  const int slen = getRangeLength(&dwc->matchSeqRange);

  double numDots = qlen/1e6*slen; /* total number of dots (millions) */
  
  if (dc->selfComp) 
    numDots /= 2;
  
  if (dc->blastMode == BLXMODE_BLASTN && !(dc->watsonOnly || dc->crickOnly)) 
    numDots *= 2;
  
  if (dc->blastMode == BLXMODE_BLASTX) 
    numDots *= 3;
  
  min = (int)(numDots/speed/60);
  sec = (int)(numDots/speed) - min*60;
  
  printf("%d vs. %d residues => %.2f million dots. ", qlen, slen, numDots);
  
  if (min+sec >= 2) {
    printf("(Takes ");
    
    if (min)
      printf("%d:%.2d minutes", min, sec);
    else 
      printf("%d seconds", sec);
    
    printf(" on an SGI MIPS R10000)");
  }
  
  printf("\n");
  fflush(stdout);

  /* Initialize lookup tables for faster execution */
  int *sIndex = allocMemoryToList(&dataList, slen * sizeof(int)) ; /* match sequence as binary values (i.e. amino-acid IDs 0 -> 23) */
  int **scoreVec = NULL;                                           /* array of precalculated scores for qseq residues */
  
  /* Find the offset of the current display range within the full range of the bit of reference sequence we have */
  const int qOffset = dc->refSeqStrand == BLXSTRAND_REVERSE ? dc->refSeqFullRange.max - dwc->refSeqRange.max : dwc->refSeqRange.min - dc->refSeqFullRange.min;
  
  /* Convert from nucleotides to peptides, if applicable */
  const int resFactor = (dc->blastMode == BLXMODE_BLASTX ? dc->numFrames : 1);
  const int pepQSeqLen = qlen / resFactor;
  const int pepQSeqOffset = qOffset / resFactor;
  const int vecLen = (dc->blastMode == BLXMODE_BLASTN ? 6 : 25);

  createScoreVec(dwc, vecLen, pepQSeqLen, &dataList, &scoreVec);
  getMatchSeqBinaryVals(dwc, slen, getTranslationTable(dc->matchSeqType, BLXSTRAND_FORWARD), sIndex);
  
  /* Allocate some vectors for use in averaging the values for whole rows at a time. */
  int *zero = allocMemoryToList(&dataList, pepQSeqLen * sizeof(int));
  int *sum1 = allocMemoryToList(&dataList, pepQSeqLen * sizeof(int));
  int *sum2 = allocMemoryToList(&dataList, pepQSeqLen * sizeof(int));

  for (i = 0; i < pepQSeqLen; ++i) 
    {
      zero[i] = 0;
      sum1[i] = 0;
      sum2[i] = 0;
    }

  if (dc->blastMode == BLXMODE_BLASTX)
    {
      for (frame = 0; frame < dc->numFrames; ++frame)
        {
          /* Re-populate the score vector for this reading frame */
          populateScoreVec(dwc, vecLen, pepQSeqLen, frame, pepQSeqOffset, atob_0, scoreVec);
        
          for (sIdx = 0; sIdx < slen; ++sIdx)
            {   
              /* Set oldsum to the previous row. (newsum will be overwritten, but we re-use the
               * same two vectors here to save having to keep allocating memory) */
              oldsum = (sIdx & 1) ? sum2 : sum1;
              newsum = (sIdx & 1) ? sum1 : sum2;
              
              /* We delete the last value that was slidingWinSize away (or zero if not at slidingWinSize yet) */
              if (sIdx >= properties->slidingWinSize) 
                {
                  delrow = scoreVec[sIndex[sIdx - properties->slidingWinSize]];
                }
              else
                {
                  delrow = zero;
                }
              
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
              
              qmax = pepQSeqLen;
              
              for ( ; qIdx < qmax ; ++qIdx) 
                {
                  ++newsum;
                  *newsum = *oldsum + *addrow - *delrow;
                  ++oldsum;
                  ++addrow;
                  ++delrow;
                  
                  if (*newsum > 0 && sIdx >= properties->slidingWinSize) 
                    {
                      dotposq = (qIdx-win2)/dwc->zoomFactor;
                      dotposs = (sIdx-win2)/dwc->zoomFactor;
                      
                      /* Only fill half the submatrix */
                      const int qPosLocal = qIdx-win2 - dotposq*dwc->zoomFactor;  /* query position in local submatrix (of one pixel) */
                      const int sPosLocal = sIdx-win2 - dotposs*dwc->zoomFactor;  /* subject position in local submatrix (of one pixel) */
                      
                      if (sPosLocal >= qPosLocal)
                        {
                          dotpos = properties->imageWidth*dotposs + dotposq;
                          
                          const int pixelmapLen = properties->imageWidth * properties->imageHeight;
                          
                          if (dotpos < 0 || dotpos >= pixelmapLen) 
                            {
                              g_critical ( "Pixel out of bounds (%d) in blastx: %d\n",
                                          pixelmapLen-1, dotpos);
                            }
                          else
                            {
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
    }

  if (dc->blastMode == BLXMODE_BLASTP || (dc->blastMode == BLXMODE_BLASTN && !dc->crickOnly)) 
    {
      populateScoreVec(dwc, vecLen, pepQSeqLen, 1, 0, getTranslationTable(dc->refSeqType, BLXSTRAND_FORWARD), scoreVec);

      for (sIdx = 0; sIdx < slen; ++sIdx)
        { 
          /* Set oldsum to the previous row. (newsum will be overwritten, but we re-use the
           * same two vectors here to save having to keep allocating memory) */
          oldsum = (sIdx & 1) ? sum2 : sum1;
          newsum = (sIdx & 1) ? sum1 : sum2;
          
          if (sIdx >= properties->slidingWinSize) 
            {
              delrow = scoreVec[sIndex[sIdx-properties->slidingWinSize]];
            }
          else
            {
              delrow = zero;
            }

          addrow = scoreVec[sIndex[sIdx]];
          *newsum = *addrow++;

          
          qmax = min(properties->slidingWinSize, qlen);
          for (qIdx = 1; qIdx < qmax ; ++qIdx)
            *++newsum = *oldsum++ + *addrow++;
          
          qmax = (dc->selfComp ? sIdx + 1 : qlen);
          
          for ( ; qIdx < qmax ; ++qIdx) 
            {
              *++newsum = *oldsum++ + *addrow++ - *delrow++;
              
              if (*newsum > 0 && sIdx >= properties->slidingWinSize) 
                {
                  dotposq = (qIdx-win2)/dwc->zoomFactor;
                  dotposs = (sIdx-win2)/dwc->zoomFactor;
                  
                  /* Only fill half the submatrix */
                  const int qPosLocal = qIdx-win2 - dotposq*dwc->zoomFactor;  /* query position in local submatrix (of one pixel) */
                  const int sPosLocal = sIdx-win2 - dotposs*dwc->zoomFactor;  /* subject position in local submatrix (of one pixel) */
                  
                  if (sPosLocal >= qPosLocal)
                    {
                      dotpos = properties->imageWidth*dotposs + dotposq;
                      
                      const int pixelmapLen = properties->imageWidth * properties->imageHeight;

                      if (dotpos < 0 || dotpos > pixelmapLen-1) 
                        {
                          g_critical ( "Pixel out of bounds (%d) in blastp/blastn-forw: %d\n", 
                                      pixelmapLen-1, dotpos);
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
  
  if (dc->blastMode == BLXMODE_BLASTN && !dc->watsonOnly) 
    {
      populateScoreVec(dwc, vecLen, pepQSeqLen, 1, 0, getTranslationTable(dc->refSeqType, BLXSTRAND_REVERSE), scoreVec);

      for (i = 0; i<qlen; i++) 
        {
          sum1[i] = 0;
          sum2[i] = 0;
        }
      
      for (sIdx = slen-1; sIdx >= 0; --sIdx)
        { 
          /* Set oldsum to the previous row. (newsum will be overwritten, but we re-use the
           * same two vectors here to save having to keep allocating memory) */
          oldsum = (sIdx & 1) ? sum2 : sum1;
          newsum = (sIdx & 1) ? sum1 : sum2;
          
          if (sIdx < slen-properties->slidingWinSize) 
            {
              delrow = scoreVec[sIndex[sIdx + properties->slidingWinSize]];
            }
          else
            {
              delrow = zero;
            }
          
          addrow = scoreVec[sIndex[sIdx]];
          *newsum = *addrow++;

          qmax = min(properties->slidingWinSize, qlen);
          
          for (qIdx = 1; qIdx < qmax ; ++qIdx)
            {
              *++newsum = *oldsum++ + *addrow++;
            }
          
          qmax = (dc->selfComp ? sIdx + 1 : qlen);
          
          for ( ; qIdx < qmax ; ++qIdx) 
            {
              *++newsum = *oldsum++ + *addrow++ - *delrow++ ;
              
              if (*newsum > 0 && sIdx <= slen-properties->slidingWinSize) 
                {
                  dotposq = (qIdx - win2)/dwc->zoomFactor;
                  
                  dotposs = (sIdx + win2)/dwc->zoomFactor;
                  
                  /* Only fill half the submatrix */
                  const int qPosLocal = qIdx - win2 - dotposq*dwc->zoomFactor;  /* query position in local submatrix (of one pixel) */
                  
                  /* Set the origin (0,0) to the bottom left corner of submatrix
                   Ugly but correct. Zoom = pixels/submatrix */
                  const int sPosLocal = dwc->zoomFactor-1 - (sIdx + win2 - dotposs * dwc->zoomFactor);  /* subject position in local submatrix (of one pixel) */
                  
                  if (sPosLocal >= qPosLocal)
                    {
                      dotpos = properties->imageWidth*dotposs + dotposq;
                      const int pixelmapLen = properties->imageWidth * properties->imageHeight;

                      if (dotpos < 0 || dotpos >= pixelmapLen) 
                        {
                          g_critical ( "Pixel out of bounds (%d) in blastn-rev: %d\n",
                                      pixelmapLen-1, dotpos);
                        }
                      else 
                        {
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
  
  if (dc->selfComp && dc->displayMirror) 
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
  
  destroyListData(dataList);
  
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

  /* Read in the format, zoom, and image size */
  unsigned char format;
  gboolean ok = fread(&format, sizeof(unsigned char), 1, loadFile) == 1;
  
  ok &= fread(&dwc->zoomFactor, sizeof(gdouble), 1, loadFile) == 1;
  ok &= fread(&properties->imageWidth, sizeof(int), 1, loadFile) == 1;
  ok &= fread(&properties->imageHeight, sizeof(int), 1, loadFile) == 1;
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Error reading file '%s'", loadFileName);
      return;
    }
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(int));
  reversebytes(&properties->imageHeight, sizeof(int));
#endif
  
  properties->lineLen = properties->imageWidth;
  
  int dotstart = 0;
  
  if (format == 1) 
    {
      dotstart = 13;
      
      /* Don't actually know these variables for sure - guess the most common */
      properties->pixelFac = 50;
      properties->slidingWinSize = 25;
    }
  else if (format == 2) 
    {
      ok &= fread(&properties->pixelFac, sizeof(int), 1, loadFile) == 1;
      ok &= fread(&properties->slidingWinSize, sizeof(int), 1, loadFile) == 1;

      /* Read in the matrix name (read the length of the name first) */
      int matrixNameLen;
      ok &= fread(&matrixNameLen, sizeof(int), 1, loadFile) == 1;
      
#ifdef ALPHA
      if (ok)
        {
          reversebytes(&properties->pixelFac, sizeof(int));
          reversebytes(&properties->slidingWinSize, sizeof(int));
          reversebytes(&matrixNameLen, sizeof(int));
        }
#endif

      ok &= matrixNameLen <= MAX_MATRIX_NAME_LENGTH;
      char matrixName[MAX_MATRIX_NAME_LENGTH + 1];
      ok &= fread(&matrixName, sizeof(char), matrixNameLen, loadFile) == matrixNameLen;

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
                  int matrixVal;
                  ok &= fread(&matrixVal, sizeof(int), 1, loadFile) == 1; 
#ifdef ALPHA
                  reversebytes(&matrixVal, 4);
#endif
                  dc->matrix[i][j] = matrixVal;
                }
            }
              
          dotstart = matrixNameLen + 2329;
        }
    }
  else 
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Unknown dotter file format version: %d", format);
      return;
    }
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Error reading file '%s'", loadFileName);
      return;
    }
  
  fseek(loadFile, 0, SEEK_END);
  int n = ftell(loadFile);

  const int pixelmapLen = properties->imageHeight*properties->imageWidth;

  if (n - dotstart != properties->imageWidth*properties->imageHeight)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_READING_FILE, "Wrong number of pixels in %s: %d. Expected %d * %-d = %d\n", 
            loadFileName, n, properties->imageWidth, properties->imageHeight, pixelmapLen);
      return;
    }
  
  /* Allocate memory for the pixmap */
  initPixmap(&properties->pixelmap, properties->imageWidth, properties->imageHeight);
  
  fseek(loadFile, dotstart, SEEK_SET);
  
  int i = 0;
  for (i = 0; i < pixelmapLen; i++)
    {
      unsigned char pixelVal;
      ok &= fread(&pixelVal, sizeof(unsigned char), 1, loadFile) == 1; 
#ifdef ALPHA
      reversebytes(&pixelVal, 4);
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
}


/* Utility to ask the user for a file to save to. Returns the file name, which should be freed with
 * g_free. */
static const char* getSaveFileName(GtkWidget *widget, const char *currentName, const char *defaultPath)
{
  const char *filename = NULL;
  
  GtkWindow *window = GTK_WINDOW(gtk_widget_get_toplevel(widget));

  GtkWidget *dialog = gtk_file_chooser_dialog_new ("Save File",
                                                   window,
                                                   GTK_FILE_CHOOSER_ACTION_SAVE,
                                                   GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                                   GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
                                                   NULL);
  
  gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (dialog), TRUE);
  
  if (defaultPath)
    {
      gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog), defaultPath);
    }
  
  if (currentName)
    {
      gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER (dialog), currentName);
    }
  else
    {
      gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER (dialog), ".dotter");
    }
  
  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
    {
      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    }
  
  gtk_widget_destroy (dialog);
  return filename;
}


/* Utility to ask the user for a file to load. Returns the file name, which should be freed with
 * g_free. */
//static const char* getLoadFileName(GtkWidget *widget)
//{
//  const char *filename = NULL;
//  
//  GtkWindow *window = GTK_WINDOW(gtk_widget_get_toplevel(widget));
//  
//  GtkWidget *dialog = gtk_file_chooser_dialog_new ("Open File",
//                                                   window,
//                                                   GTK_FILE_CHOOSER_ACTION_OPEN,
//                                                   GTK_STOCK_CANCEL, 
//                                                   GTK_RESPONSE_CANCEL,
//                                                   GTK_STOCK_OPEN, 
//                                                   GTK_RESPONSE_ACCEPT,
//                                                   NULL);
//  
//  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
//    {
//      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
//    }
//  
//  gtk_widget_destroy (dialog);
//  return filename;
//}


/* SAVEPLOT writes a 1 byte fileformat first, various params and then the dot-matrix
 */
void savePlot(GtkWidget *dotplot, const char *saveFileName, GError **error)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  int 
  i, j,
  MNlen, MNlenSave,
  mtx;
  unsigned char format = 2;

  gboolean batch = (saveFileName ? TRUE : FALSE); /* whether this is part of a batch process */

  
  /* Open the file. Use the given file name (i.e. if we're in batch mode) or ask the user to 
   * select a file. */
  static const char *fileName = NULL;
  fileName = batch ? saveFileName : getSaveFileName(dotplot, fileName, NULL);
  
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
  
  MNlen = MNlenSave = strlen(dc->matrixName);
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(int));
  reversebytes(&properties->imageHeight, sizeof(int));
  reversebytes(&pixelFac, sizeof(int));
  reversebytes(&slidingWinSize, sizeof(int));
  reversebytes(&MNlenSave, sizeof(int));
#endif
  
  gboolean ok = fwrite(&format, sizeof(unsigned char), 1, saveFile) == 1;
  ok &= fwrite(&dwc->zoomFactor, sizeof(gdouble), 1, saveFile) == 1;
  ok &= fwrite(&properties->imageWidth, sizeof(int), 1, saveFile) == 1;
  ok &= fwrite(&properties->imageHeight, sizeof(int), 1, saveFile) == 1;
  ok &= fwrite(&properties->pixelFac,  sizeof(int), 1, saveFile) == 1; /* New feature of format 2  */
  ok &= fwrite(&properties->slidingWinSize, sizeof(int), 1, saveFile) == 1; /* New feature of format 2  */
  ok &= fwrite(&MNlenSave, sizeof(int), 1, saveFile) == 1; /* New feature of format 2  */
  ok &= fwrite(dc->matrixName, sizeof(char) * MNlen, 1, saveFile) == 1; /* New feature of format 2  */
  for (i = 0; i < 24; i++)
    for (j = 0; j < 24; j++) {
      mtx = dc->matrix[i][j];
#ifdef ALPHA
      reversebytes(&mtx, sizeof(int));
#endif
      ok &= fwrite(&mtx, sizeof(int), 1, saveFile) == 1; /* New feature of format 2  */
    }
  
#ifdef ALPHA
  reversebytes(&dwc->zoomFactor, sizeof(gdouble));
  reversebytes(&properties->imageWidth, sizeof(int));
  reversebytes(&properties->imageHeight, sizeof(int));
  reversebytes(&pixelFac, sizeof(int));
  reversebytes(&slidingWinSize, sizeof(int));
#endif
  
  const int imgSize = properties->imageWidth * properties->imageHeight;
  
  for (i = 0; i < imgSize; i++)
    {
      unsigned char pixelVal = properties->pixelmap[i];
#ifdef ALPHA
      reversebytes(&pixelVal, sizeof(unsigned char));
#endif
      ok &= fwrite(&pixelVal, sizeof(unsigned char), 1, saveFile) == 1; 
    }
  
  
//  ok &= fwrite(properties->pixelmap, sizeof(unsigned char), imgSize, saveFile) == imgSize;
  
  if (!ok)
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_SAVING_FILE, "Error writing data to file '%s'.\n", fileName);
    }
  
  fclose(saveFile);
  saveFile = 0;
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


static void calculateDotplotBorders(GtkWidget *dotplot, DotplotProperties *properties)
{
  DEBUG_ENTER("calculateDotplotBorders");

  DotterContext *dc = properties->dotterWinCtx->dotterCtx;

  properties->plotRect.x = DEFAULT_X_PADDING + dc->scaleWidth;
  properties->plotRect.y = DEFAULT_Y_PADDING + dc->scaleHeight;
  properties->plotRect.width = properties->image->width;
  properties->plotRect.height = properties->image->height;
  
  const int totalWidth = properties->plotRect.x + properties->plotRect.width + DEFAULT_X_PADDING + SCALE_LINE_WIDTH;
  const int totalHeight = properties->plotRect.y + properties->plotRect.height + DEFAULT_Y_PADDING + SCALE_LINE_WIDTH;
  
  gtk_layout_set_size(GTK_LAYOUT(dotplot), totalWidth, totalHeight);
  gtk_widget_set_size_request(dotplot, totalWidth, totalHeight);
  
  widgetClearCachedDrawable(dotplot, NULL);
  gtk_widget_queue_draw(dotplot);
  
  DEBUG_EXIT("calculateDotplotBorders returning ");
}


/***********************************************************
 *                        Drawing                          *
 ***********************************************************/

/* Recalculate and redraw everything */
void redrawDotplot(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  recalculateDotplotBorders(dotplot, properties);
  
  calculateDotterExonViewBorders(properties->exonView1, properties->imageWidth);
  calculateDotterExonViewBorders(properties->exonView2, properties->imageWidth);
  
  gtk_widget_queue_draw(dotplot);
}


/* Refresh the view */
void refreshDotplot(GtkWidget *dotplot)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  gtk_widget_queue_draw(dotplot);
  gtk_widget_queue_draw(properties->exonView1);
  gtk_widget_queue_draw(properties->exonView2);
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
      DotterContext *dc = properties->dotterWinCtx->dotterCtx;
      
      GdkColor *gridColor = getGdkColor(DOTCOLOR_GRID, dc->defaultColors, FALSE, dc->usePrintColors);
      GdkGC *gc = gdk_gc_new(drawable);
      gdk_gc_set_foreground(gc, gridColor);
      
      int xEnd = horizontal ? xStart : xStart + properties->plotRect.width;
      int yEnd = horizontal ? yStart + properties->plotRect.height : yStart;
    
      gdk_draw_line(drawable, gc, xStart, yStart, xEnd, yEnd);
    }  
}


/* Utility used by drawScaleMarkers to draw labels on tick marks */
static void drawTickmarkLabel(GtkWidget *dotplot, GdkDrawable *drawable, GdkGC *gc, const int baseNum, int x, int y, const gboolean horizontal)
{
  char *displayText = convertIntToString(baseNum);
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
  const int staticBorder = horizontal ? properties->plotRect.y - SCALE_LINE_WIDTH : properties->plotRect.x - SCALE_LINE_WIDTH;

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
      const int baseNum = firstMarkCoord + (direction * basesFromFirstMark);
      const gboolean isMajorTick = (baseNum % basesPerMark == 0);

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
	  drawTickmarkLabel(dotplot, drawable, gc, baseNum, x1, y1, horizontal);
        }
    }
}


/* Draw the outer rectangle of the dot plot, with a scale along the left and top */
static void drawDotterScale(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  GdkGC *gc = gdk_gc_new(drawable);
  
  /* Draw the outline rectangle of the dot plot */
  gdk_gc_set_line_attributes(gc, SCALE_LINE_WIDTH, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  
  gdk_draw_rectangle(drawable, gc, FALSE, 
                     properties->plotRect.x - SCALE_LINE_WIDTH, properties->plotRect.y - SCALE_LINE_WIDTH, 
                     properties->plotRect.width + SCALE_LINE_WIDTH, properties->plotRect.height + SCALE_LINE_WIDTH);

  drawScaleMarkers(dotplot, drawable, gc, &properties->dotterWinCtx->refSeqRange, properties, TRUE);
  drawScaleMarkers(dotplot, drawable, gc, &properties->dotterWinCtx->matchSeqRange, properties, FALSE);
}


/* Draw the crosshair that indicates the position of the currently-selected coords */
static void drawCrosshair(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  
  if (properties->crosshairOn)
    {
      DotterWindowContext *dwc = properties->dotterWinCtx;
      DotterContext *dc = dwc->dotterCtx;

      /* Set the line color for the crosshair */
      GdkColor *lineColor = getGdkColor(DOTCOLOR_CROSSHAIR, dc->defaultColors, FALSE, dc->usePrintColors);
      GdkGC *gc = gdk_gc_new(drawable);
      gdk_gc_set_foreground(gc, lineColor);
      
      int x = UNSET_INT, y = UNSET_INT;
      getPosFromCoords(dotplot, &x, &y);

      /* Calculate the total size of the dotplot, including padding and labels etc */
      const int totalWidth = properties->plotRect.x + properties->plotRect.width + DEFAULT_X_PADDING + SCALE_LINE_WIDTH;
      const int totalHeight = properties->plotRect.y + properties->plotRect.height + DEFAULT_Y_PADDING + SCALE_LINE_WIDTH;
      
      /* Draw the horizontal line (y position is at the match sequence coord position). x coords
       * depend on whether it's across the whole widget or just the dot-plot rectangle */
      int x1 = properties->crosshairFullscreen ? 0: properties->plotRect.x;
      int width = properties->crosshairFullscreen ? totalWidth : properties->plotRect.width;
      gdk_draw_line(drawable, gc, x1, y, x1 + width, y);

      /* Draw the vertical line (x position is at the ref sequence coord position - inverted if the display is reversed) */
      int y1 = properties->crosshairFullscreen ? 0 : properties->plotRect.y;
      int height = properties->crosshairFullscreen ? totalHeight : properties->plotRect.height;
      gdk_draw_line(drawable, gc, x, y1, x, y1 + height);
      
      if (properties->crosshairCoordsOn)
        {
          /* Draw the coord text */
          char *displayText = blxprintf("%d, %d", dwc->refCoord, dwc->matchCoord);
          
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
    }
  
  DEBUG_EXIT("drawImage returning ");
}


/* Utility to cut off anything before the ':' in an MSP name. Returns the full name if 
 * there is no colon. Returns a pointer into the original name, so the result should not 
 * be free'd */
static const char* getShortMspName(const MSP const *msp)
{
  char *mspName = strchr(mspGetSName(msp), ':');  
  const char *result = mspName ? ++mspName : mspGetSName(msp);
  return result;
}


static void setHspPixmapStrength(int strength, int sx, int sy, int ex, int ey, DotplotProperties *properties)
{
  const int inc = (sx < ex) ? 1 : -1;
  const int len = (ex - sx) * inc + 1;
  const int width = properties->image->width;
  const int height = properties->image->height;
  const int pixelmapLen = width * height;
  
  unsigned char dotValue;
  
  if (strength < NUM_COLORS)
    dotValue = (unsigned char)strength;
  else 
    dotValue = NUM_COLORS - 1;
  
  int i = 0;
  for (i = 0; i < len; i++)
    {
      const int x = sx + i * inc;
      const int y = sy + i;
      
      if (x >= 0 && x < width && sy >= 0 && sy < height) 
	{
          const int dotpos = width * (y) + x;
          
          if (dotpos < 0 || dotpos > pixelmapLen-1) 
            {
              g_critical("Pixel out of bounds (0-%d) in setHspPixmapStrength: %d. Crash imminent.", pixelmapLen - 1, dotpos);
            }
          else if (dotValue > properties->hspPixmap[dotpos]) 
            {
              properties->hspPixmap[dotpos] = dotValue;
            }
	}
    }
}


static gboolean illegalSubmatrix(int sx, int sy, int shift, DotplotProperties *properties)
{
  gboolean isIllegal = FALSE;
  
  int dotposq, dotposs, ql, sl;
  
  DotterWindowContext *dwc = properties->dotterWinCtx;
  dotposq = (sx+shift) / dwc->zoomFactor;
  dotposs = (sy+shift) / dwc->zoomFactor;
  
  ql = (sx+shift) - dotposq * dwc->zoomFactor;
  sl = (sy+shift) - dotposs * dwc->zoomFactor;
  
  if (sl >= ql) 
    isIllegal = FALSE;		/* legal */
  else 
    isIllegal = TRUE;		/* illegal */
  
  return isIllegal;
}


static gboolean illegalSubmatrixRev(int sx, int sy, int shift, DotplotProperties *properties)
{
  gboolean isIllegal = FALSE;

  int dotposq, dotposs, ql, sl;
  
  DotterWindowContext *dwc = properties->dotterWinCtx;
  dotposq = (sx-shift) / dwc->zoomFactor;
  dotposs = (sy+shift) / dwc->zoomFactor;
  
  ql = (sx-shift) - dotposq * dwc->zoomFactor;
  
  /* Set the origin (0,0) to the bottom left corner of submatrix
   Zoom = pixels/submatrix */
  sl = dwc->zoomFactor - 1 - (sy - dotposs * dwc->zoomFactor);
  
  if (sl >= ql) 
    isIllegal = FALSE;		/* legal */
  else 
    isIllegal = TRUE;		/* illegal */
  
  return isIllegal;
}


/* Utility to get the screen coords for the start and end coords of the given MSP */
static void getMspScreenCoords(const MSP const *msp, DotplotProperties *properties, int *sx, int *ex, int *sy, int *ey)
{
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  const int matchlen = mspGetSRangeLen(msp) / dwc->zoomFactor;
  
  if (mspGetQStart(msp) < mspGetQEnd(msp)) 
    {
      *sx = ceil((double)(mspGetQStart(msp) - (dwc->refSeqRange.min - 1)) / dc->numFrames) - 1;
      *sy = mspGetSStart(msp) - dwc->matchSeqRange.min;
      
      /* Check if we're in an illegal part of submatrix. We only want to compress in legal submatrix parts */
      int shift = 0;
      while(illegalSubmatrix(*sx, *sy, shift, properties)) 
        {
          shift++;
        }
      
      *sx = (*sx + shift) / dwc->zoomFactor -shift;
      *sy = (*sy + shift) / dwc->zoomFactor -shift;
      
      *ex = *sx + (matchlen-1);
    }
  else 
    {
      *sx = ceil((float)(mspGetQStart(msp) - (dwc->refSeqRange.min - 1)) / dc->numFrames) -1;
      *sy = mspGetSStart(msp) - dwc->matchSeqRange.min;
      
      /* Check if we're in an illegal part of submatrix. We only want to compress in legal submatrix parts */
      int shift = 0;
      while(illegalSubmatrixRev(*sx, *sy, shift, properties)) 
        {
          shift++;
        }
      
      *sx = (*sx - shift) / dwc->zoomFactor + shift;
      *sy = (*sy + shift) / dwc->zoomFactor - shift;
      
      *ex = *sx - (matchlen-1);
    }
  
  //      if (reversedScale) 
  //        {
  //          *sx = RightBorder-LeftBorder - *sx;
  //          *ex = RightBorder-LeftBorder - *ex;
  //        }
  
  *ey = *sy + matchlen - 1;
}


/* Draw HSPs, if they are enabled */
static void drawHsps(GtkWidget *dotplot, GdkDrawable *drawable)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);

  if (properties->hspMode != DOTTER_HSPS_LINE && properties->hspMode != DOTTER_HSPS_FUNC)
    {
      /* Greyscale hsps are drawn as a pixmap via the drawImage function, so do nothing here */
      return;
    }
  
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  GdkGC *gc = gdk_gc_new(drawable);
  
  /* Loop through all MSPs */
  MSP *msp = dc->mspList;
  int sx, ex, sy, ey; /* start/end coordinates [0..n] */
  
  for (msp = dc->mspList; msp;  msp = msp->next)
    {    
      const char *mspName = getShortMspName(msp);
      if (strcmp(mspName, dc->matchSeqName))
        {
          /* Not an MSP from our sequence */
          continue;
        }
      
      getMspScreenCoords(msp, properties, &sx, &ex, &sy, &ey);
      
      const int x = properties->plotRect.x;
      const int y = properties->plotRect.y;
      
      if (properties->hspMode == DOTTER_HSPS_LINE) 
        {
          /* Draw in red */
          GdkColor hspLineColor;
          gdk_color_parse("red", &hspLineColor);
          gdk_colormap_alloc_color(gdk_colormap_get_system(), &hspLineColor, TRUE, TRUE);
          gdk_gc_set_foreground(gc, &hspLineColor);

          gdk_draw_line(drawable, gc, x + sx, y + sy, x + ex, y + ey);
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
          
          gdk_draw_line(drawable, gc, x + sx, y + sy, x + ex, y + ey);
        }
    }
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
  else
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


/* Get the x,y position of the currently selected coords */
static void getPosFromCoords(GtkWidget *dotplot, int *x, int *y)
{
  DotplotProperties *properties = dotplotGetProperties(dotplot);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  const gdouble hScaleFactor = getScaleFactor(properties, TRUE);
  const gdouble vScaleFactor = getScaleFactor(properties, FALSE);

  if (dc->hozScaleRev)
    {
      *x = properties->plotRect.x + (dwc->refSeqRange.max - dwc->refCoord) / hScaleFactor;
    }
  else
    {
      *x = properties->plotRect.x + (dwc->refCoord - dwc->refSeqRange.min) / hScaleFactor;
    }

  if (dc->vertScaleRev)
    {
      *y = properties->plotRect.y + (dwc->matchSeqRange.max - dwc->matchCoord) / vScaleFactor;
    }
  else
    {
      *y = properties->plotRect.y + (dwc->matchCoord - dwc->matchSeqRange.min) / vScaleFactor;
    }
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



