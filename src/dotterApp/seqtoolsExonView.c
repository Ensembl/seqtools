/*  File: seqtoolsExonView.c
 *  Author: Gemma Barson, 2009-12-24
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
 * Description: See seqtoolsExonView.h
 *----------------------------------------------------------------------------
 */

#include <dotterApp/seqtoolsExonView.h>
#include <dotterApp/dotter_.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <stdlib.h>

typedef struct _ExonViewProperties
  {
    GtkWidget *parent;		      /* The parent widget that this view belongs to */
    GtkCallback refreshFunc;	      /* Callback function to call on the parent when we require a refresh */
    DotterContext *dc;		      /* Dotter session context */
    DotterWindowContext *dwc;	      /* Dotter window context */
  
    BlxStrand strand;                 /* Which strand of the sequence this view displays exons for */
    gboolean horizontal;              /* Whether these exons are for the horizontal or vertical sequence  */
    const IntRange const *qRange;     /* the range of ref seq coords the exon view displays, in nucleotide coords */
    
    gboolean bumped;		      /* Whether the exon view is expanded (bumped) or compressed */
    int yPad;			      /* y padding */
    
    GdkRectangle exonViewRect;	      /* The drawing area for the exon view */
    int exonHeight;                   /* the height of an individual exon */ 
    
    gboolean showCrosshair;           /* Flag that indicates whether we should draw the crosshair */
  } ExonViewProperties;


/* Utility struct to pass around data required for drawing exons */
typedef struct _DrawData
  {
    GtkWidget *parent;
    GdkDrawable *drawable;
    GdkGC *gc;
    DotterContext *dc;

    GdkRectangle *exonViewRect;
    const IntRange const *qRange;
    
    int yPad;			      /* y padding */
    int y;			      /* y position to draw this exon at (constant if view compressed; increases if view bumped) */
    int height;			      /* height of exon box */

    const BlxStrand strand;           /* which strand of the sequence these exons are on */
    gboolean horizontal;              /* whether these exons are for the horizontal or vertical sequence */
    gboolean bumped;		      /* whether exon view is expaned or compressed */
    gboolean normalOnly;	      /* Only draw "normal" MSPs, i.e. thse that are not selected and are not in a group */
  } DrawData;



/* Local function declarations */
static ExonViewProperties*	exonViewGetProperties(GtkWidget *exonView);


/***********************************************************
 *                       Utility functions                 *
 ***********************************************************/

/* Calls the given function (passed as the data pointer) on the given widget 
 * if it is an exon view, or, if it is a container, calls the function on any
 * children/grandchildren/etc that are exon views */
//static void callFuncOnAllChildExonViews(GtkWidget *widget, gpointer data)
//{
//  const gchar *name = gtk_widget_get_name(widget);
//  if (strcmp(name, SEQTOOLS_EXON_VIEW_NAME) == 0)
//    {
//      GtkCallback func = (GtkCallback)data;
//      func(widget, NULL);
//    }
//  else if (GTK_IS_CONTAINER(widget))
//    {
//      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllChildExonViews, data);
//    }
//}


/* Draw an exon */
static void drawExon(const MSP const *msp, 
                     DrawData *data, 
                     const BlxSequence *blxSeq, 
                     const gboolean isSelected, 
                     const gint x, 
                     const gint y,
                     const gint widthIn,
                     const gint height)
{
  /* Clip to the drawing area (or just beyond, so we don't get end lines where the exon doesn't really end) */
  gint xStart = x;
  gint xEnd = x + widthIn;
  gint width = widthIn;
  
  const gint xMin = data->exonViewRect->x;
  const gint xMax = data->exonViewRect->x + data->exonViewRect->width;
  
  if (xStart <= xMax && xEnd >= xMin)
    {
      if (xStart < xMin)
        {
          xStart = xMin - 1;
          width = xEnd - xStart;
        }

      if (xEnd > xMax)
        {
          xEnd = xMax + 1;
          width = xEnd - xStart;
        }

      /* Draw the fill rectangle */
      const GdkColor *fillColor = mspGetColor(msp, data->dc->defaultColors, blxSeq, isSelected, data->dc->usePrintColors, TRUE, DOTCOLOR_EXON_FILL, DOTCOLOR_EXON_LINE, DOTCOLOR_CDS_FILL, DOTCOLOR_CDS_LINE, DOTCOLOR_UTR_FILL, DOTCOLOR_UTR_LINE);
      gdk_gc_set_foreground(data->gc, fillColor);
      gdk_draw_rectangle(data->drawable, data->gc, TRUE, xStart, y, width, height);
      
      /* Draw outline (exon box outline always the same (unselected) color; only intron lines change when selected) */
      const GdkColor *lineColor = mspGetColor(msp, data->dc->defaultColors, blxSeq, isSelected, data->dc->usePrintColors, FALSE, DOTCOLOR_EXON_FILL, DOTCOLOR_EXON_LINE, DOTCOLOR_CDS_FILL, DOTCOLOR_CDS_LINE, DOTCOLOR_UTR_FILL, DOTCOLOR_UTR_LINE);
      gdk_gc_set_foreground(data->gc, lineColor);
      gdk_draw_rectangle(data->drawable, data->gc, FALSE, xStart, y, width, height);
    }
}
  

/* Utility to actually draw the line for an intron. Clips it if necessary, maintaining the same
 * angle for the line */
static void drawIntronLine(DrawData *data, const gint x1, const gint y1, const gint x2, const gint y2, GdkRectangle *clipRect)
{
  /* Only draw anything if at least part of the line is within range. We are only ever called with
   * y values that are in range so don't bother checking them. */
  const gint xMax = clipRect->x + clipRect->width;
  
  if (x1 <= xMax && x2 >= clipRect->x)
    {
      int xStart = x1;
      int xEnd = x2;
      int yStart = y1;
      int yEnd = y2;
      
      /* Clip the start/end x values if out of range */
      if (xStart < clipRect->x)
        {
          const int origWidth = abs(xEnd - xStart);
        
          xStart = clipRect->x;

          const int newWidth = abs(xEnd - xStart);
          const int newHeight = roundNearest((double)(yEnd - yStart) * (double)newWidth / (double)origWidth); /* negative if yend < ystart */
          
          yStart = yEnd - newHeight;
        }

      if (xEnd > xMax)
        {
          const int origWidth = abs(xEnd - xStart);

          xEnd = xMax;

          const int newWidth = abs(xEnd - xStart);
          const int newHeight = roundNearest((double)(yEnd - yStart) * (double)newWidth / (double)origWidth);
          
          yEnd = yStart + newHeight;
        }
        
      gdk_draw_line(data->drawable, data->gc, xStart, yStart, xEnd, yEnd);
    }
}


/* Draw an intron */
static void drawIntron(const MSP const *msp, 
                       DrawData *data, 
                       const BlxSequence *blxSeq, 
                       const gboolean isSelected, 
                       const gint x, 
                       const gint y, 
                       const gint width, 
                       const gint height)
{
  const GdkColor *lineColor = mspGetColor(msp, data->dc->defaultColors, blxSeq, isSelected, data->dc->usePrintColors, FALSE, DOTCOLOR_EXON_FILL, DOTCOLOR_EXON_LINE, DOTCOLOR_CDS_FILL, DOTCOLOR_CDS_LINE, DOTCOLOR_UTR_FILL, DOTCOLOR_UTR_LINE);
  gdk_gc_set_foreground(data->gc, lineColor);

  int yTop = y;
  int yBottom = y + roundNearest((double)height / 2.0);
  
  /* Draw the first section, from the given x to the mid point, sloping up */
  int xStart = x;
  int xEnd = x + roundNearest((double)width / 2.0);
  drawIntronLine(data, xStart, yBottom, xEnd, yTop, data->exonViewRect);

  /* Draw the second section, from the mid point to the end, sloping down */
  xStart  = xEnd;
  xEnd = x + width;
  drawIntronLine(data, xStart, yTop, xEnd, yBottom, data->exonViewRect);
}


/* Draw the given exon/intron, if it is in range. Returns true if it was drawn */
static gboolean drawExonIntron(const MSP *msp, DrawData *data, const gboolean isSelected, const BlxSequence *blxSeq)
{
  gboolean drawn = FALSE;
  
  if (rangesOverlap(&msp->qRange, data->qRange))
    {
      drawn = TRUE;

      /* The grid pos for coords gives the left edge of the coord, so draw to max + 1 to be inclusive */
      const int qStart = msp->qRange.min;
      const int qEnd = msp->qRange.max + 1;
      
      const gint x1 = convertBaseIdxToRectPos(qStart, data->exonViewRect, data->qRange, data->dc->hozScaleRev, FALSE);
      const gint x2 = convertBaseIdxToRectPos(qEnd, data->exonViewRect, data->qRange, data->dc->hozScaleRev, FALSE);
      
      gint x = min(x1, x2);
      gint width = abs(x1 - x2);
      gint y = data->y + data->yPad;
      gint height = data->height;
      
      if (mspIsExon(msp))
	{
	  drawExon(msp, data, blxSeq, isSelected, x, y, width, height);
	}
      else if (mspIsIntron(msp))
	{
	  drawIntron(msp, data, blxSeq, isSelected, x, y, width, height);
	}
    }
    
  return drawn;
}


/* Returns true if the given MSP should be shown in this exon view */
static gboolean showMspInExonView(const MSP *msp, DrawData *drawData)
{
  /* Check it's an exon or intron */
  gboolean showMsp = mspIsExon(msp) || mspIsIntron(msp);
  
  /* Check it's in a visible layer */
  showMsp &= mspLayerIsVisible(msp);
  
  /* Check it's the correct strand */
  showMsp &= (mspGetRefStrand(msp) == drawData->strand);

  /* Check it's for the correct sequence */
  if (drawData->horizontal)
    showMsp &= stringsEqual(msp->qname, drawData->dc->refSeqName, FALSE);
  else
    showMsp &= stringsEqual(msp->qname, drawData->dc->matchSeqName, FALSE);
  
  return showMsp;
}

/* Draw the msps in the given sequence, if they are exons/introns. Use the color
 * specified in the user data */
static void drawExonIntronItem(gpointer listItemData, gpointer data)
{
  const BlxSequence *seq = (const BlxSequence*)listItemData;
  DrawData *drawData = (DrawData*)data;

  const gboolean isSelected = FALSE; /* selections not implemented in dotter yet */
  gboolean seqDrawn = FALSE;
  
  if (!drawData->normalOnly || !isSelected)
    {
      /* Loop through all msps in this sequence */
      GList *mspListItem = seq->mspList;
  
      for ( ; mspListItem; mspListItem = mspListItem->next)
	{
	  MSP *msp = (MSP*)(mspListItem->data);
      
	  if (showMspInExonView(msp, drawData))
	    {
	      seqDrawn |= drawExonIntron(msp, drawData, isSelected, seq);
	    }
	}
    }
  
  /* If the view is bumped, increase the y-coord for the next sequence */
  if (seqDrawn && drawData->bumped)
    {
      drawData->y += drawData->height + (2 * drawData->yPad) ;
    }
}


static void drawCrosshair(GtkWidget *exonView, GdkDrawable *drawable, GdkGC *gc)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  
  if (properties->showCrosshair)
    {
      DotterContext *dc = properties->dc;
      DotterWindowContext *dwc = properties->dwc;
      
      int x = UNSET_INT;
      const gdouble scaleFactor = dwc->zoomFactor * getResFactor(dc, TRUE);
      
      if (dc->hozScaleRev)
        {
          x = properties->exonViewRect.x + (properties->qRange->max - dwc->refCoord) / scaleFactor;
        }
      else
        {
          x = properties->exonViewRect.x + (dwc->refCoord - properties->qRange->min) / scaleFactor;
        }
      
      GdkColor *color = getGdkColor(DOTCOLOR_CROSSHAIR, dc->defaultColors, FALSE, dc->usePrintColors);
      gdk_gc_set_foreground(gc, color);
      gdk_draw_line(drawable, gc, x, properties->exonViewRect.y, x, properties->exonViewRect.y + properties->exonViewRect.height);
    }
}


/* Draw the exon view */
static void drawExonView(GtkWidget *exonView, GdkDrawable *drawable)
{
  DEBUG_ENTER("drawExonView");
  
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  DotterContext *dc = properties->dc;
  
  /* Set a clip rectangle for drawing the exons and introns (because they are drawn "over the
   * edges" to make sure intron lines have the correct slope etc.) */
  GdkGC *gc = gdk_gc_new(drawable);
  gdk_gc_set_clip_origin(gc, 0, 0);
  gdk_gc_set_clip_rectangle(gc, &properties->exonViewRect);
  
  /* Draw the exons and introns. Since we could have a lot of them in the loop, extract all the
   * info we need now and pass it around so we don't have to look for this stuff each time. */
  
  DrawData drawData = {
    properties->parent,
    drawable,
    gc,
    properties->dc,
    
    &properties->exonViewRect,
    properties->qRange,
    
    properties->yPad,
    properties->exonViewRect.y,
    properties->exonHeight,

    properties->strand,
    properties->horizontal,
    properties->bumped,
    FALSE
  };
  
  /* If the view is compressed (i.e. exons will overlap each other), then
   * only draw "normal" MSPs the first time round, and draw grouped/selected
   * MSPs afterwards, so that they appear on top. If the view is bumped, 
   * we can draw them all in a single loop, because they will not overlap. */
  drawData.normalOnly = !properties->bumped;
  
  /* Loop through all sequences, drawing all msps that are exons/introns */
  GList *seqList = dc->seqList;
  g_list_foreach(seqList, drawExonIntronItem, &drawData);

  if (!properties->bumped)
    {
      /* we don't currently support selected sequences in dotter */
    }

  /* Set the height based on the height of the exons that were actually drawn */
//  const int newHeight = drawData.y - properties->exonViewRect.y + drawData.yPad;
//  gtk_layout_set_size(GTK_LAYOUT(exonView), exonView->allocation.width, newHeight);

  g_object_unref(gc);
  DEBUG_EXIT("drawExonView returning ");
}


void calculateDotterExonViewHeight(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  DotterContext *dc = properties->dc;

  /* Calculate the height based on how many exon lines will actually be drawn */
  int numExons = 0;
  int maxExons = properties->bumped ? UNSET_INT : 1; /* unset means no limit */
  
  /* Loop through all sequences */
  GList *seqItem = dc->seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      /* Loop through all msps */
      const BlxSequence *seq = (BlxSequence*)(seqItem->data);
      GList *mspItem = seq->mspList;
      
      for ( ; mspItem; mspItem = mspItem->next)
	{
	  const MSP *msp = (const MSP*)(mspItem->data);
	  
	  if ((mspIsExon(msp) || mspIsIntron(msp)) && 
              mspGetRefStrand(msp) == properties->strand &&
	      rangesOverlap(&msp->qRange, properties->qRange))
            {
              ++numExons;
              break; /* break inner loop and move to next sequence */
            }
	}
      
      /* Break after we've found the maximum number of lines, if a max is specified */
      if (maxExons != UNSET_INT && numExons >= maxExons)
	{
	  break;
	}
    }
  
  properties->exonViewRect.height = (numExons * (properties->exonHeight + properties->yPad)) + (2 * properties->yPad);
  
  gtk_widget_set_size_request(exonView, -1, properties->exonViewRect.height);
}


void calculateDotterExonViewBorders(GtkWidget *exonView, const int width)
{
  DEBUG_ENTER("calculateDotterExonViewBorders(width=%d)", width);

  ExonViewProperties *properties = exonViewGetProperties(exonView);
  
  /* Calculate the area where the exon view will be drawn */
  properties->exonViewRect.x = properties->dc->scaleWidth; /* use same border as dotplot */
  properties->exonViewRect.y = 0;
  properties->exonViewRect.width = width;
  properties->exonViewRect.height = properties->exonHeight + (2 * properties->yPad);

  gtk_layout_set_size(GTK_LAYOUT(exonView), properties->exonViewRect.x + properties->exonViewRect.width, properties->exonViewRect.y + properties->exonViewRect.height);
  gtk_widget_set_size_request(exonView, properties->exonViewRect.x + properties->exonViewRect.width, properties->exonViewRect.y + properties->exonViewRect.height);
  
  widgetClearCachedDrawable(exonView, NULL);
  gtk_widget_queue_draw(exonView);
  
  DEBUG_EXIT("calculateDotterExonViewBorders returning ");
}

/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static ExonViewProperties* exonViewGetProperties(GtkWidget *exonView)
{
  return exonView ? (ExonViewProperties*)(g_object_get_data(G_OBJECT(exonView), "ExonViewProperties")) : NULL;
}

static void onDestroyExonView(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(exonView), "ExonViewProperties", NULL);
    }
}

static void exonViewCreateProperties(GtkWidget *exonView, 
				     GtkWidget *parent, 
				     GtkCallback refreshFunc,
                                     const gboolean horizontal,
				     const BlxStrand strand,
				     DotterContext *dc,
				     DotterWindowContext *dwc,
				     const int width,
				     const IntRange const *qRange,
                                     const gboolean showCrosshair)
{
  if (exonView)
    {
      ExonViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->parent	      = parent;
      properties->refreshFunc	      = refreshFunc;
      properties->dc		      = dc;
      properties->dwc		      = dwc;

      properties->strand	      = strand;
      properties->horizontal          = horizontal;
      properties->qRange	      = qRange;
      
      properties->bumped	      = FALSE;
      properties->yPad		      =	DEFAULT_EXON_YPAD;
      
      properties->exonViewRect.x      = 0;
      properties->exonViewRect.y      = DEFAULT_EXON_YPAD;
      properties->exonViewRect.width  = width;
      properties->exonViewRect.height = DEFAULT_EXON_HEIGHT;

      properties->exonHeight          = DEFAULT_EXON_HEIGHT;
      properties->showCrosshair       = showCrosshair;
      
      gtk_widget_set_size_request(exonView, 0, DEFAULT_EXON_HEIGHT + (2 * DEFAULT_EXON_YPAD));

      g_object_set_data(G_OBJECT(exonView), "ExonViewProperties", properties);
      g_signal_connect(G_OBJECT(exonView), "destroy", G_CALLBACK(onDestroyExonView), NULL);
    }
}

gboolean exonViewGetBumped(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  return properties->bumped;
}

/* Set whether the view is bumped or not */
void exonViewSetBumped(GtkWidget *exonView, const gboolean bumped)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  properties->bumped = bumped;

  if (bumped)
    {
      properties->yPad = DEFAULT_EXON_YPAD_BUMPED;
      properties->exonHeight = DEFAULT_EXON_HEIGHT_BUMPED;
    }
  else 
    {
      properties->yPad = DEFAULT_EXON_YPAD;
      properties->exonHeight = DEFAULT_EXON_HEIGHT;
    }
  
  calculateDotterExonViewHeight(exonView);

  /* Refresh the parent */
  if (properties->refreshFunc)
    properties->refreshFunc(properties->parent, NULL);
}


/* Toggle whether the view is bumped or not */
void exonViewToggleBumped(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  exonViewSetBumped(exonView, !properties->bumped);
}

void exonViewSetShowCrosshair(GtkWidget *exonView, const gboolean showCrosshair)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  properties->showCrosshair = showCrosshair;
}


/***********************************************************
 *                       Events                            *
 ***********************************************************/

static gboolean onExposeExonView(GtkWidget *exonView, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *drawable = widgetGetDrawable(exonView);
  
  if (!drawable)
    {
      /* Create a pixmap and draw the exon view onto it */
      drawable = createBlankPixmap(exonView);
      drawExonView(exonView, drawable);
    }
  
  if (drawable)
    {  
      /* Push the pixmap onto the screen */
      GdkDrawable *window = GTK_LAYOUT(exonView)->bin_window;

      GdkGC *gc = gdk_gc_new(window);
      gdk_draw_drawable(window, gc, drawable, 0, 0, 0, 0, -1, -1);

      /* Draw the crosshair over the top */
      drawCrosshair(exonView, window, gc);
    }


  return TRUE;
}

static void onSizeAllocateExonView(GtkWidget *exonView, GtkAllocation *allocation, gpointer data)
{
//  calculateDotterExonViewBorders(exonView);
}


static gboolean onButtonPressExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  return handled;
}


static gboolean onButtonReleaseExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  return TRUE;
}


static gboolean onMouseMoveExonView(GtkWidget *exonView, GdkEventMotion *event, gpointer data)
{
  return TRUE;
}


/***********************************************************
 *                       Initialisation                    *
 ***********************************************************/

/* Create the part of the view that will show the exons. Returns the outermost container of the
 * exon and sets the exonViewOut output arg to the actual layout widget for the exon view. */
GtkWidget *createDotterExonView(GtkWidget *parent, 
				GtkCallback refreshFunc,
                                const gboolean horizontal,
				const BlxStrand strand, 
				DotterWindowContext *dwc,
				const int width,
				const IntRange const *qRange,
                                const gboolean showCrosshair,
                                GtkWidget **exonViewOut)
{
  DEBUG_ENTER("createDotterExonView(width=%d, qRange=%d %d)", width, qRange->min, qRange->max);

  DotterContext *dc = dwc->dotterCtx;
  
  GtkWidget *exonView = gtk_layout_new(NULL, NULL);
  gtk_widget_set_name(exonView, SEQTOOLS_EXON_VIEW_NAME);

  if (exonViewOut)
    *exonViewOut = exonView;
  
  /* Connect signals */
  gtk_widget_add_events(exonView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(exonView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(exonView, GDK_POINTER_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(exonView),	"expose-event",		G_CALLBACK(onExposeExonView),	      NULL);
  g_signal_connect(G_OBJECT(exonView),	"size-allocate",	G_CALLBACK(onSizeAllocateExonView),   NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-press-event",   G_CALLBACK(onButtonPressExonView),    NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-release-event", G_CALLBACK(onButtonReleaseExonView),  NULL);
  g_signal_connect(G_OBJECT(exonView),	"motion-notify-event",  G_CALLBACK(onMouseMoveExonView),      NULL);

  exonViewCreateProperties(exonView, parent, refreshFunc, horizontal, strand, dc, dwc, width, qRange, showCrosshair);
  calculateDotterExonViewBorders(exonView, width);
  
  DEBUG_EXIT("createDotterExonView returning ");
  return exonView;
}
