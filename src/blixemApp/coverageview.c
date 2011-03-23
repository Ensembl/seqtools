/*  File: coverageview.c
 *  Author: Gemma Barson, 2011-03-21
 *  Copyright (c) 2011 Genome Research Ltd
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
 * Description: This view shows a plot of the number of alignments at each
 *              coordinate in the reference sequence.
 *----------------------------------------------------------------------------
 */

#include "blixemApp/coverageview.h"
#include "blixemApp/bigpicture.h"
#include "blixemApp/blixem_.h"
#include <gtk/gtk.h>


#define DEFAULT_COVERAGE_VIEW_Y_PADDING		10	  /* this provides space between the drawing area and the edge of the widget */
#define DEFAULT_NUM_V_CELLS			2	  /* number of vertical cells to show on the grid */
#define MIN_LINE_WIDTH				0.5	  /* this provides space between the drawing area and the edge of the widget */
#define COVERAGE_VIEW_NAME                      "CoverageView"


typedef struct _CoverageViewProperties
  {
    GtkWidget *blxWindow;   /* The main blixem window */

    int viewYPadding;	     /* The y padding around the view rect */
    int numVCells;	     /* The number of cells to show vertically */
    gdouble rangePerCell;    /* The range of depth values shown per grid cell on the plot */
    
    GdkRectangle viewRect;   /* The rectangle we draw in */
    GdkRectangle displayRect; /* The total display area */
    GdkRectangle highlightRect; /* The area that is highlighted (which indicates the detail-view range) */
  } CoverageViewProperties;




/***********************************************************
 *                       Properties                        *
 ***********************************************************/

CoverageViewProperties* coverageViewGetProperties(GtkWidget *widget)
{
  return widget ? (CoverageViewProperties*)(g_object_get_data(G_OBJECT(widget), "CoverageViewProperties")) : NULL;
}

static void onDestroyCoverageView(GtkWidget *widget)
{
  CoverageViewProperties *properties = coverageViewGetProperties(widget);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", NULL);
    }
}

static void coverageViewCreateProperties(GtkWidget *widget, 
                                         GtkWidget *blxWindow,
					 BlxViewContext *bc)
{
  if (widget)
    { 
      CoverageViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->blxWindow = blxWindow;
      properties->viewYPadding = DEFAULT_COVERAGE_VIEW_Y_PADDING;
      properties->numVCells = DEFAULT_NUM_V_CELLS;
      properties->rangePerCell = 0;
      
      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCoverageView), NULL); 
    }
}


/* This function should be called whenever the coverage depth data has changed */
void updateCoverageDepth(GtkWidget *coverageView, BlxViewContext *bc)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  properties->rangePerCell = (gdouble)bc->maxDepth / properties->numVCells;
  coverageViewRecalculate(coverageView);
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Clear the cached drawable and re-draw the coverage view */
void coverageViewRedraw(GtkWidget *coverageView)
{
  widgetClearCachedDrawable(coverageView, NULL);
  gtk_widget_queue_draw(coverageView);
}


/* Recalculate the size of the coverage view widget and redraw */
void coverageViewRecalculate(GtkWidget *coverageView)
{
  calculateCoverageViewBorders(coverageView);
  coverageViewRedraw(coverageView);
}


/* Draw a single bar of the coverage bar graph */
static void drawCoverageBar(const double x1,
			    const double x2,
			    double y1,
			    const double y2,
			    cairo_t *cr)
{
  double height = y2 - y1;
  double width = x2 - x1;

  if (height <= MIN_LINE_WIDTH)
    {
      height = MIN_LINE_WIDTH;
      y1 = y2 - height;
    }
  
  if (width <= 0)
    return;
  
  cairo_rectangle(cr, x1 - MIN_LINE_WIDTH, y1, width + (2 * MIN_LINE_WIDTH), height);
  cairo_fill(cr);
}


/* Draw the actual coverage data as a bar chart */
static void drawCoveragePlot(GtkWidget *coverageView, GdkDrawable *drawable)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  
  if (!bc || bc->maxDepth <= 0)
    return;
  
  const GdkColor *color = getGdkColor(BLXCOLOR_COVERAGE_PLOT, bc->defaultColors, FALSE, bc->usePrintColors);
  cairo_t *cr = gdk_cairo_create(drawable);
  gdk_cairo_set_source_color(cr, color);
  
  const double pixelsPerVal = (double)properties->viewRect.height / (double)bc->maxDepth;
  const int bottomBorder = properties->viewRect.y + properties->viewRect.height;
  
  /* Loop through each coord in the display range */
  const IntRange const *displayRange = bigPictureGetDisplayRange(bigPicture);
  
  double startX = -1.0;
  double prevX = -1.0;
  double prevY = -1.0;
  int coord = displayRange->min;
  
  for ( ; coord <= displayRange->max; ++coord)
    {
      /* Get the x position for this coord (always pass displayRev as false because
       * display coords are already inverted if the display is reversed). */
      const double x = convertBaseIdxToRectPos(coord, &properties->viewRect, displayRange, TRUE, FALSE, TRUE);
      
      /* Convert the display coord to a zero-based coord in the full ref seq
       * display range, for indexing the depth array. Note that we need to
       * un-invert the display coord if the display is reversed. */
      const int idx = invertCoord(coord, &bc->fullDisplayRange, bc->displayRev);
      const int depth = bc->depthArray[idx - bc->fullDisplayRange.min];

      /* Calculate the y position based on the depth */
      const double height = (pixelsPerVal * (double)depth);
      const double y = (double)bottomBorder - height;

      /* First time round, don't draw the line - just get the starting x and y.
       * Also, don't draw the line if y is the same as prevY; drawing is quite
       * expensive, so if we have multiple y values that are the same we remember
       * the first x coord at this y value and draw from there all in one go. */
      if (prevX == -1.0 && prevY == -1.0)
        {
          startX = x;
        }
      else if (y != prevY || coord == displayRange->max)
        {
          /* If we had multiple positions where y was the same, draw a horizontal
           * line at that y position. */
          if (prevX != startX)
	    drawCoverageBar(startX, prevX, prevY, bottomBorder, cr);
	    
          /* Now draw the sloped line from the previous y to the new y. */
	  drawCoverageBar(prevX, x, y, bottomBorder, cr);

          /* Reset the starting point */
          startX = x;
        }

      prevX = x;
      prevY = y;
    }
}


/* Main function for drawing the coverage view */
static void drawCoverageView(GtkWidget *coverageView, GdkDrawable *drawable)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);

  GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);

  drawHighlightBox(drawable,
		   &properties->highlightRect, 
		   bpProperties->highlightBoxMinWidth,
		   highlightBoxColor,
                   HIGHLIGHT_BOX_DRAW_FUNC);
  
  drawVerticalGridLines(&properties->viewRect,
			&properties->highlightRect, 
			properties->viewYPadding, 
			bc,
			bpProperties, 
			drawable);
  
  drawHorizontalGridLines(coverageView, bigPicture, &properties->viewRect, bc, bpProperties, drawable,
			  properties->numVCells, properties->rangePerCell, (gdouble)bc->maxDepth, "");
  
  drawCoveragePlot(coverageView, drawable);
}


/* Calculate the borders of the highlight box (the shaded area that indicates the
 * detail-view range). (This is just a convenience way to call calculateHighlightBoxBorders
 * from an external function.) */
void calculateCoverageViewHighlightBoxBorders(GtkWidget *coverageView)
{  
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);

  calculateHighlightBoxBorders(&properties->viewRect, &properties->highlightRect, bigPicture, 0);
}


/* Calculate the borders of the view */
void calculateCoverageViewBorders(GtkWidget *coverageView)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
  
  /* Get some info about the size of the layout */
  guint layoutWidth, layoutHeight;
  gtk_layout_get_size(GTK_LAYOUT(coverageView), &layoutWidth, &layoutHeight);
  
  /* Calculate the height based on the number of cells */
  if (properties->rangePerCell > 0)
    properties->numVCells = (gdouble)bc->maxDepth / properties->rangePerCell;
  else
    properties->numVCells = DEFAULT_NUM_V_CELLS;
  
  const int height = properties->numVCells * bigPictureGetCellHeight(bigPicture);
  
  properties->displayRect.x = 0;
  properties->displayRect.y = 0;
  properties->displayRect.width = coverageView->allocation.width;
  properties->displayRect.height = height + 2 * (bpProperties->highlightBoxYPad + DEFAULT_COVERAGE_VIEW_Y_PADDING);
  
  /* Get the boundaries of the drawing area */
  properties->viewRect.x = roundNearest(bpProperties->charWidth * (gdouble)bpProperties->leftBorderChars);
  properties->viewRect.y = bpProperties->highlightBoxYPad + DEFAULT_COVERAGE_VIEW_Y_PADDING;
  properties->viewRect.width = properties->displayRect.width - properties->viewRect.x;
  properties->viewRect.height = height;
  
  /* Get the boundaries of the highlight box */
  calculateHighlightBoxBorders(&properties->viewRect, &properties->highlightRect, bigPicture, 0);
  
  /* Set the size request to our desired height. We want a fixed heigh but don't set the
   * width, because we want the user to be able to resize that. */
  gtk_widget_set_size_request(coverageView, 0, properties->displayRect.height);
}



/***********************************************************
 *                         Events                          *
 ***********************************************************/

/* Expose handler. */
static gboolean onExposeCoverageView(GtkWidget *coverageView, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(coverageView)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(coverageView);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(coverageView);
          drawCoverageView(coverageView, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc2 = gdk_gc_new(window);
          gdk_draw_drawable(window, gc2, bitmap, 0, 0, 0, 0, -1, -1);
          
          /* Draw the preview box on top, if it is set */
          CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
          GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
          drawPreviewBox(bigPicture, window, &properties->viewRect, &properties->highlightRect, HIGHLIGHT_BOX_DRAW_FUNC);
        }
      else
	{
	  g_warning("Failed to draw coverageView [%p] - could not create bitmap.\n", coverageView);
	}
    }
  
  return TRUE;
}


static void onSizeAllocateCoverageView(GtkWidget *coverageView, GtkAllocation *allocation, gpointer data)\
{
  calculateCoverageViewBorders(coverageView);
}


static gboolean onButtonPressCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  return FALSE;
}


static gboolean onButtonReleaseCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  return FALSE;
}


static gboolean onMouseMoveCoverageView(GtkWidget *coverageView, GdkEventMotion *event, gpointer data)
{
  return FALSE;
}

/***********************************************************
 *                     Initialisation                      *
 ***********************************************************/

GtkWidget* createCoverageView(GtkWidget *blxWindow, BlxViewContext *bc)
{
  GtkWidget *coverageView = gtk_layout_new(NULL, NULL);

  /* Style properties */
  gtk_widget_set_redraw_on_allocate(coverageView, FALSE);
  gtk_widget_set_name(coverageView, COVERAGE_VIEW_NAME);

  /* Connect signals */
  gtk_widget_add_events(coverageView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(coverageView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(coverageView, GDK_POINTER_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(coverageView), "expose-event",          G_CALLBACK(onExposeCoverageView), NULL);  
  g_signal_connect(G_OBJECT(coverageView), "size-allocate",	    G_CALLBACK(onSizeAllocateCoverageView),           NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-press-event",    G_CALLBACK(onButtonPressCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-release-event",  G_CALLBACK(onButtonReleaseCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "motion-notify-event",   G_CALLBACK(onMouseMoveCoverageView),	      NULL);
  
  /* Set required data in the coverageView. */
  coverageViewCreateProperties(coverageView, blxWindow, bc);
  
  
  return coverageView;
}

