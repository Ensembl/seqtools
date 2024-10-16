/*  File: coverageview.c
 *  Author: Gemma Barson, 2011-03-21
 *  Copyright [2018-2024] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

#include "blixemApp/coverageview.hpp"
#include "blixemApp/blxcontext.hpp"
#include "blixemApp/blixem_.hpp"
#include "blixemApp/blxpanel.hpp"
#include <gtk/gtk.h>
#include <math.h>

/* We should remove references to the big picture and blxwindow here. The blxwindow is only used
 * to get the big picture */
#include "blixemApp/bigpicture.hpp"
#include "blixemApp/blxwindow.hpp"


#define DEFAULT_COVERAGE_VIEW_Y_PADDING		3	  /* this provides space between the drawing area and the edge of the widget */
#define DEFAULT_NUM_V_CELLS			4	  /* number of vertical cells to show on the grid */
#define MIN_LINE_WIDTH				0.5	  /* this provides space between the drawing area and the edge of the widget */
#define COVERAGE_VIEW_NAME                      "CoverageView"



/***********************************************************
 *                    Class member functions               *
 ***********************************************************/

CoverageViewProperties::CoverageViewProperties(GtkWidget *widget_in,
                                               GtkWidget *blxWindow_in,
                                               BlxContext *bc_in)
{
  m_widget = widget_in;
  m_bc = bc_in;
  m_panel = NULL;

  m_blxWindow = blxWindow_in;

  m_viewYPadding = DEFAULT_COVERAGE_VIEW_Y_PADDING;
  m_numVCells = DEFAULT_NUM_V_CELLS;
  m_rangePerCell = 0;

  if (bc_in)
    m_maxDepth = &bc_in->maxDepth;
}

/* Get the coverage view widget */
GtkWidget* CoverageViewProperties::widget()
{
  return m_widget;
}

/* Return the range of ref-seq coords that the coverage view currently displays */
const IntRange* CoverageViewProperties::displayRange()
{
  const IntRange *result = NULL;

  // The range is the same as the parent panel
  if (m_panel)
    result = &m_panel->displayRange;

  return result;
}

/* Return the range of ref-seq coords that should be highlighted in the highlight box, if any
 * (returns NULL if not applicable) */
const IntRange* CoverageViewProperties::highlightRange()
{
  const IntRange *result = NULL;

  // The range is the same as the parent panel
  if (m_panel)
    result = m_panel->highlightRange();

  return result;
}

/* Return the position of the left border. */
double CoverageViewProperties::contentXPos() const
{
  double result = 0.0;

  // Use the left border position of the parent panel
  if (m_panel)
    result = m_panel->contentXPos();

  return result;
}

/* Return the position of the left border. */
double CoverageViewProperties::contentWidth() const
{
  double result = 0.0;

  // Use the left border position of the parent panel
  if (m_panel)
    result = m_panel->contentWidth();

  return result;
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

CoverageViewProperties* coverageViewGetProperties(GtkWidget *widget)
{
  CoverageViewProperties *properties = NULL;

  if (widget)
    properties = (CoverageViewProperties*)(g_object_get_data(G_OBJECT(widget), "CoverageViewProperties"));

  return properties;
}

static void onDestroyCoverageView(GtkWidget *widget)
{
  CoverageViewProperties *properties = coverageViewGetProperties(widget);

  if (properties)
    {
      delete properties;
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", NULL);
    }
}

static CoverageViewProperties* coverageViewCreateProperties(GtkWidget *widget,
                                                            GtkWidget *blxWindow,
                                                            BlxContext *bc)
{
  CoverageViewProperties *properties = NULL;

  if (widget)
    {
      properties = new CoverageViewProperties(widget, blxWindow, bc);

      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCoverageView), NULL);
    }

  return properties;
}


/* Set the pointer to the parent panel that this coverage view belongs to */
void CoverageViewProperties::setPanel(BlxPanel *panel)
{
  m_panel = panel;
}


/* This function should be called whenever the coverage depth data has changed */
void CoverageViewProperties::updateDepth()
{
  /* Set up a list of 'nice' values to round to for displaying labels */
  static GSList *roundValues = NULL;

  if (!roundValues)
    {
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(200));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2500));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(250000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2500000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(250000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500000000));
    }

  if (!m_rangePerCell)
    {
      /* First time round, calculate the range per cell, aiming for
       * around 5 cells. (If we enter this function again, it's because the
       * user has manually entered the range per cell so we just need to calculate
       * the relevant number of cells) */
      m_numVCells = 5;
      m_rangePerCell = ceil((gdouble)*m_maxDepth / (gdouble)m_numVCells);

      /* Round the result and recalculate the number of cells */
      m_rangePerCell = roundUpToValueFromList(m_rangePerCell, roundValues, NULL);

      if (m_rangePerCell < 1)
        m_rangePerCell = 1;
    }

  m_numVCells = (gdouble)*m_maxDepth / m_rangePerCell;

  recalculate();
}

double CoverageViewProperties::depthPerCell()
{
  return m_rangePerCell;
}

gboolean CoverageViewProperties::setDepthPerCell(const double depthPerCell_in)
{
  if (depthPerCell_in <= 0.0)
    return FALSE;

  m_rangePerCell = depthPerCell_in;
  updateDepth();

  return TRUE;
}

/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Clear the cached drawable and re-draw the coverage view */
void CoverageViewProperties::redraw()
{
  widgetClearCachedDrawable(m_widget, NULL);
  gtk_widget_queue_draw(m_widget);
}


/* Recalculate the size of the coverage view widget and redraw */
void CoverageViewProperties::recalculate()
{
  calculateBorders();
  redraw();
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


/* Utility to get the max depth that the coverage view shows, based on
 * the number of cells and the range per cell. Note that this returns the
 * max lable value, i.e. the value of the top gridline; the real max
 * depth may be slightly greater than this, and may extend above the top
 * gridline (the height of the widget is made big enough to accommodate this). */
int CoverageViewProperties::maxLabeledDepth()
{
  /* to do: ideally we would round numcells to the nearest int rather than
   * truncating, but there is a bug with that where the horizontal grid lines
   * are sometimes not drawn with the correct labels */
  int numCells = (int)(m_numVCells);
  const int result = m_rangePerCell * numCells;
  return result;
}


double CoverageViewProperties::charWidth() const
{
  double result = 0.0;

  if (m_bc)
    result = m_bc->charWidth();

  return result;
}


/* Draw the actual coverage data as a bar chart */
void CoverageViewProperties::drawPlot(GdkDrawable *drawable)
{
  const IntRange *dispRange = displayRange();
  g_return_if_fail(m_bc && m_bc->maxDepth > 0 && dispRange);

  cairo_t *cr = gdk_cairo_create(drawable);

  const GdkColor *color = getGdkColor(BLXCOLOR_COVERAGE_PLOT, m_bc->defaultColors, FALSE, m_bc->usePrintColors);
  gdk_cairo_set_source_color(cr, color);

  const double pixelsPerVal = (double)m_viewRect.height / (double)*m_maxDepth;
  const int bottomBorder = m_viewRect.y + m_viewRect.height;

  /* Loop through each coord in the display range */

  double startX = -1.0;
  double prevX = -1.0;
  double prevY = -1.0;
  int coord = dispRange->min();

  for ( ; coord <= dispRange->max(); ++coord)
    {
      /* Get the x position for this coord (always pass displayRev as false because
       * display coords are already inverted if the display is reversed). */
      const double x = convertBaseIdxToRectPos(coord, &m_viewRect, dispRange, TRUE, FALSE, TRUE);
      const int depth = m_bc->getDepth(coord);

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
      else if (y != prevY || coord == dispRange->max())
        {
          /* If we had multiple positions where y was the same, draw a horizontal
           * line at that y position. If there was only one position at the previous y value then
           * this will draw a single column in the bar chart (i.e. startX==prevX). */
          drawCoverageBar(startX, x, prevY, bottomBorder, cr);

          /* If it's the last coord, also draw the current column, because there won't be another
           * loop to take care of this */
          if (coord == dispRange->max())
            {
              const int endX = convertBaseIdxToRectPos(coord + 1, &m_viewRect, dispRange, TRUE, FALSE, TRUE);
              drawCoverageBar(x, endX, y, bottomBorder, cr);
            }

          /* Reset the starting point */
          startX = x;
        }

      prevX = x;
      prevY = y;
    }

  cairo_destroy(cr);
}


/* Main function for drawing the coverage view */
void CoverageViewProperties::draw(GdkDrawable *drawable)
{
  /* We should move the gridlines functions to BlxPanel so that we can remove references to
   * bigpicture here */
  GtkWidget *bigPicture = blxWindowGetBigPicture(m_blxWindow);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);

  drawVerticalGridLines(&m_viewRect, &m_highlightRect,
			m_viewYPadding, m_bc, bpProperties, drawable);

  drawHorizontalGridLines(m_widget, bigPicture, &m_viewRect, m_bc, bpProperties, drawable,
			  (int)(m_numVCells), m_rangePerCell, (gdouble)*m_maxDepth, TRUE, "");

  drawPlot(drawable);
}


/* Calculate the borders of the highlight box (the shaded area that indicates the
 * selection range). (This is just a convenience way to call calculateHighlightBoxBorders
 * from an external function.) */
void CoverageViewProperties::calculateHighlightBoxBorders()
{
  if (m_bc)
    {
      m_bc->highlightBoxCalcBorders(&m_displayRect, &m_highlightRect,
                                    displayRange(), highlightRange(),
                                    0);
    }
}


/* Calculate the borders of the view */
void CoverageViewProperties::calculateBorders()
{
  /* We should move cell height to BlxPanel so that we can get rid of references to big picture here */
  GtkWidget *bigPicture = blxWindowGetBigPicture(m_blxWindow);

  /* Calculate the height based on the number of cells */
  const int height = ceil(m_numVCells * (double)bigPictureGetCellHeight(bigPicture));
  const int gridHeight = (int)m_numVCells * bigPictureGetCellHeight(bigPicture);

  m_displayRect.x = roundNearest(contentXPos());
  m_displayRect.y = height - gridHeight;

  m_viewRect.x = m_displayRect.x;
  m_viewRect.y = m_displayRect.y + HIGHLIGHT_BOX_Y_PAD + DEFAULT_COVERAGE_VIEW_Y_PADDING;

  m_displayRect.width = contentWidth();
  m_displayRect.height = height + 2 * (HIGHLIGHT_BOX_Y_PAD + DEFAULT_COVERAGE_VIEW_Y_PADDING);

  /* If the width isn't set (e.g. we don't have a parent panel) then use the allocation width */
  if (m_displayRect.width == 0.0)
    m_displayRect.width = m_widget->allocation.width - m_viewRect.x;

  m_viewRect.width = m_displayRect.width;
  m_viewRect.height = gridHeight;

  /* Get the boundaries of the highlight box */
  calculateHighlightBoxBorders();

  /* Set the size request to our desired height. We want a fixed heigh but don't set the
   * width, because we want the user to be able to resize that. */
  gtk_widget_set_size_request(m_widget, 0, m_displayRect.height);
}


/* Prepare the coverage view for printing (draws the transient hightlight box
 * onto the cached drawable). */
void CoverageViewProperties::prepareForPrinting()
{
  GdkDrawable *drawable = widgetGetDrawable(m_widget);

  if (drawable)
    {
      GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, m_bc->defaultColors, FALSE, m_bc->usePrintColors);
      drawHighlightBox(drawable, &m_highlightRect, HIGHLIGHT_BOX_MIN_WIDTH, highlightBoxColor);
    }
}


/***********************************************************
 *                         Events                          *
 ***********************************************************/

/* Expose handler. */
static gboolean onExposeCoverageView(GtkWidget *coverageView, GdkEventExpose *event, gpointer data)
{
  gboolean result = TRUE;
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  if (properties)
    result = properties->expose(event, data);

  return result;
}


gboolean CoverageViewProperties::expose(GdkEventExpose *event, gpointer data)
{
  gboolean result = TRUE;

  GdkDrawable *window = GTK_LAYOUT(m_widget)->bin_window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(m_widget);

      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(m_widget);
          draw(bitmap);
        }

      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);

          /* Draw the highlight box on top of it */
          GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, m_bc->defaultColors, FALSE, m_bc->usePrintColors);
          drawHighlightBox(window, &m_highlightRect, HIGHLIGHT_BOX_MIN_WIDTH, highlightBoxColor);

          /* Draw the preview box too, if set */
          if (m_panel)
            m_panel->drawPreviewBox(window, &m_viewRect, &m_highlightRect);
        }
      else
	{
	  g_warning("Failed to draw coverageView [%p] - could not create bitmap.\n", m_widget);
	}
    }

  return result;
}


static void onSizeAllocateCoverageView(GtkWidget *coverageView, GtkAllocation *allocation, gpointer data)\
{
  CoverageViewProperties *coverageViewP = coverageViewGetProperties(coverageView);

  if (coverageViewP)
    coverageViewP->calculateBorders();
}


static gboolean onButtonPressCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  if (properties)
    handled = properties->buttonPress(event, data);

  return handled;
}


gboolean CoverageViewProperties::buttonPress(GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->button == 2 ||
      (event->button == 1 && !handled &&
       (event->type == GDK_2BUTTON_PRESS ||
        clickedInRect(event, &m_highlightRect, HIGHLIGHT_BOX_MIN_WIDTH))))
    {
      /* Draw the preview box (draw it on the other big picture components as well) */
      int x = event->x;

      if (event->button == 1 && event->type == GDK_BUTTON_PRESS)
        x = m_highlightRect.x + m_highlightRect.width / 2;

      if (m_panel)
        m_panel->startPreviewBox(event->x, TRUE, x - event->x);

      handled = TRUE;
    }

  return handled;
}


static gboolean onButtonReleaseCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  if (properties)
    handled = properties->buttonRelease(event, data);

  return handled;
}

gboolean CoverageViewProperties::buttonRelease(GdkEventButton *event, gpointer data)
{
  if (event->button == 1 || event->button == 2) /* left or middle button */
    {
      if (m_panel)
        m_panel->finishPreviewBox(event->x, &m_viewRect, &m_highlightRect);
    }

  return TRUE;
}



/* Implement custom scrolling for horizontal mouse wheel movements over the coverageView.
 * This scrolls the position of the highlight box, i.e. it scrolls the display
 * range in the detail view. */
static gboolean onScrollCoverageView(GtkWidget *coverageView, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  if (properties)
    handled = properties->scroll(event, data);

  return handled;
}


gboolean CoverageViewProperties::scroll(GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;

  /* We should make the scroll functions here virtual functions that can be called from the base
   * class BlxPanel. Then we can get rid of references to bigpicture */

  switch (event->direction)
  {
    case GDK_SCROLL_LEFT:
    {
      scrollBigPictureLeftStep(blxWindowGetBigPicture(m_blxWindow));
      handled = TRUE;
      break;
    }

    case GDK_SCROLL_RIGHT:
    {
      scrollBigPictureRightStep(blxWindowGetBigPicture(m_blxWindow));
      handled = TRUE;
      break;
    }

    default:
    {
      handled = FALSE;
      break;
    }
  };

  return handled;
}


gboolean CoverageViewProperties::mouseMove(GdkEventMotion *event, gpointer data)
{
  if ((event->state & GDK_BUTTON1_MASK) || /* left or middle button */
       (event->state & GDK_BUTTON2_MASK))
    {
      /* Draw a preview box at the mouse pointer location */
      if (m_panel)
        m_panel->startPreviewBox(event->x, FALSE, 0);
    }

  return TRUE;
}

static gboolean onMouseMoveCoverageView(GtkWidget *coverageView, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  if (properties)
    handled = properties->mouseMove(event, data);

  return handled;
}



/***********************************************************
 *                     Initialisation                      *
 ***********************************************************/

CoverageViewProperties* createCoverageView(GtkWidget *blxWindow,
                                           BlxContext *bc)
{
  GtkWidget *coverageView = gtk_layout_new(NULL, NULL);

  /* Style properties */
  gtk_widget_set_redraw_on_allocate(coverageView, FALSE);
  gtk_widget_set_name(coverageView, COVERAGE_VIEW_NAME);

  /* Connect signals */
  gtk_widget_add_events(coverageView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(coverageView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(coverageView, GDK_POINTER_MOTION_MASK);

  g_signal_connect(G_OBJECT(coverageView), "expose-event",          G_CALLBACK(onExposeCoverageView),                 NULL);
  g_signal_connect(G_OBJECT(coverageView), "size-allocate",	    G_CALLBACK(onSizeAllocateCoverageView),           NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-press-event",    G_CALLBACK(onButtonPressCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-release-event",  G_CALLBACK(onButtonReleaseCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "motion-notify-event",   G_CALLBACK(onMouseMoveCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "scroll-event",	    G_CALLBACK(onScrollCoverageView),                 NULL);

  /* Set required data in the coverageView. */
  CoverageViewProperties *cvProperties = coverageViewCreateProperties(coverageView, blxWindow, bc);

  return cvProperties;
}
