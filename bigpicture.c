/*
 *  bigpicture.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/bigpicture.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/exonview.h>
#include <SeqTools/utilities.h>
#include <math.h>


#define DEFAULT_PREVIEW_BOX_LINE_WIDTH  1
#define DEFAULT_GRID_NUM_HEADER_LINES   1	  /* the default number of lines of text in the grid header */
#define DEFAULT_GRID_HEADER_Y_PAD	0	  /* the default y padding around the grid header */
#define DEFAULT_GRID_CELL_WIDTH		100	  /* the default cell width of the grids */
#define DEFAULT_GRID_NUM_HOZ_CELLS	5	  /* the default number of cells to show horizontally in the grids */
#define DEFAULT_PERCENT_ID_PER_CELL	20	  /* the default %ID per vertical cell to show in the grids */
#define DEFAULT_GRID_PERCENT_ID_MAX	100	  /* default maximum %ID to show on the scale */
#define MIN_NUM_V_CELLS			1	  /* minimum number of vertical cells to show in the grid */


/* Local function declarations */
static GridHeaderProperties*	    gridHeaderGetProperties(GtkWidget *gridHeader);
static IntRange*		    bigPictureGetFullRange(GtkWidget *bigPicture);
static int			    bigPictureGetInitialZoom(GtkWidget *bigPicture);

static void                         drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc);

/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

void bigPictureRedrawAll(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  widgetClearCachedDrawable(properties->header);
  widgetClearCachedDrawable(properties->fwdStrandGrid);
  widgetClearCachedDrawable(properties->revStrandGrid);
  widgetClearCachedDrawable(properties->fwdExonView);
  widgetClearCachedDrawable(properties->revExonView);
  
  gtk_widget_queue_draw(bigPicture);
}


/* Utility to get the height and (approx) width of the given widget's font */
static void getFontCharSize(GtkWidget *widget, int *charWidth, int *charHeight)
{
  PangoContext *context = gtk_widget_get_pango_context(widget);
  PangoFontMetrics *metrics = pango_context_get_metrics (context,
							 widget->style->font_desc,
							 pango_context_get_language (context));
  *charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  *charWidth = pango_font_metrics_get_approximate_digit_width(metrics) / PANGO_SCALE;
  pango_font_metrics_unref(metrics);
}


/* Calculate how many pixels wide a base is */
gdouble pixelsPerBase(const gint displayWidth, const IntRange const *displayRange)
{
  gdouble displayLen = (gdouble)(displayRange->max - displayRange->min) + 1;
  return ((gdouble)displayWidth / displayLen);
}


/* Convert a base index to an x coord within the given rectangle. Returns the number of pixels 
 * from the left edge (including the start of the rectangle) to where the base lies. displayRev 
 * should be passed as true if the display is toggled (i.e. low values on the right and high 
 * values on the left). */
gint convertBaseIdxToGridPos(const gint baseIdx, 
			     const GdkRectangle const *rect, 
			     const IntRange const *displayRange)
{
  gint result = UNSET_INT;
  
  gdouble numBasesFromEdge = (gdouble)(baseIdx - displayRange->min); /* 0-based index from edge */
  if (numBasesFromEdge < 0)
    {
      numBasesFromEdge = 0;
    }
  
  gint pixelsFromEdge = (int)(numBasesFromEdge * pixelsPerBase(rect->width, displayRange));
  
  result = rect->x + pixelsFromEdge;
  
  if (result > rect->x + rect->width)
    {
      result = rect->x + rect->width;
    }
  
  return result;
}


/* Function to round the given value to the nearest "nice" value, from the given
 * list of values to round by. Returns the value it rounded to the nearest of. */
static int roundToValueFromList(const int inputVal, GSList *roundValues, int *roundedTo)
{
  /* Decide what amount to round to the nearest number of, out of a list of possible
   * "nice" values. Find the nearest value in this list to our result. The list 
   * shouldn't be long, so this doesn't worry about efficiency */
  GSList *listItem = roundValues;
  int roundTo = UNSET_INT;
  int smallestDiff = inputVal - roundTo;
  
  for ( ; listItem; listItem = listItem->next)
    {
      int val = GPOINTER_TO_INT(listItem->data);
      int thisDiff = inputVal - val;
      
      if (roundTo == UNSET_INT || (val > 0 && abs(thisDiff) < smallestDiff))
	{
	  roundTo = val;
	  smallestDiff = abs(thisDiff);
	}
    }
  
  /* Round the input to the nearest multiple of 'roundTo'. */
  int result = roundNearest((double)inputVal / (double)roundTo) * roundTo;
  
  if (roundedTo)
    {
      *roundedTo = roundTo;
    }
  
  return result;
}


/* This function calculates the cell size and number of cells for the big picture grids.
 * It should be called whenever the big picture is resized or its display range changes. */
static void calculateBigPictureCellSize(GtkWidget *bigPicture)
{
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  GtkWidget *header = bigPictureGetGridHeader(bigPicture);
  GridHeaderProperties *headerProperties = gridHeaderGetProperties(header);
  
  /* Calculate the number of bases per cell and round it to a "nice" value from our stored list */
  const int displayStart = convertDisplayIdxToDnaIdx(properties->displayRange.min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int displayEnd = convertDisplayIdxToDnaIdx(properties->displayRange.max, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int displayWidth = abs(displayEnd - displayStart); /* use abs because can be negative if display reversed */
  
  const int defaultBasesPerCell = ceil((double)(displayWidth) / DEFAULT_GRID_NUM_HOZ_CELLS);
  properties->basesPerCell = roundToValueFromList(defaultBasesPerCell, properties->roundValues, &properties->roundTo);

  /* Adjust the cell width (and hence number of cells) proportionally wrt the rounded number of bases per cell */
  const double defaultCellWidth = (double)(headerProperties->headerRect.width) / DEFAULT_GRID_NUM_HOZ_CELLS;
  const double actualCellWidth = (defaultCellWidth * properties->basesPerCell) / defaultBasesPerCell;
  properties->numHCells = ceil((double)headerProperties->headerRect.width / actualCellWidth);
}


static void drawVerticalGridLineHeaders(GtkWidget *header, 
					GtkWidget *bigPicture, 
                                        GdkDrawable *drawable,
					GdkGC *gc, 
					const GdkColor const *textColor, 
					const GdkColor const *lineColor)
{
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  GridHeaderProperties *headerProperties = gridHeaderGetProperties(header);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);

  const int direction = bc->displayRev ? -1 : 1; /* to subtract instead of add when display reversed */
  
  /* Get the first base index (in terms of the nucleotide coords) and round it to a nice round
   * number. We'll offset all of the gridlines by the distance between this and the real start coord. */
  const int realFirstBaseIdx = convertDisplayIdxToDnaIdx(bpProperties->displayRange.min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int firstBaseIdx = roundToValue(realFirstBaseIdx, bpProperties->roundTo);
  
  /* Calculate the top and bottom heights for the lines. */
  const gint bottomBorder = headerProperties->headerRect.y + headerProperties->headerRect.height;
  const gint topBorder = bottomBorder - headerProperties->markerHeight;

  const int minX = headerProperties->headerRect.x;
  const int maxX = headerProperties->headerRect.x + headerProperties->headerRect.width;

  /* Loop through each grid cell */
  gint hCell = 0;
  for ( ; hCell <= bpProperties->numHCells; ++hCell)
    {
      /* Draw the label, showing which base index is at this x coord */
      int numBasesFromLeft = bpProperties->basesPerCell * hCell;
      int baseIdx = firstBaseIdx + (numBasesFromLeft * direction);

      const int displayIdx = convertDnaIdxToDisplayIdx(baseIdx, bc->seqType, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
      const int x = convertBaseIdxToGridPos(displayIdx, &headerProperties->headerRect, &bpProperties->displayRange);
      
      if (x > minX && x < maxX)
	{
	  gdk_gc_set_foreground(gc, textColor);
	  gchar text[numDigitsInInt(baseIdx) + 1];
	  sprintf(text, "%d", baseIdx);

	  PangoLayout *layout = gtk_widget_create_pango_layout(header, text);
	  gdk_draw_layout(drawable, gc, x, 0, layout);
	  g_object_unref(layout);
	  
	  /* Draw a marker line (small cosmetic touch to better join the label to the grid below) */
	  gdk_gc_set_foreground(gc, lineColor);
	  gdk_draw_line (drawable, gc, x, topBorder, x, bottomBorder);
	}
    }
}


/* Refresh the header - clears and redraws its bitmap */
static void redrawBigPictureGridHeader(GtkWidget *header)
{
  /* Check that the header is shown on screen */
  if (GTK_LAYOUT(header)->bin_window)
    {
      /* Create a new bitmap to draw on to and set it in the widget. (This automatically
       * deletes the old one, if there is one.) */
      GdkDrawable *bitmap = gdk_pixmap_new(GTK_LAYOUT(header)->bin_window, header->allocation.width, header->allocation.height, -1);
      gdk_drawable_set_colormap(bitmap, gdk_colormap_get_system());
      widgetSetDrawable(header, bitmap);

      /* Clear the bitmap to the background color */
      GdkGC *gc = gdk_gc_new(bitmap);
      GtkStyle *style = gtk_widget_get_style(header);
      GdkColor *bgColor = &style->bg[GTK_STATE_NORMAL];
      gdk_gc_set_foreground(gc, bgColor);
      gdk_draw_rectangle(bitmap, gc, TRUE, 0, 0, header->allocation.width, header->allocation.height);

      /* Draw the header */
      drawBigPictureGridHeader(header, bitmap, gc);
    }
}



/* Draw the big picture header onto the given drawable */
static void drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  
  /* Set the drawing properties */
  gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);
  
  /* Draw the grid headers */
  drawVerticalGridLineHeaders(header, 
			      properties->bigPicture, 
                              drawable,
			      gc,
			      getGdkColor(bc, BLXCOL_GRID_TEXT, FALSE), 
			      getGdkColor(bc, BLXCOL_GRID_LINE, FALSE));
  
  g_object_unref(gc);
}


/* Recalculate the borders for the grid header. Should be called after a resize */
void calculateGridHeaderBorders(GtkWidget *header)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Calculate the size of the grid header (zero height if it does not have one) */
  properties->headerRect.x = bigPictureProperties->charWidth * bigPictureProperties->leftBorderChars;
  properties->headerRect.y = 0;
  properties->headerRect.width = header->allocation.width - properties->headerRect.x;
  properties->headerRect.height = properties->refButton->allocation.height + (properties->headerYPad * 2);
  
  properties->markerHeight = properties->headerRect.height - (bigPictureProperties->charHeight * properties->numHeaderLines);
  if (properties->markerHeight < 0)
    properties->markerHeight = 0;
  
  gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, properties->headerRect.height);
  gtk_widget_set_size_request(header, 0, properties->headerRect.height);
}


/* Add the given child to the big picture. This is separated out into its
 * own function in case we want to change the container type in the future. */
static void addChildToBigPicture(GtkWidget *container, GtkWidget *child, gboolean expand)
{
  if (GTK_IS_BOX(container))
    {
      gtk_box_pack_start(GTK_BOX(container), child, expand, FALSE, 0);
    }
}


/* This function removes the grids from the big picture and re-adds them in the
 * correct order according to the displayRev flag. It should be called every
 * time the strands are toggled. It assumes the two grids are both already in the 
 * bigPicture container, and that the properties have been set for all 3 widgets. */
void refreshGridOrder(GtkWidget *bigPicture)
{
  GtkWidget *fwdStrandGrid = bigPictureGetFwdGrid(bigPicture);
  GtkWidget *revStrandGrid = bigPictureGetRevGrid(bigPicture);
  GtkWidget *fwdExonView = bigPictureGetFwdExonView(bigPicture);
  GtkWidget *revExonView = bigPictureGetRevExonView(bigPicture);
  
  /* Increase the reference count to make sure the widgets aren't destroyed when we remove them. */
  g_object_ref(fwdStrandGrid);
  g_object_ref(revStrandGrid);
  g_object_ref(fwdExonView);
  g_object_ref(revExonView);
  
  /* Remove them */
  gtk_container_remove(GTK_CONTAINER(bigPicture), fwdStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bigPicture), revStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bigPicture), fwdExonView);
  gtk_container_remove(GTK_CONTAINER(bigPicture), revExonView);
  
  /* Add them back, with the forward-strand grid at the top and the reverse-strand grid
   * at the bottom, or vice versa if the strands are toggled. */
  if (bigPictureGetDisplayRev(bigPicture))
    {
      addChildToBigPicture(bigPicture, revStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, revExonView, FALSE);
      addChildToBigPicture(bigPicture, fwdExonView, FALSE);
      addChildToBigPicture(bigPicture, fwdStrandGrid, TRUE);
    }
  else
    {
      addChildToBigPicture(bigPicture, fwdStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, fwdExonView, FALSE);
      addChildToBigPicture(bigPicture, revExonView, FALSE);
      addChildToBigPicture(bigPicture, revStrandGrid, TRUE);
    }
  
  /* Decrease the ref count again */
  g_object_unref(fwdStrandGrid);
  g_object_unref(revStrandGrid);
  g_object_unref(fwdExonView);
  g_object_unref(revExonView);
  
  /* Must show all child widgets because some of them may not have been in this parent before.
   * (Just calling gtk_widget_show on the individual trees doesn't seem to work.)
   * However, we then need to re-hide any that may have been previously hidden by the user. */
  gtk_widget_show_all(bigPicture);
  gtk_container_foreach(GTK_CONTAINER(bigPicture), hideUserHiddenWidget, NULL);  
}


/* Adjust upper and/or lower end of the given range as necessary to
 * keep it within the given bounds, keeping the same width if possible. */
static void boundRange(IntRange *range, IntRange *bounds)
{
  int width = range->max - range->min;

  if (range->max > bounds->max)
  {
    range->max = bounds->max;
    range->min = range->max - width;
  }
  
  if (range->min < bounds->min)
  {
    range->min = bounds->min;
    range->max = range->min + width;
    
    if (range->max > bounds->max)
      {
	range->max = bounds->max;
      }
  }   
}


/* Set the display range for the big picture, based on the given width (i.e. number of
 * bases wide). Keeps the display centred on the same range that is shown in the detail view.
 * If recalcHighlightBox is true, the highlight box borders are recalculated. */
static void setBigPictureDisplayWidth(GtkWidget *bigPicture, int width, const gboolean recalcHighlightBox)
{
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);

  IntRange *displayRange = bigPictureGetDisplayRange(bigPicture);
  IntRange *detailViewRange = detailViewGetDisplayRange(detailView);
  IntRange *fullRange = bigPictureGetFullRange(bigPicture);
  
  int detailViewWidth = detailViewRange->max - detailViewRange->min;
  int maxWidth = fullRange->max - fullRange->min;
  gboolean done = FALSE;

  if (width < detailViewWidth)
    {
      /* Don't display less than the detail view range */
      displayRange->min = detailViewRange->min;
      displayRange->max = detailViewRange->max;
      width = displayRange->max - displayRange->min;
      done = TRUE;
    }

  if (width >= maxWidth)
    {
      /* Don't display more than the full range of the reference sequence */
      displayRange->min = fullRange->min;
      displayRange->max = fullRange->max;
      done = TRUE;
    }

  if (!done)
    {
      /* Keep the view centred on what's visible in the detail view */
      const int newcentre = getRangeCentre(detailViewRange);

      /* Try to display an equal amount either side of the centre */
      const int offset = roundNearest((double)width / 2.0);

      displayRange->min = newcentre - offset;
      displayRange->max = displayRange->min + width;
      boundRange(displayRange, fullRange);
    }
  
  /* Recalculate the grid cell size */
  calculateBigPictureCellSize(bigPicture);

  /* Recalculate the exon view height, because it may have changed with more/less
   * exons being scrolled into view */
  calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
  calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
  
  /* Since we're keeping the highlight box centred, it should stay in the same place
   * if we're just scrolling. We therefore only need to recalculate its position if
   * its size has changed. */
  if (recalcHighlightBox)
    {
      callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders);
    }

  bigPictureRedrawAll(bigPicture);
}


/* This function should be called every time the detail view display range has 
 * changed. It makes sure the big picture remains centred on the highlight box:
 * it scrolls the big picture display range if necessary to keep the highlight box
 * (i.e. the range that is displayed in the detail-view) centred. If recalcHighlightBox
 * is true, the highlight box has changed size, and its boundaries need to be recalculated. */
void refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean recalcHighlightBox)
{
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);
  IntRange *displayRange = bigPictureGetDisplayRange(bigPicture);
  IntRange *detailViewRange = detailViewGetDisplayRange(detailView);
  
  if (displayRange->min == UNSET_INT && displayRange->max == UNSET_INT) /* check both in case UNSET_INT is a real coord for one end! */
    {
      /* This is the first time we've refreshed the detail view range. Set the
       * initial big picture range width to be a ratio of the detail view range width. */
      const int width = (detailViewRange->max - detailViewRange->min) * bigPictureGetInitialZoom(bigPicture);
      setBigPictureDisplayWidth(bigPicture, width, recalcHighlightBox);
    }
  else
    {
      /* Call set-width (but just use the existing width) to force the required updates. */
      const int width = displayRange->max - displayRange->min;
      setBigPictureDisplayWidth(bigPicture, width, recalcHighlightBox);
    }

  /* Recalculate the exon view height, because it may have changed with more/less
   * exons being scrolled into view */
  calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
  calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the 
 * left edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
int getLeftCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int leftCoord = centreCoord - roundNearest((double)width / 2.0);
  
  if (outerRect)
    {
      if (leftCoord < outerRect->x) 
	leftCoord = outerRect->x;
      else
	{
	  int leftCoordMax = outerRect->x + outerRect->width - width;
	  
	  if (leftCoord > leftCoordMax) 
	    {
	      leftCoord = leftCoordMax;
	    }
	}
    }
  
  return leftCoord;
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the 
 * right edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
int getRightCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int rightCoord = centreCoord + roundNearest((double)width / 2.0);
  
  if (outerRect)
    {
      if (rightCoord > outerRect->x + outerRect->width)
	rightCoord = outerRect->x + outerRect->width;
      else
	{
	  int rightCoordMin = outerRect->x + width;
	  
	  if (rightCoord < rightCoordMin) 
	    {
	      rightCoord = rightCoordMin;
	    }
	}
    }
  
  return rightCoord;
}


/* Zoom the big picture in/out. Just halves or doubles the current display range. */
void zoomBigPicture(GtkWidget *bigPicture, const gboolean zoomIn)
{
  IntRange *displayRange = bigPictureGetDisplayRange(bigPicture);
  int newWidth = UNSET_INT;
  
  if (zoomIn)
    {
      newWidth = roundNearest(((double)(displayRange->max - displayRange->min)) / 2.0);
    }
  else
    {
      newWidth = (displayRange->max - displayRange->min) * 2;
    }

  setBigPictureDisplayWidth(bigPicture, newWidth, TRUE);
}


/* Zoom the big picture out to view the whole reference sequence */
void zoomWholeBigPicture(GtkWidget *bigPicture)
{
  IntRange *displayRange = bigPictureGetDisplayRange(bigPicture);
  IntRange *fullRange = bigPictureGetFullRange(bigPicture);
  
  /* Check we're not already showing the whole range */
  if (displayRange->min != fullRange->min || displayRange->max != fullRange->max)
    {
      setBigPictureDisplayWidth(bigPicture, fullRange->max - fullRange->min, TRUE);
    }
}


/* Calculate the number of cells to show vertically in the grids, and update the
 * grid's percent-ID range to be a full number of cells. */
void calculateNumVCells(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  const int idRangeLen = properties->percentIdRange.max - properties->percentIdRange.min;
  
  properties->numVCells = ceil((double)idRangeLen / (double)properties->idPerCell); 
  
  if (properties->numVCells < MIN_NUM_V_CELLS)
    {
      properties->numVCells = MIN_NUM_V_CELLS;
    }
  
  properties->percentIdRange.min = properties->percentIdRange.max - (properties->numVCells * properties->idPerCell);
}


/* update function to be called when the percent ID range or percent-ID-per-cell
 * values have been changed */
static void updateOnPercentIdChanged(GtkWidget *bigPicture)
{
  calculateNumVCells(bigPicture);
  
  callFuncOnAllBigPictureGrids(bigPicture, calculateGridBorders);
  callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders);
  
  bigPictureRedrawAll(bigPicture);
}

/***********************************************************
 *			    Events			   *
 ***********************************************************/

static void onZoomInBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomBigPicture(bigPicture, TRUE);
}


static void onZoomOutBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomBigPicture(bigPicture, FALSE);
}


static void onZoomWholeBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomWholeBigPicture(bigPicture);
}


static gboolean onExposeGridHeader(GtkWidget *header, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(header)->bin_window;
  
  if (window)
    {
      /* Just push the stored bitmap onto the screen */
      GdkDrawable *bitmap = widgetGetDrawable(header);

      if (!bitmap)
        {
          /* If the cache is empty, create it now. */
	  redrawBigPictureGridHeader(header);
	  bitmap = widgetGetDrawable(header);
        }
      
      if (bitmap)
        {
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
        }
    }
  
  return TRUE;
}



/***********************************************************
 *                       Properties                        *
 ***********************************************************/

BigPictureProperties* bigPictureGetProperties(GtkWidget *bigPicture)
{
  return bigPicture ? (BigPictureProperties*)(g_object_get_data(G_OBJECT(bigPicture), "BigPictureProperties")) : NULL;
}

BlxViewContext* bigPictureGetContext(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetContext(blxWindow);
}

static GridHeaderProperties* gridHeaderGetProperties(GtkWidget *gridHeader)
{
  return gridHeader ? (GridHeaderProperties*)(g_object_get_data(G_OBJECT(gridHeader), "GridHeaderProperties")) : NULL;
}


static void onDestroyBigPicture(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (properties)
    {
      if (properties->roundValues)
	{
	  g_slist_free(properties->roundValues);
	  properties->roundValues = NULL;
	}
      
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", NULL);
    }
}


static void bigPictureCreateProperties(GtkWidget *bigPicture, 
				       GtkWidget *blxWindow, 
				       GtkWidget *header, 
				       GtkWidget *fwdStrandGrid,
				       GtkWidget *revStrandGrid,
				       GtkWidget *fwdExonView,
				       GtkWidget *revExonView,
				       int previewBoxCentre,
				       const int initialZoom,
				       const int lowestId)
{
  if (bigPicture)
    { 
      BigPictureProperties *properties = g_malloc(sizeof *properties);
      
      properties->blxWindow = blxWindow;
      properties->header = header;
      properties->fwdStrandGrid = fwdStrandGrid;
      properties->revStrandGrid = revStrandGrid;
      properties->fwdExonView = fwdExonView;
      properties->revExonView = revExonView;
      
      properties->numHCells = UNSET_INT;
      properties->basesPerCell = UNSET_INT;
      properties->roundTo = 25;
      
      properties->numVCells = UNSET_INT;
      properties->idPerCell = DEFAULT_PERCENT_ID_PER_CELL;
      properties->percentIdRange.min = lowestId;
      properties->percentIdRange.max = DEFAULT_GRID_PERCENT_ID_MAX;

      properties->previewBoxCentre = previewBoxCentre;
      properties->leftBorderChars = numDigitsInInt(DEFAULT_GRID_PERCENT_ID_MAX) + 3; /* Extra fudge factor because char width is approx */
      properties->highlightBoxLineWidth = DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH;
      properties->previewBoxLineWidth = DEFAULT_PREVIEW_BOX_LINE_WIDTH;
      properties->initialZoom = initialZoom;
      
      /* These will be initialized when the detail view size is first allocated,
       * because the big picture range is a ratio of the detail view range. */
      properties->displayRange.min = UNSET_INT;
      properties->displayRange.max = UNSET_INT;
      
      /* Calculate the font size */
      getFontCharSize(bigPicture, &properties->charWidth, &properties->charHeight);
      
      /* Create the list of "nice round values" to round the grid header values to.
       * Create the list in reverse order (i.e. highest values first). */
      properties->roundValues = NULL;
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(25));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(50));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(100));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(125));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(150));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(250));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(500));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(1000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(2500));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(5000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(10000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(25000));
      
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", properties);
      g_signal_connect(G_OBJECT(bigPicture), "destroy", G_CALLBACK(onDestroyBigPicture), NULL); 
    }
}


static void onDestroyGridHeader(GtkWidget *bigPicture)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(bigPicture);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(bigPicture), "GridHeaderProperties", NULL);
    }
}

static void gridHeaderCreateProperties(GtkWidget *gridHeader, GtkWidget *bigPicture, GtkWidget *refButton)
{
  if (gridHeader)
    {
      GridHeaderProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->refButton = refButton;
      properties->numHeaderLines = DEFAULT_GRID_NUM_HEADER_LINES;
      properties->headerYPad = DEFAULT_GRID_HEADER_Y_PAD;
      
      g_object_set_data(G_OBJECT(gridHeader), "GridHeaderProperties", properties);
      g_signal_connect(G_OBJECT(gridHeader), "destroy", G_CALLBACK(onDestroyGridHeader), NULL);
    }
}




GtkWidget* bigPictureGetBlxWindow(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->blxWindow : NULL;
}

GtkWidget* bigPictureGetDetailView(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetDetailView(blxWindow);
}


GtkWidget* bigPictureGetFwdGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->fwdStrandGrid : NULL;
}

GtkWidget* bigPictureGetRevGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->revStrandGrid : NULL;
}

/* Get the active grid (forward strand grid by default, reverse strand grid if display toggled) */
GtkWidget* bigPictureGetActiveGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->revStrandGrid : properties->fwdStrandGrid;
}

/* Get the in-active grid (reverse strand grid by default, forward strand grid if display toggled) */
GtkWidget* bigPictureGetInactiveGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->fwdStrandGrid : properties->revStrandGrid;
}

GtkWidget* bigPictureGetFwdExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->fwdExonView : NULL;
}

GtkWidget* bigPictureGetRevExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->revExonView : NULL;
}

GtkWidget* bigPictureGetActiveExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->revExonView : properties->fwdExonView;
}

GtkWidget* bigPictureGetInactiveExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->fwdExonView : properties->revExonView;
}

gboolean bigPictureGetDisplayRev(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindow ? blxWindowGetDisplayRev(blxWindow) : FALSE;
}

IntRange* bigPictureGetDisplayRange(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->displayRange;
}

static int bigPictureGetInitialZoom(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->initialZoom;
}

static IntRange* bigPictureGetFullRange(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetFullRange(blxWindow);
}

GtkWidget* bigPictureGetGridHeader(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->header;
}

BlxSeqType bigPictureGetSeqType(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetSeqType(blxWindow);
}

int bigPictureGetNumFrames(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetNumFrames(blxWindow);
}

int bigPictureGetIdPerCell(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->idPerCell;
}

void bigPictureSetIdPerCell(GtkWidget *bigPicture, const int idPerCell)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  properties->idPerCell = idPerCell;
  updateOnPercentIdChanged(bigPicture);
}

int bigPictureGetNumVCells(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->numVCells;
}

IntRange* bigPictureGetPercentIdRange(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->percentIdRange;
}

void bigPictureSetMaxPercentId(GtkWidget *bigPicture, const int newValue)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  properties->percentIdRange.max = newValue;
  updateOnPercentIdChanged(bigPicture);
}

void bigPictureSetMinPercentId(GtkWidget *bigPicture, const int newValue)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  properties->percentIdRange.min = newValue;
  updateOnPercentIdChanged(bigPicture);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* Create a tool button for the big picture header */
static GtkWidget* createButton(GtkWidget *container, char *label, char *tooltip, GtkSignalFunc callback_func, gpointer data)
{
  GtkWidget *eventBox = gtk_event_box_new();
  gtk_box_pack_start(GTK_BOX(container), eventBox, FALSE, FALSE, 0);
  
  GtkWidget *button = gtk_button_new_with_label(label);
  gtk_container_add(GTK_CONTAINER(eventBox), button);
  
  /* gtk_widget_set_tooltip_text(button, tooltip); */ /* not in pre-GTK-2.12 */
  
  g_signal_connect(GTK_OBJECT(button), "clicked", callback_func, data);
  
  return button;
}


/* Create the header for the big picture section. This is a custom header
 * that contains some tool buttons and the header labels for the grid. We
 * squash them all into the same horizontal section to save space, and also
 * so that we can hide any of the grids but still keep the header visible. */
static GtkWidget *createBigPictureGridHeader(GtkWidget *bigPicture)
{
  GtkWidget *header = gtk_layout_new(NULL, NULL);
  
  /* Create the header buttons. Put them in an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_layout_put(GTK_LAYOUT(header), hbox, 0, 0);
  
  /* Store a ref to the first button so we can find out the height of the
   * buttons when we come to calculate the grid borders */
  GtkWidget *refButton = createButton(hbox, "Zoom in", "Zoom in (Ctrl-=)", (GtkSignalFunc)onZoomInBigPicture, bigPicture);
  createButton(hbox, "Zoom out", "Zoom out (Ctrl--)", (GtkSignalFunc)onZoomOutBigPicture, bigPicture);
  createButton(hbox, "Whole", "Zoom out to whole width (Shift-Ctrl--_", (GtkSignalFunc)onZoomWholeBigPicture, bigPicture);
  
  /* Create the header properties */
  gridHeaderCreateProperties(header, bigPicture, refButton);
  
  /* Conect signals */
  g_signal_connect(G_OBJECT(header), "expose-event", G_CALLBACK(onExposeGridHeader), NULL);
  
  return header;
}


GtkWidget* createBigPicture(GtkWidget *blxWindow, 
			    GtkWidget *container,
			    GtkWidget **fwdStrandGrid, 
			    GtkWidget **revStrandGrid,
			    const int initialZoom,
			    const int lowestId)
{
  /* Create the main big picture widget, which will contain all of the 
   * individual big-picture grids, plus a header. */
  GtkWidget *bigPicture = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(container), bigPicture, FALSE, TRUE, 0);
  
  /* Our big picture needs to have a header plus 3 panes:
   * 1. the top grid (fwd strand);
   * 2. the exon view; and
   * 3. bottom grid (reverse strand) */
  
  GtkWidget *header = createBigPictureGridHeader(bigPicture);
  addChildToBigPicture(bigPicture, header, FALSE);

  *fwdStrandGrid = createBigPictureGrid(bigPicture, BLXSTRAND_FORWARD);
  *revStrandGrid = createBigPictureGrid(bigPicture, BLXSTRAND_REVERSE);

  GtkWidget *fwdExonView = createExonView(bigPicture, BLXSTRAND_FORWARD);
  GtkWidget *revExonView = createExonView(bigPicture, BLXSTRAND_REVERSE);
  
  /* By default, make the forward strand the top grid */
  addChildToBigPicture(bigPicture, *fwdStrandGrid, FALSE);
  addChildToBigPicture(bigPicture, fwdExonView, FALSE);
  addChildToBigPicture(bigPicture, revExonView, FALSE);
  addChildToBigPicture(bigPicture, *revStrandGrid, TRUE);

  /* Set the big picture properties */
  bigPictureCreateProperties(bigPicture, 
			     blxWindow,
			     header, 
			     *fwdStrandGrid,
			     *revStrandGrid,
			     fwdExonView,
			     revExonView,
			     UNSET_INT,
			     initialZoom,
			     lowestId);
  
  return bigPicture;
}



