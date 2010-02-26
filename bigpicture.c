/*
 *  bigpicture.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/bigpicture.h>
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/exonview.h>
#include <SeqTools/utilities.h>
#include <math.h>


#define DEFAULT_PREVIEW_BOX_LINE_WIDTH  1
#define DEFAULT_GRID_NUM_HEADER_LINES   1
#define DEFAULT_GRID_HEADER_Y_PAD	0
#define DEFAULT_GRID_CELL_WIDTH		100
#define DEFAULT_GRID_NUM_HOZ_CELLS	5


/* Local function declarations */
static GridHeaderProperties*	    gridHeaderGetProperties(GtkWidget *gridHeader);
static int			    bigPictureGetNumHCells(GtkWidget *bigPicture);
static int			    bigPictureGetCellWidth(GtkWidget *bigPicture);
static IntRange*		    bigPictureGetFullRange(GtkWidget *bigPicture);

static void                         drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc);

/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

/* Utility to get the height and (approx) width of the given widget's font */
static void getFontCharSize(GtkWidget *widget, int *charWidth, int *charHeight)
{
  PangoContext *context = gtk_widget_get_pango_context(widget);
  PangoFontMetrics *metrics = pango_context_get_metrics (context,
							 widget->style->font_desc,
							 pango_context_get_language (context));
  *charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  *charWidth = pango_font_metrics_get_approximate_char_width(metrics) / PANGO_SCALE;
  pango_font_metrics_unref(metrics);
}


/* Calculate how many pixels wide a base is */
gdouble pixelsPerBase(const gint displayWidth, const IntRange const *displayRange)
{
  gdouble displayLen = (gdouble)(displayRange->max - displayRange->min) + 1;
  return ((gdouble)displayWidth / displayLen);
}


/* Convert a base index to an x coord within the given rectangle. Returns the number of pixels 
 * from the left edge (including the start of the rectangle) to where the base lies. rightToLeft 
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


static void drawVerticalGridLineHeaders(GtkWidget *header, 
					GtkWidget *bigPicture, 
                                        GdkDrawable *drawable,
					GdkGC *gc, 
					const GdkColor const *textColour, 
					const GdkColor const *lineColour)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  const gint bottomBorder = properties->headerRect.y + properties->headerRect.height;
  const gint topBorder = bottomBorder - properties->markerHeight;

  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(mainWindow);
  const IntRange const *refSeqRange = mainWindowGetRefSeqRange(mainWindow);
  const BlxSeqType seqType = mainWindowGetSeqType(mainWindow);
  const int numFrames = bigPictureGetNumReadingFrames(bigPicture);

  const IntRange const *displayRange = bigPictureGetDisplayRange(bigPicture);
  int numCells = bigPictureGetNumHCells(bigPicture);
  int cellWidth = bigPictureGetCellWidth(bigPicture);
  
  gint basesPerCell = ceil((double)(displayRange->max - displayRange->min) / numCells);
  
  /* Omit the first and last lines */
  gint hCell = 1;
  for ( ; hCell < numCells; ++hCell)
    {
      gint x = properties->headerRect.x + (gint)((gdouble)hCell * cellWidth);
      
      /* Draw the label, showing which base index is at this x coord */
      int numBasesFromLeft = basesPerCell * hCell;
      int baseIdx = rightToLeft ? displayRange->max - numBasesFromLeft : displayRange->min + numBasesFromLeft;
      
      /* Convert the display coord to a ref seq coord */
      baseIdx = convertDisplayIdxToDnaIdx(baseIdx, seqType, 1, 1, numFrames, rightToLeft, refSeqRange);
      
      gdk_gc_set_foreground(gc, textColour);
      gchar text[numDigitsInInt(baseIdx) + 1];
      sprintf(text, "%d", baseIdx);

      PangoLayout *layout = gtk_widget_create_pango_layout(header, text);
      gdk_draw_layout(drawable, gc, x, 0, layout);
      g_object_unref(layout);
      
      /* Draw the marker line */
      gdk_gc_set_foreground(gc, lineColour);
      gdk_draw_line (drawable, gc, x, topBorder, x, bottomBorder);
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

      /* Clear the bitmap to the background colour */
      GdkGC *gc = gdk_gc_new(bitmap);
      GtkStyle *style = gtk_widget_get_style(header);
      GdkColor *bgColour = &style->bg[GTK_STATE_NORMAL];
      gdk_gc_set_foreground(gc, bgColour);
      gdk_draw_rectangle(bitmap, gc, TRUE, 0, 0, header->allocation.width, header->allocation.height);

      /* Draw the header */
      drawBigPictureGridHeader(header, bitmap, gc);
    }
}



/* Draw the big picture header onto the given drawable */
static void drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Set the drawing properties */
  gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);
  
  /* Draw the grid headers */
  drawVerticalGridLineHeaders(header, 
			      properties->bigPicture, 
                              drawable,
			      gc,
			      &bigPictureProperties->gridTextColour, 
			      &bigPictureProperties->gridLineColour);
  
  g_object_unref(gc);
}


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
 * correct order according to the strandsToggled flag. It should be called every
 * time the strands are toggled. It assumes the two grids are both already in the 
 * bigPicture container, and that the properties have been set for all 3 widgets. */
void refreshGridOrder(GtkWidget *bigPicture)
{
  GtkWidget *fwdStrandGrid = bigPictureGetFwdGrid(bigPicture);
  GtkWidget *revStrandGrid = bigPictureGetRevGrid(bigPicture);
  GtkWidget *exonView = bigPictureGetExonView(bigPicture);
  
  /* Increase the reference count to make sure the widgets aren't destroyed when we remove them. */
  g_object_ref(fwdStrandGrid);
  g_object_ref(exonView);
  g_object_ref(revStrandGrid);
  
  /* Remove them */
  gtk_container_remove(GTK_CONTAINER(bigPicture), fwdStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bigPicture), exonView);
  gtk_container_remove(GTK_CONTAINER(bigPicture), revStrandGrid);
  
  /* Add them back, with the forward-strand grid at the top and the reverse-strand grid
   * at the bottom, or vice versa if the strands are toggled. */
  if (bigPictureGetStrandsToggled(bigPicture))
    {
      addChildToBigPicture(bigPicture, revStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, exonView, FALSE);
      addChildToBigPicture(bigPicture, fwdStrandGrid, TRUE);
    }
  else
    {
      addChildToBigPicture(bigPicture, fwdStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, exonView, FALSE);
      addChildToBigPicture(bigPicture, revStrandGrid, TRUE);
    }
  
  /* Decrease the ref count again */
  g_object_unref(fwdStrandGrid);
  g_object_unref(exonView);
  g_object_unref(revStrandGrid);
  
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
      int newcentre = getRangeCentre(detailViewRange);

      /* Try to display an equal amount either side of the centre */
      int offset = roundNearest((double)width / 2.0);

      displayRange->min = newcentre - offset;
      displayRange->max = newcentre + offset;
      boundRange(displayRange, fullRange);
    }
  
  /* Since we're keeping the highlight box centred, it should stay in the same place
   * if we're just scrolling. We therefore only need to recalculate its position if
   * its size has changed. */
  if (recalcHighlightBox)
    {
      callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders);
    }

  /* Recreate the grids and grid header */
  callFuncOnAllBigPictureGrids(bigPicture, widgetClearCachedDrawable);
  widgetClearCachedDrawable(bigPictureGetGridHeader(bigPicture));
  gtk_widget_queue_draw(bigPicture);
}


/* This function makes sure the big picture remains centred on the highlight box:
 * it scrolls the big picture display range if necessary to keep the highlight box
 * (i.e. the range that is displayed in the detail-view) centred. If recalcHighlightBox
 * is true, the highlight box has changed size, and its boundaries need to be recalculated. */
void refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean recalcHighlightBox)
{
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);
  IntRange *displayRange = bigPictureGetDisplayRange(bigPicture);
  IntRange *detailViewRange = detailViewGetDisplayRange(detailView);
  
  if (getRangeCentre(displayRange) != getRangeCentre(detailViewRange))
    {
      int width = displayRange->max - displayRange->min;
      setBigPictureDisplayWidth(bigPicture, width, recalcHighlightBox);
    }
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

static GridHeaderProperties* gridHeaderGetProperties(GtkWidget *gridHeader)
{
  return gridHeader ? (GridHeaderProperties*)(g_object_get_data(G_OBJECT(gridHeader), "GridHeaderProperties")) : NULL;
}


static void onDestroyBigPicture(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", NULL);
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


GtkWidget* bigPictureGetMainWindow(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->mainWindow : NULL;
}

GtkWidget* bigPictureGetDetailView(GtkWidget *bigPicture)
{
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetDetailView(mainWindow);
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
  return bigPictureGetStrandsToggled(bigPicture) ? properties->revStrandGrid : properties->fwdStrandGrid;
}

/* Get the in-active grid (reverse strand grid by default, forward strand grid if display toggled) */
GtkWidget* bigPictureGetInactiveGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetStrandsToggled(bigPicture) ? properties->fwdStrandGrid : properties->revStrandGrid;
}

GtkWidget* bigPictureGetExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->exonView : NULL;
}

gboolean bigPictureGetStrandsToggled(GtkWidget *bigPicture)
{
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindow ? mainWindowGetStrandsToggled(mainWindow) : FALSE;
}

IntRange* bigPictureGetDisplayRange(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->displayRange;
}

static IntRange* bigPictureGetFullRange(GtkWidget *bigPicture)
{
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetFullRange(mainWindow);
}

static int bigPictureGetNumHCells(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->numHCells;
}

int bigPictureGetCellWidth(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->cellWidth;
}

GtkWidget* bigPictureGetGridHeader(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->header;
}

BlxSeqType bigPictureGetSeqType(GtkWidget *bigPicture)
{
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetSeqType(mainWindow);
}

int bigPictureGetNumReadingFrames(GtkWidget *bigPicture)
{
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetNumReadingFrames(mainWindow);
}

static void bigPictureCreateProperties(GtkWidget *bigPicture, 
				       GtkWidget *mainWindow, 
				       GtkWidget *header, 
				       GtkWidget *fwdStrandGrid,
				       GtkWidget *revStrandGrid,
				       GtkWidget *exonView,
				       const IntRange const *initDisplayRange, 
				       int cellWidth, 
				       int numHCells, 
				       int previewBoxCentre)
{
  if (bigPicture)
    { 
      BigPictureProperties *properties = g_malloc(sizeof *properties);
      
      properties->mainWindow = mainWindow;
      properties->header = header;
      properties->fwdStrandGrid = fwdStrandGrid;
      properties->revStrandGrid = revStrandGrid;
      properties->exonView = exonView;
      properties->displayRange.min = initDisplayRange->min;
      properties->displayRange.max = initDisplayRange->max;
      properties->cellWidth = cellWidth;
      properties->numHCells = numHCells;
      properties->previewBoxCentre = previewBoxCentre;
      properties->leftBorderChars = numDigitsInInt(DEFAULT_GRID_PERCENT_ID_MAX) + 3; /* Extra fudge factor because char width is approx */
      properties->highlightBoxLineWidth = DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH;
      properties->previewBoxLineWidth = DEFAULT_PREVIEW_BOX_LINE_WIDTH;
      
      properties->gridLineColour = getGdkColor(GDK_YELLOW);
      properties->gridTextColour = getGdkColor(GDK_BLACK);
      properties->highlightBoxColour = getGdkColor(GDK_BLUE);
      properties->previewBoxColour = getGdkColor(GDK_DARK_GREY);
      properties->mspLineHighlightColour = getGdkColor(GDK_CYAN);
      properties->mspLineColour = getGdkColor(GDK_BLACK);
      
      /* Calculate the font size */
      getFontCharSize(bigPicture, &properties->charWidth, &properties->charHeight);
      
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", properties);
      g_signal_connect(G_OBJECT(bigPicture), "destroy", G_CALLBACK(onDestroyBigPicture), NULL); 
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
  
  gtk_widget_set_tooltip_text(button, tooltip);
  
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
  GtkWidget *refButton = createButton(hbox, "+", "Zoom in\tCtrl =", (GtkSignalFunc)onZoomInBigPicture, bigPicture);
  createButton(hbox, "-", "Zoom out\tCtrl -", (GtkSignalFunc)onZoomOutBigPicture, bigPicture);
  createButton(hbox, "Whole", "Zoom out to whole width\tCtrl-W", (GtkSignalFunc)onZoomWholeBigPicture, bigPicture);
  
  /* Create the header properties */
  gridHeaderCreateProperties(header, bigPicture, refButton);
  
  /* Conect signals */
  g_signal_connect(G_OBJECT(header), "expose-event", G_CALLBACK(onExposeGridHeader), NULL);
  
  return header;
}


GtkWidget* createBigPicture(GtkWidget *mainWindow, 
			    GtkWidget *container,
			    GtkWidget **fwdStrandGrid, 
			    GtkWidget **revStrandGrid, 
			    const IntRange const *initDisplayRange)
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

  *fwdStrandGrid = createBigPictureGrid(bigPicture, FORWARD_STRAND);
  *revStrandGrid = createBigPictureGrid(bigPicture, REVERSE_STRAND);

  GtkWidget *exonView = createExonView(bigPicture);
  
  /* By default, make the forward strand the top grid */
  addChildToBigPicture(bigPicture, *fwdStrandGrid, FALSE);
  addChildToBigPicture(bigPicture, exonView, FALSE);
  addChildToBigPicture(bigPicture, *revStrandGrid, TRUE);

  /* Set the big picture properties */
  bigPictureCreateProperties(bigPicture, 
			     mainWindow,
			     header, 
			     *fwdStrandGrid,
			     *revStrandGrid,
			     exonView,
			     initDisplayRange, 
			     DEFAULT_GRID_CELL_WIDTH, 
			     DEFAULT_GRID_NUM_HOZ_CELLS, 
			     UNSET_INT);
  
  
  return bigPicture;
}



