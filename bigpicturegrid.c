/*
 *  bigpicturegrid.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <SeqTools/bigpicturegrid.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/bigpicturemspline.h>
#include <SeqTools/blxviewMainWindow.h>
#include <math.h>
#include <string.h>


#define BIG_PICTURE_MSP_LINE_NAME	"BigPictureMspLine"
#define DEFAULT_GRID_NUM_VERT_CELLS	5
#define DEFAULT_GRID_PERCENT_ID_MIN	0
#define DEFAULT_MSP_LINE_HEIGHT		5
#define DEFAULT_GRID_Y_PADDING		5
#define DEFAULT_GRID_CELL_Y_PADDING	-2
#define DEFAULT_HIGHLIGHT_BOX_Y_PAD	2

/* Local function declarations */
static GtkAdjustment*	    gridGetAdjustment(GtkWidget *grid);
static GtkAdjustment*	    gridGetDetailView(GtkWidget *grid);
static IntRange*	    gridGetDisplayRange(GtkWidget *grid);


/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

/* Calls the given function (passed as the data pointer) on the given widget 
 * if it is a grid in the big picture view, or, if it is a container, 
 * calls the function on all children/grandchildren/etc that are grids */
void callFuncOnAllBigPictureGrids(GtkWidget *widget, gpointer data)
{
  const gchar *name = gtk_widget_get_name(widget);
  if (strcmp(name, BIG_PICTURE_GRID_NAME) == 0)
    {
      GtkCallback func = (GtkCallback)data;
      func(widget, NULL);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllBigPictureGrids, data);
    }
}


/* Get the height of the grid cells */
static gint gridGetCellHeight(GtkWidget *grid)
{
  /* Base the cell height on the font height */
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  return bigPictureProperties->charHeight + (2 * properties->cellYPadding);
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the 
 * left edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
static int getLeftCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int leftCoord = centreCoord - (width / 2);
  
  if (outerRect)
    {
      if (leftCoord < outerRect->x) 
	leftCoord = outerRect->x;
      else
	{
	  int leftCoordMax = outerRect->x + outerRect->width - width;
	  if (leftCoord > leftCoordMax) 
	    leftCoord = leftCoordMax;
	}
    }
  
  return leftCoord;
}


/* Get the x coord of the centre of the preview box within the given grid. 
 * Returns UNSET_INT if the preview box is not displayed. */
static int gridGetPreviewBoxCentre(GtkWidget *grid)
{
  int result = UNSET_INT;
  
  /* This data lives in the detail view because it is common to all grids */
  GridProperties *gridProperties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(gridProperties->bigPicture);
  if (bigPictureProperties)
    result = bigPictureProperties->previewBoxCentre;
  
  return result;
}


/* Calculate how many pixels wide a base is */
static gdouble pixelsPerBase(const gint gridWidth, const IntRange const *displayRange)
{
  gdouble displayLen = (gdouble)(displayRange->max - displayRange->min);
  return ((gdouble)gridWidth / displayLen);
}


/* Convert a base index to a grid position x coord. Returns the number of pixels from the
 * left edge of the widget to where the base lies. */
gint convertBaseIdxToGridPos(const gint baseIdx, const GdkRectangle const *gridRect, const IntRange const *displayRange)
{
  return gridRect->x + (gint)((gdouble)baseIdx * pixelsPerBase(gridRect->width, displayRange));
}


/* Convert an x coord on the given grid to a base index */
static gint convertGridPosToBaseIdx(const gint gridPos, const GdkRectangle const *gridRect,  const IntRange const *displayRange)
{
  return round((gdouble)(gridPos - gridRect->x) / pixelsPerBase(gridRect->width, displayRange));
}


///* Convert a y coord on the given grid to an ID% value */
//static gint convertGridPosToValue(GtkWidget *grid, const gint gridPos)
//{
//  GridProperties *properties = gridGetProperties(grid);
//  IntRange *valRange = &properties->percentIdRange;
//  return ((gridPos - properties->gridRect.y) / (properties->gridRect.height * (valRange->max - valRange->min)));
//}


/* Convert an ID% value to the y coord on the given grid */
gint convertValueToGridPos(GtkWidget *grid, const gint value)
{
  /* The top line of the grid is drawn one cell height down from the top of the grid border */
  GridProperties *properties = gridGetProperties(grid);
  
  IntRange *valRange = &properties->percentIdRange;
  gdouble percent = (gdouble)(valRange->max - value) / (gdouble)(valRange->max - valRange->min);
  
  /* Make sure we do the multiplication before truncating to int */
  gint result = properties->gridRect.y + (gint)((gdouble)properties->gridRect.height * percent);
  return result;
}


/* Set the x coord of the centre of the preview box within the given grid.
 * Setting it to UNSET_INT means the preview box will not be displayed. */
static void gridSetPreviewBoxCentre(GtkWidget *grid, int previewBoxCentre)
{
  /* This data lives in the detail view because it is common to all grids */
  GridProperties *gridProperties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(gridProperties->bigPicture);
  if (bigPictureProperties)
    bigPictureProperties->previewBoxCentre = previewBoxCentre;
}


/* Draw the highlight box. */
static void drawHighlightBox(GtkWidget *widget, GdkGC *gc, const GdkRectangle const *rect,
			     const gint lineWidth, GdkColor *lineColor)
{
  gdk_gc_set_rgb_fg_color(gc, lineColor);
  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  gdk_draw_rectangle(GTK_LAYOUT(widget)->bin_window, gc, FALSE, 
		     rect->x, rect->y, rect->width, rect->height);
}


/* Draw the preview box, if applicable */
static void drawPreviewBox(GtkWidget *grid, GdkGC *gc, 
			   const GdkRectangle const *highlightRect,
			   const gint lineWidth, GdkColor *lineColor)
{
  /* If a preview position is set, draw a preview of the highlight rect 
   * centred at that position */
  int previewBoxCentre = gridGetPreviewBoxCentre(grid);
  if (previewBoxCentre != UNSET_INT)
    {
      GridProperties *properties = gridGetProperties(grid);
      
      /* Find the x coord for the left edge of the preview box. */
      int x = getLeftCoordFromCentre(previewBoxCentre, highlightRect->width, &properties->gridRect);
      
      /* The other dimensions of the preview box are the same as the current highlight box. */
      GdkRectangle previewRect = {x, highlightRect->y, highlightRect->width, highlightRect->height};
      drawHighlightBox(grid, gc, &previewRect, lineWidth, lineColor);
    }
}


/* Draw the vertical gridlines in the big picture view */
static void drawVerticalGridLines(GtkWidget *grid, GdkGC *gc, 
				  const gint numCells, const gdouble cellWidth, 
				  const GdkColor const *lineColor)
{
  GridProperties *properties = gridGetProperties(grid);
  
  const gint topBorder = properties->highlightRect.y - properties->gridYPadding;
  const gint bottomBorder = properties->gridRect.y + properties->gridRect.height;
  
  /* Omit the first and last lines */
  gint hCell = 1;
  for ( ; hCell < numCells; ++hCell)
    {
      gint x = properties->gridRect.x + (gint)((gdouble)hCell * cellWidth);
      gdk_gc_set_rgb_fg_color(gc, lineColor);
      gdk_draw_line (GTK_LAYOUT(grid)->bin_window, gc, x, topBorder, x, bottomBorder);
    }
}


/* Draw the horizontal grid lines for the big picture view */
static void drawHorizontalGridLines(GtkWidget *grid, GdkGC *gc, 
				    const gint numCells, 
				    const gint rangePerCell, const gint maxVal, 
				    const GdkColor const *textColor,
				    const GdkColor const *lineColor)
{
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  const gint rightBorder = properties->gridRect.x + properties->gridRect.width;
  
  gint vCell = 0;
  for ( ; vCell <= numCells; ++vCell)
    {
      gint y = properties->gridRect.y + (gint)((gdouble)vCell * gridGetCellHeight(grid));
      
      /* Label this gridline with the %ID */
      gint percent = maxVal - (rangePerCell * vCell);
      gdk_gc_set_rgb_fg_color(gc, textColor);
      
      char text[bigPictureProperties->leftBorderChars + 1];
      sprintf(text, "%d%%", percent);
      
      PangoLayout *layout = gtk_widget_create_pango_layout(grid, text);
      gdk_draw_layout(GTK_LAYOUT(grid)->bin_window, gc, 0, y - gridGetCellHeight(grid)/2, layout);
      g_object_unref(layout);
      
      /* Draw the gridline */
      gdk_gc_set_rgb_fg_color(gc, lineColor);
      gdk_draw_line (GTK_LAYOUT(grid)->bin_window, gc, properties->gridRect.x, y, rightBorder, y);
    }
}


/* Calculates the size and position of an MSP line in the given grid. */
void calculateMspLineDimensions(GtkWidget *grid, const MSP const *msp, 
				int *x, int *y, int *width, int *height)
{
  GridProperties *gridProperties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(gridProperties->bigPicture);

  /* Find the coordinates of the start and end base in this match sequence */
  *x = convertBaseIdxToGridPos(msp->qstart, &gridProperties->gridRect, &bigPictureProperties->displayRange);
  int xEnd = convertBaseIdxToGridPos(msp->qend + 1, &gridProperties->gridRect, &bigPictureProperties->displayRange);
  *width = xEnd - *x;
  
  /* Find where in the y axis we should draw the line, based on the %ID value */
  *y = convertValueToGridPos(grid, msp->id);
  *height = gridProperties->mspLineHeight;
}


/* Draw a line on the grid to represent the MSP in the given tree path */
static gboolean drawMspLine(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *grid = GTK_WIDGET(data);
  GridProperties *gridProperties = gridGetProperties(grid);
  GtkWidget *tree = gridProperties->tree;

  /* Set some colours, and our drawing context */
  GdkGC *gc = gdk_gc_new(grid->window);
  gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);

  GdkColor normalColor, highlightColor;
  setGdkColorBlack(&normalColor);
  setGdkColorCyan(&highlightColor);
  
  /* Set the colour depending on the selection status */
  if (treePathIsSelected(GTK_TREE_VIEW(tree), path, model))
    gdk_gc_set_rgb_fg_color(gc, &highlightColor);
  else
    gdk_gc_set_rgb_fg_color(gc, &normalColor);
 
  /* Calculate where it should go */
  const MSP* msp = treeGetMsp(model, iter);
  int x, y, width, height;
  calculateMspLineDimensions(grid, msp, &x, &y, &width, &height);
  
  /* Draw it */
  gdk_draw_rectangle(GTK_LAYOUT(grid)->bin_window, gc, TRUE, x, y, width, height);
  
  /* Clean up */
  g_object_unref(gc);
  
  return FALSE;
}


/* Draw the MSP lines for the given grid */
static void drawMspLines(GtkWidget *grid)
{
  /* Get the corresponding tree */
  GridProperties *properties = gridGetProperties(grid);
  
  /* Loop through each row in the base (i.e. unfiltered) model for this tree. */
  if (properties && properties->tree)
    {
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(properties->tree));
      gtk_tree_model_foreach(model, drawMspLine, grid);
    }
}


/* Draw a big picture grid */
static void drawBigPictureGrid(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Set the drawing properties */
  GdkGC *gc = gdk_gc_new(grid->window);
  
  /* Calculate some factors for scaling */
  const gint percentPerCell = 
    (properties->percentIdRange.max - properties->percentIdRange.min) / 
    (properties->numVCells > 0 ? properties->numVCells : 1);
  
  /* Draw the grid */
  drawVerticalGridLines(grid, gc, bigPictureProperties->numHCells, bigPictureProperties->cellWidth,
			&bigPictureProperties->gridLineColor);
  
  drawHorizontalGridLines(grid, gc, properties->numVCells, percentPerCell, 
			  properties->percentIdRange.max, &bigPictureProperties->gridTextColor, 
			  &bigPictureProperties->gridLineColor);
  
  /* Draw lines corresponding to the MSPs */
  drawMspLines(grid);
  
  /* Draw the highlight box */
  drawHighlightBox(grid, gc, &properties->highlightRect, bigPictureProperties->highlightBoxLineWidth,
		   &bigPictureProperties->highlightColor);
  
  drawPreviewBox(grid, gc, &properties->highlightRect, bigPictureProperties->previewBoxLineWidth, 
		 &bigPictureProperties->previewColor);
  
  g_object_unref(gc);
}


void calculateHighlightBoxBorders(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  
  /* Calculate how many pixels from the left edge of the widget to the first base in the range. Truncating
   * the double to an int after the multiplication means we can be up to 1 pixel out, but this should be fine. */
  GtkAdjustment *adjustment = gridGetAdjustment(grid);
  if (adjustment)
    {
      IntRange *displayRange = gridGetDisplayRange(grid);

      properties->highlightRect.x = convertBaseIdxToGridPos(adjustment->value, &properties->gridRect, displayRange);
      properties->highlightRect.y = properties->gridRect.y - properties->highlightBoxYPad;
      
      properties->highlightRect.width = (gint)((gdouble)adjustment->page_size * pixelsPerBase(properties->gridRect.width, displayRange));
      properties->highlightRect.height = properties->gridRect.height + properties->mspLineHeight + (2 * properties->highlightBoxYPad);
    }
}


/* Calculate the borders for the big picture view */
void calculateGridBorders(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Get some info about the size of the layout and the font used */
  guint layoutWidth, layoutHeight;
  gtk_layout_get_size(GTK_LAYOUT(grid), &layoutWidth, &layoutHeight);
  
  properties->displayRect.x = 0;
  properties->displayRect.y = 0;
  properties->displayRect.width = grid->allocation.width;
  
  /* Get the boundaries of the grid */
  properties->gridRect.x = bigPictureProperties->charWidth * bigPictureProperties->leftBorderChars;
  properties->gridRect.y = properties->highlightBoxYPad + properties->gridYPadding;
  properties->gridRect.width = properties->displayRect.width - properties->gridRect.x;
  properties->gridRect.height = gridGetCellHeight(grid) * properties->numVCells;
  
  /* Get the boundaries of the highlight box */
  calculateHighlightBoxBorders(grid);
  
  /* Get the total display height required. Set the layout size to fit. */
  properties->displayRect.height = properties->highlightRect.height + (2 * properties->gridYPadding);
  gtk_layout_set_size(GTK_LAYOUT(grid), properties->displayRect.width, properties->displayRect.height);
  
  /* Set the size request to our desired height. We want a fixed heigh but don't set the
   * width, because we want the user to be able to resize that. */
  gtk_widget_set_size_request(grid, 0, properties->displayRect.height);
  
  /* Update the cell width dynamically with the grid width */
  bigPictureProperties->cellWidth = properties->gridRect.width / bigPictureProperties->numHCells;
}


/***********************************************************
 *			    Events			   *
 ***********************************************************/


static gboolean onExposeGrid(GtkWidget *grid, GdkEventExpose *event, gpointer data)
{
  drawBigPictureGrid(grid);
  return TRUE;
}


/* Re-calculate size info about the given grid following a resize */
static void onSizeAllocateBigPictureGrid(GtkWidget *grid, GtkAllocation *allocation, gpointer data)
{
  /* Recalculate the borders for the grids and the header */
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  calculateGridHeaderBorders(bigPictureProperties->header);
  callFuncOnAllBigPictureGrids(properties->bigPicture, calculateGridBorders);
  
  /* Update the child msp lines */
//  callFuncOnAllMspLines(grid, configureMspLine);
}


static gboolean onButtonPressGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  if (event->button == 2) /* middle button */
    {
      /* Draw a preview box at the clicked position */
      gridSetPreviewBoxCentre(grid, event->x);
      
      /* Force immediate update so that it happens before the button-release event */
      gtk_widget_queue_draw(grid);
      gdk_window_process_updates(grid->window, TRUE);
    }
  
  return FALSE;
}


static gboolean onButtonReleaseGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  if (event->button == 2) /* middle button */
    {
      GridProperties *properties = gridGetProperties(grid);
      BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
      
      /* Find the base index where the preview box starts */
      int x = getLeftCoordFromCentre(event->x, properties->highlightRect.width, &properties->gridRect);
      int baseIdx = convertGridPosToBaseIdx(x, &properties->gridRect, &bigPictureProperties->displayRange);
      
      /* Clear the preview box */
      gridSetPreviewBoxCentre(grid, UNSET_INT);
      
      /* Update the detail view's scroll pos to start at the start base */
      GtkAdjustment *adjustment = properties->tree ? treeGetAdjustment(properties->tree) : NULL;
      if (adjustment)
	setDetailViewScrollPos(gridGetDetailView(grid), baseIdx);
      
      gtk_widget_queue_draw(grid);
    }
  
  return TRUE;
}


static gboolean onMouseMoveGrid(GtkWidget *grid, GdkEventMotion *event, gpointer data)
{
  if (event->state == GDK_BUTTON2_MASK) /* middle button */
    {
      /* Draw a preview box at the mouse pointer location */
      gridSetPreviewBoxCentre(grid, event->x);
      gtk_widget_queue_draw(grid);
      gdk_window_process_updates(grid->window, TRUE);
    }
  
  return TRUE;
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

GridProperties* gridGetProperties(GtkWidget *widget)
{
  return widget ? (GridProperties*)(g_object_get_data(G_OBJECT(widget), "GridProperties")) : NULL;
}


static void onDestroyGrid(GtkWidget *widget)
{
  GridProperties *properties = gridGetProperties(widget);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "GridProperties", NULL);
    }
}


static void gridCreateProperties(GtkWidget *widget, 
				 GtkWidget *tree, 
				 GtkWidget *bigPicture,
				 gulong exposeHandlerId,
				 Strand strand)
{
  if (widget)
    { 
      /* Work out how many characters we need to display in the left border for
       * the vertical grid axis label.  This is the number of digits in the 
       * largest label value plus one (for a '%' character)
       * Add a fudge factor to give more space to allow for the fact that 
       * the calculated char width is approximate and may not give enough space */
      GridProperties *properties = g_malloc(sizeof *properties);
      
      properties->tree = tree;
      properties->bigPicture = bigPicture;
      properties->strand = strand;
      
      properties->gridYPadding = DEFAULT_GRID_Y_PADDING;
      properties->cellYPadding = DEFAULT_GRID_CELL_Y_PADDING;
      properties->highlightBoxYPad = DEFAULT_HIGHLIGHT_BOX_Y_PAD + DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH;
      
      properties->numVCells = DEFAULT_GRID_NUM_VERT_CELLS;
      properties->percentIdRange.min = DEFAULT_GRID_PERCENT_ID_MIN;
      properties->percentIdRange.max = DEFAULT_GRID_PERCENT_ID_MAX;
      properties->mspLineHeight = DEFAULT_MSP_LINE_HEIGHT;
      properties->exposeHandlerId = exposeHandlerId;
      properties->ignoreSizeAlloc = FALSE;
      
      g_object_set_data(G_OBJECT(widget), "GridProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyGrid), NULL); 
    }
}


//static GtkWidget* gridGetMainWindow(GtkWidget *grid)
//{
//  GridProperties *properties = grid ? gridGetProperties(grid) : NULL;
//  return properties ? bigPictureGetMainWindow(properties->bigPicture) : NULL;
//}

Strand gridGetStrand(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties->strand;
}

GtkWidget* gridGetBigPicture(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties->bigPicture;
}

static IntRange* gridGetDisplayRange(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetDisplayRange(bigPicture);
}

static GtkAdjustment* gridGetAdjustment(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  return detailViewGetAdjustment(detailView);
}

static GtkAdjustment* gridGetDetailView(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetDetailView(mainWindow);
}


/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

GtkWidget* createBigPictureGrid(GtkWidget *bigPicture,
				gboolean hasHeaders, 
				Strand strand)
{
  /* Create a layout area for the big picture */
  GtkWidget *grid = gtk_layout_new(NULL, NULL);
  gtk_widget_set_redraw_on_allocate(grid, FALSE);
  gtk_widget_set_name(grid, BIG_PICTURE_GRID_NAME);
  
  /* Connect signals */
  gtk_widget_add_events(grid, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(grid, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(grid, GDK_POINTER_MOTION_MASK);
  
  gulong exposeHandlerId = g_signal_connect(G_OBJECT(grid), "expose-event", G_CALLBACK(onExposeGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "size-allocate", G_CALLBACK(onSizeAllocateBigPictureGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "button-press-event", G_CALLBACK(onButtonPressGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "button-release-event", G_CALLBACK(onButtonReleaseGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "motion-notify-event", G_CALLBACK(onMouseMoveGrid), NULL);
  
  /* Set required data in the grid. We can't set the tree yet because it hasn't been created yet. */
  gridCreateProperties(grid, NULL, bigPicture, exposeHandlerId, strand);
  
  return grid;
}


