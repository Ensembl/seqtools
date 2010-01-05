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
#include <SeqTools/blxviewMainWindow.h>
#include <math.h>
#include <string.h>


#define BIG_PICTURE_MSP_LINE_NAME	"BigPictureMspLine"
#define DEFAULT_GRID_NUM_VERT_CELLS	5
#define DEFAULT_GRID_PERCENT_ID_MIN	0
#define DEFAULT_MSP_LINE_HEIGHT		3
#define DEFAULT_GRID_Y_PADDING		5
#define DEFAULT_GRID_CELL_Y_PADDING	-2
#define DEFAULT_HIGHLIGHT_BOX_Y_PAD	2


/* Local function declarations */
static GtkAdjustment*	    gridGetAdjustment(GtkWidget *grid);
static IntRange*	    gridGetDisplayRange(GtkWidget *grid);
static gboolean		    gridGetStrandsToggled(GtkWidget *grid);
static GdkColor*	    gridGetMspLineHighlightColour(GtkWidget *grid);
static GdkColor*	    gridGetMspLineColour(GtkWidget *grid);
static GtkWidget*	    gridGetTree(GtkWidget *grid);
static GtkWidget*	    gridGetDetailView(GtkWidget *grid);

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


/* Convert an x coord on the given grid to a base index */
static gint convertGridPosToBaseIdx(const gint gridPos, 
				    const GdkRectangle const *gridRect,  
				    const IntRange const *displayRange,
				    gboolean rightToLeft)
{
  gint result = UNSET_INT;
  
  gdouble basesFromEdge = round(((gdouble)gridPos - (gdouble)gridRect->x) / pixelsPerBase(gridRect->width, displayRange));
  
  if (rightToLeft)
    {
      result = displayRange->max - basesFromEdge;
    }
  else
    {
      result = displayRange->min + basesFromEdge;
    }
  
  return result;
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
static void drawHighlightBox(GtkWidget *widget, 
			     GdkGC *gc, 
			     const GdkRectangle const *rect,
			     const gint lineWidth, 
			     GdkColor *lineColour)
{
  gdk_gc_set_foreground(gc, lineColour);
  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  gdk_draw_rectangle(GTK_LAYOUT(widget)->bin_window, gc, FALSE, 
		     rect->x, rect->y, rect->width, rect->height);
}


/* Draw the preview box, if applicable */
static void drawPreviewBox(GtkWidget *grid, 
			   GdkGC *gc, 
			   const GdkRectangle const *highlightRect,
			   const gint lineWidth, 
			   GdkColor *lineColour)
{
  /* If a preview position is set, draw a preview of the highlight rect 
   * centred at that position */
  int previewBoxCentre = gridGetPreviewBoxCentre(grid);
  if (previewBoxCentre != UNSET_INT)
    {
      GridProperties *properties = gridGetProperties(grid);
      IntRange *displayRange = gridGetDisplayRange(grid);
      gboolean rightToLeft = gridGetStrandsToggled(grid);

      /* Find the x coord for the left edge of the preview box. Convert it to the base
       * index and back again so that we get it rounded to the position of the nearest base. */
      int xLeft = getLeftCoordFromCentre(previewBoxCentre, highlightRect->width, &properties->gridRect);
      int baseIdx = convertGridPosToBaseIdx(xLeft, &properties->gridRect, displayRange, rightToLeft);
      int xRounded = convertBaseIdxToGridPos(baseIdx, &properties->gridRect, displayRange, rightToLeft);
      
//      /* Find the x coord for the left edge of the preview box. */
//      int x = getLeftCoordFromCentre(previewBoxCentre, highlightRect->width, &properties->gridRect);
      
      /* The other dimensions of the preview box are the same as the current highlight box. */
      GdkRectangle previewRect = {xRounded, highlightRect->y, highlightRect->width, highlightRect->height};
      drawHighlightBox(grid, gc, &previewRect, lineWidth, lineColour);
    }
}


/* Draw the vertical gridlines in the big picture view */
static void drawVerticalGridLines(GtkWidget *grid, GdkGC *gc, 
				  const gint numCells, const gdouble cellWidth, 
				  const GdkColor const *lineColour)
{
  GridProperties *properties = gridGetProperties(grid);
  
  const gint topBorder = properties->highlightRect.y - properties->gridYPadding;
  const gint bottomBorder = properties->gridRect.y + properties->gridRect.height;
  
  /* Omit the first and last lines */
  gint hCell = 1;
  for ( ; hCell < numCells; ++hCell)
    {
      gint x = properties->gridRect.x + (gint)((gdouble)hCell * cellWidth);
      gdk_gc_set_foreground(gc, lineColour);
      gdk_draw_line (GTK_LAYOUT(grid)->bin_window, gc, x, topBorder, x, bottomBorder);
    }
}


/* Draw the horizontal grid lines for the big picture view */
static void drawHorizontalGridLines(GtkWidget *grid, GdkGC *gc, 
				    const gint numCells, 
				    const gint rangePerCell, const gint maxVal, 
				    const GdkColor const *textColour,
				    const GdkColor const *lineColour)
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
      gdk_gc_set_foreground(gc, textColour);
      
      char text[bigPictureProperties->leftBorderChars + 1];
      sprintf(text, "%d%%", percent);
      
      PangoLayout *layout = gtk_widget_create_pango_layout(grid, text);
      gdk_draw_layout(GTK_LAYOUT(grid)->bin_window, gc, 0, y - gridGetCellHeight(grid)/2, layout);
      g_object_unref(layout);
      
      /* Draw the gridline */
      gdk_gc_set_foreground(gc, lineColour);
      gdk_draw_line (GTK_LAYOUT(grid)->bin_window, gc, properties->gridRect.x, y, rightBorder, y);
    }
}


/* Calculates the size and position of an MSP line in the given grid. Return
 * args can be null if not required. */
void calculateMspLineDimensions(GtkWidget *grid, const MSP const *msp, 
				int *x, int *y, int *width, int *height)
{
  GridProperties *gridProperties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(gridProperties->bigPicture);

  gboolean rightToLeft = bigPictureGetStrandsToggled(gridProperties->bigPicture);
  
  /* Find the coordinates of the start and end base in this match sequence */
  int qSeqMin, qSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, NULL, NULL);
  int x1 = convertBaseIdxToGridPos(qSeqMin, &gridProperties->gridRect, &bigPictureProperties->displayRange, rightToLeft);
  int x2 = convertBaseIdxToGridPos(qSeqMax + 1, &gridProperties->gridRect, &bigPictureProperties->displayRange, rightToLeft);

  /* We'll start drawing at the lowest x coord, so set x to the min of x1 and x2 */
  if (x)
    *x = min(x1, x2);
  
  if (width)
    *width = abs(x1 - x2);
  
  /* Find where in the y axis we should draw the line, based on the %ID value */
  if (y)
    *y = convertValueToGridPos(grid, msp->id);
  
  if (height)
    *height = gridProperties->mspLineHeight;
}


/* Returns true if the given msp is displayed in the grid, i.e. is not
 * an intron or something */
static gboolean mspShownInGrid(const MSP const *msp)
{
  return !mspIsFake(msp) && !mspIsIntron(msp) && !mspIsExon(msp);
}


/* Draw a line on the given grid to represent the given match */
static void drawMspLine(GtkWidget *grid, GdkColor *colour, const MSP const *msp)
{
  /* Ignore "fake" MSPs and introns. */
  if (mspShownInGrid(msp))
    {
      GdkGC *gc = gdk_gc_new(grid->window);
      gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);
      gdk_gc_set_foreground(gc, colour);
      
      /* Calculate where it should go */
      int x, y, width, height;
      calculateMspLineDimensions(grid, msp, &x, &y, &width, &height);
      
      /* Draw it */
      gdk_draw_rectangle(GTK_LAYOUT(grid)->bin_window, gc, TRUE, x, y, width, height);
      
      g_object_unref(gc);
    }
}


/* Wrapper function around drawMspLine that only draws unselected MSPs. */
static gboolean drawUnselectedMspLine(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *grid = GTK_WIDGET(data);
  GtkWidget *tree = gridGetTree(grid);

  if (!treePathIsSelected(GTK_TREE_VIEW(tree), path, model))
    {
      const MSP* msp = treeGetMsp(model, iter);
      drawMspLine(grid, gridGetMspLineColour(grid), msp);
    }
  
  return FALSE;
}


/* Wrapper function around drawMspLine that only draws selected MSPs. */
static gboolean drawSelectedMspLine(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *grid = GTK_WIDGET(data);
  GtkWidget *tree = gridGetTree(grid);
  
  if (treePathIsSelected(GTK_TREE_VIEW(tree), path, model))
    {
      const MSP* msp = treeGetMsp(model, iter);
      drawMspLine(grid, gridGetMspLineHighlightColour(grid), msp);
    }
  
  return FALSE;
}


/* Draw a line for each MSP in the given grid */
static void drawMspLines(GtkWidget *grid)
{
  GtkWidget *tree = gridGetTree(grid);

  if (tree)
    {
      /* We'll loop through every row in the base (i.e. unfiltered) data. */
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      
      /* Loop twice, first drawing unselected MSPs and then selected MSPs, so
       * that selected MSPs are always drawn on top. */
      gtk_tree_model_foreach(model, drawUnselectedMspLine, grid);
      gtk_tree_model_foreach(model, drawSelectedMspLine, grid);
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
  drawVerticalGridLines(grid, 
			gc, 
			bigPictureProperties->numHCells, 
			bigPictureProperties->cellWidth,
			&bigPictureProperties->gridLineColour);
  
  drawHorizontalGridLines(grid, 
			  gc, 
			  properties->numVCells, 
			  percentPerCell, 
			  properties->percentIdRange.max, 
			  &bigPictureProperties->gridTextColour, 
			  &bigPictureProperties->gridLineColour);
  
  /* Draw lines corresponding to the MSPs */
  drawMspLines(grid);
  
  /* Draw the highlight box */
  drawHighlightBox(grid, 
		   gc,
		   &properties->highlightRect, 
		   bigPictureProperties->highlightBoxLineWidth,
		   &bigPictureProperties->highlightBoxColour);
  
  drawPreviewBox(grid,
		 gc,
		 &properties->highlightRect, 
		 bigPictureProperties->previewBoxLineWidth, 
		 &bigPictureProperties->previewBoxColour);
  
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
      gboolean rightToLeft = gridGetStrandsToggled(grid);
      int firstBaseIdx = rightToLeft ? adjustment->value + adjustment->page_size + 1 : adjustment->value + 1;

      properties->highlightRect.x = convertBaseIdxToGridPos(firstBaseIdx, &properties->gridRect, displayRange, rightToLeft);
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


/* Show a preview box centred on the given x coord */
void showPreviewBox(GtkWidget *grid, const int x)
{
  /* Draw a preview box at the clicked position */
  gridSetPreviewBoxCentre(grid, x);
  
  /* Force immediate update so that it happens before the button-release event */
  gtk_widget_queue_draw(grid);
  gdk_window_process_updates(grid->window, TRUE);
}


/* Scroll to the position of the preview box, and clear the preview */
void acceptAndClearPreviewBox(GtkWidget *grid, const int xCentre)
{
  GridProperties *properties = gridGetProperties(grid);
  IntRange *displayRange = gridGetDisplayRange(grid);
  gboolean rightToLeft = gridGetStrandsToggled(grid);
  
  /* Find the base index where the new scroll range will start. This is the leftmost
   * edge of the preview box if numbers increase in the normal left-to-right direction, 
   * or the rightmost edge if the display is reversed. */
  int x = rightToLeft
    ? getRightCoordFromCentre(xCentre, properties->highlightRect.width, &properties->gridRect)
    : getLeftCoordFromCentre(xCentre, properties->highlightRect.width, &properties->gridRect);
  
  int baseIdx = convertGridPosToBaseIdx(x, &properties->gridRect, displayRange, rightToLeft);
  
  /* Clear the preview box */
  gridSetPreviewBoxCentre(grid, UNSET_INT);
  
  /* Update the detail view's scroll pos to start at the start base */
  GtkAdjustment *adjustment = properties->tree ? treeGetAdjustment(properties->tree) : NULL;
  
  if (adjustment)
    {
      setDetailViewScrollPos(adjustment, baseIdx - 1); /* adjust for 0-indexing */
    }
  
  gtk_widget_queue_draw(grid);
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
  calculateGridBorders(grid);
}


/* Mark the given row as selected if its MSP line contains the given coords,
 * Returns true if it was selected. */
static gboolean selectRowIfContainsCoords(GtkWidget *grid, 
					  GtkWidget *tree, 
					  GtkTreeModel *model, 
					  GtkTreeIter *iter,
					  const int x,
					  const int y,
					  gboolean deselectOthers)
{
  gboolean wasSelected = FALSE;
  
  const MSP *msp = treeGetMsp(model, iter);
  
  if (mspShownInGrid(msp))
    {
      int mspX, mspY, mspWidth, mspHeight;
      calculateMspLineDimensions(grid, msp, &mspX, &mspY, &mspWidth, &mspHeight);
      
      if (x >= mspX && x <= mspX + mspWidth && y >= mspY && y <= mspY + mspHeight)
	{
	  /* It's a hit */
	  wasSelected = TRUE;
	  
	  if (deselectOthers)
	    deselectAllSiblingTrees(tree, TRUE);
	  
	  selectRow(GTK_TREE_VIEW(tree), model, iter);
	}
    }
  
  return wasSelected;
}


/* Loop through all the msp lines for this grid and mark them as selected
 * if they contain the coords of the mouse press */
static void selectClickedMspLines(GtkWidget *grid, GdkEventButton *event)
{
  /* The msp info lives in the tree */
  GtkWidget *tree = gridGetTree(grid);

  if (tree)
    {
      /* We'll loop through every row in the base (i.e. unfiltered) data */
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      
      /* If multiple overlap, just select the first one we find for now. */
      gboolean done = FALSE;
      
      GtkTreeIter iter;
      gboolean validIter = gtk_tree_model_get_iter_first(model, &iter);
      
      while (validIter && !done)
	{
	  done = selectRowIfContainsCoords(grid, tree, model, &iter, event->x, event->y, TRUE);
	  validIter = gtk_tree_model_iter_next(model, &iter);
	}
    }
}


static gboolean onButtonPressGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->button == 1) /* left button */
    {
      /* If we clicked on top of an msp line, select that msp */
      selectClickedMspLines(grid, event);
      handled = TRUE;
    }
  else if (event->button == 2) /* middle button */
    {
      showPreviewBox(grid, event->x);
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  if (event->button == 2) /* middle button */
    {
      acceptAndClearPreviewBox(grid, event->x);
    }
  
  return TRUE;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the grid.
 * This scrolls the position of the highlight box, i.e. it scrolls the display
 * range in the detail view. */
static gboolean onScrollGrid(GtkWidget *grid, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
    {
      case GDK_SCROLL_LEFT:
	{
	  scrollDetailViewLeftStep(gridGetDetailView(grid));
	  handled = TRUE;
	  break;
	}
	
      case GDK_SCROLL_RIGHT:
	{
	  scrollDetailViewRightStep(gridGetDetailView(grid));
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


static gboolean onMouseMoveGrid(GtkWidget *grid, GdkEventMotion *event, gpointer data)
{
  if (event->state == GDK_BUTTON2_MASK) /* middle button */
    {
      /* Draw a preview box at the mouse pointer location */
      showPreviewBox(grid, event->x);
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

static GtkWidget* gridGetTree(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties->tree;
}

static IntRange* gridGetDisplayRange(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetDisplayRange(bigPicture);
}

static gboolean gridGetStrandsToggled(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetStrandsToggled(bigPicture);
}

static GtkAdjustment* gridGetAdjustment(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  return detailViewGetAdjustment(detailView);
}

static GtkWidget* gridGetDetailView(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  GtkWidget *mainWindow = bigPictureGetMainWindow(bigPicture);
  return mainWindowGetDetailView(mainWindow);
}

static GdkColor *gridGetMspLineColour(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->mspLineColour;
}

static GdkColor *gridGetMspLineHighlightColour(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->mspLineHighlightColour;
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

GtkWidget* createBigPictureGrid(GtkWidget *bigPicture, Strand strand)
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
  g_signal_connect(G_OBJECT(grid), "size-allocate",	    G_CALLBACK(onSizeAllocateBigPictureGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "button-press-event",    G_CALLBACK(onButtonPressGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "button-release-event",  G_CALLBACK(onButtonReleaseGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "motion-notify-event",   G_CALLBACK(onMouseMoveGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "scroll-event",	    G_CALLBACK(onScrollGrid),		      NULL);
  
  /* Set required data in the grid. We can't set the tree yet because it hasn't been created yet. */
  gridCreateProperties(grid, NULL, bigPicture, exposeHandlerId, strand);
  
  return grid;
}


