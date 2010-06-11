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
#include <SeqTools/blxwindow.h>
#include <SeqTools/utilities.h>
#include <math.h>
#include <string.h>


#define BIG_PICTURE_MSP_LINE_NAME	"BigPictureMspLine"
#define DEFAULT_MSP_LINE_HEIGHT		3	  /* the height of the MSP lines in the grid */
#define DEFAULT_GRID_Y_PADDING		5	  /* this provides space between the grid and the edge of the widget */
#define DEFAULT_GRID_CELL_Y_PADDING	-2	  /* this controls the vertical space between the labels on the y axis */
#define DEFAULT_HIGHLIGHT_BOX_Y_PAD	2	  /* this provides space between highlight box and the top/bottom of the grid */
#define MIN_MSP_LINE_WIDTH		1	  /* used to make sure that MSP lines never shrink to nothing */

typedef struct _DrawGridData
{
  GtkWidget *grid;
  GdkDrawable *drawable;
  GdkGC *gc;
  GdkColor *color;
  GdkColor *shadowColor;
} DrawGridData;

  

/* Local function declarations */
static BlxViewContext*	    gridGetContext(GtkWidget *grid);
static GtkAdjustment*	    gridGetAdjustment(GtkWidget *grid);
static IntRange*	    gridGetDisplayRange(GtkWidget *grid);
static GdkColor*	    gridGetMspLineColor(GtkWidget *grid, const gboolean selected);
static GtkWidget*	    gridGetDetailView(GtkWidget *grid);
static GtkWidget*	    gridGetBlxWindow(GtkWidget *grid);
static void                 drawBigPictureGrid(GtkWidget *grid, GdkDrawable *drawable);

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


static int gridGetNumVCells(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetNumVCells(bigPicture);
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
				    const IntRange const *displayRange)
{
  gint result = UNSET_INT;
  
  int distFromEdge = (int)((gdouble)gridPos - (gdouble)gridRect->x);
  int basesFromEdge = distFromEdge / pixelsPerBase(gridRect->width, displayRange);
  result = displayRange->min + basesFromEdge;
  
  return result;
}


/* Convert an ID% value to the y coord on the given grid */
gint convertValueToGridPos(GtkWidget *grid, const gint value)
{
  /* The top line of the grid is drawn one cell height down from the top of the grid border */
  GridProperties *properties = gridGetProperties(grid);
  
  IntRange *valRange = bigPictureGetPercentIdRange(properties->bigPicture);
  gdouble percent = (gdouble)(valRange->max - value) / (gdouble)(valRange->max - valRange->min);
  
  /* Make sure we do the multiplication on doubles before rounding to int */
  gint result = properties->gridRect.y + roundNearest((gdouble)properties->gridRect.height * percent);
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
static void drawHighlightBox(GtkWidget *grid, 
			     GdkDrawable *drawable,
			     GdkGC *gc, 
			     const GdkRectangle const *rect,
			     const gint lineWidth, 
			     GdkColor *lineColor)
{
  gdk_gc_set_foreground(gc, lineColor);
  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  gdk_draw_rectangle(drawable, gc, FALSE, 
		     rect->x, rect->y, rect->width, rect->height);
}


/* Draw the preview box. This is drawn directly onto the grid window (not onto the
 * cached bitmap that the rest of the grid is drawn onto). Only does anything if the
 * preview box centre is set. */
static void drawPreviewBox(GtkWidget *grid, GdkDrawable *drawable, GdkGC *gc)
{
  int previewBoxCentre = gridGetPreviewBoxCentre(grid);
  
  if (previewBoxCentre == UNSET_INT)
    {
      return;
    }
  
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  IntRange *displayRange = gridGetDisplayRange(grid);

  /* Find the x coord for the left edge of the preview box (or the right edge, if
   * the display is right-to-left). */
  int x = getLeftCoordFromCentre(previewBoxCentre, properties->highlightRect.width, &properties->gridRect);
  
  /* Convert it to the base index and back again so that we get it rounded to the position of
   * the nearest base. */
  int baseIdx = convertGridPosToBaseIdx(x, &properties->gridRect, displayRange);
  int xRounded = convertBaseIdxToGridPos(baseIdx, &properties->gridRect, displayRange);
  
  /* The other dimensions of the preview box are the same as the current highlight box. */
  GdkRectangle previewRect = {xRounded, properties->highlightRect.y, properties->highlightRect.width, properties->highlightRect.height};

  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  GdkColor *previewBoxColor = getGdkColor(bc, BLXCOL_PREVIEW_BOX, FALSE);
  drawHighlightBox(grid, drawable, gc, &previewRect, bigPictureProperties->previewBoxLineWidth, previewBoxColor);
}


/* Draw the vertical gridlines in the big picture view */
static void drawVerticalGridLines(GtkWidget *grid, 
                                  GdkDrawable *drawable,
                                  GdkGC *gc,
				  const GdkColor const *lineColor)
{
  BlxViewContext *bc = gridGetContext(grid);
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
  
  const int direction = bc->displayRev ? -1 : 1; /* to subtract instead of add when display reversed */
  
  /* Get the first base index (in terms of the nucleotide coords) and round it to a nice round
   * number. We'll offset all of the gridlines by the distance between this and the real start coord. */
  const int realFirstBaseIdx = convertDisplayIdxToDnaIdx(bpProperties->displayRange.min,bc-> seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int firstBaseIdx = roundToValue(realFirstBaseIdx, bpProperties->roundTo);
  
  /* Calculate the top and bottom heights for the lines. */
  const gint topBorder = properties->highlightRect.y - properties->gridYPadding;
  const gint bottomBorder = properties->gridRect.y + properties->gridRect.height;
  
  const int minX = properties->gridRect.x;
  const int maxX = properties->gridRect.x + properties->gridRect.width;

  gint hCell = 0;
  for ( ; hCell <= bpProperties->numHCells; ++hCell)
    {
      /* Get the base index for this grid line and calc its x coord */
      int numBasesFromLeft = bpProperties->basesPerCell * hCell;
      int baseIdx = firstBaseIdx + (numBasesFromLeft * direction);

      const int displayIdx = convertDnaIdxToDisplayIdx(baseIdx, bc->seqType, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
      const int x = convertBaseIdxToGridPos(displayIdx, &properties->gridRect, &bpProperties->displayRange);

      if (x > minX && x < maxX)
	{
	  gdk_gc_set_foreground(gc, lineColor);
	  gdk_draw_line (drawable, gc, x, topBorder, x, bottomBorder);
	}
    }
}


/* Draw the horizontal grid lines for the big picture view */
static void drawHorizontalGridLines(GtkWidget *grid,
                                    GdkDrawable *drawable,
                                    GdkGC *gc,
				    const gint numCells, 
				    const gint rangePerCell, 
				    const gint maxVal, 
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
      gdk_gc_set_foreground(gc, textColor);
      
      char text[bigPictureProperties->leftBorderChars + 1];
      sprintf(text, "%d%%", percent);
      
      PangoLayout *layout = gtk_widget_create_pango_layout(grid, text);
      gdk_draw_layout(drawable, gc, 0, y - gridGetCellHeight(grid)/2, layout);
      g_object_unref(layout);
      
      /* Draw the gridline */
      gdk_gc_set_foreground(gc, lineColor);
      gdk_draw_line (drawable, gc, properties->gridRect.x, y, rightBorder, y);
    }
}


/* Calculates the size and position of an MSP line in the given grid. Return
 * args can be null if not required. */
void calculateMspLineDimensions(GtkWidget *grid, 
				const MSP const *msp, 
				int *x, 
				int *y, 
				int *width, 
				int *height)
{
  BlxViewContext *bc = gridGetContext(grid);
  GridProperties *gridProperties = gridGetProperties(grid);

  const IntRange const *displayRange = bigPictureGetDisplayRange(gridProperties->bigPicture);
  const int frame = mspGetRefFrame(msp, bc->seqType);

  /* Find the coordinates of the start and end base in this match sequence, and convert to display coords */
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qstart, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qend, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  
  /* Convert the coords to grid positions. The grid positions we use are for the left edge
   * of the coord, so to draw the end coord inclusively we need to increase the max coord by 1 */
  const int minCoord = min(coord1, coord2);
  const int maxCoord = max(coord1, coord2) + 1;
  
  int xMin = convertBaseIdxToGridPos(minCoord, &gridProperties->gridRect, displayRange);
  int xMax = convertBaseIdxToGridPos(maxCoord, &gridProperties->gridRect, displayRange);

  if (x)
    {
      *x = xMin;
    }
    
  if (width)
    {
      *width = max((xMax - xMin), MIN_MSP_LINE_WIDTH);
    }
  
  /* Find where in the y axis we should draw the line, based on the %ID value */
  if (y)
    {
      *y = convertValueToGridPos(grid, msp->id);
    }
    
  if (height)
    {
      *height = gridProperties->mspLineHeight;
    }
}


/* Returns true if the given msp is displayed in the given grid, i.e. is the 
 * correct strand, is not an intron or exon, and is not out of range */
static gboolean mspShownInGrid(const MSP const *msp, GtkWidget *grid)
{
  gboolean result = FALSE;
  
  if (mspIsBlastMatch(msp) && mspGetRefStrand(msp) == gridGetStrand(grid))
    {
      /* Convert the msp's dna coords to display coords */
      BlxViewContext *bc = gridGetContext(grid);
      const IntRange const *displayRange = gridGetDisplayRange(grid);
      
      const int mspStart = convertDnaIdxToDisplayIdx(msp->qstart, bc->seqType, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
      const int mspEnd = convertDnaIdxToDisplayIdx(msp->qend, bc->seqType, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
      const int qMin = min(mspStart, mspEnd);
      const int qMax = max(mspStart, mspEnd);
      
      if (qMin <= displayRange->max && qMax >= displayRange->min)
	{
	  result = TRUE;
	}
    }
  
  return result;
}


/* Draw a line on the given grid to represent the given match */
static void drawMspLine(const MSP const *msp, DrawGridData *drawData)
{
  /* Ignore introns. */
  if (mspShownInGrid(msp, drawData->grid))
    {
      gdk_gc_set_subwindow(drawData->gc, GDK_INCLUDE_INFERIORS);
      gdk_gc_set_foreground(drawData->gc, drawData->color);
      
      /* Calculate where it should go */
      int x, y, width, height;
      calculateMspLineDimensions(drawData->grid, msp, &x, &y, &width, &height);
      
      /* Draw a block rectangle */
      gdk_draw_rectangle(drawData->drawable, drawData->gc, TRUE, x, y, width, height);
      
      /* Draw a drop-shadow, to make sure it is always visible even in the paler highlight colors */
      const int shadowHt = 1;
      gdk_gc_set_foreground(drawData->gc, drawData->shadowColor);
      gdk_gc_set_line_attributes(drawData->gc, shadowHt, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      int yBottom = y + height - shadowHt;
      gdk_draw_line(drawData->drawable, drawData->gc, x, yBottom, x + width, yBottom);
      
      /* only draw the right-hand-side if the width is great enough to fit it and still be visible */
      if (width > shadowHt + 1)
      {
	int xRight = x + width - shadowHt;
	gdk_draw_line(drawData->drawable, drawData->gc, xRight, y, xRight, y + height);
      }
    }
}


/* Draw the MSPs for the given sequence in the given color. */
static void drawSequenceMspLines(gpointer listItemData, gpointer data)
{
  const BlxSequence *seq = (const BlxSequence*)listItemData;
  DrawGridData *drawData = (DrawGridData*)data;  
  GList *mspListItem = seq->mspList;

  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      if (mspShownInGrid(msp, drawData->grid))
	{
	  drawMspLine(msp, drawData);
	}
    }
}


/* Draw MSP lines that are in groups. Use the highlight color of the group if one is set */
static void drawGroupedMspLines(gpointer listItemData, gpointer data)
{
  SequenceGroup *seqGroup = (SequenceGroup*)listItemData;
  DrawGridData *drawData = (DrawGridData*)data;
  GdkColor *origColor = drawData->color;
  GdkColor *origShadowColor = drawData->shadowColor;

  /* Use the group's highlight color, if this sequence should be highlighted */
  if (seqGroup->highlighted)
    {
      /* Use the group's highlight color */
      drawData->color = &seqGroup->highlightColor;
      GdkColor shadowColor;
      getDropShadowColor(drawData->color, &shadowColor);
      drawData->shadowColor = &shadowColor;
      
      g_list_foreach(seqGroup->seqList, drawSequenceMspLines, drawData);

      /* Set the draw data back to its original values */
      drawData->color = origColor;
      drawData->shadowColor = origShadowColor;
    }
  else
    {
      g_list_foreach(seqGroup->seqList, drawSequenceMspLines, drawData);
    }
  
}


/* Draw a line for each MSP in the given grid */
static void drawMspLines(GtkWidget *grid, GdkDrawable *drawable, GdkGC *gc)
{
  BlxViewContext *bc = gridGetContext(grid);

  DrawGridData drawData = {
    grid, 
    drawable, 
    gc, 
    gridGetMspLineColor(grid, FALSE), 
    gridGetMspLineColor(grid, FALSE)
  };
  
  /* Draw all MSPs for this grid */ 
  g_list_foreach(bc->matchSeqs, drawSequenceMspLines, &drawData);  

  /* Now draw MSPs that are in groups (to do: it would be good to do this in reverse
   * Sort Order, so that those ordered first get drawn last and therefore appear on top) */
  g_list_foreach(bc->sequenceGroups, drawGroupedMspLines, &drawData);
  
  /* Finally, draw selected sequences. These will appear on top of everything else. */
  drawData.color = gridGetMspLineColor(grid, TRUE);
  GdkColor shadowColor;
  getDropShadowColor(drawData.color, &shadowColor);
  drawData.shadowColor = &shadowColor;
  g_list_foreach(bc->selectedSeqs, drawSequenceMspLines, &drawData);
}


/* Main function that does the drawing for the grid. Drawing is done onto a pixmap
 * which is then stored in the grid properties */
static void drawBigPictureGrid(GtkWidget *grid, GdkDrawable *drawable)
{
  /* Calculate some factors for scaling */
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);

  const gint percentPerCell = bigPictureGetIdPerCell(properties->bigPicture);
  const gint numVCells = gridGetNumVCells(grid);

  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  GdkColor *gridLineColor = getGdkColor(bc, BLXCOL_GRID_LINE, FALSE);
  GdkColor *gridTextColor = getGdkColor(bc, BLXCOL_GRID_TEXT, FALSE);
  GdkColor *highlightBoxColor = getGdkColor(bc, BLXCOL_HIGHLIGHT_BOX, FALSE);
  
  /* Draw the grid lines */
  GdkGC *gc = gdk_gc_new(drawable);
  drawVerticalGridLines(grid, drawable, gc, gridLineColor);
  drawHorizontalGridLines(grid, drawable, gc, numVCells, percentPerCell, bigPictureProperties->percentIdRange.max, gridTextColor, gridLineColor);
  
  /* Draw lines corresponding to the MSPs */
  drawMspLines(grid, drawable, gc);
  
  /* Draw the highlight box */
  drawHighlightBox(grid, 
		   drawable,
		   gc,
		   &properties->highlightRect, 
		   bigPictureProperties->highlightBoxLineWidth,
		   highlightBoxColor);
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
      int firstBaseIdx = adjustment->value;

      properties->highlightRect.x = convertBaseIdxToGridPos(firstBaseIdx, &properties->gridRect, displayRange);
      properties->highlightRect.y = properties->gridRect.y - properties->highlightBoxYPad;
      
      properties->highlightRect.width = roundNearest((gdouble)adjustment->page_size * pixelsPerBase(properties->gridRect.width, displayRange));
      properties->highlightRect.height = properties->gridRect.height + properties->mspLineHeight + (2 * properties->highlightBoxYPad);
    }
}


/* Calculate the borders for the big picture view */
void calculateGridBorders(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  int numVCells = gridGetNumVCells(grid);
  
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
  properties->gridRect.height = gridGetCellHeight(grid) * numVCells;
  
  /* Get the boundaries of the highlight box */
  calculateHighlightBoxBorders(grid);
  
  /* Get the total display height required. Set the layout size to fit. */
  properties->displayRect.height = properties->highlightRect.height + (2 * properties->gridYPadding);
  gtk_layout_set_size(GTK_LAYOUT(grid), properties->displayRect.width, properties->displayRect.height);
  
  /* Set the size request to our desired height. We want a fixed heigh but don't set the
   * width, because we want the user to be able to resize that. */
  gtk_widget_set_size_request(grid, 0, properties->displayRect.height);
}


/* Show a preview box centred on the given x coord */
void showPreviewBox(GtkWidget *grid, const int x)
{
  /* Set the centre position of the preview box. When this is set the box will
   * be drawn when the grid is exposed. */
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
  
  /* Find the base index where the new scroll range will start. This is the leftmost
   * edge of the preview box if numbers increase in the normal left-to-right direction, 
   * or the rightmost edge if the display is reversed. */
  int x = getLeftCoordFromCentre(xCentre, properties->highlightRect.width, &properties->gridRect);
  int baseIdx = convertGridPosToBaseIdx(x, &properties->gridRect, displayRange);
  
  /* Clear the preview box */
  gridSetPreviewBoxCentre(grid, UNSET_INT);
  
  /* Update the detail view's scroll pos to start at the new base. The base index is in terms of
   * the display coords, so the coord's sequence type is whatever the displayed sequence type is */
  GtkWidget *detailView = gridGetDetailView(grid);
  const BlxSeqType seqType = detailViewGetSeqType(detailView);  /* displayed seq type */
  setDetailViewStartIdx(detailView, baseIdx, seqType);

  gtk_widget_queue_draw(grid);
}


/***********************************************************
 *			    Events			   *
 ***********************************************************/

/* Expose handler for the grid. If the grid has a bitmap cached from the previous
 * call, this function just pushes that back to screen and then (optionally) draws
 * the preview box directly to the window over the top of it. If the cached bitmap
 * is null, this function first calls drawBigPictureGrid to create and cache
 * the bitmap first. */
static gboolean onExposeGrid(GtkWidget *grid, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(grid)->bin_window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(grid);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(grid);
          drawBigPictureGrid(grid, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc2 = gdk_gc_new(window);
          gdk_draw_drawable(window, gc2, bitmap, 0, 0, 0, 0, -1, -1);

          /* Draw the preview box on top, if it is set */
          drawPreviewBox(grid, window, gc2);
        }
      else
	{
	  g_warning("Failed to draw grid [%p] - could not create bitmap.\n", grid);
	}
    }
  
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


/* Mark the given MSP's sequence as selected if this MSP line contains the given coords.
 * Returns true if it was selected. */
static gboolean selectMspIfContainsCoords(GtkWidget *grid, 
					  const MSP *msp,
					  const int x,
					  const int y,
					  gboolean deselectOthers)
{
  gboolean wasSelected = FALSE;
  
  if (mspShownInGrid(msp, grid))
    {
      int mspX, mspY, mspWidth, mspHeight;
      calculateMspLineDimensions(grid, msp, &mspX, &mspY, &mspWidth, &mspHeight);
      
      if (x >= mspX && x <= mspX + mspWidth && y >= mspY && y <= mspY + mspHeight)
	{
	  /* It's a hit. Select this sequence. */
	  GtkWidget *blxWindow = gridGetBlxWindow(grid);
	  
	  if (deselectOthers)
	    {
	      blxWindowDeselectAllSeqs(blxWindow);
	    }
	  
	  blxWindowSelectSeq(blxWindow, msp->sSequence);
	  
	  /* Update the selected strand */
	  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
	  detailViewSetSelectedStrand(detailView, gridGetStrand(grid));

	  /* Scroll the detail view trees to bring the new selection into view */
	  callFuncOnAllDetailViewTrees(detailView, treeScrollSelectionIntoView);

	  wasSelected = TRUE;
	}
    }
  
  return wasSelected;
}


/* Loop through all the msp lines for this grid and mark them as selected
 * if they contain the coords of the mouse press */
static void selectClickedMspLines(GtkWidget *grid, GdkEventButton *event, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  /* Loop through all the MSPs until we find one under the click coords */
  GtkWidget *blxWindow = gridGetBlxWindow(grid);
  const MSP *msp = blxWindowGetMspList(blxWindow);
  const gboolean deselectOthers = !ctrlModifier && !shiftModifier; /* whether to deselect all others first */
  
  for ( ; msp; msp = msp->next)
    {
      const gboolean found = selectMspIfContainsCoords(grid, msp, event->x, event->y, deselectOthers);

      if (found)
	break;
    }
}


static gboolean onButtonPressGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->button == 1) /* left button */
    {
      /* If we clicked on top of an msp line, select that msp */
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
      const gboolean shiftModifier = ((event->state & modifiers) == GDK_SHIFT_MASK);

      selectClickedMspLines(grid, event, ctrlModifier, shiftModifier);
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
  if (event->state & GDK_BUTTON2_MASK) /* middle button */
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

static BlxViewContext* gridGetContext(GtkWidget *grid)
{
  GtkWidget *blxWindow = gridGetBlxWindow(grid);
  return blxWindowGetContext(blxWindow);
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
				 GtkWidget *bigPicture,
				 gulong exposeHandlerId,
				 BlxStrand strand)
{
  if (widget)
    { 
      /* Work out how many characters we need to display in the left border for
       * the vertical grid axis label.  This is the number of digits in the 
       * largest label value plus one (for a '%' character)
       * Add a fudge factor to give more space to allow for the fact that 
       * the calculated char width is approximate and may not give enough space */
      GridProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->strand = strand;
      
      properties->gridYPadding = DEFAULT_GRID_Y_PADDING;
      properties->cellYPadding = DEFAULT_GRID_CELL_Y_PADDING;
      properties->highlightBoxYPad = DEFAULT_HIGHLIGHT_BOX_Y_PAD + DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH;
      
      properties->mspLineHeight = DEFAULT_MSP_LINE_HEIGHT;
      properties->exposeHandlerId = exposeHandlerId;
      properties->ignoreSizeAlloc = FALSE;
      
      g_object_set_data(G_OBJECT(widget), "GridProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyGrid), NULL); 
    }
}


static GtkWidget* gridGetBlxWindow(GtkWidget *grid)
{
  GridProperties *properties = grid ? gridGetProperties(grid) : NULL;
  return properties ? bigPictureGetBlxWindow(properties->bigPicture) : NULL;
}


BlxStrand gridGetStrand(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties->strand;
}

GtkWidget* gridGetBigPicture(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties ? properties->bigPicture : NULL;
}

static IntRange* gridGetDisplayRange(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetDisplayRange(bigPicture);
}

static GtkAdjustment* gridGetAdjustment(GtkWidget *grid)
{
  GtkWidget *detailView = gridGetDetailView(grid);
  return detailViewGetAdjustment(detailView);
}

static GtkWidget* gridGetDetailView(GtkWidget *grid)
{
  GtkWidget *blxWindow = gridGetBlxWindow(grid);
  return blxWindowGetDetailView(blxWindow);
}

static GdkColor *gridGetMspLineColor(GtkWidget *grid, const gboolean selected)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  return getGdkColor(bc, BLXCOL_MSP_LINE, selected);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

GtkWidget* createBigPictureGrid(GtkWidget *bigPicture, BlxStrand strand)
{
  /* Create a layout area for the grid */
  GtkWidget *grid = gtk_layout_new(NULL, NULL);
  
  /* Grid style properties */
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
  gridCreateProperties(grid, bigPicture, exposeHandlerId, strand);
  
  return grid;
}


