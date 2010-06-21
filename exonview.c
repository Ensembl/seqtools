/*
 *  exonview.c
 *  blixem
 *
 *  Created by Gemma Barson on 24/12/2009.
 *
 */

#include <SeqTools/exonview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/bigpicturegrid.h>
#include <SeqTools/detailview.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/utilities.h>

#define DEFAULT_EXON_HEIGHT			 10
#define	DEFAULT_EXON_YPAD			 7

typedef struct _ExonViewProperties
  {
    GtkWidget *bigPicture;	      /* The big picture that this view belongs to */
    BlxStrand currentStrand;	      /* Which strand of the ref seq this view displays exons for */
    
    gboolean expanded;		      /* Whether the exon view is expanded or compressed */
    
    int yPad;			      /* y padding */
    
    GdkRectangle exonViewRect;	      /* The drawing area for the exon view */
  } ExonViewProperties;


/* Utility struct to pass around data required for drawing exons */
typedef struct _DrawData
  {
    GdkDrawable *drawable;
    GdkGC *gc;
    GdkColor *cdsFillColor;
    GdkColor *cdsFillColorSelected;
    GdkColor *cdsLineColor;
    GdkColor *cdsLineColorSelected;
    GdkColor *utrFillColor;
    GdkColor *utrFillColorSelected;
    GdkColor *utrLineColor;
    GdkColor *utrLineColorSelected;
    GdkRectangle *exonViewRect;
    GtkWidget *blxWindow;
    BlxViewContext *bc;
    const BlxStrand strand;
    const IntRange const *displayRange;
    const IntRange const *refSeqRange;
    const gboolean displayRev;
    const int numFrames;
    const BlxSeqType seqType;
    gboolean expanded;		      /* whether exon view is expaned or compressed */
    gboolean normalOnly;	      /* Only draw "normal" MSPs, i.e. thse that are not selected and are not in a group */
    int y;			      /* y position to draw this exon at (constant if view compressed; increases if view bumped) */
    int height;			      /* height of exon box */
  } DrawData;



/* Local function declarations */
static GtkWidget*		exonViewGetBigPicture(GtkWidget *exonView);
static GtkWidget*		exonViewGetBlxWindow(GtkWidget *exonView);
static ExonViewProperties*	exonViewGetProperties(GtkWidget *exonView);
static GtkWidget*		exonViewGetTopGrid(GtkWidget *exonView);


/***********************************************************
 *                       Utility functions                 *
 ***********************************************************/

static GdkColor* getFillColor(const MSP const *msp, const gboolean isSelected, DrawData *data, const BlxSequence *blxSeq)
{
  if (mspIsCds(msp, blxSeq))
    return (isSelected ? data->cdsFillColorSelected : data->cdsFillColor);
  else /* UTR and undefined exon types */
    return (isSelected ? data->utrFillColorSelected : data->utrFillColor);
}

static GdkColor* getLineColor(const MSP const *msp, const gboolean isSelected, DrawData *data, const BlxSequence *blxSeq)
{
  if (mspIsCds(msp, blxSeq))
    return (isSelected ? data->cdsLineColorSelected : data->cdsLineColor);
  else 
    return (isSelected ? data->utrLineColorSelected : data->utrLineColor);
}


/* Draw an exon */
static void drawExon(const MSP const *msp, GdkDrawable *drawable, DrawData *data, const BlxSequence *blxSeq, const gboolean isSelected, int x, int y, int width, int height)
{
  /* Draw the fill rectangle */
  gdk_gc_set_foreground(data->gc, getFillColor(msp, isSelected, data, blxSeq));
  gdk_draw_rectangle(drawable, data->gc, TRUE, x, y, width, height);
  
  /* Draw outline (exon box outline always the same (unselected) color; only intron lines change when selected) */
  gdk_gc_set_foreground(data->gc, getLineColor(msp, FALSE, data, blxSeq));
  gdk_draw_rectangle(drawable, data->gc, FALSE, x, y, width, height);
  
}


/* Draw an intron */
static void drawIntron(const MSP const *msp, GdkDrawable *drawable, DrawData *data, const BlxSequence *blxSeq, const gboolean isSelected, int x, int y, int width, int height)
{
  gdk_gc_set_foreground(data->gc, getLineColor(msp, isSelected, data, blxSeq));

  int xMid = x + roundNearest((double)width / 2.0);
  int xEnd = x + width;
  int yMid = y + roundNearest((double)height / 2.0);
  
  gdk_draw_line(drawable, data->gc, x, yMid, xMid, y);
  gdk_draw_line(drawable, data->gc, xMid, y, xEnd, yMid);
}


/* Draw the given exon/intron, if it is in range. Returns true if it was drawn */
static gboolean drawExonIntron(const MSP *msp, DrawData *data, const gboolean isSelected, const BlxSequence *blxSeq)
{
  gboolean drawn = FALSE;
  
  const int frame = mspGetRefFrame(msp, data->seqType);

  /* Find the coordinates of the start and end base in this msp, converting to display coords */
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, data->seqType, frame, data->numFrames, data->displayRev, data->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, data->seqType, frame, data->numFrames, data->displayRev, data->refSeqRange, NULL);
  
  if (valueWithinRange(coord1, data->displayRange) || valueWithinRange(coord2, data->displayRange))
    {
      drawn = TRUE;

      /* The grid pos gives the left edge of the coord, so to be inclusive we draw to the max coord + 1 */
      const int minCoord = min(coord1, coord2);
      const int maxCoord = max(coord1, coord2) + 1;
      
      int xMin = convertBaseIdxToGridPos(minCoord, data->exonViewRect, data->displayRange);
      int xMax = convertBaseIdxToGridPos(maxCoord, data->exonViewRect, data->displayRange);
      
      int x = xMin;
      int width = xMax - xMin;
      int y = data->y;
      int height = data->height;
      
      if (mspIsExon(msp))
	{
	  drawExon(msp, data->drawable, data, blxSeq, isSelected, x, y, width, height);
	}
      else if (mspIsIntron(msp))
	{
	  drawIntron(msp, data->drawable, data, blxSeq, isSelected, x, y, width, height);
	}
    }
  
  return drawn;
}


/* Returns true if the given MSP should be shown in this exon view */
static gboolean showMspInExonView(const MSP *msp, DrawData *drawData)
{
  /* Check it's an exon or intron */
  gboolean showMsp = mspIsExon(msp) || mspIsIntron(msp);
  
  /* Check it's the correct strand */
  showMsp &= (mspGetRefStrand(msp) == drawData->strand);
  
  return showMsp;
}

/* Draw the msps in the given sequence, if they are exons/introns. Use the color
 * specified in the user data */
static void drawExonIntronItem(gpointer listItemData, gpointer data)
{
  const BlxSequence *seq = (const BlxSequence*)listItemData;
  DrawData *drawData = (DrawData*)data;

  const gboolean isSelected = blxWindowIsSeqSelected(drawData->blxWindow, seq);
  SequenceGroup *group = blxWindowGetSequenceGroup(drawData->blxWindow, seq);
  gboolean seqDrawn = FALSE;
  
  if (!drawData->normalOnly || (!isSelected && !group))
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
  
  /* If the view is expanded, increase the y-coord for the next sequence */
  if (seqDrawn && drawData->expanded)
    {
      drawData->y += drawData->height + DEFAULT_EXON_YPAD;
    }
}


/* Draw the exon view */
static void drawExonView(GtkWidget *exonView, GdkDrawable *drawable)
{
  GtkWidget *blxWindow = exonViewGetBlxWindow(exonView);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  GdkGC *gc = gdk_gc_new(drawable);
  
  DrawData drawData = {
    drawable,
    gc,
    getGdkColor(bc, BLXCOL_EXON_FILL_CDS, FALSE),
    getGdkColor(bc, BLXCOL_EXON_FILL_CDS, TRUE),
    getGdkColor(bc, BLXCOL_EXON_LINE_CDS, FALSE),
    getGdkColor(bc, BLXCOL_EXON_LINE_CDS, TRUE),
    getGdkColor(bc, BLXCOL_EXON_FILL_UTR, FALSE),
    getGdkColor(bc, BLXCOL_EXON_FILL_UTR, TRUE),
    getGdkColor(bc, BLXCOL_EXON_LINE_UTR, FALSE),
    getGdkColor(bc, BLXCOL_EXON_LINE_UTR, TRUE),
    &properties->exonViewRect,
    blxWindow,
    blxContext,
    properties->currentStrand,
    bigPictureGetDisplayRange(properties->bigPicture),
    &blxContext->refSeqRange,
    blxContext->displayRev,
    blxContext->numFrames,
    blxContext->seqType,
    properties->expanded,
    FALSE,
    properties->exonViewRect.y,
    properties->exonViewRect.height
  };
  
  /* If the view is compressed (i.e. exons will overlap each other), then
   * only draw "normal" MSPs the first time round, and draw grouped/selected
   * MSPs afterwards, so that they appear on top. If the view is expanded, 
   * we can draw them all in a single loop, because they will not overlap. */
  drawData.normalOnly = !properties->expanded;
  
  /* Loop through all sequences, drawing all msps that are exons/introns */
  GList *seqList = blxWindowGetAllMatchSeqs(blxWindow);
  g_list_foreach(seqList, drawExonIntronItem, &drawData);

  if (!properties->expanded)
    {
      drawData.normalOnly = FALSE;
  
//      /* Draw all msps that are in groups */
//      GList *groupItem = blxContext->sequenceGroups;
//      for ( ; groupItem; groupItem = groupItem->next)
//	{
//	  SequenceGroup *group = (SequenceGroup*)(groupItem->data);
//	  drawData.color = group->highlighted ? &group->highlightColor : drawData.exonColor;
//	  g_list_foreach(group->seqList, drawExonIntronItem, &drawData);
//	}
      
      /* Draw all selected msps */
      g_list_foreach(blxContext->selectedSeqs, drawExonIntronItem, &drawData);
      
      /* Increment the y value when finished, because we calculate the view height based on this */
      drawData.y += drawData.height + DEFAULT_EXON_YPAD;
    }

  /* Set the height based on the height of the exons that were actually drawn */
  const int newHeight = drawData.y - properties->exonViewRect.y + DEFAULT_EXON_YPAD;
  gtk_layout_set_size(GTK_LAYOUT(exonView), exonView->allocation.width, newHeight);
//  gtk_widget_set_size_request(exonView, -1, newHeight);
  
  g_object_unref(gc);
}


void calculateExonViewHeight(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  const IntRange const *displayRange = bigPictureGetDisplayRange(properties->bigPicture);

  /* Calculate the height based on how many exon lines will actually be drawn */
  int numExons = 0;
  int maxExons = properties->expanded ? UNSET_INT : 1; /* unset means no limit */
  
  /* Loop through all sequences */
  GList *seqItem = blxWindowGetAllMatchSeqs(bigPictureGetBlxWindow(properties->bigPicture));
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      /* Loop through all msps */
      const BlxSequence *seq = (BlxSequence*)(seqItem->data);
      GList *mspItem = seq->mspList;
      
      for ( ; mspItem; mspItem = mspItem->next)
	{
	  const MSP *msp = (const MSP*)(mspItem->data);
	  
	  if ((mspIsExon(msp) || mspIsIntron(msp)) && mspGetRefStrand(msp) == properties->currentStrand)
	    {
	      const int frame = mspGetRefFrame(msp, bc->seqType);
	      const int startCoord = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
	      const int endCoord = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
	      
	      if (valueWithinRange(startCoord, displayRange) || valueWithinRange(endCoord, displayRange))
		{
		  ++numExons;
		  break; /* break inner loop and move to next sequence */
		}
	    }
	}
      
      /* Break after we've found the maximum number of lines, if a max is specified */
      if (maxExons != UNSET_INT && numExons >= maxExons)
	{
	  break;
	}
    }
  
  const int newHeight = (numExons * (DEFAULT_EXON_HEIGHT + DEFAULT_EXON_YPAD)) + (2 * DEFAULT_EXON_YPAD);
  
  gtk_widget_set_size_request(exonView, -1, newHeight);
}


void calculateExonViewBorders(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Calculate the size of the exon view */
  properties->exonViewRect.x = bigPictureProperties->charWidth * bigPictureProperties->leftBorderChars;
  properties->exonViewRect.width = exonView->allocation.width - properties->exonViewRect.x;
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
				     GtkWidget *bigPicture, 
				     const BlxStrand currentStrand)
{
  if (exonView)
    {
      ExonViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture	      = bigPicture;
      properties->currentStrand	      = currentStrand;
      
      properties->expanded	      = FALSE;
      properties->yPad		      =	DEFAULT_EXON_YPAD;
      
      properties->exonViewRect.x      = 0;
      properties->exonViewRect.y      = DEFAULT_EXON_YPAD;
      properties->exonViewRect.width  = 0;
      properties->exonViewRect.height = DEFAULT_EXON_HEIGHT;
      
      gtk_widget_set_size_request(exonView, 0, DEFAULT_EXON_HEIGHT + (2 * DEFAULT_EXON_YPAD));

      g_object_set_data(G_OBJECT(exonView), "ExonViewProperties", properties);
      g_signal_connect(G_OBJECT(exonView), "destroy", G_CALLBACK(onDestroyExonView), NULL);
    }
}

static GtkWidget* exonViewGetBigPicture(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  return properties->bigPicture;
}

static GtkWidget* exonViewGetBlxWindow(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  return bigPictureGetBlxWindow(bigPicture);
}

static GtkWidget* exonViewGetTopGrid(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  
  GtkWidget *topGrid = NULL;
  
  if (bigPictureGetDisplayRev(bigPicture))
    {
      topGrid = bigPictureGetRevGrid(bigPicture);
    }
  else
    {
      topGrid = bigPictureGetFwdGrid(bigPicture);
    }
      
  return topGrid;
}

gboolean exonViewGetExpanded(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  return properties->expanded;
}

/* Set whether the view is expanded or not */
void exonViewSetExpanded(GtkWidget *exonView, const gboolean expanded)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  properties->expanded = expanded;

  calculateExonViewHeight(exonView);
  bigPictureRedrawAll(properties->bigPicture);
}


/* Toggle whether the view is expanded or not */
void exonViewToggleExpanded(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  exonViewSetExpanded(exonView, !properties->expanded);
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
      GdkGC *gc = gdk_gc_new(GTK_LAYOUT(exonView)->bin_window);
      gdk_draw_drawable(GTK_LAYOUT(exonView)->bin_window, gc, drawable, 0, 0, 0, 0, -1, -1);
    }
  
  return TRUE;
}

static void onSizeAllocateExonView(GtkWidget *exonView, GtkAllocation *allocation, gpointer data)
{
  calculateExonViewBorders(exonView);
}


static gboolean onButtonPressExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->button == 2) /* middle button */
    {
      GtkWidget *grid = exonViewGetTopGrid(exonView);
      showPreviewBox(grid, event->x);
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  if (event->button == 2) /* middle button */
    {
      GtkWidget *grid = exonViewGetTopGrid(exonView);
      acceptAndClearPreviewBox(grid, event->x);
    }
  
  return TRUE;
}


static gboolean onMouseMoveExonView(GtkWidget *exonView, GdkEventMotion *event, gpointer data)
{
  if (event->state == GDK_BUTTON2_MASK) /* middle button */
    {
      /* Draw a preview box at the mouse pointer location */
      GtkWidget *grid = exonViewGetTopGrid(exonView);
      showPreviewBox(grid, event->x);
    }
  
  return TRUE;
}


/***********************************************************
 *                       Initialisation                    *
 ***********************************************************/

/* Create the part of the view that will show the exons */
GtkWidget *createExonView(GtkWidget *bigPicture, const BlxStrand currentStrand)
{
  GtkWidget *exonView = gtk_layout_new(NULL, NULL);
  
  /* Connect signals */
  gtk_widget_add_events(exonView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(exonView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(exonView, GDK_POINTER_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(exonView),	"expose-event",		G_CALLBACK(onExposeExonView),	      NULL);
  g_signal_connect(G_OBJECT(exonView),	"size-allocate",	G_CALLBACK(onSizeAllocateExonView),   NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-press-event",   G_CALLBACK(onButtonPressExonView),    NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-release-event", G_CALLBACK(onButtonReleaseExonView),  NULL);
  g_signal_connect(G_OBJECT(exonView),	"motion-notify-event",  G_CALLBACK(onMouseMoveExonView),      NULL);

  exonViewCreateProperties(exonView, bigPicture, currentStrand);

  return exonView;
}
