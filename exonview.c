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

#define EXON_VIEW_DEFAULT_COMPRESSED_HEIGHT      10
#define EXON_VIEW_DEFAULT_EXPANDED_HEIGHT        20
#define	EXON_VIEW_DEFAULT_Y_PADDING		 7

typedef struct _ExonViewProperties
  {
    GtkWidget *bigPicture;	      /* The big picture that this view belongs to */
    Strand currentStrand;		      /* Which strand of the ref seq this view displays exons for */
    
    gboolean expanded;		      /* Whether the exon view is expanded or compressed */
    
    int compressedHeight;	      /* The height of the view when it is compressed */
    int expandedHeight;		      /* The height of the view when it is expanded */
    int yPad;			      /* y padding */
    
    GdkRectangle exonViewRect;	      /* The drawing area for the exon view */
    GdkColor exonColour;	      /* The colour to draw normal exons */
    GdkColor exonColourSelected;      /* The colour to draw selected exons */
  } ExonViewProperties;


/* Utility struct to pass around data required for drawing exons */
typedef struct _DrawData
  {
    GdkDrawable *drawable;
    GdkGC *gc;
    GdkColor *colour;
    GdkRectangle *exonViewRect;
    GtkWidget *blxWindow;
    const Strand strand;
    const IntRange const *displayRange;
    const IntRange const *refSeqRange;
    const gboolean displayRev;
    const int numFrames;
    const BlxSeqType seqType;
  } DrawData;



/* Local function declarations */
static GtkWidget*		exonViewGetBigPicture(GtkWidget *exonView);
static GtkWidget*		exonViewGetBlxWindow(GtkWidget *exonView);
static ExonViewProperties*	exonViewGetProperties(GtkWidget *exonView);
static GtkWidget*		exonViewGetTopGrid(GtkWidget *exonView);


/***********************************************************
 *                       Utility functions                 *
 ***********************************************************/

/* Draw an exon */
static void drawExon(GdkDrawable *drawable, GdkGC *gc, int x, int y, int width, int height)
{
  gdk_draw_rectangle(drawable, gc, FALSE, x, y, width, height);
  
}


/* Draw an intron */
static void drawIntron(GdkDrawable *drawable, GdkGC *gc, int x, int y, int width, int height)
{
  int xMid = x + roundNearest((double)width / 2.0);
  int xEnd = x + width;
  int yMid = y + roundNearest((double)height / 2.0);
  
  gdk_draw_line(drawable, gc, x, yMid, xMid, y);
  gdk_draw_line(drawable, gc, xMid, y, xEnd, yMid);
}


/* Draw the given exon/intron */
static void drawExonIntron(const MSP *msp, DrawData *data)
{
  const int frame = mspGetRefFrame(msp, data->seqType);

  /* Find the coordinates of the start and end base in this msp, converting to display coords */
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qstart, data->seqType, frame, data->numFrames, data->displayRev, data->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qend, data->seqType, frame, data->numFrames, data->displayRev, data->refSeqRange, NULL);
  
  /* The grid pos gives the left edge of the coord, so to be inclusive we draw to the max coord + 1 */
  const int minCoord = min(coord1, coord2);
  const int maxCoord = max(coord1, coord2) + 1;
  
  int xMin = convertBaseIdxToGridPos(minCoord, data->exonViewRect, data->displayRange);
  int xMax = convertBaseIdxToGridPos(maxCoord, data->exonViewRect, data->displayRange);
  
  int x = xMin;
  int width = xMax - xMin;
  
  int y = data->exonViewRect->y;
  int height = data->exonViewRect->height;
  
  gdk_gc_set_foreground(data->gc, data->colour);

  if (mspIsExon(msp))
    {
      drawExon(data->drawable, data->gc, x, y, width, height);
    }
  else if (mspIsIntron(msp))
    {
      drawIntron(data->drawable, data->gc, x, y, width, height);
    }
}


/* Draw the msps in the given sequence, if they are exons/introns. Use the colour
 * specified in the user data */
static void drawExonIntronItem(gpointer listItemData, gpointer data)
{
  const char *seqName = (const char*)listItemData;
  DrawData *drawData = (DrawData*)data;
  
  /* Loop through all msps in this sequence */
  GList *mspListItem = blxWindowGetSequenceMsps(drawData->blxWindow, seqName);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if ((mspIsExon(msp) || mspIsIntron(msp)) && mspGetRefStrand(msp) == drawData->strand)
	{
	  drawExonIntron(msp, drawData);
	}
    }
}


/* Draw the exon view */
static void drawExonView(GtkWidget *exonView)
{
  GtkWidget *blxWindow = exonViewGetBlxWindow(exonView);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  const MSP *msp = blxContext->mspList;

  GdkDrawable *drawable = widgetGetDrawable(exonView);
  GdkGC *gc = gdk_gc_new(drawable);
  
  DrawData drawData = {
    drawable,
    gc,
    &properties->exonColour,
    &properties->exonViewRect,
    blxWindow,
    properties->currentStrand,
    bigPictureGetDisplayRange(properties->bigPicture),
    &blxContext->refSeqRange,
    blxContext->displayRev,
    blxContext->numFrames,
    blxContext->seqType
  };
  
  /* Loop through all msps drawing only unselected ones */
  for ( ; msp; msp = msp->next)
    {
      if ((mspIsExon(msp) || mspIsIntron(msp)) && mspGetRefStrand(msp) == properties->currentStrand)
	{
	  if (!blxWindowIsSeqSelected(blxWindow, msp->sname));
	    {
	      drawExonIntron(msp, &drawData);
	    }
	}
    }
  
  /* Draw grouped msps */
  GList *groupItem = blxContext->sequenceGroups;
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      drawData.colour = group->highlighted ? &group->highlightColour : &properties->exonColour;
      g_list_foreach(group->seqNameList, drawExonIntronItem, &drawData);
    }
  
  /* Draw all selected msps */
  drawData.colour = &properties->exonColourSelected;
  g_list_foreach(blxContext->selectedSeqs, drawExonIntronItem, &drawData);
  
  g_object_unref(gc);
}


void calculateExonViewBorders(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Calculate the size of the exon view */
  properties->exonViewRect.x = bigPictureProperties->charWidth * bigPictureProperties->leftBorderChars;
  properties->exonViewRect.width = exonView->allocation.width - properties->exonViewRect.x;
  
  gtk_layout_set_size(GTK_LAYOUT(exonView), exonView->allocation.width, properties->exonViewRect.height);
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
				     const Strand currentStrand)
{
  if (exonView)
    {
      ExonViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture	      = bigPicture;
      properties->currentStrand	      = currentStrand;
      
      properties->expanded	      = FALSE;      
      properties->yPad		      =	EXON_VIEW_DEFAULT_Y_PADDING;
      properties->compressedHeight    =	EXON_VIEW_DEFAULT_COMPRESSED_HEIGHT + (2 * EXON_VIEW_DEFAULT_Y_PADDING);
      properties->expandedHeight      = EXON_VIEW_DEFAULT_EXPANDED_HEIGHT + (2 * EXON_VIEW_DEFAULT_Y_PADDING);
      
      properties->exonViewRect.x      = 0;
      properties->exonViewRect.y      = EXON_VIEW_DEFAULT_Y_PADDING;
      properties->exonViewRect.width  = 0;
      properties->exonViewRect.height = EXON_VIEW_DEFAULT_COMPRESSED_HEIGHT;
      
      properties->exonColour = getGdkColor(GDK_BLUE);
      properties->exonColourSelected = getGdkColor(GDK_CYAN);
      
      gtk_widget_set_size_request(exonView, 0, properties->compressedHeight);

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


/***********************************************************
 *                       Events                            *
 ***********************************************************/

static gboolean onExposeExonView(GtkWidget *exonView, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *drawable = gdk_pixmap_new(GTK_LAYOUT(exonView)->bin_window, exonView->allocation.width, exonView->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(exonView, drawable);

  GdkGC *gc = gdk_gc_new(GTK_LAYOUT(exonView)->bin_window);
  GtkStyle *style = gtk_widget_get_style(exonView);
  GdkColor *bgColour = &style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColour);
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, exonView->allocation.width, exonView->allocation.height);

  /* Draw the exon view onto the pixmap */
  drawExonView(exonView);
  
  /* Push the pixmap onto the screen */
  gdk_draw_drawable(GTK_LAYOUT(exonView)->bin_window, gc, drawable, 0, 0, 0, 0, -1, -1);
  
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
GtkWidget *createExonView(GtkWidget *bigPicture, const Strand currentStrand)
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

