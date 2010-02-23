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
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/utilities.h>

#define EXON_VIEW_DEFAULT_COMPRESSED_HEIGHT      10
#define EXON_VIEW_DEFAULT_EXPANDED_HEIGHT        20
#define	EXON_VIEW_DEFAULT_Y_PADDING		 7

typedef struct _ExonViewProperties
  {
    GtkWidget *bigPicture;	      /* The big picture that this view belongs to */
    
    gboolean expanded;		      /* Whether the exon view is expanded or compressed */
    
    int compressedHeight;	      /* The height of the view when it is compressed */
    int expandedHeight;		      /* The height of the view when it is expanded */
    int yPad;			      /* y padding */
    
    GdkRectangle exonViewRect;	      /* The drawing area for the exon view */
    GdkColor exonColour;	      /* The colour to draw normal exons */
    GdkColor exonColourSelected;      /* The colour to draw selected exons */
  } ExonViewProperties;


/* Local function declarations */
static GtkWidget*		exonViewGetBigPicture(GtkWidget *exonView);
static GtkWidget*		exonViewGetMainWindow(GtkWidget *exonView);
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
static void drawExonIntron(const MSP *msp, GtkWidget *exonView, const gboolean isSelected)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  IntRange *displayRange = bigPictureGetDisplayRange(properties->bigPicture);
  gboolean rightToLeft = bigPictureGetStrandsToggled(properties->bigPicture);
  
  int x1 = convertBaseIdxToGridPos(msp->displayStart, &properties->exonViewRect, displayRange, rightToLeft);
  int x2 = convertBaseIdxToGridPos(msp->displayEnd, &properties->exonViewRect, displayRange, rightToLeft);
  int x = min(x1, x2);
  int width = abs(x2 - x1);
  
  int y = properties->exonViewRect.y;
  int height = properties->exonViewRect.height;
  
  GdkDrawable *drawable = widgetGetDrawable(exonView);
  GdkGC *gc = gdk_gc_new(drawable);
  
//  /* to do: If this msp is in a group, use the group's colour. Otherwise use the standard exon colour. */
//  SequenceGroup *group = mainWindowGetSequenceGroup(bigPictureGetMainWindow(properties->bigPicture), msp->sname);
//  if (group)
//    {
//      gdk_gc_set_foreground(gc, &group->highlightColour);
//    }
//  else
    {
      gdk_gc_set_foreground(gc, isSelected ? &properties->exonColourSelected : &properties->exonColour);
    }

  if (mspIsExon(msp))
    {
      drawExon(drawable, gc, x, y, width, height);
    }
  else if (mspIsIntron(msp))
    {
      drawIntron(drawable, gc, x, y, width, height);
    }
    
  g_object_unref(gc);
}


/* Draw the msp in the given row if it is an exon/intron and is unselected */
static gboolean drawUnselectedExonIntron(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  /* One row can contain multiple MSPs. Draw them all. */
  const GList *mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
	{
	  GtkWidget *exonView = GTK_WIDGET(data);
	  const gboolean isSelected = mainWindowIsSeqSelected(exonViewGetMainWindow(exonView), msp->sname);
	  
	  if (!isSelected)
	    {
	      drawExonIntron(msp, exonView, isSelected);
	    }
	}
    }

  return FALSE;
}


/* Draw the msp in the given row if it is an exon/intron and is selected */
static gboolean drawSelectedExonIntron(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  /* One row can contain multiple MSPs. Draw them all. */
  const GList *mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
	{
	  GtkWidget *exonView = GTK_WIDGET(data);
	  const gboolean isSelected = mainWindowIsSeqSelected(exonViewGetMainWindow(exonView), msp->sname);
	  
	  if (isSelected)
	    {
	      drawExonIntron(msp, exonView, isSelected);
	    }
	}
    }
  
  return FALSE;
}


/* Draw all of the exons/introns in the given tree */
static void drawExonsIntronsForTree(GtkWidget *exonView, GtkWidget *tree)
{
  /* Loop through all of the (unfiltered) rows in the tree. Loop twice, first
   * drawing unselected msp then selected ones, so that the selected ones
   * are drawn on top */
  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
  
  gtk_tree_model_foreach(model, drawUnselectedExonIntron, exonView);
  gtk_tree_model_foreach(model, drawSelectedExonIntron, exonView);
}


/* Draw the exon view */
static void drawExonView(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);

  int numFrames = detailViewGetNumReadingFrames(detailView);
  int frame = 1;
  
  for ( ; frame <= numFrames; ++frame)
    {
      GtkWidget *tree = detailViewGetTree(detailView, FORWARD_STRAND, frame);
      drawExonsIntronsForTree(exonView, tree);

      tree = detailViewGetTree(detailView, REVERSE_STRAND, frame);
      drawExonsIntronsForTree(exonView, tree);
    }
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

static void exonViewCreateProperties(GtkWidget *exonView, GtkWidget *bigPicture)
{
  if (exonView)
    {
      ExonViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture	      = bigPicture;
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

static GtkWidget* exonViewGetMainWindow(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  return bigPictureGetMainWindow(bigPicture);
}

static GtkWidget* exonViewGetTopGrid(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  
  GtkWidget *topGrid = NULL;
  
  if (bigPictureGetStrandsToggled(bigPicture))
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
GtkWidget *createExonView(GtkWidget *bigPicture)
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

  exonViewCreateProperties(exonView, bigPicture);

  return exonView;
}

