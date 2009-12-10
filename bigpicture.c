/*
 *  bigpicture.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/bigpicture.h>
#include <SeqTools/blxviewMainWindow.h>
#include <math.h>


#define DEFAULT_PREVIEW_BOX_LINE_WIDTH  1
#define DEFAULT_GRID_NUM_HEADER_LINES   1
#define DEFAULT_GRID_HEADER_Y_PAD	0
#define DEFAULT_GRID_CELL_WIDTH		100
#define DEFAULT_GRID_NUM_HOZ_CELLS	5


/* Local function declarations */
static GridHeaderProperties*	    gridHeaderGetProperties(GtkWidget *gridHeader);

/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

/* Utility to calculate how many digits are in an integer */
static int numDigitsInInt(int val)
{
  int count = 0;
  while (val > 0)
    {
      ++count;
      val /= 10;
    }
  
  return count;
}


/* Utility functions to specify some standard colours */
void setGdkColorBlack(GdkColor *color)
{
  color->red = 0;
  color->green = 0;
  color->blue = 0;
}
void setGdkColorYellow(GdkColor *color)
{
  color->red = 65535;
  color->green = 65535;
  color->blue = 0;
}
void setGdkColorBlue(GdkColor *color)
{
  color->red = 0;
  color->green = 0;
  color->blue = 65535;
}
void setGdkColorCyan(GdkColor *color)
{
  color->red = 0;
  color->green = 65535;
  color->blue = 65535;
}


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


static void drawVerticalGridLineHeaders(GtkWidget *header, GdkGC *gc, 
					const gint numCells, const gdouble cellWidth, 
					const gint rangePerCell, const GdkColor const *textColor, const GdkColor const *lineColor)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  const gint bottomBorder = properties->headerRect.y + properties->headerRect.height;
  const gint topBorder = bottomBorder - properties->markerHeight;
  
  /* Omit the first and last lines */
  gint hCell = 1;
  for ( ; hCell < numCells; ++hCell)
    {
      gint x = properties->headerRect.x + (gint)((gdouble)hCell * cellWidth);
      
      /* Draw the label */
      gdk_gc_set_rgb_fg_color(gc, textColor);
      gchar text[20];
      sprintf(text, "%d", rangePerCell * hCell);
      PangoLayout *layout = gtk_widget_create_pango_layout(header, text);
      gdk_draw_layout(GTK_LAYOUT(header)->bin_window, gc, x, 0, layout);
      g_object_unref(layout);
      
      /* Draw the marker line */
      gdk_gc_set_rgb_fg_color(gc, lineColor);
      gdk_draw_line (GTK_LAYOUT(header)->bin_window, gc, x, topBorder, x, bottomBorder);
    }
}


/* Draw a big picture header */
static void drawBigPictureGridHeader(GtkWidget *header)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Set the drawing properties */
  GdkGC *gc = gdk_gc_new(header->window);
  gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);
  
  /* Calculate some factors for scaling */
  const gint displayLen = bigPictureProperties->displayRange.max - bigPictureProperties->displayRange.min;
  const gint basesPerCell = ceil((double)displayLen / bigPictureProperties->numHCells);
  
  /* Draw the grid */
  drawVerticalGridLineHeaders(header, gc, bigPictureProperties->numHCells, 
			      bigPictureProperties->cellWidth, basesPerCell,
			      &bigPictureProperties->gridTextColor, &bigPictureProperties->gridLineColor);
  
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
  
  /* Increase the reference count to make sure the widgets aren't destroyed when we remove them. */
  g_object_ref(fwdStrandGrid);
  g_object_ref(revStrandGrid);
  
  /* Remove them */
  gtk_container_remove(GTK_CONTAINER(bigPicture), fwdStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bigPicture), revStrandGrid);
  
  /* Add them back, with the forward-strand grid at the top, or vice versa
   * if the strands are toggled. */
  if (bigPictureGetStrandsToggled(bigPicture))
    {
      addChildToBigPicture(bigPicture, revStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, fwdStrandGrid, TRUE);
    }
  else
    {
      addChildToBigPicture(bigPicture, fwdStrandGrid, FALSE);
      addChildToBigPicture(bigPicture, revStrandGrid, TRUE);
    }
  
  /* Decrease the ref count again */
  g_object_unref(fwdStrandGrid);
  g_object_unref(revStrandGrid);
}


/***********************************************************
 *			    Events			   *
 ***********************************************************/

static void onZoomInBigPicture(GtkButton *button, gpointer data)
{
  printf("onZoomInBigPicture\n");
}

static void onZoomOutBigPicture(GtkButton *button, gpointer data)
{
  printf("onZoomOutBigPicture\n");
}

static void onZoomWholeBigPicture(GtkButton *button, gpointer data)
{
  printf("onZoomWholeBigPicture\n");
}


static gboolean onExposeGridHeader(GtkWidget *header, GdkEventExpose *event, gpointer data)
{
  drawBigPictureGridHeader(header);
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


static void bigPictureCreateProperties(GtkWidget *bigPicture, 
				       GtkWidget *mainWindow, 
				       GtkWidget *header, 
				       GtkWidget *fwdStrandGrid,
				       GtkWidget *revStrandGrid,
				       IntRange *displayRange, 
				       IntRange *fullRange, 
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
      properties->displayRange.min = displayRange->min;
      properties->displayRange.max = displayRange->max;
      properties->fullRange.min = fullRange->min;
      properties->fullRange.max = fullRange->max;
      properties->cellWidth = cellWidth;
      properties->numHCells = numHCells;
      properties->previewBoxCentre = previewBoxCentre;
      properties->leftBorderChars = numDigitsInInt(DEFAULT_GRID_PERCENT_ID_MAX) + 3; /* Extra fudge factor because char width is approx */
      properties->highlightBoxLineWidth = DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH;
      properties->previewBoxLineWidth = DEFAULT_PREVIEW_BOX_LINE_WIDTH;
      
      /* Calculate the font size */
      getFontCharSize(bigPicture, &properties->charWidth, &properties->charHeight);
      
      /* Set the drawing colours */
      setGdkColorYellow(&properties->gridLineColor);
      setGdkColorBlack(&properties->gridTextColor);
      setGdkColorBlue(&properties->highlightColor);
      setGdkColorBlack(&properties->previewColor);
      
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
static GtkWidget* createButton(GtkWidget *container, char *label, GtkSignalFunc callback_func)
{
  GtkWidget *eventBox = gtk_event_box_new();
  gtk_box_pack_start(GTK_BOX(container), eventBox, FALSE, FALSE, 0);
  
  GtkWidget *button = gtk_button_new_with_label(label);
  gtk_container_add(GTK_CONTAINER(eventBox), button);
  
  g_signal_connect(GTK_OBJECT(button), "clicked", callback_func, NULL);
  
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
  GtkWidget *refButton = createButton(hbox, "+", (GtkSignalFunc)onZoomInBigPicture);
  createButton(hbox, "-", (GtkSignalFunc)onZoomOutBigPicture);
  createButton(hbox, "Whole", (GtkSignalFunc)onZoomWholeBigPicture);
  
  /* Create the header properties */
  gridHeaderCreateProperties(header, bigPicture, refButton);
  
  /* Conect signals */
  g_signal_connect(G_OBJECT(header), "expose-event", G_CALLBACK(onExposeGridHeader), NULL);
  
  return header;
}


GtkWidget* createBigPicture(GtkWidget *mainWindow, 
			    GtkWidget **header,
			    GtkWidget **fwdStrandGrid, 
			    GtkWidget **revStrandGrid, 
			    IntRange *displayRange, 
			    IntRange *fullRange)
{
  /* Create the main big picture widget, which will contain all of the 
   * individual big-picture grids, plus a header. */
  GtkWidget *bigPicture = gtk_vbox_new(FALSE, 0);
  gtk_paned_pack1(GTK_PANED(mainWindow), bigPicture, FALSE, TRUE);
  
  /* Our big picture needs to have a header plus 3 panes:
   * 1. the top grid (fwd strand);
   * 2. the exon view; and
   * 3. bottom grid (reverse strand) */
  
  *header = createBigPictureGridHeader(bigPicture);
  addChildToBigPicture(bigPicture, *header, FALSE);
  messout("Grid header created [%x]", *header);

  /* By default, make the forward strand the top grid */
  *fwdStrandGrid = createBigPictureGrid(bigPicture, TRUE, FORWARD_STRAND);
  *revStrandGrid = createBigPictureGrid(bigPicture, FALSE, REVERSE_STRAND);
  
  messout("Grid1 created [%x]", *fwdStrandGrid);
  messout("Grid2 created [%x]", *revStrandGrid);
  
  addChildToBigPicture(bigPicture, *fwdStrandGrid, FALSE);
  addChildToBigPicture(bigPicture, *revStrandGrid, TRUE);

  /* Set the big picture properties */
  bigPictureCreateProperties(bigPicture, 
			     mainWindow,
			     *header, 
			     *fwdStrandGrid,
			     *revStrandGrid,
			     displayRange, 
			     fullRange, 
			     DEFAULT_GRID_CELL_WIDTH, 
			     DEFAULT_GRID_NUM_HOZ_CELLS, 
			     UNSET_INT);
  
  
  return bigPicture;
}
