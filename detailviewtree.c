/*
 *  detailviewtree.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/detailviewtree.h>
#include <SeqTools/detailview.h>
#include <SeqTools/bigpicturegrid.h>
#include <SeqTools/bigpicturemspline.h>
#include <SeqTools/sequencecellrenderer.h>
#include <string.h>

#define DETAIL_VIEW_TREE_NAME		"DetailViewTreeName"


/***************************************************************
 *                       Sequence column                       *
 * A specially-formatted column for displaying match sequences *
 ***************************************************************/

typedef struct _SeqColProperties
  {
    GtkWidget *tree;
  } SeqColProperties;

static SeqColProperties* seqColGetProperties(GtkTreeViewColumn *column)
{
  return column ? (SeqColProperties*)(g_object_get_data(G_OBJECT(column), "SeqColProperties")) : NULL;
}

static void onDestroySeqCol(GtkTreeViewColumn *column)
{
  SeqColProperties *properties = seqColGetProperties(column);
  if (properties)
    {
      free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(column), "SeqColProperties", NULL);
    }
}

static void seqColCreateProperties(GtkTreeViewColumn *column, GtkWidget *tree)
{
  if (column)
    {
      SeqColProperties *properties = malloc(sizeof *properties);
      properties->tree = tree;
      g_object_set_data(G_OBJECT(column), "SeqColProperties", properties);
      g_signal_connect(G_OBJECT(column), "destroy", G_CALLBACK(onDestroySeqCol), NULL); 
    }
}


/***********************************************************
 *                Tree - utility functions                 *
 ***********************************************************/

GtkAdjustment *treeGetAdjustment(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(properties->detailView);
  return detailViewProperties ? detailViewProperties->adjustment : NULL;
}

GtkCellRenderer *treeGetRenderer(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->renderer : NULL;
}

static GtkWidget *treeGetDetailView(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->detailView : NULL;
}


/* Return the current strand that this tree is displaying (takes into account toggle status) */
Strand treeGetStrand(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return gridGetStrand(properties->grid);
}

int treeGetNumReadingFrames(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetNumReadingFrames(properties->detailView) : UNSET_INT;
}

int treeGetSelectedBaseIdx(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetSelectedBaseIdx(properties->detailView) : UNSET_INT;
}

static void treeSetSelectedBaseIdx(GtkWidget *tree, const int selectedBaseIdx)
{
  TreeProperties *properties = treeGetProperties(tree);
  if (properties)
    {
      DetailViewProperties *detailViewProperties = detailViewGetProperties(properties->detailView);
      
      /* Only update if things have changed */
      if (detailViewProperties && detailViewProperties->selectedBaseIdx != selectedBaseIdx)
	{
	  detailViewProperties->selectedBaseIdx = selectedBaseIdx;
	  
	  /* Redraw all trees in the detail view */
	  gtk_widget_queue_draw(properties->detailView);
	}
    }
}

char* treeGetRefSeq(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetRefSeq(properties->detailView) : NULL;
}

IntRange* treeGetDisplayRange(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetDisplayRange(properties->detailView) : NULL;
}

/* Calls the given function on the given widget if it is a tree in the detail-view, or, if it
 * is a container, calls the function on all children/grandchildren/etc that are dettail-view trees */
void callFuncOnAllDetailViewTrees(GtkWidget *widget, gpointer data)
{
  GtkCallback func = (GtkCallback)data;

  const gchar *name = gtk_widget_get_name(widget);
  if (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0)
    {
      func(widget, NULL);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllDetailViewTrees, func);
    }
}


/* Deselect all rows in the given tree */
void deselectAllRows(GtkWidget *tree)
{
  GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
  gtk_tree_selection_unselect_all(selection);
}


/* Deselect all rows in the given widget (if it is a tree) or in all of its child 
 * trees (if it is a container), unless it is the tree that is passed as the data pointer,
 * in which case ignore it. */
void deselectAllRowsNotInTree(GtkWidget *widget, gpointer data)
{
  if (widget != data)
    {
      const gchar *name = gtk_widget_get_name(widget);
      if (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0)
	{
	  deselectAllRows(widget);
	}
      else if (GTK_IS_CONTAINER(widget))
	{
	  gtk_container_foreach(GTK_CONTAINER(widget), deselectAllRowsNotInTree, data);
	}
    }
}


/* Call the deselect-all function on all trees that are in the same detail view 
 * as the given tree, but NOT on the given tree itself. */
void deselectAllSiblingTrees(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  gtk_container_foreach(GTK_CONTAINER(detailView), deselectAllRowsNotInTree, tree);
}


/* Return the msp in a given tree row */
const MSP* treeGetMsp(GtkTreeView *tree, GtkTreeModel *model, GtkTreeIter *iter)
{
  GValue val = {0};
  gtk_tree_model_get_value(model, iter, MSP_COL, &val);
  return (MSP*)g_value_get_pointer(&val);
}


/* For the given tree view, return the data model used for the visible view
 * (which may be a filtered tree model). */
GtkTreeModel* treeGetVisibleDataModel(GtkTreeView *tree)
{
  /* Just return the model stored directly in the tree view. This will be the filter
   * if there is one, otherwise it will be the original data model. */
  return tree ? gtk_tree_view_get_model(tree) : NULL;
}


/* For the given tree view, return the base data model i.e. all data in the tree
 * without any filtering. Returns the same as treeGetVisibleDataModel if the tree does
 * not have a filter. */
GtkTreeModel* treeGetBaseDataModel(GtkTreeView *tree)
{
  GtkTreeModel *result = NULL;
  
  GtkTreeView *treeView = tree;
  if (treeView)
    {
      GtkTreeModel *model = gtk_tree_view_get_model(treeView);
      if (GTK_IS_TREE_MODEL_FILTER(model))
	{
	  result = gtk_tree_model_filter_get_model(GTK_TREE_MODEL_FILTER(model));
	}
      else
	{
	  result = model;
	}
    }
  
  return result;
}


/* Get the path in the tree view's filtered model that corresponds to the given
 * path in the base (unfiltered) model */
static GtkTreePath *treeConvertBasePathToVisiblePath(GtkTreeView *tree, GtkTreePath *basePath)
{
  GtkTreePath *result = NULL;
  
  if (tree && basePath)
    {
      /* Convert the child path to the equivalent path in the filtered model. The
       * result may be null if the child row does not appear in the filtered model. */
      GtkTreeModel *filter = treeGetVisibleDataModel(tree);
      if (GTK_IS_TREE_MODEL_FILTER(filter))
	{
	  result = gtk_tree_model_filter_convert_child_path_to_path(GTK_TREE_MODEL_FILTER(filter), basePath);
	}
      else
	{
	  result = basePath;
	}
    }
  
  return result;
}
								    

/* Return true if the given path is selected in the given tree view. The given 
 * path must be in the given model, but this can be in either the base model
 * or the filtered model of the tree - we will check and do any conversion necessary. */
gboolean treePathIsSelected(GtkTreeView *tree, GtkTreePath *path, GtkTreeModel *model)
{
  gboolean result = FALSE;
  
  GtkTreeModel *baseModel = treeGetBaseDataModel(tree);
  GtkTreePath *visiblePath = (model == baseModel) ? treeConvertBasePathToVisiblePath(tree, path) : path;
  
  if (visiblePath)
    {
      GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
      result = gtk_tree_selection_path_is_selected(selection, visiblePath);
    }
  
  return result;
}


/* Refilter the data for the given tree */
void refilterTree(GtkWidget *tree, gpointer data)
{
  GtkTreeModelFilter *filter = GTK_TREE_MODEL_FILTER(gtk_tree_view_get_model(GTK_TREE_VIEW(tree)));
  gtk_tree_model_filter_refilter(filter);
}


/* Queue a redraw for the given tree, and the grid that corresponds to it */
void refreshTreeAndGrid(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  if (properties->grid)
    gtk_widget_queue_draw(properties->grid);
  
  gtk_widget_queue_draw(tree);
}


/* Set the font for the detail view cell renderer. Should be called after
 * the font size is changed by zooming in/out of the detail view. */
gint setCellRendererFont(GtkWidget *tree, GtkCellRenderer *renderer, const char *fontName, const int fontSize)
{
  gint newSize = fontSize * PANGO_SCALE;
  
  /* Update the widget */
  PangoFontDescription *font_desc = pango_font_description_copy (tree->style->font_desc);
  pango_font_description_set_size(font_desc, newSize);
  pango_font_description_set_family(font_desc, fontName);
  gtk_widget_modify_font(tree, font_desc);
  
  /* Calculate the row height from the font metrics */
  PangoContext *context = gtk_widget_get_pango_context (tree);
  PangoFontMetrics *metrics = pango_context_get_metrics (context,
							 font_desc,
							 pango_context_get_language (context));
  pango_font_description_free (font_desc);
  
  gint row_height = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  pango_font_metrics_unref (metrics);
  
  /* Subtract the gaps from top and bottom so that we render over them, giving the impression of no gaps. */
  gint vertical_separator;
  gtk_widget_style_get (tree, "vertical-separator", &vertical_separator, NULL);
  row_height -= vertical_separator * 2;
  
  /* Update the cell renderer */
  gtk_cell_renderer_set_fixed_size(renderer, 0, row_height);
  
  return row_height;
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


/* Get the character index at the given x coordinate in the given tree view column. */
static int getCharIndexAtTreeCoord(GtkWidget *tree, GtkTreeViewColumn *col, const int x)
{
  int result = UNSET_INT;
  
  /* Get the cell dimensions */
  int colWidth, startPos;
  GtkCellRenderer *renderer = treeGetRenderer(tree);
  gtk_tree_view_column_cell_get_position(col, renderer, &startPos, &colWidth);
  
  /* Allow for any padding */
  int xPad = renderer->xpad;
  gint horiz_separator;
  gtk_widget_style_get (tree, "horizontal-separator", &horiz_separator, NULL);
  
  
  /* Check it's in range */
  if (x > xPad + horiz_separator && x < colWidth - xPad)
    {
      int charWidth, charHeight;
      getFontCharSize(tree, &charWidth, &charHeight);
      
      result = ((x - xPad - horiz_separator - 1) / charWidth) + 1;
    }
  
  return result;
}


/* In the detail view, get the base index at the given coords */
static int getBaseIndexAtTreeCoords(GtkWidget *tree, const int x, const int y)
{
  int baseIdx = UNSET_INT;
  
  GtkTreePath *path = NULL;
  GtkTreeViewColumn *col = NULL;
  int cell_x, cell_y; /* coords relative to cell */
  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), x, y, &path, &col, &cell_x, &cell_y);
  
  /* See if we clicked in the sequence column */
  if (col == gtk_tree_view_get_column(GTK_TREE_VIEW(tree), MSP_COL))
    {
      /* Get the base index at the clicked position */
      int charIdx = getCharIndexAtTreeCoord(tree, col, cell_x);
      GtkAdjustment *adjustment = treeGetAdjustment(tree);
      baseIdx = charIdx + adjustment->value; /* Adjustment is 0-based, display range is 1-based */
    }
  
  return baseIdx;
}


/* Returns true if the given row in the given tree model should be visible. */
gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  gboolean bDisplay = FALSE;
  if (!bDisplay)
    {
      /* Find the MSP in this row */
      GValue val = {0};
      gtk_tree_model_get_value(model, iter, MSP_COL, &val);
      MSP *msp = g_value_get_pointer(&val);
      
      /* Show this row if any part of the MSP's range is inside the displayed range */
      GtkAdjustment *adjustment = treeGetAdjustment(tree);
      if (adjustment)
	{
	  int displayStart = adjustment->value;
	  int displayEnd = displayStart + adjustment->page_size;
	  bDisplay = !(msp->qstart > displayEnd || msp->qend < displayStart);
	}
    }
  
  return bDisplay;
}

/***********************************************************
 *                       Properties                        *
 ***********************************************************/

TreeProperties* treeGetProperties(GtkWidget *widget)
{
  return widget ? (TreeProperties*)(g_object_get_data(G_OBJECT(widget), "TreeProperties")) : NULL;
}

static void onDestroyTree(GtkWidget *widget)
{
  TreeProperties *properties = treeGetProperties(widget);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "TreeProperties", NULL);
    }
}

static void treeCreateProperties(GtkWidget *widget, GtkWidget *grid, GtkWidget *detailView, GtkCellRenderer *renderer)
{
  if (widget)
    { 
      TreeProperties *properties = g_malloc(sizeof *properties);
      
      properties->grid = grid;
      properties->detailView = detailView;
      properties->renderer = renderer;
      
      g_object_set_data(G_OBJECT(widget), "TreeProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyTree), NULL); 
    }
}


/***********************************************************
 *                       Tree events                       *
 ***********************************************************/

static gboolean onExposeDetailViewTree(GtkWidget *tree, GdkEventExpose *event, gpointer data)
{
  return FALSE;
}


static gboolean onButtonPressTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Left button: Let the default handler select the row */
  if (event->button == 1)
    {
      /* First, deselect anything in any other trees than this one. Then let the
       * default handler select the row(s) in this tree. */
      deselectAllSiblingTrees(tree);
      
      handled = FALSE;
    }
  
  /* Right button: show context menu (to do) */
  if (event->button == 3)
    {
      handled = TRUE;
    }
  
  /* Middle button: scrolls */
  if (event->button == 2)
    {
      treeSetSelectedBaseIdx(tree, getBaseIndexAtTreeCoords(tree, event->x, event->y));
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Left button: Let the default handler select the row */
  if (event->button == 1)
    {
      handled = FALSE;
    }

  /* Right button: show context menu (to do) */
  if (event->button == 3)
    {
      handled = TRUE;
    }
  
  /* Middle button: scrolls */
  if (event->button == 2)
    {
      /* Move the scrollbar so that the currently-selected base index is at the centre */
      GtkWidget *detailView = GTK_WIDGET(data);
      DetailViewProperties *properties = detailViewGetProperties(detailView);
      
      if (properties->selectedBaseIdx != UNSET_INT)
	{
	  GtkAdjustment *adjustment = treeGetAdjustment(tree);
	  if (adjustment)
	    {
	      int scrollStart = properties->selectedBaseIdx - (adjustment->page_size / 2);
	      setDetailViewScrollPos(adjustment, scrollStart);
	    }
	}
      
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onMouseMoveTree(GtkWidget *tree, GdkEventMotion *event, gpointer data)
{
  /* Middle mouse button down */
  if (event->state == GDK_BUTTON2_MASK)
    {
      treeSetSelectedBaseIdx(tree, getBaseIndexAtTreeCoords(tree, event->x, event->y));
    }
  
  return TRUE;
}


static void onDragBeginTree(GtkWidget *widget, GdkDragContext *event, gpointer data)
{
}


static void onDragEndTree(GtkWidget *widget, GdkDragContext *event, gpointer data)
{
}


static gboolean onDragMotionTree(GtkWidget *widget, GdkDragContext *event, gint x, gint y, guint time, gpointer data)
{
  return FALSE;
}


static gboolean onEnterTree(GtkWidget *tree, GdkEventCrossing *event, gpointer data)
{
  /* Return true to stop the default handler re-drawing when the focus changes */
  return TRUE;
}


static gboolean onLeaveTree(GtkWidget *tree, GdkEventCrossing *event, gpointer data)
{
  /* Return true to stop the default handler re-drawing when the focus changes */
  return TRUE;
}


static void onSelectionChangedTree(GtkWidget *widget, gpointer data)
{
  /* Redraw the corresponding grid */
  GtkWidget *tree = GTK_WIDGET(data);
  TreeProperties *properties = treeGetProperties(tree);
  
  if (properties->grid)
    gtk_widget_queue_draw(properties->grid);
}


/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

static GtkTreeViewColumn* initColumn(GtkWidget *tree, GtkCellRenderer *renderer, char *colName, 
				     char *rendererProperty, int colNum, gboolean expand, const int width)
{
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(colName, renderer, rendererProperty, colNum, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
  
  //  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
  //  gtk_tree_view_column_set_min_width(column, 0);
  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_fixed_width(column, width);
  gtk_tree_view_column_set_resizable(column, TRUE);
  
  /* If expand is true, this column gobbles up any extra width in the widget
   * when the window or any other columns are resized. */
  if (expand)
    gtk_tree_view_column_set_expand(column, TRUE);
  
  return column;
}


static MSP* createRefSeqMsp(GtkWidget *tree, gboolean fwd)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  MSP *msp = g_malloc(sizeof(MSP));

  msp->next = NULL;
  msp->type = BLX_MSP_INVALID;
  msp->score = -1;
  msp->id = -1;
  msp->qname = "Reference";
  msp->qframe[0] = '(';
  msp->qframe[1] = fwd ? '+' : '-';
  msp->qframe[2] = '0';
  msp->qframe[3] = ')';
  msp->qstart = detailViewProperties->displayRange.min;
  msp->qend = detailViewProperties->displayRange.max;
  msp->sname = "Reference";
  msp->sframe[0] = '(';
  msp->sframe[1] = fwd ? '+' : '-';
  msp->sframe[2] = '0';
  msp->sframe[3] = ')';
  msp->slength = detailViewProperties->displayRange.max - detailViewProperties->displayRange.min;
  msp->sstart = detailViewProperties->displayRange.min;
  msp->send = detailViewProperties->displayRange.max;
  msp->sseq = detailViewProperties->refSeq;
  
  return msp;
}


static void addMspToTreeStore(GtkListStore *store, const MSP *msp)
{
  GtkTreeIter iter;
  gtk_list_store_append(store, &iter);
  gtk_list_store_set(store, &iter,
		     S_NAME_COL, msp->sname,
		     SCORE_COL, msp->score,
		     ID_COL, msp->id,
		     START_COL, msp->sstart,
		     MSP_COL, msp,
		     END_COL, msp->send,
		     -1);
}


/* Add the MSPs in the given list to the given list store */
static void addMspsToTreeStore(GtkListStore *store, const MSP const *mspList)
{
  const MSP *msp = mspList;
  for ( ; msp; msp = msp->next)
    {
      addMspToTreeStore(store, msp);
    }  
}


static void addTreeColumns(GtkWidget *tree, GtkCellRenderer *renderer)
{
  /* Add the columns */
  initColumn(tree,  renderer,  "Name",     "text", S_NAME_COL, FALSE, 90);
  initColumn(tree,  renderer,  "Score",    "text", SCORE_COL,  FALSE, 30);
  initColumn(tree,  renderer,  "%ID",      "text", ID_COL,	  FALSE, 30);
  initColumn(tree,  renderer,  "Start",    "text", START_COL,  FALSE, 40);
  GtkTreeViewColumn *seqCol = initColumn(tree,  renderer,  "Sequence", "msp",  MSP_COL,  TRUE, 40);
  initColumn(tree,  renderer,  "End",      "text", END_COL,	  FALSE, 30);
  
  /* Set it to use a fixed width font. This also sets the row height based on the font size.
   * For now we just use the default size from the theme. */
  int origSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE);
  setCellRendererFont(tree, renderer, FIXED_WIDTH_FONT, origSize);
  
  seqColCreateProperties(seqCol, tree);
}


/* Create the data store for a detail view tree */
static void createTreeDataStore(GtkWidget *tree,
				const MSP const *mspList, 
				GtkCellRenderer *renderer)
{
  /* Create the data store for the tree view */
  GtkListStore *store = gtk_list_store_new(N_COLUMNS,
					   G_TYPE_STRING, 
					   G_TYPE_INT, 
					   G_TYPE_INT, 
					   G_TYPE_INT, 
					   G_TYPE_POINTER, 
					   G_TYPE_INT);
  
  /* Add a 'fake' msp for the ref sequence */
  addMspToTreeStore(store, createRefSeqMsp(tree, TRUE));

  /* Add our data to the store */
  addMspsToTreeStore(store, mspList);
  
  /* Create a filtered version of the data store to only show rows that are within the display range. */
  GtkTreeModel *filter = gtk_tree_model_filter_new(GTK_TREE_MODEL(store), NULL);
  g_object_unref (G_OBJECT (store));
  
  gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(filter), 
					 (GtkTreeModelFilterVisibleFunc)isTreeRowVisible, 
					 tree, 
					 NULL);
  
  /* Add the filtered store to the tree view */
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), GTK_TREE_MODEL(filter));
  g_object_unref (G_OBJECT (filter));
}


/* Create a widget in the given grid for every MSP in the given tree view's
 * data store. */
//static void createMspLineWidgets(GtkWidget *grid, GtkWidget *treeView)
//{
//  /* Get the data model for the detail view, which contains info about the MSPs.
//   * The model stored in the tree is filtered - we want to get the original model
//   * so that we know it contains all of the MSPs. */
//  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(treeView));
//  
//  /* Create a corresponding widget in the big picture for each row in the tree view */
//  gtk_tree_model_foreach(model, createMspLineWidget, grid);
//}


static void setTreeStyle(GtkTreeView *tree, gboolean hasHeaders)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_VERTICAL);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_MULTIPLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, hasHeaders);

  /* Set the expander size to 0 so that we can have tiny rows (otherwise the min is 12pt) */
  char parseString[500];
  sprintf(parseString, "style \"packedTree\"\n"
	  "{\n"
	  "GtkTreeView::expander-size = 0\n"
	  "GtkTreeView::vertical-separator = 2\n"
	  "}"
	  "widget \"*%s*\" style \"packedTree\"", DETAIL_VIEW_TREE_NAME);
  gtk_rc_parse_string(parseString);
}


GtkWidget* createDetailViewTree(GtkWidget *grid, 
				GtkWidget *detailView, 
				gboolean hasHeaders, 
				const MSP const *mspList)
{
  /* Create a tree view for the list of match sequences */
  GtkWidget *tree = gtk_tree_view_new();
  setTreeStyle(GTK_TREE_VIEW(tree), hasHeaders);
  
  /* Put it in a scrolled window for vertical scrolling only (hoz scrolling will be via our
   * custom adjustment). Always display the scrollbars because we assume the column widths 
   * are the same for all trees and they won't be if one shows a scrollbar and another doesn't. */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
  gtk_container_add(GTK_CONTAINER(scrollWin), tree);
  
  /* Create a custom cell renderer to render the match sequences in this tree */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();

  /* The tree needs to know which grid and renderer it corresponds to, and vice versa */
  treeCreateProperties(tree, grid, detailView, renderer);
  gridGetProperties(grid)->tree = tree;
  SEQUENCE_CELL_RENDERER(renderer)->tree = tree;
  
  /* The detailView needs one of the trees as a reference. Give it the first one we make. */
  static int firstTree = TRUE;
  if (firstTree)
    {
      DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
      detailViewProperties->refTree = tree;
      firstTree = FALSE;
    }
  
  /* Add the columns and data, and create widgets on the grid to represent the msps */
  addTreeColumns(tree, renderer);
  createTreeDataStore(tree, mspList, renderer);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);
  g_signal_connect(G_OBJECT(tree), "button-press-event", G_CALLBACK(onButtonPressTree), detailView);
  g_signal_connect(G_OBJECT(tree),  "button-release-event", G_CALLBACK(onButtonReleaseTree), detailView);
  g_signal_connect(G_OBJECT(tree), "motion-notify-event", G_CALLBACK(onMouseMoveTree), detailView);
  g_signal_connect(G_OBJECT(tree), "enter-notify-event", G_CALLBACK(onEnterTree), NULL);
  g_signal_connect(G_OBJECT(tree), "leave-notify-event", G_CALLBACK(onLeaveTree), NULL);
  g_signal_connect(G_OBJECT(tree), "drag-begin", G_CALLBACK(onDragBeginTree), NULL);
  g_signal_connect(G_OBJECT(tree), "drag-end", G_CALLBACK(onDragEndTree), NULL);
  g_signal_connect(G_OBJECT(tree), "drag-motion", G_CALLBACK(onDragMotionTree), NULL);
  g_signal_connect(G_OBJECT(tree), "expose-event", G_CALLBACK(onExposeDetailViewTree), NULL);
  
  GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
  g_signal_connect(G_OBJECT(selection), "changed", G_CALLBACK(onSelectionChangedTree), tree);
  
  return scrollWin;
}