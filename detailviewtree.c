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
#include <SeqTools/blixem_.h>
#include <string.h>

#define DETAIL_VIEW_TREE_NAME		"DetailViewTreeName"
#define NO_SUBJECT_SELECTED_TEXT	"<no subject selected>"
#define FONT_INCREMENT_SIZE		1
#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20


/* Local function declarations */
static void		updateFeedbackBoxForAllTrees(GtkWidget *tree);
static void		onSelectionChangedTree(GObject *selection, gpointer data);
static gint		setCellRendererFont(GtkWidget *tree, GtkCellRenderer *renderer, const int fontSize, const char *fontFamily);
static GtkTreePath*	treeConvertBasePathToVisiblePath(GtkTreeView *tree, GtkTreePath *basePath);


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

GdkColor getGdkColor(gulong colour)
{
  GdkColor result = {colour};
  return result;
}


/* Returns true if the given MSP is a "fake" MSP (i.e. one that that is used just
 * to display the reference sequence and does not represent a real match.) */
gboolean mspIsFake(const MSP const *msp)
{
  /* This is not a great way to identify fake msp's but it'll do for now. */
  return (msp && msp->score == 0);
}

gboolean mspIsExon(const MSP const *msp)
{
  return (msp && msp->score == -1);
}

gboolean mspIsIntron(const MSP const *msp)
{
  return (msp && msp->score == -2);
}

gboolean mspIsBlastMatch(const MSP const *msp)
{
  return (msp && msp->score > 0);
}


static void assertTree(GtkWidget *tree)
{
  if (!tree)
    messcrash("Tree is null", tree);
  
  if (!GTK_IS_WIDGET(tree))
    messcrash("Tree is not a valid widget [%x]", tree);
    
  if (!GTK_IS_TREE_VIEW(tree))
    messcrash("Tree is not a valid tree view [%x]", tree);

  if (strcmp(gtk_widget_get_name(tree), DETAIL_VIEW_TREE_NAME) != 0)
    messcrash("Tree is not a valid detail-view tree [%x]", tree);

  if (!treeGetProperties(tree))
    messcrash("Tree properties not set [widget=%x]", tree);
}

GtkAdjustment *treeGetAdjustment(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(properties->detailView);
  return detailViewProperties ? detailViewProperties->adjustment : NULL;
}

GtkWidget *treeGetGrid(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->grid : NULL;
}

GtkCellRenderer *treeGetRenderer(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->renderer : NULL;
}

GtkWidget *treeGetDetailView(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->detailView : NULL;
}


Strand treeGetStrand(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return gridGetStrand(properties->grid);
}


gboolean treeGetStrandsToggled(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetStrandsToggled(detailView);
}

BlxBlastMode treeGetBlastMode(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetBlastMode(detailView);
}

int treeGetNumReadingFrames(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetNumReadingFrames(properties->detailView) : UNSET_INT;
}

int treeGetSelectedBaseIdx(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetSelectedBaseIdx(properties->detailView) : UNSET_INT;
}

static void treeSetSelectedBaseIdx(GtkWidget *tree, const int selectedBaseIdx)
{
  assertTree(tree);

  TreeProperties *properties = treeGetProperties(tree);
  if (properties)
    {
      DetailViewProperties *detailViewProperties = detailViewGetProperties(properties->detailView);
      
      /* Only update if things have changed */
      if (detailViewProperties && detailViewProperties->selectedBaseIdx != selectedBaseIdx)
	{
	  detailViewProperties->selectedBaseIdx = selectedBaseIdx;
	  
	  /* Update the feedback box */
	  updateFeedbackBoxForAllTrees(tree);

	  /* Re-render all trees in the detail view */
	  gtk_widget_queue_draw(properties->detailView);
	}
    }
}

char* treeGetRefSeq(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetRefSeq(detailView);
}

IntRange* treeGetDisplayRange(GtkWidget *tree)
{
  assertTree(tree);
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


/* Selects the given row iterator in the given tree. The iter must be in the given
 * model, but the model can be either the base or filtered model for this tree
 * (this function converts accordingly) */
void selectRow(GtkTreeView *tree, GtkTreeModel *model, GtkTreeIter *iter)
{
  GtkTreeSelection *selection = gtk_tree_view_get_selection(tree);
  GtkTreePath *path = gtk_tree_model_get_path(model, iter);

  /* Convert to the visible (filtered) model if necessary */
  GtkTreeModel *visibleModel = treeGetVisibleDataModel(tree);
  GtkTreePath *visiblePath = (model == visibleModel) ? path : treeConvertBasePathToVisiblePath(tree, path);
  
  if (visiblePath)
    gtk_tree_selection_select_path(selection, visiblePath);
}


/* Deselect all rows in the given tree */
void deselectAllRows(GtkWidget *tree)
{
  assertTree(tree);
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
 * as the given tree (but NOT on the given tree itself unless includeCurrent is true). */
void deselectAllSiblingTrees(GtkWidget *tree, gboolean includeCurrent)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  
  if (includeCurrent)
    {
      /* Just call the recursive function */
      callFuncOnAllDetailViewTrees(detailView, deselectAllRows);
    }
  else
    {
      /* Need to pass extra data (i.e. the tree to omit), so do this call manually */
      gtk_container_foreach(GTK_CONTAINER(detailView), deselectAllRowsNotInTree, tree);
    }
}


/* Return the msp in a given tree row */
const MSP* treeGetMsp(GtkTreeModel *model, GtkTreeIter *iter)
{
  MSP *msp = NULL;
  gtk_tree_model_get(model, iter, MSP_COL, &msp, -1);
  return msp;
}


/* Return the msp in a given tree row */
static GtkWidget* treeGetFeedbackBox(GtkWidget *tree)
{
  assertTree(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(treeGetDetailView(tree));
  return detailViewProperties->feedbackBox;
}

GdkColor* treeGetRefSeqColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->refSeqColour;
}

GdkColor* treeGetRefSeqSelectedColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->refSeqSelectedColour;
}

GdkColor* treeGetMatchColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->matchColour;
}

GdkColor* treeGetMatchSelectedColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->matchSelectedColour;
}

GdkColor* treeGetMismatchColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->mismatchColour;
}

GdkColor* treeGetMismatchSelectedColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->mismatchSelectedColour;
}

GdkColor* treeGetExonColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->exonColour;
}

GdkColor* treeGetExonSelectedColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->exonSelectedColour;
}

GdkColor* treeGetGapColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->gapColour;
}

GdkColor* treeGetGapSelectedColour(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return &properties->gapSelectedColour;
}

const char* treeGetFontFamily(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetFontFamily(detailView);
}


/* For the given tree view, return the data model used for the visible view
 * (which may be a filtered tree model). */
GtkTreeModel* treeGetVisibleDataModel(GtkTreeView *tree)
{
  /* Just return the model stored directly in the tree view. This will be the filter
   * if there is one, otherwise it will be the original data model. */
  assertTree(GTK_WIDGET(tree));
  return gtk_tree_view_get_model(tree);
}

/* For the given tree view, return the base data model i.e. all data in the tree
 * without any filtering. Returns the same as treeGetVisibleDataModel if the tree does
 * not have a filter. */
GtkTreeModel* treeGetBaseDataModel(GtkTreeView *tree)
{
  assertTree(GTK_WIDGET(tree));

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


/* Get the path in the tree view's filtered (visible) model that corresponds to the given
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
	

///* Get the path in the tree view's base (unfiltered) model that corresponds to the given
// * path in the visible (filtered) model */
//static GtkTreePath *treeConvertVisiblePathToBasePath(GtkTreeView *tree, GtkTreePath *visiblePath)
//{
//  GtkTreePath *result = NULL;
//  
//  if (tree && visiblePath)
//    {
//      /* Convert the path to the equivalent path in the unfiltered (child) model. */
//      GtkTreeModel *filter = treeGetVisibleDataModel(tree);
//      if (GTK_IS_TREE_MODEL_FILTER(filter))
//	{
//	  result = gtk_tree_model_filter_convert_path_to_child_path(GTK_TREE_MODEL_FILTER(filter), visiblePath);
//	}
//      else
//	{
//	  result = visiblePath;
//	}
//    }
//  
//  return result;
//}


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
  assertTree(tree);
  GtkTreeModelFilter *filter = GTK_TREE_MODEL_FILTER(gtk_tree_view_get_model(GTK_TREE_VIEW(tree)));
  gtk_tree_model_filter_refilter(filter);
}


/* Queue a redraw for the given tree, and the grid that corresponds to it */
void refreshTreeAndGrid(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  if (properties->grid)
    gtk_widget_queue_draw(properties->grid);
  
//  gtk_widget_queue_draw(tree); /* re-rendering seems to happen before this... */
}


/* Increase the font size in the detail view trees (i.e. effectively zoom in) */
void treeIncrementFontSize(GtkWidget *tree, gpointer data)
{
  int newSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE) + FONT_INCREMENT_SIZE;
  
  if (newSize <= MAX_FONT_SIZE)
    {
      GtkCellRenderer *renderer = treeGetRenderer(tree);
      setCellRendererFont(tree, renderer, newSize, treeGetFontFamily(tree));
    }
}


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
void treeDecrementFontSize(GtkWidget *tree, gpointer data)
{
  int newSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE) - FONT_INCREMENT_SIZE;
  
  if (newSize >= MIN_FONT_SIZE)
    {
      GtkCellRenderer *renderer = treeGetRenderer(tree);
      setCellRendererFont(tree, renderer, newSize, treeGetFontFamily(tree));
    }
}


/* Set the font for the detail view cell renderer. Should be called after
 * the font size is changed by zooming in/out of the detail view. */
static gint setCellRendererFont(GtkWidget *tree, GtkCellRenderer *renderer, const int fontSize, const char *fontFamily)
{
  /* Update the widget with the new font size. */
  PangoFontDescription *fontDesc = pango_font_description_copy(tree->style->font_desc);
  
  if (fontFamily)
    {
      pango_font_description_set_family(fontDesc, fontFamily);
    }
  
  pango_font_description_set_size  (fontDesc, fontSize * PANGO_SCALE);
  gtk_widget_modify_font(tree, fontDesc);
  
  /* Calculate the row height from the font metrics */
  PangoContext *context = gtk_widget_get_pango_context (tree);
  PangoFontMetrics *metrics = pango_context_get_metrics (context, fontDesc, pango_context_get_language(context));
  
  gint charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  gint charWidth = pango_font_metrics_get_approximate_char_width(metrics) / PANGO_SCALE;

  pango_font_description_free(fontDesc);
  pango_font_metrics_unref(metrics);

  /* Cache these results in the cell renderer */
  SEQUENCE_CELL_RENDERER(renderer)->charHeight = charHeight;
  SEQUENCE_CELL_RENDERER(renderer)->charWidth = charWidth;
  
  /* Set the row height. Subtract the size of the vertical separators . This makes the
   * row smaller than it is, but we will render it at normal size, so we end up drawing over
   * the gaps, hence giving the impression that there are no gaps. */
  gint vertical_separator;
  gtk_widget_style_get (tree, "vertical-separator", &vertical_separator, NULL);
  gint rowHeight = charHeight - vertical_separator * 2;
  gtk_cell_renderer_set_fixed_size(renderer, 0, rowHeight);
  
  return rowHeight;
}


/* Get the character index at the given x coordinate in the given tree view column. */
static int getCharIndexAtTreeCoord(GtkWidget *tree, GtkTreeViewColumn *col, const int x)
{
  int result = UNSET_INT;
  
  /* Get the cell dimensions */
  int colWidth, startPos;
  GtkCellRenderer *renderer = treeGetRenderer(tree);
  gtk_tree_view_column_cell_get_position(col, renderer, &startPos, &colWidth);
  
  /* Find the direction of the display. If strands are toggled, the reference
   * sequence is shown right-to-left (i.e. lowest value on the right). */
  gboolean rightToLeft = treeGetStrandsToggled(tree);
  
  /* Check that the given coord is within the cell's display area */
  if (x > renderer->xpad && x < colWidth)
    {
      gint charWidth = SEQUENCE_CELL_RENDERER(renderer)->charWidth;

      /* Calculate the distance from the left edge (or right edge, if display reversed). */
      int distFromEdge = rightToLeft ? colWidth - (x - renderer->xpad * 2) : x - renderer->xpad * 2;
      result = (int)((double)(distFromEdge) / (double)charWidth);
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
      baseIdx = charIdx + adjustment->value + 1; /* Adjustment is 0-based, display range is 1-based */
    }
  
  return baseIdx;
}


/* Get the text displayed in the user feedback box based on the given row's sequence name (if a
 * row iterater and model are given), and also the currently-selected base index (if there is one). 
 * The string returned by this function must be free'd with g_free. */
static char* getFeedbackText(GtkWidget *tree, GtkTreeModel *model, GtkTreeIter *iter)
{
  /* The info we need to find... */
  int qIdx = treeGetSelectedBaseIdx(tree);
  int sIdx = UNSET_INT;
  char *sname = NULL;

  /* Make sure we have enough space for all the bits to go in the string */
  int msgLen = numDigitsInInt(qIdx);
  
  /* If a row is given, set the sequence name to that row's sequence */
  if (iter != NULL && model != NULL)
    {
      const MSP *msp = treeGetMsp(model, iter);
      if (msp)
	{
	  sname = msp->sname;
	  msgLen += strlen(sname);
      
	  /* If a ref sequence base is selected, find the corresponding base in the match sequence */
	  sIdx = UNSET_INT;
	  if (msp && qIdx != UNSET_INT)
	    {
	      /* If the base doesn't exist in the match sequence, just get the "nearest" value */
	      gapCoord(msp, 
		       qIdx, 
		       treeGetNumReadingFrames(tree), 
		       treeGetStrand(tree), 
		       treeGetStrandsToggled(tree), 
		       &sIdx);
	      
	      msgLen += numDigitsInInt(sIdx);
	    }
	}
    }

  if (!sname)
    {
      msgLen += strlen(NO_SUBJECT_SELECTED_TEXT);
    }

  msgLen += 10; /* for format text */
  char *messageText = g_malloc(sizeof(char) * msgLen);
  
  /* Create the message text. */
  if (sname != NULL && qIdx != UNSET_INT)
    {
      sprintf(messageText, "%d   %s: %d", qIdx, sname, sIdx);		/* display all */
    }
  else if (sname != NULL)
    {
      sprintf(messageText, "%s", sname);				/* just the name */
    }
  else if (qIdx != UNSET_INT)
    {
      sprintf(messageText, "%d   %s", qIdx, NO_SUBJECT_SELECTED_TEXT);  /* just the index */
    }
  else
    {
      sprintf(messageText, " ");
    }
  
  return messageText;
}


/* Updates the feedback box with info about any currently-selected row/base in the
 * given tree.  This currently assumes single-row selections, but could be extended
 * in the future to display, say, summary information about multiple rows. Returns true
 * if the given tree has rows selected. */
static gboolean updateFeedbackBoxForTree(GtkWidget *tree)
{
  gboolean done = FALSE;
  char *messageText = NULL;
  
  GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
  GtkTreeModel *model = NULL;
  GList *selectedRows = gtk_tree_selection_get_selected_rows(selection, &model);

  if (g_list_length(selectedRows) > 0)
    {
      done = TRUE;

      /* We currently only handle single-selection, so only process the first selected row */
      GtkTreePath *path = (GtkTreePath*)(selectedRows->data);

      GtkTreeIter iter;
      gtk_tree_model_get_iter(model, &iter, path);
      
      messageText = getFeedbackText(tree, model, &iter);
    }
  else
    {
      /* No rows selected. Just see if a base index is selected. */
      messageText = getFeedbackText(tree, NULL, NULL);
    }

  GtkWidget *feedbackBox = treeGetFeedbackBox(tree);
  gtk_entry_set_text(GTK_ENTRY(feedbackBox), messageText);
  gtk_widget_queue_draw(feedbackBox);

  return done;
}


/* Updates the feedback box with info about the currently-selected row/base. Takes
 * into account selections in any of the trees. This should be called every time 
 * the row selection or base selection is changed. 
 * TNote that there shouldn't be rows selected in different trees at the same time,
 * so this function doesn't really deal with that situation - what will happen is
 * that the contents of the feedback box will be overwritten by the last tree the
 * update function is called for. */
static void updateFeedbackBoxForAllTrees(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);


  /* Loop through all of the trees. Stop if we find one that has selected rows. */
  int numFrames = detailViewGetNumReadingFrames(detailView);
  gboolean done = FALSE;

  int frame = 1;
  for ( ; frame <= numFrames && !done; ++frame)
    {
      /* Forward strand tree for this reading frame */
      GtkWidget *curTree = detailViewGetFrameTree(detailView, TRUE, frame);
      done = updateFeedbackBoxForTree(curTree);

      if (!done)
	{
	  /* Reverse strand tree for this reading frame */
	  GtkWidget *curTree = detailViewGetFrameTree(detailView, FALSE, frame);
	  done = updateFeedbackBoxForTree(curTree);
	}
    }
}


/* Returns true if the given row in the given tree model should be visible. */
gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  gboolean bDisplay = FALSE;
  if (!bDisplay)
    {
      /* Find the MSP in this row */
      const MSP *msp = NULL;
      gtk_tree_model_get(model, iter, MSP_COL, &msp, -1);
      
      if (msp)
	{
	  /* Show this row if any part of the MSP's range is inside the displayed range */
	  GtkAdjustment *adjustment = treeGetAdjustment(tree);
	  if (adjustment)
	    {
	      int displayStart = adjustment->value + 1;
	      int displayEnd = displayStart + adjustment->page_size;
	      
	      int qSeqMin, qSeqMax;
	      getMspRangeExtents(msp, &qSeqMin, &qSeqMax, NULL, NULL);
	      
	      bDisplay = !(qSeqMin > displayEnd || qSeqMax < displayStart);
	    }
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


static void treeCreateProperties(GtkWidget *widget, 
				 GtkWidget *grid, 
				 GtkWidget *detailView, 
				 GtkCellRenderer *renderer)
{
  if (widget)
    { 
      TreeProperties *properties = g_malloc(sizeof *properties);
      
      properties->grid = grid;
      properties->detailView = detailView;
      properties->renderer = renderer;
      
      properties->refSeqColour = getGdkColor(GDK_YELLOW);
      properties->refSeqSelectedColour = getGdkColor(GDK_DARK_YELLOW);
      properties->matchColour = getGdkColor(GDK_CYAN);
      properties->matchSelectedColour = getGdkColor(GDK_DARK_CYAN);
      properties->mismatchColour = getGdkColor(GDK_GREY);
      properties->mismatchSelectedColour = getGdkColor(GDK_DARK_GREY);
      properties->exonColour = getGdkColor(GDK_YELLOW);
      properties->exonSelectedColour = getGdkColor(GDK_DARK_YELLOW);
      properties->gapColour = getGdkColor(GDK_GREY);
      properties->gapSelectedColour = getGdkColor(GDK_DARK_GREY);
      
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
  
  switch (event->button)
    {
    
    case 1:
      {
	/* Left button: select row */
	/* First, deselect anything in any other trees than this one. Then let the
	 * default handler select the row. */
	deselectAllSiblingTrees(tree, FALSE);
	
	/* Let the default handler select the row */
	handled = FALSE;
	break;
      }

    case 2:
      {
	/* Middle button: scroll to centre on clicked base */
	/* Select the base index at the clicked coords */
	treeSetSelectedBaseIdx(tree, getBaseIndexAtTreeCoords(tree, event->x, event->y));
	handled = TRUE;
	break;
      }
      
    case 3:
      {
	/* Right button: show context menu (to do) */
	handled = TRUE;
	break;
      }
      
    default:
      break;
    };
  
  return handled;
}


static gboolean onButtonReleaseTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Left button: let default handler select the row */
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


/* Implement custom scrolling for horizontal mouse wheel movements over the tree.
 * This scrolls our custom horizontal scrollbar for the whole Detail View. Leaves
 * vertical scrolling to the default handler. */
static gboolean onScrollTree(GtkWidget *tree, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
    {
      case GDK_SCROLL_LEFT:
	{
	  scrollDetailViewLeftStep(treeGetDetailView(tree));
	  handled = TRUE;
	  break;
	}
	
      case GDK_SCROLL_RIGHT:
	{
	  scrollDetailViewRightStep(treeGetDetailView(tree));
	  handled = TRUE;
	  break;
	}

      default:
	{
	  /* Default handler can handle vertical scrolling */
	  handled = FALSE;
	  break;
	}
    };
  
  return handled;
}


static gboolean onMouseMoveTree(GtkWidget *tree, GdkEventMotion *event, gpointer data)
{
  if (event->state == GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base
       * (but don't re-centre on the selected base until the mouse button is released). */
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


static void onSelectionChangedTree(GObject *selection, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);

  /* Update the feedback box to tell the user which sequence is selected. */
  updateFeedbackBoxForAllTrees(tree);
  
  /* Redraw the corresponding grid */
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


/* calcID: caculated percent identity of an MSP
 * 
 * There seems to be a general problem with this routine for protein
 * alignments, the existing code certainly does not do the right thing.
 * I have fixed this routine for gapped sequence alignments but not for
 * protein stuff at all.
 * 
 * To be honest I think this routine is a _waste_ of time, the alignment
 * programs that feed data to blixem produce an identity anyway so why
 * not use that...why reinvent the wheel......
 * 
 * */
static void calcID(MSP *msp, GtkWidget *tree)
{
  static int id, i, len ;
  char *qseq ;
  
  BOOL sForward = strchr(msp->sframe, '+') ? TRUE : FALSE;
  BlxBlastMode blastMode = treeGetBlastMode(tree);

  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  if (msp->sseq /* to do: is this required? && msp->sseq != padseq */)
    {
      /* Note that getqseq() will reverse complement if qstart > qend, this means that
       * where there is no gaps array the comparison is trivial as coordinates can be
       * ignored and the two sequences just whipped through. */
      char *refSeq = treeGetRefSeq(tree);
      
      if (!(qseq = getqseq(msp->qstart, msp->qend, refSeq)))
	{
	  messout ( "calcID failed: Don't have genomic sequence %d - %d\n",
		   msp->qstart, msp->qend);
	  msp->id = 0;
	  
	  return;
	}
      
      
      /* NOTE, tblastn/x are not implemented for gaps yet. */
      if (!(msp->gaps) || arrayMax(msp->gaps) == 0)
	{
	  /* Ungapped alignments. */
	  
	  if (blastMode == BLXMODE_TBLASTN)
	    {
	      len = msp->qend - msp->qstart + 1;
	      
	      for (i=0, id=0; i < len; i++)
		if (freeupper(msp->sseq[i]) == freeupper(qseq[i]))
		  id++;
	    }
	  else if (blastMode == BLXMODE_TBLASTX)
	    {
	      len = abs(msp->qend - msp->qstart + 1)/3;
	      
	      for (i=0, id=0; i < len; i++)
		if (freeupper(msp->sseq[i]) == freeupper(qseq[i]))
		  id++;
	    }
	  else						    /* blastn, blastp & blastx */
	    {
	      len = sSeqMax - sSeqMin + 1 ;
	      
	      for (i=0, id=0; i< len; i++)
                {
                  int sIndex = sForward ? sSeqMin + i - 1 : sSeqMax - i - 1;
		  if (freeupper(msp->sseq[sIndex]) == freeupper(qseq[i]))
		    id++;
                }
	    }
	}
      else
	{
	  /* Gapped alignments. */
	  
	  /* To do tblastn and tblastx is not imposssible but would like to work from
	   * examples to get it right.... */
	  if (blastMode == BLXMODE_TBLASTN)
	    {
	      printf("not implemented yet\n") ;
	    }
	  else if (blastMode == BLXMODE_TBLASTX)
	    {
	      printf("not implemented yet\n") ;
	    }
	  else
	    {
	      /* blastn and blastp remain simple but blastx is more complex since the query
               * coords are nucleic not protein. */
	      
              BOOL qForward = (strchr(msp->qframe, '+')) ? TRUE : FALSE ;
              Array gaps = msp->gaps;
	      int numFrames = treeGetNumReadingFrames(tree);
	      
	      len = 0 ;
              int i ;
	      for (i = 0, id = 0 ; i < arrayMax(gaps) ; i++)
		{
		  SMapMap *m = arrp(gaps, i, SMapMap) ;
		  
                  int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
		  getSMapMapRangeExtents(m, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
		  
                  len += sRangeMax - sRangeMin + 1;
                  
                  /* Note that qseq has been cut down to just the section relating to this msp.
		   * We need to translate the first coord in the range (which is in terms of the full
		   * reference sequence) into coords in the cut-down ref sequence. */
                  int q_start = qForward ? (qRangeMin - qSeqMin) / numFrames : (qSeqMax - qRangeMax) / numFrames;
		  
		  /* We can index sseq directly (but we need to adjust by 1 for zero-indexing). We'll loop forwards
		   * through sseq if we have the forward strand or backwards if we have the reverse strand,
		   * so start from the lower or upper end accordingly. */
                  int s_start = sForward ? sRangeMin - 1 : sRangeMax - 1 ;
		  
                  int j = s_start, k = q_start ;
		  while ((sForward && j < sRangeMax) || (!sForward && j >= sRangeMin - 1))
		    {
		      
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
		      char *sub, *query ;
		      char *dummy ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
		      
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
		      /* for debug..... */
		      sub = msp->sseq + j ;
		      query = qseq + k ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
		      
		      
		      if (freeupper(msp->sseq[j]) == freeupper(qseq[k]))
			id++ ;
		      
                      
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
		      else
			dummy = sub ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
                      
		      
                      /* Move to the next base. The qseq is always forward, but we might have to
                       * traverse the s sequence in reverse. */
                      ++k ;
                      if (sForward) ++j ;
                      else --j ;
		    }
		}
	    }
	}
      
      msp->id = (int)((float)100*id/len + .5);
      
      messfree(qseq);
    }
  else
    msp->id = 0 ;
  
  return ;
}


static MSP* createRefSeqMsp(GtkWidget *tree, gboolean fwd, char *refSeq, IntRange *displayRange)
{
  MSP *msp = g_malloc(sizeof(MSP));

  msp->next = NULL;
  msp->type = BLX_MSP_INVALID;
  msp->score = 0;
  msp->id = -1;
  msp->qname = REFERENCE_SEQUENCE_NAME;
  msp->qframe[0] = '(';
  msp->qframe[1] = fwd ? '+' : '-';
  msp->qframe[2] = '0';
  msp->qframe[3] = ')';
  msp->qstart = displayRange->min;
  msp->qend = displayRange->max;
  msp->sname = REFERENCE_SEQUENCE_NAME;
  msp->sframe[0] = '(';
  msp->sframe[1] = fwd ? '+' : '-';
  msp->sframe[2] = '0';
  msp->sframe[3] = ')';
  msp->slength = displayRange->max - displayRange->min;
  msp->sstart = displayRange->min;
  msp->send = displayRange->max;
  msp->sseq = refSeq;

  return msp;
}


void addMspToTreeModel(GtkTreeModel *model, MSP *msp, GtkWidget *tree)
{
  /* We don't display introns in the detail view */
  if (!mspIsIntron(msp))
    {
      /* Calculate the id */
      calcID(msp, tree);
      
      /* Add it to the tree's data store */
      GtkListStore *store = GTK_LIST_STORE(model);
      
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
}


/* Cell data function for the "start" column. This displays the start coord of the match 
 * sequence in normal left-to-right display, but the end coord if the display is reversed */
static void cellDataFunctionStartCol(GtkTreeViewColumn *column,
				     GtkCellRenderer *renderer, 
				     GtkTreeModel *model, 
				     GtkTreeIter *iter, 
				     gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  gboolean rightToLeft = treeGetStrandsToggled(tree);
  
  /* Get the msp for this row */
  const MSP *msp = treeGetMsp(model, iter);
  
  /* We want to display the start coord, unless the display is reversed, in which case display the end */
  int coord = rightToLeft ? msp->send : msp->sstart;

  char displayText[numDigitsInInt(coord) + 1];
  sprintf(displayText, "%d", coord);
  g_object_set(renderer, "text", displayText, NULL);
}


/* Cell data function for the "end" column. This displays the end coord of the match 
 * sequence in normal left-to-right display, but the start coord if the display is reversed */
static void cellDataFunctionEndCol(GtkTreeViewColumn *column, 
				   GtkCellRenderer *renderer, 
				   GtkTreeModel *model, 
				   GtkTreeIter *iter, 
				   gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  gboolean rightToLeft = treeGetStrandsToggled(tree);
  
  /* Get the msp for this row */
  const MSP *msp = treeGetMsp(model, iter);
  
  /* We want to display the end coord, unless the display is reversed, in which case display the start */
  int coord = rightToLeft ? msp->sstart : msp->send;
  
  char displayText[numDigitsInInt(coord) + 1];
  sprintf(displayText, "%d", coord);
  g_object_set(renderer, "text", displayText, NULL);
}


static void addTreeColumns(GtkWidget *tree, GtkCellRenderer *renderer, const char *fontFamily)
{
  /* Add the columns */
  initColumn(tree,  renderer,  "Name",     "text", S_NAME_COL, FALSE, 90);
  initColumn(tree,  renderer,  "Score",    "text", SCORE_COL,  FALSE, 30);
  initColumn(tree,  renderer,  "%ID",      "text", ID_COL,	  FALSE, 30);
  GtkTreeViewColumn *startCol = initColumn(tree,  renderer,  "Start",    "text", START_COL,  FALSE, 40);
  GtkTreeViewColumn *seqCol   = initColumn(tree,  renderer,  "Sequence", "msp",  MSP_COL,  TRUE, 40);
  GtkTreeViewColumn *endCol   = initColumn(tree,  renderer,  "End",      "text", END_COL,	  FALSE, 30);

  /* Set custom data functions for start/end cols (so that coords flip when display is reversed) */
  gtk_tree_view_column_set_cell_data_func(startCol, renderer, cellDataFunctionStartCol, tree, NULL);
  gtk_tree_view_column_set_cell_data_func(endCol, renderer, cellDataFunctionEndCol, tree, NULL);
  
  /* Set it to use a fixed width font. This also sets the row height based on the font size.
   * For now we just use the default size from the theme. */
  int origSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE);
  setCellRendererFont(tree, renderer, origSize, fontFamily);
  
  seqColCreateProperties(seqCol, tree);
}


/* Create the base data store for a detail view tree */
void treeCreateBaseDataModel(GtkWidget *tree, gpointer data)
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
  char *refSeq = treeGetRefSeq(tree);
  
  MSP *refSeqMsp = createRefSeqMsp(tree,
				   treeGetStrand(tree) == FORWARD_STRAND, 
				   refSeq, 
				   treeGetDisplayRange(tree));
  
  addMspToTreeModel(GTK_TREE_MODEL(store), refSeqMsp, tree);
		    
  
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), GTK_TREE_MODEL(store));
  g_object_unref(G_OBJECT(store));
}


/* Create a filtered version of the data store to only show rows that are within the display range. */
void treeCreateFilteredDataModel(GtkWidget *tree, gpointer data)
{
  GtkTreeModel *baseModel = gtk_tree_view_get_model(GTK_TREE_VIEW(tree));
  GtkTreeModel *filter = gtk_tree_model_filter_new(GTK_TREE_MODEL(baseModel), NULL);
  
  gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(filter), 
					 (GtkTreeModelFilterVisibleFunc)isTreeRowVisible, 
					 tree, 
					 NULL);
  
  /* Add the filtered store to the tree view */
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), GTK_TREE_MODEL(filter));
  g_object_unref (G_OBJECT (filter));
}


static void setTreeStyle(GtkTreeView *tree)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_VERTICAL);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_SINGLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, FALSE);

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
				GList **treeList,
				const char *fontFamily)
{
  /* Create a tree view for the list of match sequences */
  GtkWidget *tree = gtk_tree_view_new();
  setTreeStyle(GTK_TREE_VIEW(tree));
  
  /* Put it in a scrolled window for vertical scrolling only (hoz scrolling will be via our
   * custom adjustment). Always display the scrollbars because we assume the column widths 
   * are the same for all trees and they won't be if one shows a scrollbar and another doesn't. */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
  gtk_container_add(GTK_CONTAINER(scrollWin), tree);
  g_object_ref(scrollWin); /* so we can remove/re-add without worrying about it getting destroyed */
  
  /* Add the tree to the given list. We add its overall container, so we can treat the thing as a whole. */
  *treeList = g_list_append(*treeList, scrollWin);

  /* Create a custom cell renderer to render the match sequences in this tree */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();

  /* The tree needs to know which grid and renderer it corresponds to, and vice versa */
  treeCreateProperties(tree, grid, detailView, renderer);
  gridGetProperties(grid)->tree = tree;
  SEQUENCE_CELL_RENDERER(renderer)->tree = tree;

  /* Add the columns */
  addTreeColumns(tree, renderer, fontFamily);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);
  
  g_signal_connect(G_OBJECT(tree), "button-press-event",    G_CALLBACK(onButtonPressTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "button-release-event",  G_CALLBACK(onButtonReleaseTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "motion-notify-event",   G_CALLBACK(onMouseMoveTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "scroll-event",	    G_CALLBACK(onScrollTree),		detailView);
  g_signal_connect(G_OBJECT(tree), "enter-notify-event",    G_CALLBACK(onEnterTree),		NULL);
  g_signal_connect(G_OBJECT(tree), "leave-notify-event",    G_CALLBACK(onLeaveTree),		NULL);
  g_signal_connect(G_OBJECT(tree), "drag-begin",	    G_CALLBACK(onDragBeginTree),	NULL);
  g_signal_connect(G_OBJECT(tree), "drag-end",		    G_CALLBACK(onDragEndTree),		NULL);
  g_signal_connect(G_OBJECT(tree), "drag-motion",	    G_CALLBACK(onDragMotionTree),	NULL);
  g_signal_connect(G_OBJECT(tree), "expose-event",	    G_CALLBACK(onExposeDetailViewTree), NULL);
  
  GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
  g_signal_connect(G_OBJECT(selection), "changed", G_CALLBACK(onSelectionChangedTree), tree);
  
  messout("Created detail-view tree [%x] in container [%x]", tree, scrollWin);
  
  GTK_WIDGET(tree);
  return scrollWin;
}