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
#include <SeqTools/sequencecellrenderer.h>
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/utilities.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>

enum {SORT_BY_NAME, SORT_BY_ID, SORT_BY_SCORE, SORT_BY_POS} SortType;


/* Local function declarations */
static void		onSelectionChangedTree(GObject *selection, gpointer data);
static gint		sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data);
static int		calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree);
static gboolean		isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data);

/***********************************************************
 *                Tree - utility functions                 *
 ***********************************************************/

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
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetRenderer(detailView);
}

GtkWidget *treeGetDetailView(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->detailView : NULL;
}

GtkWidget *treeGetMainWindow(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetMainWindow(detailView);
}

GHashTable *treeGetSeqTable(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties->seqTable;
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

int treeGetFrame(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->readingFrame : UNSET_INT;
}

BlxSeqType treeGetSeqType(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetSeqType(detailView);
}

int treeGetSelectedBaseIdx(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetSelectedBaseIdx(properties->detailView) : UNSET_INT;
}

/* Sets the currently-selected base index (and the base number within the reading frame
 * for protein matches). If allowScroll is TRUE, the display range will be scrolled if
 * necessary to keep the selected base index in view. */
static void treeSetSelectedBaseIdx(GtkWidget *tree, const int selectedBaseIdx, const int frame, const int baseNum, const gboolean allowScroll)
{
  assertTree(tree);

  TreeProperties *properties = treeGetProperties(tree);
  if (properties)
    {
      DetailViewProperties *detailViewProperties = detailViewGetProperties(properties->detailView);
      
      /* Only update if things have changed */
      if (detailViewProperties->selectedBaseIdx != selectedBaseIdx ||
	  detailViewProperties->selectedBaseNum != baseNum ||
	  detailViewProperties->selectedFrame != frame)
	{
	  detailViewSetSelectedBaseIdx(properties->detailView, selectedBaseIdx, frame, baseNum, allowScroll);
	}
    }
}

/* Get the reference sequence (always the forward strand and always the DNA seq) */
char* treeGetRefSeq(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetRefSeq(detailView);
}

/* Returns the currently-displayed range in the tree view */
IntRange* treeGetDisplayRange(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetDisplayRange(properties->detailView) : NULL;
}

/* Returns the full display range */
IntRange* treeGetFullRange(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetFullRange(detailView);
}

/* Returns the range of the reference sequence (as a DNA sequence, regardless of
 * whether it is shown as a peptide sequence in reality) */
IntRange* treeGetRefSeqRange(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetRefSeqRange(detailView);
}

PangoFontDescription* treeGetFontDesc(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetFontDesc(detailView);
}

int treeGetCellXPadding(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCellXPadding(detailView);
}

int treeGetCellYPadding(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCellYPadding(detailView);
}

int treeGetCharWidth(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCharWidth(detailView);
}

int treeGetCharHeight(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCharHeight(detailView);
}

GdkColor* treeGetRefSeqColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetRefSeqColour(detailView, selected);
}

GdkColor* treeGetMatchColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetMatchColour(detailView, selected);
}

GdkColor* treeGetMismatchColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetMismatchColour(detailView, selected);
}

GdkColor* treeGetConsColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetConsColour(detailView, selected);
}

GdkColor* treeGetExonColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetExonColour(detailView, selected);
}

GdkColor* treeGetGapColour(GtkWidget *tree, const gboolean selected)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetGapColour(detailView, selected);
}

GdkColor* treeGetExonBoundaryColour(GtkWidget *tree, const gboolean isStart)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetExonBoundaryColour(detailView, isStart);
}

int treeGetExonBoundaryWidth(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetExonBoundaryWidth(detailView);
}

GdkLineStyle treeGetExonBoundaryStyle(GtkWidget *tree, const gboolean isStart)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetExonBoundaryStyle(detailView, isStart);
}


/* Call the given function on all trees in the detail view */
void callFuncOnAllDetailViewTrees(GtkWidget *detailView, gpointer data)
{
  int numReadingFrames = detailViewGetNumReadingFrames(detailView);
  GtkCallback func = (GtkCallback)data;
  
  /* Call the function on the forward strand tree and reverse strand tree
   * for each frame. */
  int frame = 1;
  for ( ; frame <= numReadingFrames; ++frame)
    {
      GtkWidget *fwdTree = detailViewGetTree(detailView, FORWARD_STRAND, frame);
      GtkWidget *revTree = detailViewGetTree(detailView, REVERSE_STRAND, frame);
      
      if (fwdTree)
	{
	  func(fwdTree, NULL);
	}
      
      if (revTree)
	{
	  func(revTree, NULL);
	}
    }
}


/* Add the list of MSPs in the hash table under the given key to the tree. */
static void addMspListToTreeRow(gpointer key, gpointer value, gpointer data)
{
  GList *mspGList = (GList*)value;
  
  if (g_list_length(mspGList) > 0)
    {
      /* Get the first msp in the list */
      MSP *msp = (MSP*)(mspGList->data);
      
      /* Add a row to the tree store */
      GtkListStore *store = GTK_LIST_STORE(data);
      GtkTreeIter iter;
      gtk_list_store_append(store, &iter);
      
      /* Populate the row. If there are multiple MSPs, most of the column info is not applicable */
      if (g_list_length(mspGList) == 1)
	{
	  gtk_list_store_set(store, &iter,
			     S_NAME_COL, msp->sname,
			     SCORE_COL, msp->score,
			     ID_COL, msp->id,
			     START_COL, msp->sstart,
			     SEQUENCE_COL, mspGList,
			     END_COL, msp->send,
			     -1);
	}
      else
	{
	  gtk_list_store_set(store, &iter,
			     S_NAME_COL, msp->sname,
			     SCORE_COL, NULL,
			     ID_COL, NULL,
			     START_COL, NULL,
			     SEQUENCE_COL, mspGList,
			     END_COL, NULL,
			     -1);
	}
    }
}


void addSequencesToTree(GtkWidget *tree, gpointer data)
{
  GtkListStore *store = gtk_list_store_new(N_COLUMNS,
					   G_TYPE_STRING, 
					   G_TYPE_INT, 
					   G_TYPE_INT, 
					   G_TYPE_INT, 
					   G_TYPE_POINTER, 
					   G_TYPE_INT);
  
  /* Set the sort functions for each column */
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), S_NAME_COL, sortColumnCompareFunc, (gpointer)SORT_BY_NAME, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), ID_COL, sortColumnCompareFunc, (gpointer)SORT_BY_ID, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), SCORE_COL, sortColumnCompareFunc, (gpointer)SORT_BY_SCORE, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), SEQUENCE_COL, sortColumnCompareFunc, (gpointer)SORT_BY_POS, NULL);

  /* Add the rows. Use the hash table that groups MSPs by sequence name and qFrame/qStrand */
  GHashTable *seqTable = treeGetSeqTable(tree);
  g_hash_table_foreach(seqTable, addMspListToTreeRow, store);

  /* Create a filtered version which will only show sequences that are in the display range */
  GtkTreeModel *filter = gtk_tree_model_filter_new(GTK_TREE_MODEL(store), NULL);
  gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(filter), (GtkTreeModelFilterVisibleFunc)isTreeRowVisible, tree, NULL);
  g_object_unref(G_OBJECT(store));
  
  /* Remember the tree store in the properties so we can switch between this and the 'unsquashed' tree model */
  TreeProperties *properties = treeGetProperties(tree);
  properties->seqTreeModel = GTK_TREE_MODEL(filter);
}


/* Return the MSP(s) in a given tree row */
GList* treeGetMsps(GtkTreeModel *model, GtkTreeIter *iter)
{
  GList *mspGList = NULL;
  gtk_tree_model_get(model, iter, SEQUENCE_COL, &mspGList, -1);
  return mspGList;
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
//static GtkTreePath *treeConvertBasePathToVisiblePath(GtkTreeView *tree, GtkTreePath *basePath)
//{
//  GtkTreePath *result = NULL;
//  
//  if (tree && basePath)
//    {
//      /* Convert the child path to the equivalent path in the filtered model. The
//       * result may be null if the child row does not appear in the filtered model. */
//      GtkTreeModel *filter = treeGetVisibleDataModel(tree);
//      if (GTK_IS_TREE_MODEL_FILTER(filter))
//	{
//	  result = gtk_tree_model_filter_convert_child_path_to_path(GTK_TREE_MODEL_FILTER(filter), basePath);
//	}
//      else
//	{
//	  result = basePath;
//	}
//    }
//  
//  return result;
//}


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
void treeUpdateFontSize(GtkWidget *tree, gpointer data)
{
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);  
  gtk_widget_modify_font(tree, fontDesc);
}


/* Get the character index at the given x coordinate in the given tree view column. */
static int getCharIndexAtTreeCoord(GtkWidget *tree, GtkTreeViewColumn *col, const int x, int *baseNum)
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
  double leftEdge = (double)renderer->xpad + renderer->xalign;
  double rightEdge = colWidth - renderer->xpad;
  
  if (x >= leftEdge && x <= rightEdge)
    {
      /* Calculate the number of bases from the left edge */
      gint charWidth = treeGetCharWidth(tree);
      
      const double numBasesFromLeftDbl = ((double)x - leftEdge) / (double)charWidth;
      const int numBasesFromLeft = (int)numBasesFromLeftDbl;
      
      /* If requested, also calculate the base number within the frame */
      if (baseNum)
	{
	  const int distFromCharStart = roundNearest((numBasesFromLeftDbl - numBasesFromLeft) * (double)charWidth);
	  const double frameWidth = charWidth / treeGetNumReadingFrames(tree);
	  *baseNum = (distFromCharStart / frameWidth) + 1;
	}
      
      if (rightToLeft)
	{
	  /* Numbers are displayed increasing from right-to-left. Subtract the 
	   * number of bases from the total number of bases displayed. */
	  IntRange *displayRange = treeGetDisplayRange(tree);
	  int displayLen = displayRange->max - displayRange->min;
	  result = displayLen - numBasesFromLeft;
	  
	  /* x could be in an empty gap at the end, so bounds check the result.
	   * Note that our index is 0-based but the display range is 1-based. */
	  if (result < 0)
	    {
	      result = UNSET_INT;
	    }
	}
      else
	{
	  /* Normal left-to-right display. */
	  result = numBasesFromLeft;
	}
    }
  
  return result;
}


/* In the detail view, get the base index at the given coords */
static int getBaseIndexAtTreeCoords(GtkWidget *tree, const int x, const int y, int *baseNum)
{
  int baseIdx = UNSET_INT;
  
  GtkTreePath *path = NULL;
  GtkTreeViewColumn *col = NULL;
  int cell_x, cell_y; /* coords relative to cell */
  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), x, y, &path, &col, &cell_x, &cell_y);
  
  /* See if we clicked in the sequence column */
  if (col == gtk_tree_view_get_column(GTK_TREE_VIEW(tree), SEQUENCE_COL))
    {
      /* Get the base index at the clicked position */
      int charIdx = getCharIndexAtTreeCoord(tree, col, cell_x, baseNum);
      
      if (charIdx != UNSET_INT)
	{
	  GtkAdjustment *adjustment = treeGetAdjustment(tree);
	  baseIdx = charIdx + adjustment->value + 1; /* Adjustment is 0-based, display range is 1-based */
	}
    }
  
  return baseIdx;
}


/* This function makes multiple MSPs in the same sequence appear in the same row in the tree. */
void treeSquashMatches(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), properties->seqTreeModel);

  /* Refilter the new data store, because it might not be up to date */
  refilterTree(tree, NULL);
}


/* This function makes all MSPs appear in their own individual rows in the tree. */
void treeUnsquashMatches(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), properties->mspTreeModel);
  
  /* Refilter the new data store, because it might not be up to date */
  refilterTree(tree, NULL);
}


/* Returns true if the matches are squashed */
gboolean treeGetMatchesSquashed(GtkWidget *tree)
{
  gboolean result = FALSE;
  
  /* Check which stored model the tree view is currently looking at */
  GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(tree));
  TreeProperties *properties = treeGetProperties(tree);
  
  if (model == properties->mspTreeModel)
    {
      result = FALSE;
    }
  else if (model == properties->seqTreeModel)
    {
      result = TRUE;
    }
  else
    {
      messerror("Unexpected tree data store [%x]. Expected either [%x] (normal data) or [%x] (condensed data)", model, properties->mspTreeModel, properties->seqTreeModel);
    }
  
  return result;
}

/***********************************************************
 *			      Updates			   *
 ***********************************************************/

/* Refresh the tree header widgets */
void refreshTreeHeaders(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  GList *header = properties->treeColumnHeaderList;
  
  for ( ; header; header = header->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)header->data;
      if (headerInfo && headerInfo->headerWidget)
	{
	  /* Set the background colour. Set it in the parent too seeing as labels don't have a window themselves. */
//	  GdkColor *bgColour = treeGetRefSeqColour(tree, FALSE); //to do: should use this but doesn't work - comes out black for some reason
	  GdkColor bgColour2;
	  gdk_color_parse("yellow", &bgColour2);
	  gtk_widget_modify_bg(headerInfo->headerWidget, GTK_STATE_NORMAL, &bgColour2);

	  GtkWidget *parent = gtk_widget_get_parent(headerInfo->headerWidget);
	  gtk_widget_modify_bg(parent, GTK_STATE_NORMAL, &bgColour2);
	  
	  /* Update the font and widget size */
	  gtk_widget_modify_font(headerInfo->headerWidget, treeGetFontDesc(tree));
	  int width = calculateColumnWidth(headerInfo, tree);
	  gtk_widget_set_size_request(headerInfo->headerWidget, width, -1);
	  gtk_widget_queue_resize(headerInfo->headerWidget);
	  
	  if (headerInfo->refreshFunc)
	    {
	      /* Refresh the contents */
	      headerInfo->refreshFunc(headerInfo->headerWidget, tree);
	    }
	}
      else
	{
	  messerror("refreshTreeHeaders: Warning - Invalid tree header info. Tree header may not refresh properly.");
	}
    }
}


/***********************************************************
 *			  Selections			   *
 ***********************************************************/

/* Deselect all rows in the given tree */
void deselectAllRows(GtkWidget *tree, gpointer data)
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
	  deselectAllRows(widget, NULL);
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


gboolean treeIsSeqSelected(GtkWidget *tree, const char *seqName)
{
  GtkWidget *mainWindow = treeGetMainWindow(tree);
  return mainWindowIsSeqSelected(mainWindow, seqName);
}


/* Select the given row if its sequence is marked as selected. */
static gboolean selectRowIfSeqSelected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GList *mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      /* Check any MSP in this row to see if it's seq is selected (they should all have the same seq name) */
      MSP *msp = (MSP*)(mspGList->data);
  
      if (treeIsSeqSelected(tree, msp->sname))
      {
	GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
	gtk_tree_selection_select_path(selection, path);
      }
    }
  
  return FALSE;
}


/* Select all rows in this tree whose sequences are marked as selected */
void selectRowsForSelectedSeqs(GtkWidget *tree, gpointer data)
{
  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
  gtk_tree_model_foreach(model, selectRowIfSeqSelected, tree);
}


/* Mark the given row's sequence as selected in the main window's list of selected seqs.
 * Does not allow the main window to re-update the tree selection */
static void markRowMspSelected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *mainWindow = treeGetMainWindow(tree);
  
  GList *mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      /* Get the sequence name for any MSP in this row (they should all have the same seq) */
      MSP *msp = (MSP*)(mspGList->data);
      
      mainWindowSelectSeq(mainWindow, msp->sname, FALSE);
      gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(tree), path, NULL, FALSE, 0.0, 0.0);
    }
}

static void onSelectionChangedTree(GObject *selection, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  gtk_tree_selection_selected_foreach(GTK_TREE_SELECTION(selection), markRowMspSelected, tree);
}


/***********************************************************
 *                   Sorting and filtering                 *
 ***********************************************************/



/* Refilter the data for the given tree */
void refilterTree(GtkWidget *tree, gpointer data)
{
  assertTree(tree);
  
  GtkTreeModelFilter *filter = GTK_TREE_MODEL_FILTER(gtk_tree_view_get_model(GTK_TREE_VIEW(tree)));
  gtk_tree_model_filter_refilter(filter);

  selectRowsForSelectedSeqs(tree, NULL);
}


/* Filter function. Returns true if the given row in the given tree model should be visible. */
static gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  gboolean bDisplay = FALSE;

  GtkWidget *tree = GTK_WIDGET(data);
  const IntRange const *displayRange = treeGetDisplayRange(tree);
  
  /* Loop through all MSPs in this row */
  GList *mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP* msp = (MSP*)(mspListItem->data);
      
      /* Don't show introns */
      if (msp && !mspIsIntron(msp))
	{
	  /* Only show this row if part of the MSP's range is inside the displayed range */
	  GtkAdjustment *adjustment = treeGetAdjustment(tree);
	  if (adjustment)
	    {
	      int qSeqMin = min(msp->displayStart, msp->displayEnd);
	      int qSeqMax = max(msp->displayStart, msp->displayEnd);
	      
	      if (!(qSeqMin > displayRange->max || qSeqMax < displayRange->min))
		{
		  bDisplay = TRUE;
		  break;
		}
	    }
	}
    }
  
  return bDisplay;
}


void treeSortByName(GtkWidget *tree, gpointer data)
{
  GtkTreeSortable *model = GTK_TREE_SORTABLE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  gtk_tree_sortable_set_sort_column_id(model, S_NAME_COL, GTK_SORT_ASCENDING);
}

void treeSortById(GtkWidget *tree, gpointer data)
{
  GtkTreeSortable *model = GTK_TREE_SORTABLE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  gtk_tree_sortable_set_sort_column_id(model, ID_COL, GTK_SORT_ASCENDING);
}

void treeSortByScore(GtkWidget *tree, gpointer data)
{
  GtkTreeSortable *model = GTK_TREE_SORTABLE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  gtk_tree_sortable_set_sort_column_id(model, SCORE_COL, GTK_SORT_ASCENDING);
}

void treeSortByPos(GtkWidget *tree, gpointer data)
{
  GtkTreeSortable *model = GTK_TREE_SORTABLE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  
  /* Sort ascending if the reference sequence is displayed in the normal left-to-right
   * direction, otherwise sort descending. */
//  gboolean rightToLeft = treeGetStrandsToggled(tree);
  gtk_tree_sortable_set_sort_column_id(model, SEQUENCE_COL, GTK_SORT_ASCENDING);
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
				 GtkCellRenderer *renderer,
				 const int frame,
				 GList *treeColumnHeaderList)
{
  if (widget)
    { 
      TreeProperties *properties = g_malloc(sizeof *properties);
      
      properties->grid = grid;
      properties->detailView = detailView;
      properties->renderer = renderer;
      properties->readingFrame = frame;
      properties->treeColumnHeaderList = treeColumnHeaderList;
      properties->seqTable = g_hash_table_new(g_str_hash, g_str_equal);
      properties->mspTreeModel = NULL;
      properties->seqTreeModel = NULL;

      g_object_set_data(G_OBJECT(widget), "TreeProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyTree), NULL); 
    }
}


/***********************************************************
 *                       Tree events                       *
 ***********************************************************/

static gboolean onExposeDetailViewTree(GtkWidget *tree, GdkEventExpose *event, gpointer data)
{
  /* Create a new drawable to draw to. Our custom cell renderer will draw to this as
   * well as the main window. (Ideally we'd just draw to the bitmap and then push
   * this to the screen, but I'm not sure if it's possible to detect when the 
   * cell renderer has finished drawing.) */
  GdkDrawable *drawable = gdk_pixmap_new(tree->window, tree->allocation.width, tree->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(tree, drawable);

  /* Draw a blank rectangle of the required widget background colour */
  GdkGC *gc = gdk_gc_new(drawable);
  //GtkStyle *style = gtk_widget_get_style(tree);
  GdkColor bgColour = getGdkColor(GDK_WHITE);
  gdk_gc_set_foreground(gc, &bgColour);
  
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, tree->allocation.width, tree->allocation.height);
  
  /* Let the default handler continue */
  return FALSE;
}


static gboolean onButtonPressTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
    {
    case 1:
      {
	/* Left button: selects a row, or shows the pfetch window if this was a double-click. */
	if (event->type == GDK_BUTTON_PRESS)
	  {
	    mainWindowDeselectAllSeqs(treeGetMainWindow(tree), FALSE);
	    deselectAllSiblingTrees(tree, FALSE);
	    
	    /* If we've clicked a different tree to the previously selected one, make this tree's 
	     * reading frame the active one. Select the first base in the triplet (or the last
	     * if the display is reversed). */
	    const int firstBaseNum = treeGetStrandsToggled(tree) ? treeGetNumReadingFrames(tree) : 1;
	    treeSetSelectedBaseIdx(tree, treeGetSelectedBaseIdx(tree), treeGetFrame(tree), firstBaseNum, FALSE);
	    
	    /* Let the default handler select the row */
	    handled = FALSE;
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
	    /* Get the selected sequence (assumes that only one is selected) */
	    GtkWidget *mainWindow = treeGetMainWindow(tree);
	    GList *selectedSeqs = mainWindowGetSelectedSeqs(mainWindow);
	    
	    if (selectedSeqs)
	      {
		const char *seqName = (const char*)selectedSeqs->data;
		blxDisplayMSP(seqName, 0, mainWindow);
	      }
	      
	    handled = TRUE;
	  }
	  
	break;
      }

    case 2:
      {
	/* Middle button: select the base index at the clicked coords */
	int baseIdx = getBaseIndexAtTreeCoords(tree, event->x, event->y, NULL);
	
	/* Select the first base in this peptide's triplet (or the last if the display is reversed) */
	const int firstBaseNum = treeGetStrandsToggled(tree) ? treeGetNumReadingFrames(tree) : 1;
	treeSetSelectedBaseIdx(tree, baseIdx, treeGetFrame(tree), firstBaseNum, TRUE);
	
	handled = TRUE;
	break;
      }
      
    case 3:
      {
	/* Right button: show context menu. */
	GtkWidget *mainmenu = mainWindowGetMainMenu(treeGetMainWindow(tree));
	gtk_menu_popup (GTK_MENU(mainmenu), NULL, NULL, NULL, NULL, event->button, event->time);
	handled = TRUE;
	break;
      }
      
    default:
      break;
    };
  
  return handled;
}


static gboolean onKeyPressTree(GtkWidget *tree, GdkEventKey *event, gpointer data)
{
  if (event->keyval == GDK_Up || event->keyval == GDK_Down)
    {
      /* Mark the currently selected MSPs as no longer selected unless selection
       * has not changed */
      GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
      GtkTreeModel *model;
      GList *listItem = gtk_tree_selection_get_selected_rows(selection, &model);
      
      if (listItem)
	{
	  GtkTreePath *path = (GtkTreePath*)(listItem->data);
	  
	  if (event->keyval == GDK_Up)
	    {
	      if (gtk_tree_path_prev(path))
		{
		  mainWindowDeselectAllSeqs(treeGetMainWindow(tree), FALSE);
		}
	    }
	  else
	    {
	      /* For some reason gtk_tree_path_next doesn't indicate if next is null, so we have to get an iter... */
	      GtkTreeIter iter;
	      gtk_tree_model_get_iter(model, &iter, path);
	      
	      if (gtk_tree_model_iter_next(model, &iter))
		{
		  mainWindowDeselectAllSeqs(treeGetMainWindow(tree), FALSE);
		}
	    }
	}
    }
  
  return FALSE;
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
  
  /* Middle button: scroll the selected base index to the centre (unless CTRL is pressed) */
  if (event->button == 2)
    {
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);

      if (!ctrlModifier)
	{
	  /* Move the scrollbar so that the currently-selected base index is at the centre */
	  GtkWidget *detailView = GTK_WIDGET(data);
	  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
	  
	  if (selectedBaseIdx != UNSET_INT)
	    {
	      /* The coord is in terms of the display coords, i.e. whatever the displayed seq type is. */
	      const BlxSeqType seqType = detailViewGetSeqType(detailView);
	      const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
	      
	      int newStart = selectedBaseIdx - (getRangeLength(displayRange) / 2);
	      setDetailViewStartIdx(detailView, newStart, seqType);
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
  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base
       * (but don't re-centre on the selected base until the mouse button is released). */
      const int selectedBaseIdx = getBaseIndexAtTreeCoords(tree, event->x, event->y, NULL);
      const int frame = treeGetFrame(tree);

      /* For protein matches, get the 1st base in the peptide (or last base if display reversed) */
      int baseNum = treeGetStrandsToggled(tree) ? treeGetNumReadingFrames(tree) : 1;
      
      treeSetSelectedBaseIdx(tree, selectedBaseIdx, frame, baseNum, TRUE);
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


/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

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
  const gboolean sForward = (strchr(msp->sframe, '+')) ? TRUE : FALSE ;
  const gboolean qForward = (strchr(msp->qframe, '+')) ? TRUE : FALSE ;

  const BlxBlastMode blastMode = treeGetBlastMode(tree);
  const int numFrames = treeGetNumReadingFrames(tree);

  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  msp->id = 0;
  
  if (msp->sseq /* to do: is this required? && msp->sseq != padseq */)
    {
      /* Note that this will reverse complement the ref seq if it is the reverse 
       * strand. This means that where there is no gaps array the comparison is trivial
       * as coordinates can be ignored and the two sequences just whipped through. */
      GtkWidget *detailView = treeGetDetailView(tree);
      char *refSeqSegment = getSequenceSegment(detailViewGetMainWindow(detailView),
					       detailViewGetRefSeq(detailView),
					       detailViewGetRefSeqRange(detailView),
					       msp->qstart, 
					       msp->qend, 
					       treeGetStrand(tree), 
					       BLXSEQ_DNA, /* msp q coords are always on the dna sequence */
					       treeGetFrame(tree),
					       treeGetNumReadingFrames(tree),
					       !qForward,
					       TRUE);

      if (!refSeqSegment)
	{
	  messout ( "calcID failed: Don't have genomic sequence %d - %d, requested for match sequence '%s' (match coords = %d - %d)\n", msp->qstart, msp->qend, msp->sname, msp->sstart, msp->send);
	  msp->id = 0;
	  return;
	}

      /* We need to find the number of characters that match out of the total number */
      int numMatchingChars = 0;
      int totalNumChars = 0;
      
      if (!(msp->gaps) || arrayMax(msp->gaps) == 0)
	{
	  /* Ungapped alignments. */
	  totalNumChars = (qSeqMax - qSeqMin + 1) / numFrames;

	  if (blastMode == BLXMODE_TBLASTN || blastMode == BLXMODE_TBLASTX)
	    {
	      int i = 0;
	      for ( ; i < totalNumChars; i++)
		{
		  if (freeupper(msp->sseq[i]) == freeupper(refSeqSegment[i]))
		    {
		      numMatchingChars++;
		    }
		}
	    }
	  else						    /* blastn, blastp & blastx */
	    {
	      int i = 0;
	      for ( ; i < totalNumChars; i++)
                {
                  int sIndex = sForward ? sSeqMin + i - 1 : sSeqMax - i - 1;
		  if (freeupper(msp->sseq[sIndex]) == freeupper(refSeqSegment[i]))
		    {
		      numMatchingChars++;
		    }
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
	      
              Array gaps = msp->gaps;
	      
              int i = 0;
	      for ( ; i < arrayMax(gaps) ; i++)
		{
		  SMapMap *m = arrp(gaps, i, SMapMap) ;
		  
                  int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
		  getSMapMapRangeExtents(m, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
		  
                  totalNumChars += sRangeMax - sRangeMin + 1;
                  
                  /* Note that refSeqSegment is just the section of the ref seq relating to this msp.
		   * We need to translate the first coord in the range (which is in terms of the full
		   * reference sequence) into coords in the cut-down ref sequence. */
                  int q_start = qForward ? (qRangeMin - qSeqMin) / numFrames : (qSeqMax - qRangeMax) / numFrames;
		  
		  /* We can index sseq directly (but we need to adjust by 1 for zero-indexing). We'll loop forwards
		   * through sseq if we have the forward strand or backwards if we have the reverse strand,
		   * so start from the lower or upper end accordingly. */
                  int s_start = sForward ? sRangeMin - 1 : sRangeMax - 1 ;
		  
                  int sIdx = s_start, qIdx = q_start ;
		  while ((sForward && sIdx < sRangeMax) || (!sForward && sIdx >= sRangeMin - 1))
		    {
		      if (freeupper(msp->sseq[sIdx]) == freeupper(refSeqSegment[qIdx]))
			numMatchingChars++ ;
		      
                      /* Move to the next base. The refSeqSegment is always forward, but we might have to
                       * traverse the s sequence in reverse. */
                      ++qIdx ;
                      if (sForward) ++sIdx ;
                      else --sIdx ;
		    }
		}
	    }
	}
      
      msp->id = (int)((100.0 * numMatchingChars / totalNumChars) + 0.5);
      
      g_free(refSeqSegment);
    }
  
  return ;
}


/* Calculate the display coordinates for the given MSP and store them in the
 * MSP struct. Reference sequence (q) coords are always in terms of the DNA
 * sequence, whereas the display may be in terms of the peptide sequence if
 * we are viewing protein matches. */
static void calcDisplayCoords(MSP *msp, GtkWidget *tree)
{
  const int frame = treeGetFrame(tree);
  const int numReadingFrames = treeGetNumReadingFrames(tree);

  msp->displayStart = convertDnaToPeptide(msp->qstart, frame, numReadingFrames, NULL); 
  msp->displayEnd = convertDnaToPeptide(msp->qend, frame, numReadingFrames, NULL);
}


/* Add a row to the given tree containing the given MSP */
static void addMspToTreeRow(MSP *msp, GtkWidget *tree)
{
  GtkListStore *store = GTK_LIST_STORE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  
  GtkTreeIter iter;
  gtk_list_store_append(store, &iter);
  
  /* The SequenceCellRenderer expects a GList of MSPs, so put our MSP in a list */
  GList *mspGList = NULL;
  mspGList = g_list_append(NULL, msp);
  
  gtk_list_store_set(store, &iter,
		     S_NAME_COL, msp->sname,
		     SCORE_COL, msp->score,
		     ID_COL, msp->id,
		     START_COL, msp->sstart,
		     SEQUENCE_COL, mspGList,
		     END_COL, msp->send,
		     -1);
}


/* Add the given msp as a row in the given tree view, and also adds it to the
 * tree's hash tables. */
void addMspToTree(GtkWidget *tree, MSP *msp)
{
  /* First calculate the ID and display coords and update them in the MSP */
  calcID(msp, tree);
  calcDisplayCoords(msp, tree);

  /* Add the MSP as an individual row */
  addMspToTreeRow(msp, tree);
  
  /* Add the MSP to the hash table that will group this tree's MSPs by sequence name */
  GHashTable *seqTable = treeGetSeqTable(tree);
  addMspToHashTable(seqTable, msp, msp->sname);
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
  
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) == 1)
    {
      MSP *msp = (MSP*)(mspGList->data);
      
      /* We want to display the start coord, unless the display is reversed, in which case display the end */
      int coord = rightToLeft ? msp->send : msp->sstart;

      char displayText[numDigitsInInt(coord) + 1];
      sprintf(displayText, "%d", coord);
      g_object_set(renderer, START_COLUMN_PROPERTY_NAME, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, START_COLUMN_PROPERTY_NAME, "", NULL);
    }
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
  
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) == 1)
    {
      MSP *msp = (MSP*)(mspGList->data);
      
      /* We want to display the end coord, unless the display is reversed, in which case display the start */
      int coord = rightToLeft ? msp->sstart : msp->send;
      
      char displayText[numDigitsInInt(coord) + 1];
      sprintf(displayText, "%d", coord);
      g_object_set(renderer, END_COLUMN_PROPERTY_NAME, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, START_COLUMN_PROPERTY_NAME, "", NULL);
    }
}


/* Cell data function for a generic integer column. */
static void cellDataFunctionIntCol(GtkTreeViewColumn *column, 
				   GtkCellRenderer *renderer, 
				   GtkTreeModel *model, 
				   GtkTreeIter *iter, 
				   gpointer data)
{
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) != 1)
    {
      g_object_set(renderer, START_COLUMN_PROPERTY_NAME, "", NULL);
    }
}


/* Create a single column in the tree. */
static void initColumn(GtkWidget *tree, 
		       GtkCellRenderer *renderer, 
		       DetailViewColumnInfo *columnInfo)
{
  /* Create the column */
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(columnInfo->title, 
								       renderer, 
								       columnInfo->propertyName, columnInfo->columnId, 
								       "data", SEQUENCE_COL, /* always set msp so each column has access to all the data */
								       NULL);

  /* Reduce the width of the end col by the scrollbar width. (This is so it matches the width
   * of the header, which does not have a scrollbar. ) */
  int width = columnInfo->width;
  if (columnInfo->columnId == END_COL)
    {
      /* Create a temp scrollbar to find the default width from the style properties. */
      GtkWidget *scrollbar = gtk_vscrollbar_new(NULL);
      
      gint sliderWidth, separatorWidth, troughBorder, stepperSpacing;
      gtk_widget_style_get(scrollbar, "slider-width", &sliderWidth, NULL);
      gtk_widget_style_get(scrollbar, "separator-width", &separatorWidth, NULL);
      gtk_widget_style_get(scrollbar, "trough-border", &troughBorder, NULL);
      gtk_widget_style_get(scrollbar, "stepper-spacing", &stepperSpacing, NULL);

      gtk_widget_destroy(scrollbar);
      
      width = width - sliderWidth - separatorWidth*2 - troughBorder*2 - stepperSpacing*2 - 4; /* to do: find out why the extra fudge factor is needed here */
    }
  
  /* Set the column properties and add the column to the tree */
  gtk_tree_view_column_set_fixed_width(column, width);
  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_resizable(column, TRUE);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
  
  /* Special treatment for specific columns */
  switch (columnInfo->columnId)
  {
    case SEQUENCE_COL:
      gtk_tree_view_column_set_expand(column, TRUE);
      break;
      
    case START_COL:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionStartCol, tree, NULL);
      break;
      
    case END_COL:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionEndCol, tree, NULL);
      break;
      
    case SCORE_COL: /* fall through */
    case ID_COL:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionIntCol, tree, NULL);
      break;
      
    default:
      break;
  }
}


/* Set the markup in the given label, so that the label displays the given sequence
 * segment, with the given base index is highlighted in the given colour. If the given
 * base index is UNSET_INT, or is out of the current display range, then the label 
 * is just given the plain, non-marked-up sequence. */
static void createMarkupForLabel(GtkLabel *label, 
				 const char *segmentToDisplay, 
				 const int selectedBaseIdx,
				 GdkColor *colour,
				 const gboolean rightToLeft,
				 const IntRange const *displayRange)
{
  int charIdx = rightToLeft ? displayRange->max - selectedBaseIdx : selectedBaseIdx - displayRange->min;
  const int segLen = strlen(segmentToDisplay);
  
  if (charIdx >= 0 && charIdx < segLen)
    {
      /* Cut the sequence into 3 sections: the text before the selected base,
       * the text after it, and the selected base itself. */
      char textBefore[charIdx];
      char textAfter[segLen - charIdx + 200];
      char textSelection[2];
      
      int i = 0;
      for ( ; i < charIdx; ++i)
	{
	  textBefore[i] = segmentToDisplay[i];
	}
      textBefore[i] = '\0';
      
      textSelection[0] = segmentToDisplay[i];
      textSelection[1] = '\0';

      int j = 0;
      ++i;
      for ( ; i < segLen; ++i, ++j)
	{
	  textAfter[j] = segmentToDisplay[i];
	}
      textAfter[j] = '\0';

      /* Create the marked-up text */
      char markupText[segLen + 200];
      sprintf(markupText, "%s<span background='#%x'>%s</span>%s", 
	      textBefore, colour->pixel, textSelection, textAfter);

      gtk_label_set_markup(label, markupText);
    }
  else
    {
      /* Selected index is out of range, so no markup required */
      gtk_label_set_markup(label, segmentToDisplay);
    }
}


/* Refresh the sequence column header. This header shows the section of reference
 * sequence for the current display range, so it needs to be refreshed after
 * scrolling, zooming etc. */
static void refreshSequenceColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *detailView = treeGetDetailView(tree);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  
  /* Find the segment of the ref sequence to display (complemented if this tree is
   * displaying the reverse strand, and reversed if the display is toggled) */
  IntRange *displayRange = treeGetDisplayRange(tree);
  IntRange *refSeqRange = mainWindowGetRefSeqRange(mainWindow);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(mainWindow);

  gchar *segmentToDisplay = getSequenceSegment(mainWindow,
					       mainWindowGetRefSeq(mainWindow),
					       refSeqRange,
					       displayRange->min, 
					       displayRange->max, 
					       treeGetStrand(tree), 
					       mainWindowGetSeqType(mainWindow),
					       treeGetFrame(tree), 
					       mainWindowGetNumReadingFrames(mainWindow),
					       rightToLeft,
					       TRUE);

  if (segmentToDisplay)
    {
      const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);

      if (selectedBaseIdx == UNSET_INT)
        {
          /* Just draw plain text */
          gtk_label_set_markup(GTK_LABEL(headerWidget), segmentToDisplay);
        }
      else
        {
          GdkColor *selectedColour = treeGetRefSeqColour(tree, TRUE);
          createMarkupForLabel(GTK_LABEL(headerWidget), segmentToDisplay, selectedBaseIdx, selectedColour, rightToLeft, displayRange);
        }
      
      g_free(segmentToDisplay);
  }
}


/* Refresh the start column header. This displays the start index of the current display range */
static void refreshStartColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  if (GTK_IS_LABEL(headerWidget))
    {
      int displayVal = getStartDnaCoord(treeGetDisplayRange(tree), 
					treeGetFrame(tree),
					treeGetSeqType(tree), 
					treeGetStrandsToggled(tree), 
					treeGetNumReadingFrames(tree));
      
      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(headerWidget), displayText);
    }
  else
    {
      messerror("refreshStartColHeader: Column header is an unexpected widget type");
    }
}


/* Refresh the start column header. This displays the end index of the current display range */
static void refreshEndColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  if (GTK_IS_LABEL(headerWidget))
    {
      int displayVal = getEndDnaCoord(treeGetDisplayRange(tree),
				      treeGetFrame(tree),
				      treeGetSeqType(tree), 
				      treeGetStrandsToggled(tree), 
				      treeGetNumReadingFrames(tree));
      
      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(headerWidget), displayText);
    }
  else
    {
      messerror("refreshStartColHeader: Column header is an unexpected widget type");
    }
}


static int calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  
  /* Sum the width of all columns that this header includes */
  GList *listItem = headerInfo->columnIds;
  int width = 0;
  
  for ( ; listItem; listItem = listItem->next)
    {
      int columnId = GPOINTER_TO_INT(listItem->data);
      width = width + getDetailViewColumnWidth(detailView, columnId);
    }
  
  return width;
}


/* Create the header widget for the given column in the tree header. The tree
 * header bar shows information about the reference sequence. */
static void createTreeColHeader(GList **headerWidgets, 
				DetailViewColumnInfo *columnInfo,
				GtkWidget *headerBar,
				GtkWidget *tree,
				const char const *refSeqName,
				const int frame,
				const Strand strand)
{
  /* Create a header for this column, if required. Create a list of other 
   * columns we wish to merge under the same header. */
  GtkWidget *headerWidget = NULL;
  GList *columnIds = NULL;
  GtkCallback refreshFunc = NULL;
  
  switch (columnInfo->columnId)
    {
      case S_NAME_COL:
	{
	  /* The header above the name column will display the reference sequence name.
	   * This header will also span the score and id columns, seeing as we don't need
	   * to show any info in those columns. */
	  int textLen = strlen(refSeqName) + numDigitsInInt(frame) + 5;
	  char displayText[textLen];
	  sprintf(displayText, "%s (%s%d)", refSeqName, (strand == FORWARD_STRAND ? "+" : "-"), frame);
	  
	  headerWidget = createLabel(displayText, 0.0, 1.0, TRUE, TRUE);
	  
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(S_NAME_COL));
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(SCORE_COL));
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(ID_COL));
	  break;
	}
	
      case SEQUENCE_COL:
	{
	  /* The sequence column header displays the reference sequence. This needs a custom
	   * refresh callback function because it needs to be updated after scrolling etc. */
	  headerWidget = createLabel(NULL, 0.0, 1.0, TRUE, TRUE);
	  gtk_label_set_use_markup(GTK_LABEL(headerWidget), TRUE);
	  refreshFunc = refreshSequenceColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
	  break;
	}
	
      case START_COL:
	{
	  /* The start column header displays the start index of the current display range */
	  headerWidget = createLabel("", 0.0, 1.0, TRUE, TRUE);
	  refreshFunc = refreshStartColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
	  break;
	}
	
      case END_COL:
	{
	  /* The end column header displays the start index of the current display range */
	  headerWidget = createLabel("", 0.0, 1.0, TRUE, TRUE);
	  refreshFunc = refreshEndColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
	  break;
	}

      case SCORE_COL: /* fall through */
      case ID_COL:    /* fall through */
      default:
	break;        /* do nothing */
    };

  
  if (headerWidget != NULL)
    { 
      /* Put the widget in an event box so that we can colour its background. */
      GtkWidget *parent = gtk_event_box_new();
      gtk_container_add(GTK_CONTAINER(parent), headerWidget);
      
      /* Put the event box into the header bar */
      gtk_box_pack_start(GTK_BOX(headerBar), parent, (columnInfo->columnId == SEQUENCE_COL), TRUE, 0);
      
      /* Create the header info */
      TreeColumnHeaderInfo *headerInfo = g_malloc(sizeof(TreeColumnHeaderInfo));
      
      headerInfo->columnIds = columnIds;
      headerInfo->headerWidget = headerWidget;
      headerInfo->refreshFunc = refreshFunc;

      *headerWidgets = g_list_append(*headerWidgets, headerInfo);
    }
}


/* Create the columns. Returns a list of header info for the column headers */
static GList* addTreeColumns(GtkWidget *tree, 
			     GtkCellRenderer *renderer, 
			     const BlxSeqType seqType,
			     GList *columnList,
			     GtkWidget *headerBar,
			     const char const *refSeqName,
			     const int frame,
			     const Strand strand)
{
  /* We'll create the headers as we create the columns */
  GList *headerWidgets = NULL;

  /* The columns are defined by the columnList from the detail view. This list contains the
   * info such as the column width and title for each column. */
  GList *column = columnList;
  
  for ( ; column; column = column->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)column->data;
      
      if (columnInfo)
	{
	  initColumn(tree, renderer, columnInfo);
	  createTreeColHeader(&headerWidgets, columnInfo, headerBar, tree, refSeqName, frame, strand);
	}
      else
	{
	  messerror("addTreeColumns: error creating column - invalid column info in detail-view column-list.");
	}
    }
  
  return headerWidgets;
}


static gint sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data)
{
  gint result = UNSET_INT;

  /* Extract the sort column and sort order */
  gint sortColumn;
  GtkSortType sortOrder;
  gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model), &sortColumn, &sortOrder);
//  gboolean ascending = (sortOrder == GTK_SORT_ASCENDING);

  /* Extract the MSP lists from the tree rows */
  GList *mspGList1 = treeGetMsps(model, iter1);
  GList *mspGList2 = treeGetMsps(model, iter2);

  /* Get the first MSP in each list. */
  MSP *msp1 = (MSP*)(mspGList1->data);
  MSP *msp2 = (MSP*)(mspGList2->data);

  const gboolean multipleMsps = g_list_length(mspGList1) > 1 || g_list_length(mspGList2) > 1;

  gint sortType = (gint)data;
  
  switch (sortType)
    {
      case SORT_BY_NAME:
	{
	  MSP *msp1 = (MSP*)(mspGList1->data);
	  MSP *msp2 = (MSP*)(mspGList2->data);
	  result = strcmp(msp1->sname, msp2->sname);
	  break;
	}
	
      case SORT_BY_SCORE:
	{
	  result = multipleMsps ? 0 : msp1->score - msp2->score;
	  break;
	}
	
      case SORT_BY_ID:
	{
	  result = multipleMsps ? 0 : msp1->id - msp2->id;
	  break;
	}
	
      case SORT_BY_POS:
	{
	  if (multipleMsps)
	    {
	      result = 0;
	    }
	  else
	    {
	      /* Use the low end of the reference sequence range */
	      int qMin1, qMin2;
	      getMspRangeExtents(msp1, &qMin1, NULL, NULL, NULL);
	      getMspRangeExtents(msp2, &qMin2, NULL, NULL, NULL);
	      result = qMin1 - qMin2;
	    }
	  
	  break;
	}

      default:
	break;
    };
  
  /* If values are the same, further sort by name, then position, then length */
  if (!result && sortType != SORT_BY_NAME)
    {
      result = strcmp(msp1->sname, msp2->sname);
    }
  
  if (!multipleMsps && !result && sortType != SORT_BY_POS)
    {
      int sMin1, sMin2;
      getMspRangeExtents(msp1, NULL, NULL, &sMin1, NULL);
      getMspRangeExtents(msp2, NULL, NULL, &sMin2, NULL);
      result = sMin1 - sMin2;
    }
  
  if (!multipleMsps && !result)
    {
      result = msp1->slength - msp2->slength;
    }
  
  return result;
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
  
  /* Set the sort functions for each column */
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), S_NAME_COL, sortColumnCompareFunc, (gpointer)SORT_BY_NAME, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), ID_COL, sortColumnCompareFunc, (gpointer)SORT_BY_ID, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), SCORE_COL, sortColumnCompareFunc, (gpointer)SORT_BY_SCORE, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), SEQUENCE_COL, sortColumnCompareFunc, (gpointer)SORT_BY_POS, NULL);
  
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

  /* Keep a reference to the model in the properties so we can switch between this and the 'squashed' model */
  TreeProperties *properties = treeGetProperties(tree);
  properties->mspTreeModel = GTK_TREE_MODEL(filter);
}


static void setTreeStyle(GtkTreeView *tree)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_VERTICAL);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_MULTIPLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, FALSE);
  
  /* Set the background colour for the rows to be the same as the widget's background colour */
  gtk_widget_modify_base(GTK_WIDGET(tree), GTK_STATE_NORMAL, GTK_WIDGET(tree)->style->bg);

  /* The default text colour when rows are selected is white. This doesn't work 
   * well against our default background colour of cyan, so use the same text colour
   * as unselected rows. */
  gtk_widget_modify_text(GTK_WIDGET(tree), GTK_STATE_SELECTED, GTK_WIDGET(tree)->style->text);
  gtk_widget_modify_text(GTK_WIDGET(tree), GTK_STATE_ACTIVE, GTK_WIDGET(tree)->style->text);
  
  /* Set the expander size to 0 so that we can have tiny rows (otherwise the min is 12pt).
   * Also set the vertical separator to 0 so that we can have the option of the smallest
   * fonts possible. (The vertical separator causes gaps between rows. We want rows to be
   * flush, so we render over this gap. However, the bigger the vertical separator, the
   * bigger the row height (and hence font size) we must have to cover the gaps.) */
  char parseString[500];
  sprintf(parseString, "style \"packedTree\"\n"
	  "{\n"
	  "GtkTreeView::expander-size	      = 0\n"
	  "GtkTreeView::vertical-separator    = 0\n"
	  "GtkTreeView::horizontal-separator  = 0\n"
	  "}"
	  "widget \"*%s*\" style \"packedTree\"", DETAIL_VIEW_TREE_NAME);
  gtk_rc_parse_string(parseString);
}


/* Create the widget that will contain the header labels. Not much to do here
 * because the labels are added when the columns are added. */
static GtkWidget *createDetailViewTreeHeader()
{
  GtkWidget *header = gtk_hbox_new(FALSE, 0);
  return header;
}


GtkWidget* createDetailViewTree(GtkWidget *grid, 
				GtkWidget *detailView, 
				GtkCellRenderer *renderer,
				GList **treeList,
				GList *columnList,
				BlxSeqType seqType,
				const char const *refSeqName,
				const int frame)
{
  /* Create a tree view for the list of match sequences */
  GtkWidget *tree = gtk_tree_view_new();
  setTreeStyle(GTK_TREE_VIEW(tree));

  /* Put it in a scrolled window for vertical scrolling only (hoz scrolling will be via our
   * custom adjustment). Always display the scrollbars because we assume the column widths 
   * are the same for all trees and they won't be if one shows a scrollbar and another doesn't. */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_widget_set_name(scrollWin, DETAIL_VIEW_TREE_CONTAINER_NAME);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
  gtk_container_add(GTK_CONTAINER(scrollWin), tree);

  /* Create a header, and put the tree and header in a vbox */
  GtkWidget *treeHeader = createDetailViewTreeHeader();
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_widget_set_name(vbox, DETAIL_VIEW_TREE_CONTAINER_NAME);
  gtk_box_pack_start(GTK_BOX(vbox), treeHeader, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), scrollWin, TRUE, TRUE, 0);
  
  /* Add the columns */
  GList *treeColumnHeaderList = addTreeColumns(tree, renderer, seqType, columnList, treeHeader, refSeqName, frame, gridGetStrand(grid));
  
  /* Set the essential tree properties */
  treeCreateProperties(tree, grid, detailView, renderer, frame, treeColumnHeaderList);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);
  g_signal_connect(G_OBJECT(tree), "button-press-event",    G_CALLBACK(onButtonPressTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "button-release-event",  G_CALLBACK(onButtonReleaseTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "key-press-event",	    G_CALLBACK(onKeyPressTree),	detailView);
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
  
  /* Add the tree's outermost container to the given tree list. Also increase its ref count so that
   * we can add/remove it from its parent (which we do to switch panes when we toggle strands)
   * without worrying about it being destroyed. */
  *treeList = g_list_append(*treeList, vbox);
  g_object_ref(vbox);
  
  return vbox;
}

