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
#include <SeqTools/blxwindow.h>
#include <SeqTools/utilities.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>


#define DETAIL_VIEW_STATUSBAR_CONTEXT     "statusBarCtx"


/* Local function declarations */
static GtkWidget*	treeGetDetailView(GtkWidget *tree);
static gboolean		onSelectionChangedTree(GObject *selection, gpointer data);
static gint		sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data);
static int		calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree);
static gboolean		isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data);
static GtkSortType	treeGetColumnSortOrder(GtkWidget *tree, const BlxColumnId columnId);
static int		scrollBarWidth();
static gboolean		onExposeRefSeqHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);
static GList*		treeGetSequenceRows(GtkWidget *tree, const BlxSequence *clickedSeq);
static BlxSequence*	treeGetSequence(GtkTreeModel *model, GtkTreeIter *iter);

/***********************************************************
 *                Tree - utility functions                 *
 ***********************************************************/

static void assertTree(GtkWidget *tree)
{
  if (!tree)
    g_error("Tree is null\n");
  
  if (!GTK_IS_WIDGET(tree))
    g_error("Tree is not a valid widget [%p]\n", tree);
    
  if (!GTK_IS_TREE_VIEW(tree))
    g_error("Tree is not a valid tree view [%p]\n", tree);

  if (strcmp(gtk_widget_get_name(tree), DETAIL_VIEW_TREE_NAME) != 0)
    g_error("Tree is not a valid detail-view tree [%p]\n", tree);

  if (!treeGetProperties(tree))
    g_error("Tree properties not set [widget=%p]\n", tree);
}


static GtkWidget *treeGetDetailView(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? properties->detailView : NULL;
}

GtkWidget *treeGetBlxWindow(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetBlxWindow(detailView);
}

BlxStrand treeGetStrand(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return gridGetStrand(properties->grid);
}

static gboolean treeGetDisplayRev(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetDisplayRev(detailView);
}

static int treeGetFrame(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties->readingFrame;
}

static gboolean treeHasSnpHeader(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties->hasSnpHeader;
}

static int treeGetSelectedBaseIdx(GtkWidget *tree)
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
	  detailViewSetSelectedBaseIdx(properties->detailView, selectedBaseIdx, frame, baseNum, allowScroll, TRUE);
	}
    }
}

/* Returns the currently-displayed range in the tree view */
static IntRange* treeGetDisplayRange(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetDisplayRange(properties->detailView) : NULL;
}

static PangoFontDescription* treeGetFontDesc(GtkWidget *tree)
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

static int treeGetCharWidth(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCharWidth(detailView);
}

static TreeColumnHeaderInfo* treeColumnGetHeaderInfo(GtkWidget *tree, BlxColumnId columnId)
{
  /* Loop through each headerinfo item until we find the one that contains this column id */
  TreeProperties *properties = treeGetProperties(tree);
  GList *headerInfoItem = properties->treeColumnHeaderList;
  TreeColumnHeaderInfo *result = NULL;
  
  for ( ; headerInfoItem && !result; headerInfoItem = headerInfoItem->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)(headerInfoItem->data);
      
      /* Loop through each column id included under this header */
      GList *columnIdItem = headerInfo->columnIds;
      
      for ( ; columnIdItem && !result; columnIdItem = columnIdItem->next)
	{
	  if (GPOINTER_TO_INT(columnIdItem->data) == columnId)
	    {
	      result = headerInfo;
	    }
	}
    }
  
  return result;
}

/***********************************************************
 *			utility functions		   *
 ***********************************************************/

/* Call the given function on all trees in the detail view */
void callFuncOnAllDetailViewTrees(GtkWidget *detailView, GtkCallback func, gpointer data)
{
  int numFrames = detailViewGetNumFrames(detailView);
  
  /* Call the function on the forward strand tree and reverse strand tree
   * for each frame. */
  int frame = 1;
  for ( ; frame <= numFrames; ++frame)
    {
      GtkWidget *fwdTree = detailViewGetTree(detailView, BLXSTRAND_FORWARD, frame);
      GtkWidget *revTree = detailViewGetTree(detailView, BLXSTRAND_REVERSE, frame);
      
      if (fwdTree)
	{
	  func(fwdTree, data);
	}
      
      if (revTree)
	{
	  func(revTree, data);
	}
    }
}


/* Add a BlxSequence to as a row in the given tree store */
static void addSequenceStructToRow(gpointer listItemData, gpointer data)
{
  BlxSequence *subjectSeq = (BlxSequence*)listItemData;
  
  if (subjectSeq && subjectSeq->mspList && g_list_length(subjectSeq->mspList) > 0)
    {
      /* Add a row to the tree store */
      GtkListStore *store = GTK_LIST_STORE(data);
      GtkTreeIter iter;
      gtk_list_store_append(store, &iter);
      
      if (g_list_length(subjectSeq->mspList) == 1)
	{
	  /* Only one MSP - get specific info about this MSP. */
	  MSP *msp = (MSP*)(subjectSeq->mspList->data);

	  gtk_list_store_set(store, &iter,
			     BLXCOL_SEQNAME, mspGetSName(msp),
			     BLXCOL_SOURCE, msp->source,
			     BLXCOL_ORGANISM, NULL,
			     BLXCOL_GENE_NAME, NULL,
			     BLXCOL_TISSUE_TYPE, NULL,
			     BLXCOL_STRAIN, NULL,
			     BLXCOL_GROUP, NULL,
			     BLXCOL_SCORE, msp->score,
			     BLXCOL_ID, msp->id,
			     BLXCOL_START, msp->sRange.min,
			     BLXCOL_SEQUENCE, subjectSeq->mspList,
			     BLXCOL_END, msp->sRange.max,
			     -1);
	}
      else
	{
	  /* Add generic info about the sequence */
	  gtk_list_store_set(store, &iter,
			     BLXCOL_SEQNAME, blxSequenceGetDisplayName(subjectSeq),
			     BLXCOL_SOURCE, NULL,
			     BLXCOL_ORGANISM, NULL,
			     BLXCOL_GENE_NAME, NULL,
			     BLXCOL_TISSUE_TYPE, NULL,
			     BLXCOL_STRAIN, NULL,
			     BLXCOL_GROUP, NULL,
			     BLXCOL_SCORE, 0.0,
			     BLXCOL_ID, 0.0,
			     BLXCOL_START, blxSequenceGetStart(subjectSeq),
			     BLXCOL_SEQUENCE, subjectSeq->mspList,
			     BLXCOL_END, blxSequenceGetEnd(subjectSeq),
			     -1);
	}
    }
}


void addSequencesToTree(GtkWidget *tree, gpointer data)
{
  GtkListStore *store = gtk_list_store_new(BLXCOL_NUM_COLUMNS, TREE_COLUMN_TYPE_LIST);
  
  /* Set the sort function for each column */
  int colNum = 0;
  for ( ; colNum < BLXCOL_NUM_COLUMNS; ++colNum)
    {
      gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), colNum,	sortColumnCompareFunc, tree, NULL);
    }

  /* Add the rows - one row per sequence. Use the list we've already compiled of all
   * sequences as BlxSequences */
  GList *seqList = blxWindowGetAllMatchSeqs(treeGetBlxWindow(tree));
  g_list_foreach(seqList, addSequenceStructToRow, store);

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
  gtk_tree_model_get(model, iter, BLXCOL_SEQUENCE, &mspGList, -1);
  return mspGList;
}


/* For the given tree view, return the data model used for the visible view
 * (which may be a filtered tree model). */
static GtkTreeModel* treeGetVisibleDataModel(GtkTreeView *tree)
{
  /* Just return the model stored directly in the tree view. This will be the filter
   * if there is one, otherwise it will be the original data model. */
  assertTree(GTK_WIDGET(tree));
  return gtk_tree_view_get_model(tree);
}

/* Get the status bar that this tree should report details about moused-over rows to.
 * This may be different to the main window status bar. */
static GtkStatusbar* treeGetStatusBar(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  return GTK_STATUSBAR(dvProperties->statusBar);
  
  //GTK_STATUSBAR(treeGetContext(tree)->statusBar);
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


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
void treeUpdateFontSize(GtkWidget *tree, gpointer data)
{
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);  
  gtk_widget_modify_font(tree, fontDesc);
}


/* This function does the work to squash/unsquash tree rows */
static void treeSetSquashMatches(GtkWidget *tree, const gboolean squash)
{
  /* Get the current model (and find which column it is currently sorted by) */
  gint sortColumn;
  GtkTreeModel *origModel = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
  gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(origModel), &sortColumn, NULL);
  
  /* Replace it with the squashed/unsquashed model as requested (if not already set to that) */
  TreeProperties *properties = treeGetProperties(tree);
  gboolean changed = FALSE;
  
  if (squash && origModel != properties->seqTreeModel)
    {
      gtk_tree_view_set_model(GTK_TREE_VIEW(tree), properties->seqTreeModel);
      changed = TRUE;
    }
  else if (!squash && origModel != properties->mspTreeModel)
    {
      gtk_tree_view_set_model(GTK_TREE_VIEW(tree), properties->mspTreeModel);
      changed = TRUE;
    }
  
  /* If we changed the model, re-sort and re-filter it */
  if (changed)
    {  
      /* We sort the base data model, not the filtered one (i.e. sort all rows, not just visible ones) */
      GtkTreeModel *newModel = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      GtkSortType sortOrder = treeGetColumnSortOrder(tree, sortColumn);
      gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(newModel), sortColumn, sortOrder);
      
      /* Re-filter, because this model may not be up to date for the current display range */
      refilterTree(tree, NULL);
    }
}


/* This function makes multiple MSPs in the same sequence appear in the same row in the tree. */
void treeSquashMatches(GtkWidget *tree, gpointer data)
{
  treeSetSquashMatches(tree, TRUE);
}


/* This function makes all MSPs appear in their own individual rows in the tree. */
void treeUnsquashMatches(GtkWidget *tree, gpointer data)
{
  treeSetSquashMatches(tree, FALSE);
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
      g_critical("Unexpected tree data store [%p]. Expected either [%p] (normal data) or [%p] (condensed data).\n", model, properties->mspTreeModel, properties->seqTreeModel);
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

  /* Loop through all widgets in the header and call refreshTextHeader. This
   * updates the font etc. if it is a type of widget that requires that. */
  gtk_container_foreach(GTK_CONTAINER(properties->treeHeader), refreshTextHeader, treeGetDetailView(tree));
  
  /* Loop through all column headers and call their individual refresh functions
   * (this sets the specific data for the column headers) */
  GList *headerItem = properties->treeColumnHeaderList;
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);
  
  for ( ; headerItem; headerItem = headerItem->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)headerItem->data;
      
      if (headerInfo && headerInfo->headerWidget)
	{
	  /* Set the background color (in the parent, seeing as the label doesn't have a window) */
	  BlxViewContext *bc = treeGetContext(tree);
	  GtkWidget *parent = gtk_widget_get_parent(headerInfo->headerWidget);
	  GdkColor *bgColor = getGdkColor(BLXCOLOR_REF_SEQ, bc->defaultColors, FALSE, bc->usePrintColors);
	  gtk_widget_modify_bg(parent, GTK_STATE_NORMAL, bgColor);

	  /* Update the font, in case its size has changed */
	  gtk_widget_modify_font(headerInfo->headerWidget, fontDesc);

	  /* Call the refresh function, if there is one */
	  if (headerInfo->refreshFunc)
	    {
	      headerInfo->refreshFunc(headerInfo->headerWidget, tree);
	    }
	}
    }
}


/* Resize tree header widgets (should be called after the column width has changed) */
static void resizeTreeHeaders(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  GList *header = properties->treeColumnHeaderList;
  
  for ( ; header; header = header->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)header->data;
      
      if (headerInfo && headerInfo->headerWidget && g_list_length(headerInfo->columnIds) > 0)
	{
	  int firstColId = GPOINTER_TO_INT(headerInfo->columnIds->data);

          /* For the sequence column, don't set the size request to the real size, or we 
           * won't be able to shrink the window. (The sequence col header will be resized
           * dynamically anyway.) */
	  if (firstColId != BLXCOL_SEQUENCE)
	    {
	      /* For other columns, we can set the size request: they're small enough
	       * that we can live without the need to shrink below their sum widths. */
	      const int width = calculateColumnWidth(headerInfo, tree);
	      gtk_widget_set_size_request(headerInfo->headerWidget, width, -1);
	    }
	}
    }
  
  refreshTreeHeaders(tree, NULL);
}


/* Resize the columns for this tree. Should be called after the column width
 * is changed manually. */
void resizeTreeColumns(GtkWidget *tree, gpointer data)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  
  GList *listItem = detailViewGetColumnList(detailView);
  int sumWidth = 0;
  
  for ( ; listItem; listItem = listItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(listItem->data);
      
      if (columnInfo->columnId != BLXCOL_SEQUENCE)
	{
	  GtkTreeViewColumn *treeColumn = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), columnInfo->columnId);
      
	  int width = columnInfo->width;
	  if (columnInfo->columnId == BLXCOL_END)
	    {
	      width -= scrollBarWidth();
	    }

	  sumWidth += width;

	  if (width > 0)
	    {
	      gtk_tree_view_column_set_visible(treeColumn, TRUE);
	      gtk_tree_view_column_set_fixed_width(treeColumn, width);
	    }
	  else
	    {
	      /* Can't have 0 width, so hide the column instead */
	      gtk_tree_view_column_set_visible(treeColumn, FALSE);
	    }
	}
    }

  /* The sequence column width is the full width minus the other columns' widths. 
   * (Can't use gtk_tree_view_column_get_width here because it is not up to date yet) */
  DetailViewColumnInfo *seqColumnInfo = detailViewGetColumnInfo(detailView, BLXCOL_SEQUENCE);
  seqColumnInfo->width = tree->allocation.width - sumWidth;
  

  updateSeqColumnSize(treeGetDetailView(tree));
  resizeTreeHeaders(tree, NULL);
}


/***********************************************************
 *			  Selections			   *
 ***********************************************************/


/* Scrolls the currently-selected row(s) into view */
void treeScrollSelectionIntoView(GtkWidget *tree, gpointer data)
{
  /* Get the last selected sequence */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  const BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
  
  /* Get the tree row(s) that contain that sequence. (If multiple, just use 1st one) */
  GList *rows = treeGetSequenceRows(tree, lastSelectedSeq);
  
  if (g_list_length(rows) > 0)
    {
      GtkTreePath *path = (GtkTreePath*)(rows->data);
      gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(tree), path, NULL, FALSE, 0.0, 0.0);
    }
}


/* We handle row selections ourselves */
static gboolean onSelectionChangedTree(GObject *selection, gpointer data)
{
  return TRUE;
}


/***********************************************************
 *                   Sorting and filtering                 *
 ***********************************************************/

/* Refilter the data for the given tree */
void refilterTree(GtkWidget *tree, gpointer data)
{
  assertTree(tree);
  
  GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(tree));
  gtk_tree_model_filter_refilter(GTK_TREE_MODEL_FILTER(model));
}


/* Determine the sort order for this column. Columns may have a different default
 * sort order (e.g. name is sorted ascending, but score is sorted descending.)
 * However, the opposite sort order is returned if the detail view's invert-sort-
 * order flag is set. */
static GtkSortType treeGetColumnSortOrder(GtkWidget *tree, const BlxColumnId columnId)
{
  GtkSortType result = GTK_SORT_ASCENDING;
  BlxViewContext *bc = treeGetContext(tree);

  switch (columnId)
    {
      case BLXCOL_SCORE: /* fall through */
      case BLXCOL_ID:
	result = GTK_SORT_DESCENDING;
	break;
	
      default: /* all others ascending */
	break;
    };

  if (blxContextGetFlag(bc, BLXFLAG_INVERT_SORT))
    {
      result = (result == GTK_SORT_ASCENDING) ? GTK_SORT_DESCENDING : GTK_SORT_ASCENDING;
    }

  return result;
}


/* Re-sort the data for the given tree */
void resortTree(GtkWidget *tree, gpointer data)
{
  /* Not sure if there's a better way to do this, but we can force a re-sort by 
   * setting the sort column to something else and then back again. */
  gint sortColumn;
  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
  gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model), &sortColumn, NULL);
  
  /* Set it to unsorted */
  GtkSortType sortOrder = treeGetColumnSortOrder(tree, sortColumn);
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID, sortOrder);
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), sortColumn, sortOrder);
}



static void mspGetVisibleRange(const MSP const *msp,
			       const BlxViewContext *bc,
			       const gboolean showUnalignedSeq,
			       const gboolean limitUnalignedBases,
			       const int numUnalignedBases,
			       IntRange *qRange_out, 
			       IntRange *sRange_out)
{
  if (showUnalignedSeq && mspIsBlastMatch(msp))
    {
      /* We need to include some/all unaligned parts of the match sequence too. 
       * First, get the full range of the match sequence. */
      int sMin = 1;
      int sMax = mspGetMatchSeqLen(msp);
      
      if (limitUnalignedBases)
        {
          /* Limit the range to the start/end of the alignment but with the given number of additional 
           * bases. Make sure we don't go outside the range the full match sequence range, though. */
          sMin = max(sMin, msp->sRange.min - numUnalignedBases);
          sMax = min(sMax, msp->sRange.max + numUnalignedBases);
        }
      
      if (sRange_out)
        {
          sRange_out->min = sMin;
          sRange_out->max = sMax;
        }
      
      if (qRange_out)
        {
          /* Find out how much the new s coords are offset from the start/end of the MSP range */
          int startOffset = msp->sRange.min - sMin;
          int endOffset = sMax - msp->sRange.max;

          /* Get the q coords. The s coords are in terms of display coords so we need the q coords in 
           * display coords too. Note that the conversion may invert the coords, so we have to work out 
           * which is the max/min again. */
          const int qFrame = mspGetRefFrame(msp, bc->seqType);
          
          const int qIdx1 = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
          const int qIdx2 = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
          
          int qMin = min(qIdx1, qIdx2);
          int qMax = max(qIdx1, qIdx2);
          
          /* Adjust the q coords by the offset. If the strands are in opposite directions, or if
           * the display is reversed (which inverts the q coords), then we need to apply the offset
           * at the opposite end of the reference sequence. */
          const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
          
          if (sameDirection != bc->displayRev)
            {
              qMin -= startOffset;
              qMax += endOffset;
            }
          else
            {
              qMin -= endOffset;
              qMax += startOffset;
            }
          
          qRange_out->min = qMin;
          qRange_out->max = qMax;
        }
    }
  else
    {
      /* Just return the coords */
      if (qRange_out)
        {
          const int qFrame = mspGetRefFrame(msp, bc->seqType);

          const int qIdx1 = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
          const int qIdx2 = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
          
          qRange_out->min = min(qIdx1, qIdx2);
          qRange_out->max = max(qIdx1, qIdx2);
        }
      
      if (sRange_out)
        {
          sRange_out->min = msp->sRange.min;
          sRange_out->max = msp->sRange.max;
        }
    }
}


/* Utility that returns true if the given MSP is currently shown in the tree with the given
 * strand/frame */
static gboolean isMspVisible(const MSP const *msp, 
			     GtkWidget *blxWindow, 
			     const BlxViewContext *bc, 
			     const BlxStrand strand, 
			     const int frame, 
			     const IntRange const *displayRange,
			     const int numUnalignedBases)
{
  /* Check if the MSP is in a visible layer */
  gboolean result = mspLayerIsVisible(msp);

  /* The tree view only displays blast matches or exons (or polyA tails, if the MSP's sequence is 
   * selected) */
  result &= mspIsBlastMatch(msp) || mspIsExon(msp) ||
		    (mspIsPolyATail(msp) && blxWindowIsSeqSelected(blxWindow, msp->sSequence));
		    
  /* Check that it is in this tree's frame and strand */
  result &= (mspGetRefStrand(msp) == strand);
  result &= (mspGetRefFrame(msp, bc->seqType) == frame);
  
  /* Check it's in the current display range */
  if (result)
    {
      IntRange qRange;
      const gboolean showUnalignedSeq = blxContextGetFlag(bc, BLXFLAG_SHOW_UNALIGNED_SEQ);
      const gboolean limitUnalignedBases = blxContextGetFlag(bc, BLXFLAG_LIMIT_UNALIGNED_BASES);
    
      mspGetVisibleRange(msp, bc, showUnalignedSeq, limitUnalignedBases, numUnalignedBases, &qRange, NULL);

      result &= (qRange.min <= displayRange->max && qRange.max >= displayRange->min);
    }
    
  return result;
}


/* Filter function. Returns true if the given row in the given tree model should be visible. */
static gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  gboolean bDisplay = FALSE;
  
  GList *mspList = treeGetMsps(model, iter);
  
  if (g_list_length(mspList) > 0)
    {
      GtkWidget *tree = GTK_WIDGET(data);
      GtkWidget *blxWindow = treeGetBlxWindow(tree);

      /* Check the first msp to see if this sequence is in a group that's hidden */
      const MSP *firstMsp = (const MSP*)(mspList->data);
      SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, firstMsp->sSequence);
      
      if (!group || !group->hidden)
	{
	  BlxViewContext *bc = treeGetContext(tree);
          GtkWidget *detailView = treeGetDetailView(tree);
          DetailViewProperties *dvProperties = detailViewGetProperties(detailView);

	  const int frame = treeGetFrame(tree);
	  const BlxStrand strand = treeGetStrand(tree);
	  const IntRange const *displayRange = treeGetDisplayRange(tree);

	  /* Show the row if any MSP in the list is an exon or blast match in the correct frame/strand
	   * and within the display range */
	  GList *mspListItem = mspList;
	
	  for ( ; mspListItem; mspListItem = mspListItem->next)
	    {
	      const MSP* msp = (const MSP*)(mspListItem->data);
	      
	      if (isMspVisible(msp, blxWindow, bc, strand, frame, displayRange, dvProperties->numUnalignedBases))
		{
		  bDisplay = TRUE;
		  break;
		}
	    }	  
	}
    }
  
  return bDisplay;
}


void treeSetSortColumn(GtkWidget *tree, gpointer data)
{
  GtkTreeSortable *model = GTK_TREE_SORTABLE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
  
  BlxColumnId sortColumn = GPOINTER_TO_INT(data);
  GtkSortType sortOrder = treeGetColumnSortOrder(tree, sortColumn);
  
  gtk_tree_sortable_set_sort_column_id(model, sortColumn, sortOrder);
}

/***********************************************************
 *                       Properties                        *
 ***********************************************************/

TreeProperties* treeGetProperties(GtkWidget *widget)
{
  return widget ? (TreeProperties*)(g_object_get_data(G_OBJECT(widget), "TreeProperties")) : NULL;
}

BlxViewContext* treeGetContext(GtkWidget *tree)
{
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  return blxWindowGetContext(blxWindow);
}

static void onDestroyTree(GtkWidget *widget)
{
  TreeProperties *properties = treeGetProperties(widget);
  
  if (properties)
    {
      if (properties->treeColumnHeaderList)
	{
	  g_list_free(properties->treeColumnHeaderList);
	  properties->treeColumnHeaderList = NULL;
	}
    
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "TreeProperties", NULL);
    }
}

static void treeCreateProperties(GtkWidget *widget, 
				 GtkWidget *grid, 
				 GtkWidget *detailView, 
				 const int frame,
				 GtkWidget *treeHeader,
				 GList *treeColumnHeaderList,
				 const gboolean hasSnpHeader)
{
  if (widget)
    { 
      TreeProperties *properties = g_malloc(sizeof *properties);
      
      properties->grid = grid;
      properties->detailView = detailView;
      properties->readingFrame = frame;
      properties->treeHeader = treeHeader;
      properties->treeColumnHeaderList = treeColumnHeaderList;
      properties->hasSnpHeader = hasSnpHeader;
      properties->mspTreeModel = NULL;
      properties->seqTreeModel = NULL;

      g_object_set_data(G_OBJECT(widget), "TreeProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyTree), NULL); 
    }
}


/***********************************************************
 *                       Tree events                       *
 ***********************************************************/

/* Handler for when vertical scrollbar for the tree has changed (either value or range) */
static void onScrollChangedTree(GtkObject *object, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  /* Remove any message from the detail-view status bar, because this shows info for the
   * currently-moused-over item, which may have changed now we've scrolled. */
   GtkStatusbar *statusBar = treeGetStatusBar(tree);
   
   if (statusBar)
     {
       guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
       gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
     }
}

/* Find the position of the left edge of the source widget wrt the dest widget and offset 
 * the event by this amount, then propagate the event to the dest widget. This allows us
 * to get the correct x coord on the parent widget for a click on the child widget. */
static void propagateEventButton(GtkWidget *srcWidget, GtkWidget *destWidget, GdkEventButton *event)
{
  int xOffset;
  gtk_widget_translate_coordinates(srcWidget, destWidget, 0, 0, &xOffset, NULL);
  
  event->x += xOffset;
  
  gtk_propagate_event(destWidget, (GdkEvent*)event);
}


/* Utility to translate the coords in the given event from coords relative to the source
 * widget to coords relative to the destination widget, and then propagate the event to the
 * dest widget. */
static void propagateEventMotion(GtkWidget *srcWidget, GtkWidget *destWidget, GdkEventMotion *event)
{
  /* Find the position of 0,0 of the source widget wrt the dest widget and offset by this amount. */
  int xOffset, yOffset;
  gtk_widget_translate_coordinates(srcWidget, destWidget, 0, 0, &xOffset, &yOffset);
  
  event->x += xOffset;
  event->y += yOffset;  
  
  gtk_propagate_event(destWidget, (GdkEvent*)event);
}


/* Draw the part of the tree header that shows the reference sequence (either as
 * nucleotide sequence, or a peptide sequence if viewing protein matches) */
static void drawRefSeqHeader(GtkWidget *headerWidget, GtkWidget *tree)
{
  GdkDrawable *drawable = createBlankPixmap(headerWidget);
  
  BlxViewContext *bc = treeGetContext(tree);
  const BlxStrand strand = treeGetStrand(tree);
  const int frame = treeGetFrame(tree);
  
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  const gboolean showSnps = treeHasSnpHeader(tree);
  
  /* Find the segment of the ref seq to display. */
  GError *error = NULL;
  
  gchar *segmentToDisplay = getSequenceSegment(bc,
					       bc->refSeq,
					       properties->displayRange.min, 
					       properties->displayRange.max, 
					       strand, 
					       bc->seqType,
					       frame, 
					       bc->displayRev,
					       bc->displayRev,	/* show backwards if display reversed */
					       TRUE,		/* always complement reverse strand */
					       TRUE,		/* always translate peptide sequences */
					       &error);
  
  if (!segmentToDisplay)
    {
      g_assert(error);
      prefixError(error, "Could not draw reference sequence header. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
  
  GdkGC *gc = gdk_gc_new(drawable);
  
  /* Find out if there are any bases in the introns that need highlighting. */
  const int qIdx1 = convertDisplayIdxToDnaIdx(properties->displayRange.min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int qIdx2 = convertDisplayIdxToDnaIdx(properties->displayRange.max, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  IntRange qRange = {min(qIdx1, qIdx2), max(qIdx1, qIdx2)};
  
  GHashTable *intronBases = getIntronBasesToHighlight(detailView, &qRange, bc->seqType, strand);

  const int incrementValue = bc->displayRev ? -1 * bc->numFrames : bc->numFrames;
  int displayIdx = properties->displayRange.min;
  int dnaIdx = qIdx1;
  
  while (displayIdx >= properties->displayRange.min && displayIdx <= properties->displayRange.max)
    {
      /* Set the background color depending on whether this base is selected or
       * is affected by a SNP */
      const gboolean displayIdxSelected = (displayIdx == properties->selectedBaseIdx);
      const char baseChar = segmentToDisplay[displayIdx - properties->displayRange.min];
      
      const int x = (displayIdx - properties->displayRange.min) * properties->charWidth;
      const int y = 0;

      drawHeaderChar(bc, properties, dnaIdx, baseChar, strand, frame, bc->seqType, displayIdxSelected, displayIdxSelected, TRUE, showSnps, FALSE, BLXCOLOR_REF_SEQ, drawable, gc, x, y, intronBases);
      
      dnaIdx += incrementValue;
      ++displayIdx;
    }
  
  /* Mark up the text to highlight the selected base, if there is one */
  PangoLayout *layout = gtk_widget_create_pango_layout(detailView, segmentToDisplay);
  pango_layout_set_font_description(layout, detailViewGetFontDesc(detailView));
  
  if (layout)
    {
      gtk_paint_layout(headerWidget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, 0, 0, layout);
      g_object_unref(layout);
    }
  
  drawColumnSeparatorLine(headerWidget, drawable, gc, bc);
  
  g_free(segmentToDisplay);
  g_object_unref(gc);
}



/* expose function to push a cached bitmap to screen */
static gboolean onExposeRefSeqHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkWindow *window = GTK_IS_LAYOUT(headerWidget) ? GTK_LAYOUT(headerWidget)->bin_window : headerWidget->window;
  
  if (window)
    {
      /* See if there's a cached drawable and, if not, create it */
      GdkDrawable *bitmap = widgetGetDrawable(headerWidget);
      
      if (!bitmap)
	{
	  /* There isn't a bitmap yet. Create it now. */
	  GtkWidget *tree = GTK_WIDGET(data);
	  drawRefSeqHeader(headerWidget, tree);
	  bitmap = widgetGetDrawable(headerWidget);
	}
      
      if (bitmap)
	{
	  /* Push the bitmap onto the window */
	  GdkGC *gc2 = gdk_gc_new(window);
	  gdk_draw_drawable(window, gc2, bitmap, 0, 0, 0, 0, -1, -1);
	}
    }
  
  return FALSE;
}



static gboolean onExposeDetailViewTree(GtkWidget *tree, GdkEventExpose *event, gpointer data)
{
  /* Create a new drawable to draw to. Our custom cell renderer will draw to this as
   * well as the widget's window. (Ideally we'd just draw to the bitmap and then push
   * this to the screen, but I'm not sure if it's possible to detect when the 
   * cell renderer has finished drawing.) */
  GdkDrawable *drawable = gdk_pixmap_new(tree->window, tree->allocation.width, tree->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(tree, drawable);

  /* Draw a blank rectangle of the required widget background color */
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *bgColor = tree->style->bg;
  gdk_gc_set_foreground(gc, bgColor);
  
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, tree->allocation.width, tree->allocation.height);
  
  /* Let the default handler continue */
  return FALSE;
}


/* Select the range of rows in the given tree from the first path to the last inclusive */
static void treeSelectRowRange(GtkWidget *blxWindow, GtkTreeModel *model, GtkTreePath *firstPath, GtkTreePath *lastPath)
{
  /* We currently only ever have a list, so this function assumes the model is a list. If 
   * the alignment lists are ever changed to be trees, this function will need updating. */
  if (GTK_IS_TREE_STORE(model))
    {
      printf("Error selecting range of rows: function not implemented (expected alignments to be in a list store, not a tree store).");
      return;
    }
  
  /* Loop through all the rows from the last selected one to the one that was clicked
   * on. The direction depends on which one appears first in the list. */
  gint *firstIdx = gtk_tree_path_get_indices(firstPath);
  gint *lastIdx = gtk_tree_path_get_indices(lastPath);
  
  gint currentIdx = *firstIdx;
  const gint incrementValue = *firstIdx < *lastIdx ? 1 : -1;
  gboolean done = FALSE;
  
  while (!done)
    {
      /* Select the sequence in the current row */
      GtkTreePath *currentPath = gtk_tree_path_new_from_indices(currentIdx, -1);
      
      GtkTreeIter currentIter;
      if (gtk_tree_model_get_iter(model, &currentIter, currentPath))
	{
	  BlxSequence *currentSeq = treeGetSequence(model, &currentIter);
	  blxWindowSelectSeq(blxWindow, currentSeq);
	}
      else
	{
	  printf("Warning: invalid iterator found when selecting a range of rows\n");
	  done = TRUE;
	}
      
      /* Get the next path, unless we're already at the last one */
      if (currentIdx == *lastIdx)
	{
	  done = TRUE;
	}
      else
	{
	  currentIdx += incrementValue;
	}
    }
}


static gboolean treeSelectRow(GtkWidget *tree, GdkEventButton *event)
{
  guint modifiers = gtk_accelerator_get_default_mod_mask();
  const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
  const gboolean shiftModifier = ((event->state & modifiers) == GDK_SHIFT_MASK);
  
  if (!ctrlModifier && !shiftModifier)
    {
      blxWindowDeselectAllSeqs(treeGetBlxWindow(tree));
    }
  
  /* If we've clicked a different tree to the previously selected one, make this tree's 
   * reading frame the active one. Select the first base in the triplet */
  treeSetSelectedBaseIdx(tree, treeGetSelectedBaseIdx(tree), treeGetFrame(tree), 1, FALSE);
  
  /* Find which row was clicked */
  GtkTreePath *clickedPath = NULL;
  GtkTreeViewColumn *clickedCol = NULL;
  int cell_x, cell_y;
  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), event->x, event->y, &clickedPath, &clickedCol, &cell_x, &cell_y);

  if (clickedPath)
    {
      GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
      GtkTreeIter clickedIter;
      gtk_tree_model_get_iter(model, &clickedIter, clickedPath);
      
      /* Get the sequence of the MSP(s) in this row */
      BlxSequence *clickedSeq = treeGetSequence(model, &clickedIter);

      if (clickedSeq)
	{
	  GtkWidget *blxWindow = treeGetBlxWindow(tree);
	  
	  if (!ctrlModifier && !shiftModifier)
	    {
	      /* No modifiers: select the sequence */
	      blxWindowSelectSeq(blxWindow, clickedSeq);
	    }
	  else if (ctrlModifier && !shiftModifier)
	    {
	      /* Ctrl pressed: toggle the selection state of the sequence */
	      blxWindowSetSeqSelected(blxWindow, clickedSeq, !blxWindowIsSeqSelected(blxWindow, clickedSeq));
	    }
	  else if (!ctrlModifier && shiftModifier)
	    {
	      /* Shift pressed: select all rows between the last-selected sequence and the clicked row */
	      const BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
	      GList *lastSelectedRows = treeGetSequenceRows(tree, lastSelectedSeq);
	      
	      if (g_list_length(lastSelectedRows) > 0) /* do nothing if no rows were previously selected */
		{
		  /* Get the tree row for the last selected sequence. (If there are multiple rows for this
		   * sequence, the desired behaviour here is a bit ambiguous, so for now just use any of them */
		  GtkTreePath *lastSelectedPath = (GtkTreePath*)(lastSelectedRows->data);
		  treeSelectRowRange(blxWindow, model, lastSelectedPath, clickedPath);
		}
	    }
	}

      GtkWidget *detailView = treeGetDetailView(tree);
      detailViewSetSelectedStrand(detailView, treeGetStrand(tree));
      detailViewRedrawAll(detailView);
    }
  
  return TRUE; /* handled */
}


static gboolean treePfetchRow(GtkWidget *tree)
{
  /* Get the selected sequence (assumes that only one is selected) */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (selectedSeqs)
    {
      const BlxSequence *clickedSeq = (const BlxSequence*)selectedSeqs->data;
      char *seqName = clickedSeq->fullName;
      
#ifdef ACEDB
      fetchAndDisplaySequence(seqName, 0, blxWindow);
#else
      fetchAndDisplaySequence(seqName, blxWindow);
#endif
    }

  return TRUE;
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
	    handled = treeSelectRow(tree, event);
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
	    handled = treePfetchRow(tree);
	  }

	break;
      }
      
    default:
      break;
    };
  
  if (!handled)
    { 
      propagateEventButton(tree, treeGetDetailView(tree), event);
    }
  
  return handled;
}


static gboolean onButtonPressTreeHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  GtkWidget *tree = GTK_WIDGET(data);

  switch (event->button)
    {
      case 1:
      {
	GtkWidget *detailView = treeGetDetailView(tree);
	
	if (event->type == GDK_BUTTON_PRESS)
	  {
	    /* Select the SNP that was clicked on.  */
	    blxWindowDeselectAllSeqs(detailViewGetBlxWindow(detailView));
	    int clickedBase = 1; /* only get in here for DNA matches; so frame/base number is always one */

	    selectClickedSnp(header, NULL, detailView, event->x, event->y, FALSE, FALSE, clickedBase); /* SNPs are always un-expanded in the DNA track */
	    
	    refreshDetailViewHeaders(detailView);
	    callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
	    BlxViewContext *bc = treeGetContext(tree);
            gboolean showSnpTrack = !blxContextGetFlag(bc, BLXFLAG_SHOW_SNP_TRACK);
            
	    blxContextSetFlag(bc, BLXFLAG_SHOW_SNP_TRACK, showSnpTrack);
            detailViewUpdateShowSnpTrack(detailView, showSnpTrack);
	  }
	
	handled = TRUE;
	break;
      }
      
    default:
      break;
    };
  
  if (!handled)
    {
      /* Propagate to the detail view, translating event coords to detail view widget coords */
      GtkWidget *detailView = treeGetDetailView(tree);
      propagateEventButton(header, detailView, event);
    }
  
  return handled;
}


static BlxSequence* treeGetSequence(GtkTreeModel *model, GtkTreeIter *iter)
{
  BlxSequence *result = NULL;
  GList *mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      const MSP *firstMsp = (const MSP*)(mspGList->data);
      result = firstMsp->sSequence;
    }
  
  return result;
}


/* Get all the tree rows that contain MSPs from the given sequence */
static GList *treeGetSequenceRows(GtkWidget *tree, const BlxSequence *clickedSeq)
{
  GList *resultList = NULL;
  
  /* Loop through all tree rows looking for any with a matching sequence */
  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
  
  GtkTreeIter iter;
  gboolean validIter = gtk_tree_model_get_iter_first(model, &iter);
  
  while (validIter)
    {
      GList *mspGList = treeGetMsps(model, &iter);
      
      if (g_list_length(mspGList) > 0)
	{
	  const MSP *firstMsp = (const MSP*)(mspGList->data);
	  if (firstMsp->sSequence == clickedSeq)
	    {
	      GtkTreePath *path = gtk_tree_model_get_path(model, &iter);
	      resultList = g_list_append(resultList, path);
	    }
	}
      
      validIter = gtk_tree_model_iter_next(model, &iter);
    }

  return resultList;
}


/* Move the current row selection up/down by one row */
gboolean treeMoveRowSelection(GtkWidget *tree, const gboolean moveUp, const gboolean shiftModifier)
{
  /* Get the last selected sequence */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
  
  /* Get the row in this tree that contain MSPs from that sequence (if there are multiple 
   * rows behaviour is a bit ambiguous; for now just act on the first one that was found) */
  GList *rows = treeGetSequenceRows(tree, lastSelectedSeq);

  if (g_list_length(rows) > 0)
    {
      GtkTreePath *path = (GtkTreePath*)(rows->data);
      GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
      
      GtkTreeIter newRowIter;
      gboolean newRowExists = FALSE;
      
      /* We have to use a path to get the previous row and an iter to get the next row (because
       * gtk_tree_path_next is cyclic, so we can't tell if we're at the bottom already). */
      if (moveUp && gtk_tree_path_prev(path))
	{
	  newRowExists = gtk_tree_model_get_iter(model, &newRowIter, path);
	}
      else if (!moveUp)
	{
	  newRowExists = gtk_tree_model_get_iter(model, &newRowIter, path);
	  newRowExists &= gtk_tree_model_iter_next(model, &newRowIter);
	}
      
      if (newRowExists)
	{
	  BlxSequence *newSeq = treeGetSequence(model, &newRowIter);

	  if (shiftModifier && blxWindowIsSeqSelected(blxWindow, newSeq))
	    {
	      /* The new row is already selected; deselect the last-selected row */
	      blxWindowDeselectSeq(blxWindow, lastSelectedSeq);
	    }
	  else if (shiftModifier)
	    {
	      /* Add new row to current selection */
	      blxWindowSelectSeq(blxWindow, newSeq);
	    }
	  else
	    {
	      /* New row becomes sole selection */
	      blxWindowDeselectAllSeqs(treeGetBlxWindow(tree));
	      blxWindowSelectSeq(blxWindow, newSeq);
	    }
	  
	  treeScrollSelectionIntoView(tree, NULL);
	}
    }
  
  return TRUE;
}


static gboolean onButtonReleaseTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  propagateEventButton(tree, treeGetDetailView(tree), event);
  return FALSE;
}


static gboolean onButtonReleaseTreeHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  propagateEventButton(tree, treeGetDetailView(tree), event);
  return FALSE;
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
  if (!(event->state & GDK_BUTTON1_MASK) &&
      !(event->state & GDK_BUTTON2_MASK) &&
      !(event->state & GDK_BUTTON3_MASK) &&
      !(event->state & GDK_BUTTON4_MASK) &&
      !(event->state & GDK_BUTTON5_MASK))
    {
      GtkStatusbar *statusBar = treeGetStatusBar(tree);
      
      /* Remove any previous message */
      guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
      gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
      
      /* Add details about the current row that is being hovered over */
      int bin_x = event->x, bin_y = event->y;
//      gtk_tree_view_convert_widget_to_bin_window_coords(GTK_TREE_VIEW(tree), event->x, event->y, &bin_x, &bin_y); //only from GTK v2.12
      
      GtkTreePath *path = NULL;
      GtkTreeViewColumn *column = NULL;
      
      if (gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), bin_x, bin_y, &path, &column, NULL, NULL))
	{
	  GtkTreeIter iter;
	  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
	  gtk_tree_model_get_iter(model, &iter, path);

	  GList *mspList = treeGetMsps(model, &iter);
	  
	  if (g_list_length(mspList) == 1)
	    {
	      const MSP const *msp = (const MSP const*)(mspList->data);
	      
	      char *displayText = mspGetSummaryInfo(msp);
	      
	      if (displayText)
		{
		  gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
		  g_free(displayText);
		}
	    }
	  else if (g_list_length(mspList) > 0)
	    {
	      gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, "<multiple sequences>");
	    }
	}
	
      return TRUE;
    }
  else
    {
      propagateEventMotion(tree, treeGetDetailView(tree), event);
      return FALSE;
    }
}

 
static gboolean onMouseMoveTreeHeader(GtkWidget *header, GdkEventMotion *event, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  /* The start of the header widget is at the start of the sequence column, so offset the
   * coords before propagating to the detail view */
  GtkWidget *detailView = treeGetDetailView(tree);
  IntRange xRange;
  detailViewGetColumnXCoords(detailView, BLXCOL_SEQUENCE, &xRange);
  event->x += xRange.min;

  propagateEventMotion(tree, treeGetDetailView(tree), event);
  
  return FALSE;
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
  /* Remove any message that was added by mousing over the tree rows */
  GtkStatusbar *statusBar = treeGetStatusBar(tree);
  guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);

  gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);

  /* Return true to stop the default handler re-drawing when the focus changes */
  return TRUE;
}


/* Add a row to the given tree containing the given MSP */
static void addMspToTreeRow(MSP *msp, GtkWidget *tree)
{
  if (tree)
    {
      GtkListStore *store = GTK_LIST_STORE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
      
      GtkTreeIter iter;
      gtk_list_store_append(store, &iter);
      
      /* The SequenceCellRenderer expects a GList of MSPs, so put our MSP in a list */
      GList *mspGList = NULL;
      mspGList = g_list_append(NULL, msp);
      
      gtk_list_store_set(store, &iter,
			 BLXCOL_SEQNAME, mspGetSName(msp),
			 BLXCOL_SOURCE, msp->source,
                         BLXCOL_ORGANISM, NULL,
                         BLXCOL_GENE_NAME, NULL,
                         BLXCOL_TISSUE_TYPE, NULL,
                         BLXCOL_STRAIN, NULL,
                         BLXCOL_GROUP, NULL,
			 BLXCOL_SCORE, msp->score,
			 BLXCOL_ID, msp->id,
			 BLXCOL_START, msp->sRange.min,
			 BLXCOL_SEQUENCE, mspGList,
			 BLXCOL_END, msp->sRange.max,
			 -1);
   }
}


/* Add the given msp as a row in the given tree view, and also adds it to the
 * tree's hash tables. */
void addMspToTree(GtkWidget *tree, MSP *msp)
{
  /* Add the MSP as an individual row */
  addMspToTreeRow(msp, tree);
}


/* Cell data function for the "name" column. This displays the sequence name, but
 * abbreviated to fit the cell, and with a symbol appended to indicate the s strand */
static void cellDataFunctionNameCol(GtkTreeViewColumn *column,
				     GtkCellRenderer *renderer, 
				     GtkTreeModel *model, 
				     GtkTreeIter *iter, 
				     gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  /* Get the MSP(s) for this row. They should all have the same sequence name. */
  GList	*mspGList = treeGetMsps(model, iter);

  if (g_list_length(mspGList) > 0)
    {
      MSP *msp = (MSP*)(mspGList->data);

      const int colWidth = gtk_tree_view_column_get_width(column);
      const int charWidth = treeGetCharWidth(tree);
      const int maxLen = (int)(colWidth / charWidth);

      if (maxLen > 2)
	{
	  /* Ignore any text before the colon (if there is one) */
	  const char *name = strchr(mspGetSName(msp), ':');
	  if (name)
	    {
	      name++; /* start from the char after the colon */
	    }
	  else
	    {
	      name = mspGetSName(msp); /* use the full name */
	    }

	  /* Abbreviate the name to fit the column */
	  char *displayName = abbreviateText(name, maxLen - 2);
	  char displayText[maxLen + 1];
	  sprintf(displayText, "%s", displayName);
	  
	  int i = strlen(displayName);
	  for ( ; i < maxLen - 1; ++i)
	    {
	      displayText[i] = ' ';
	    }

	  displayText[maxLen - 1] = getStrandAsChar(mspGetMatchStrand(msp));
	  displayText[maxLen] = 0;
	  
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
	  
	  g_free(displayName);
	}
      else
	{
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
	}
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

  /* Get the MSPs in this row */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {  
      /* We want to display the min coord if we're in the same direction as the q strand,
       * or the max coord if we're in the opposite direction (unless the display is reversed,
       * in which case it's vice versa). */
      const MSP const *msp = (const MSP const *)(mspGList->data);
      const gboolean sameDirection = (treeGetStrand(tree) == mspGetMatchStrand(msp));
      const gboolean findMin = (treeGetDisplayRev(tree) != sameDirection);
      
      const int coord = findMspListSExtent(mspGList, findMin);

      char displayText[numDigitsInInt(coord) + 1];
      sprintf(displayText, "%d", coord);
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
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

  /* Get the MSPs in this row */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {  
      /* We want to display the max coord if we're in the same direction as the q strand,
       * or the min coord if we're in the opposite direction (unless the display is reversed,
       * in which case it's vice versa). */
      const MSP const *msp = (const MSP const *)(mspGList->data);
      const gboolean sameDirection = (treeGetStrand(tree) == mspGetMatchStrand(msp));
      const gboolean findMin = (treeGetDisplayRev(tree) == sameDirection);
      
      const int coord = findMspListSExtent(mspGList, findMin);

      char displayText[numDigitsInInt(coord) + 1];
      sprintf(displayText, "%d", coord);
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Cell data function for the score column. */
static void cellDataFunctionScoreCol(GtkTreeViewColumn *column, 
				    GtkCellRenderer *renderer, 
				    GtkTreeModel *model, 
				    GtkTreeIter *iter, 
				    gpointer data)
{
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);

  if (g_list_length(mspGList) == 1)
    {
      const MSP const *msp = (const MSP const *)(mspGList->data);
      const gdouble score = msp->score;
      char displayText[numDigitsInInt((int)score) + 3]; /* +3 to include decimal point, 1 dp, and terminating nul */
      
//      sprintf(displayText, "%1.1f", score); 
      sprintf(displayText, "%d", (int)score); 
      
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Cell data function for the id column. */
static void cellDataFunctionIdCol(GtkTreeViewColumn *column, 
				  GtkCellRenderer *renderer, 
				  GtkTreeModel *model, 
				  GtkTreeIter *iter, 
				  gpointer data)
{
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);
  
    if (g_list_length(mspGList) == 1)
    {
      const MSP const *msp = (const MSP const *)(mspGList->data);
      const gdouble id = msp->id;
      char displayText[numDigitsInInt((int)id) + 3]; /* +3 to include decimal point, 1 dp, and terminating nul */
      
//      sprintf(displayText, "%1.1f", id); 
      sprintf(displayText, "%d", (int)id); 

      g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Cell data function for the Group column. */
static void cellDataFunctionGroupCol(GtkTreeViewColumn *column, 
                                     GtkCellRenderer *renderer, 
                                     GtkTreeModel *model, 
                                     GtkTreeIter *iter, 
                                     gpointer data)
{
  /* Get the MSP(s) for this row and find out they are in a group. All MSPs in a row should
   * be in the same sequence. */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      const MSP const *msp = (const MSP const*)(mspGList->data);

      GtkWidget *tree = GTK_WIDGET(data);
      GtkWidget *blxWindow = treeGetBlxWindow(tree);

      SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, msp->sSequence);
      
      if (group)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, group->groupName, NULL);
        }
    }
}


/* Cell data function for the Organism column. */
static void cellDataFunctionOrganismCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  if (g_list_length(mspGList) > 0)
    {
      const MSP const *msp = (const MSP const*)(mspGList->data);
      if (msp->sSequence && msp->sSequence->organism)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, msp->sSequence->organism->str, NULL);
        }
    }
}


/* Cell data function for the Gene Name column. */
static void cellDataFunctionGeneNameCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  if (g_list_length(mspGList) > 0)
    {
      const MSP const *msp = (const MSP const*)(mspGList->data);
      if (msp->sSequence && msp->sSequence->geneName)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, msp->sSequence->geneName->str, NULL);
        }
    }
}


/* Cell data function for the Tissue Type column. */
static void cellDataFunctionTissueTypeCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  if (g_list_length(mspGList) > 0)
    {
      const MSP const *msp = (const MSP const*)(mspGList->data);
      if (msp->sSequence && msp->sSequence->tissueType)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, msp->sSequence->tissueType->str, NULL);
        }
    }
}


/* Cell data function for the Strain column. */
static void cellDataFunctionStrainCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  if (g_list_length(mspGList) > 0)
    {
      const MSP const *msp = (const MSP const*)(mspGList->data);
      if (msp->sSequence && msp->sSequence->strain)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, msp->sSequence->strain->str, NULL);
        }
    }
}


/* Utility function to calculate the width of a vertical scrollbar */
static int scrollBarWidth()
{
  static int result = UNSET_INT;
  
  if (result == UNSET_INT)
    {
      /* Create a temp scrollbar and find the default width from the style properties. */
      GtkWidget *scrollbar = gtk_vscrollbar_new(NULL);
      
      gint sliderWidth = 0, separatorWidth = 0, troughBorder = 0, stepperSpacing = 0;
      gtk_widget_style_get(scrollbar, "slider-width", &sliderWidth, NULL);
      gtk_widget_style_get(scrollbar, "separator-width", &separatorWidth, NULL);
      gtk_widget_style_get(scrollbar, "trough-border", &troughBorder, NULL);
      gtk_widget_style_get(scrollbar, "stepper-spacing", &stepperSpacing, NULL);
      
      gtk_widget_destroy(scrollbar);

      result = sliderWidth + separatorWidth*2 + troughBorder*2 + stepperSpacing*2 + 4; /* to do: find out why the extra fudge factor is needed here */
    }
  
  return result;
}


/* Callback called when the width of the sequence column has changed. (To do: this
 * will be needed if we come to do drag-resizing of columns. At the moment it is
 * not needed because column resizing only happens when the window is resized, so it
 * is handed by the detail view's onSizeAllocate function.) */
static void onSeqColWidthChanged(GtkTreeViewColumn *column, GParamSpec *paramSpec, gpointer data)
{
//  GtkWidget *tree = GTK_WIDGET(data);
//  GtkWidget *detailView = treeGetDetailView(tree);
//  
//  DetailViewColumnInfo *columnInfo = detailViewGetColumnInfo(detailView, BLXCOL_SEQUENCE);
//  columnInfo->width = gtk_tree_view_column_get_width(column);
//  printf ("set seq col w = %d\n", columnInfo->width);
//  
//  updateSeqColumnSize(treeGetDetailView(tree));
}


/* Create a single column in the tree. */
static GtkTreeViewColumn* createTreeColumn(GtkWidget *tree, 
                                           GtkCellRenderer *renderer, 
                                           DetailViewColumnInfo *columnInfo)
{
  /* Create the column */
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(
    columnInfo->title, 
    renderer, 
    columnInfo->propertyName, columnInfo->columnId,  /* set the given property for this column */
    RENDERER_DATA_PROPERTY, BLXCOL_SEQUENCE,         /* also set the data property so all columns have access to the MSP data */
    NULL);

  /* Reduce the width of the end col by the scrollbar width. (This is so it matches the width
   * of the header, which does not have a scrollbar. ) */
  int width = columnInfo->width;
  if (columnInfo->columnId == BLXCOL_END)
    {
      width = width - scrollBarWidth();
    }
  
  if (columnInfo->columnId == BLXCOL_SEQUENCE)
  {
    g_signal_connect(G_OBJECT(column), "notify::width", G_CALLBACK(onSeqColWidthChanged), tree);
  }
  
  /* Set the column properties and add the column to the tree */
  if (width > 0)
    {
      gtk_tree_view_column_set_visible(column, TRUE);
      gtk_tree_view_column_set_fixed_width(column, width);
    }
  else
    {
      /* Can't have 0 width, so hide the column instead */
      gtk_tree_view_column_set_visible(column, FALSE);
    }
  
  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_resizable(column, TRUE);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
  
  /* Special treatment for specific columns */
  switch (columnInfo->columnId)
  {
    case BLXCOL_SEQUENCE:
      gtk_tree_view_column_set_expand(column, TRUE);
      break;
      
    case BLXCOL_SEQNAME:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionNameCol, tree, NULL);
      break;

    case BLXCOL_START:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionStartCol, tree, NULL);
      break;
      
    case BLXCOL_END:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionEndCol, tree, NULL);
      break;
      
    case BLXCOL_SCORE:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionScoreCol, tree, NULL);
      break;

    case BLXCOL_ID:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionIdCol, tree, NULL);
      break;

    case BLXCOL_GROUP:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionGroupCol, tree, NULL);
      break;

    case BLXCOL_ORGANISM:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionOrganismCol, tree, NULL);
      break;
      
    case BLXCOL_GENE_NAME:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionGeneNameCol, tree, NULL);
      break;
      
    case BLXCOL_TISSUE_TYPE:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionTissueTypeCol, tree, NULL);
      break;
      
    case BLXCOL_STRAIN:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionStrainCol, tree, NULL);
      break;
      
    default:
      break;
  }
  
  return column;
}


/* Refresh the name column header. This displays an abbreviated version of the
 * reference sequence name. It needs refreshing when the columns change size,
 * so we can re-abbreviate with more/less text as required. */
static void refreshNameColHeader(GtkWidget *headerWidget, gpointer data)
{
  if (GTK_IS_LABEL(headerWidget))
    {
      GtkWidget *tree = GTK_WIDGET(data);

      /* Update the font, in case its size has changed */
      gtk_widget_modify_font(headerWidget , treeGetFontDesc(tree));

      TreeColumnHeaderInfo *headerInfo = treeColumnGetHeaderInfo(tree, BLXCOL_SEQNAME);
      const int colWidth = calculateColumnWidth(headerInfo, tree);

      /* Abbreviate the name */
      const char *refSeqName = blxWindowGetRefSeqName(treeGetBlxWindow(tree));
      const int maxLen = (colWidth / treeGetCharWidth(tree));
      
      char stringToAppend[] = "(+0)";
      stringToAppend[1] = (treeGetStrand(tree) == BLXSTRAND_FORWARD ? '+' : '-');
      stringToAppend[2] = *convertIntToString(treeGetFrame(tree));
      const int numCharsToAppend = strlen(stringToAppend);

      gchar *displayText = NULL;
      
      if (maxLen > numCharsToAppend)
	{
	  /* Abbreviate the name and then append the strand/frame */
	  gchar *displayName = abbreviateText(refSeqName, maxLen - numCharsToAppend);
	  displayText = g_strconcat(displayName, stringToAppend, NULL);
	  g_free(displayName);
	}
      else
	{
	  /* No space to concatenate the frame and strand. Just include whatever of the name we can */
	  displayText = abbreviateText(refSeqName, maxLen);
	}

      if (displayText)
	{
	  gtk_label_set_text(GTK_LABEL(headerWidget), displayText);
	  g_free(displayText);
	}
    }
  else
    {
      g_warning("Unexpected widget type for Name column header; header may not refresh properly.\n");
    }
}


/* Refresh the start column header. This displays the start index of the current display range */
static void refreshStartColHeader(GtkWidget *headerWidget, gpointer data)
{
  if (GTK_IS_LABEL(headerWidget))
    {
      GtkWidget *tree = GTK_WIDGET(data);
      BlxViewContext *bc = treeGetContext(tree);

      /* Update the font, in case its size has changed */
      gtk_widget_modify_font(headerWidget, treeGetFontDesc(tree));

      int displayVal = getStartDnaCoord(treeGetDisplayRange(tree), 
					treeGetFrame(tree),
					bc->seqType, 
					bc->displayRev, 
					bc->numFrames,
					&bc->refSeqRange);
      
      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(headerWidget), displayText);
    }
  else
    {
      g_warning("Unexpected widget type for Start column header; header may not refresh properly.\n");
    }
}


/* Refresh the start column header. This displays the end index of the current display range */
static void refreshEndColHeader(GtkWidget *headerWidget, gpointer data)
{
  if (GTK_IS_LABEL(headerWidget))
    {
      GtkWidget *tree = GTK_WIDGET(data);
      BlxViewContext *bc = treeGetContext(tree);
    
      int displayVal = getEndDnaCoord(treeGetDisplayRange(tree),
				      treeGetFrame(tree),
				      bc->seqType, 
				      bc->displayRev, 
				      bc->numFrames,
				      &bc->refSeqRange);
      
      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(headerWidget), displayText);
    }
  else
    {
      g_warning("Unexpected widget type for End column header; header may not refresh properly.\n");
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
      width = width + detailViewGetColumnWidth(detailView, columnId);
    }
  
  return width;
}


/* Utility function to create a TreeColumnHeaderInfo struct and initialise its values */
static TreeColumnHeaderInfo* createTreeColumnHeaderInfo(GtkWidget *headerWidget, 
							GtkWidget *tree, 
							GList *columnIds, 
							GtkCallback refreshFunc)
{
  TreeColumnHeaderInfo *headerInfo = g_malloc(sizeof(TreeColumnHeaderInfo));
  
  headerInfo->headerWidget = headerWidget;
  headerInfo->tree = tree;
  headerInfo->columnIds = columnIds;
  headerInfo->refreshFunc = refreshFunc;
  
  return headerInfo;
}


/* Create the column header widget for the given column in the tree header. It gets placed in
 * a TreecolumnHeaderInfo struct along with other header properties. The tree header shows
 * information about the reference sequence. */
static TreeColumnHeaderInfo* createTreeColHeader(GList **columnHeaders, 
                                                 GtkTreeViewColumn *treeColumn,
                                                 DetailViewColumnInfo *columnInfo,
                                                 TreeColumnHeaderInfo* firstTreeCol,
                                                 GtkWidget *headerBar,
                                                 GtkWidget *tree,
                                                 GtkWidget *detailView,
                                                 const char const *refSeqName,
                                                 const int frame,
                                                 const BlxStrand strand)
{
  /* Create a header widget for this column, if required. Also create a list of other 
   * column IDs we wish to merge under the same header (i.e. the 'Name' header spans over
   * the name column as well as the score and Id columns). */
  GtkWidget *columnHeader = NULL;
  GList *columnIds = NULL;
  GtkCallback refreshFunc = NULL;
  
  switch (columnInfo->columnId)
    {
      case BLXCOL_SEQNAME:
	{
	  /* The header above the name column will display the reference sequence name.
	   * This header will also span the score and id columns, seeing as we don't need
	   * to show any info in those columns. */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE);
	  refreshFunc = refreshNameColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
	  g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}
	
      case BLXCOL_SEQUENCE:
	{
	  /* The sequence column header contains the reference sequence. */
	  columnHeader = gtk_layout_new(NULL, NULL);
	  
	  seqColHeaderSetRow(columnHeader, frame);
	  gtk_widget_set_name(columnHeader, DNA_TRACK_HEADER_NAME);
	  g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeRefSeqHeader), tree);
	  
	  refreshFunc = refreshTextHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));

          gtk_widget_add_events(columnHeader, GDK_BUTTON_PRESS_MASK);
          gtk_widget_add_events(columnHeader, GDK_BUTTON_RELEASE_MASK);
          gtk_widget_add_events(columnHeader, GDK_BUTTON2_MOTION_MASK);
          g_signal_connect(G_OBJECT(columnHeader), "button-press-event",    G_CALLBACK(onButtonPressTreeHeader), tree);
          g_signal_connect(G_OBJECT(columnHeader), "button-release-event",  G_CALLBACK(onButtonReleaseTreeHeader), tree);
          g_signal_connect(G_OBJECT(columnHeader), "motion-notify-event",   G_CALLBACK(onMouseMoveTreeHeader), tree);
          
	  break;
	}
	
      case BLXCOL_START:
	{
	  /* The start column header displays the start index of the current display range */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE);
	  refreshFunc = refreshStartColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
          g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}
	
      case BLXCOL_END:
	{
	  /* The end column header displays the start index of the current display range */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE);
	  refreshFunc = refreshEndColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
          g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}

      default:
        /* Any columns not specified above will not have their own header, so add them under the
         * first column's header instead. */
        {
          if (firstTreeCol)
            {
              firstTreeCol->columnIds = g_list_append(firstTreeCol->columnIds, GINT_TO_POINTER(columnInfo->columnId));
            }
          
          break;
        }
    };

  TreeColumnHeaderInfo *headerInfo = NULL;
  
  if (columnHeader)
    {
      /* Create a header info struct for each column, even if its contents are null */
      headerInfo = createTreeColumnHeaderInfo(columnHeader, tree, columnIds, refreshFunc);
      *columnHeaders = g_list_append(*columnHeaders, headerInfo);
    }
  
  return headerInfo;
}


///* Tooltip function for detail-view trees */
//static gboolean onQueryTooltip(GtkWidget *tree, gint x, gint y, gboolean keyboard_mode, GtkTooltip *tooltip, gpointer data)
//{
//  int binx, biny;
//  gtk_tree_view_convert_widget_to_bin_window_coords(GTK_TREE_VIEW(tree), x, y, &binx, &biny);
//  
//  GtkTreeViewColumn *column = NULL;
//  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), binx, biny, NULL, &column, NULL, NULL);
//  
//  if (column == gtk_tree_view_get_column(GTK_TREE_VIEW(tree), BLXCOL_SEQNAME))
//    {
//        gtk_tooltip_set_text(tooltip, "hello");
//
//    }
//  
//  return TRUE;
//}


/* Create the columns. Returns a list of header info for the column headers */
static GList* createTreeColumns(GtkWidget *tree, 
                                GtkWidget *detailView,
                                GtkCellRenderer *renderer, 
                                const BlxSeqType seqType,
                                GList *columnList,
                                GtkWidget *columnHeaderBar,
                                const char const *refSeqName,
                                const int frame,
                                const BlxStrand strand)
{
  /* We'll create a list of tree column headers */
  GList *treeColumns = NULL;

  /* The columns are defined by the columnList from the detail view. This list contains the
   * info such as the column width and title for each column. Loop through and create the tree
   * column header for each detail-view column. If any columns do not require their own header
   * in the tree, add them under the first column's header instead. */
  TreeColumnHeaderInfo* firstTreeCol = NULL;
  GList *column = columnList;
  
  for ( ; column; column = column->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)column->data;
      
      if (columnInfo)
	{
	  GtkTreeViewColumn *treeColumn = createTreeColumn(tree, renderer, columnInfo);
          
	  TreeColumnHeaderInfo* headerInfo = createTreeColHeader(&treeColumns, 
                                                                 treeColumn, 
                                                                 columnInfo, 
                                                                 firstTreeCol, 
                                                                 columnHeaderBar, 
                                                                 tree, 
                                                                 detailView, 
                                                                 refSeqName, 
                                                                 frame, 
                                                                 strand);
          
          if (!firstTreeCol)
            {
              firstTreeCol = headerInfo;
            }
	}
      else
	{
	  g_warning("Error creating column; invalid column info in detail-view column-list.\n");
	}
    }
  
  /* Set a tooltip to display the sequence name. To do: at the moment this shows
   * the tooltip when hovering anywhere over the tree: really we probably just
   * want to show it when hovering over the name column. */
//  gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(tree), BLXCOL_SEQNAME); /* only in GTK 2.12 and higher */
//  gtk_widget_set_has_tooltip(tree, TRUE);
//  g_signal_connect(G_OBJECT(tree), "query-tooltip", G_CALLBACK(onQueryTooltip), NULL);

  return treeColumns;
}


/* Add the header widget for each column to the tree's header bar */
static void addColumnsToTreeHeader(GtkWidget *headerBar, GList *columnList)
{
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      TreeColumnHeaderInfo *columnInfo = (TreeColumnHeaderInfo*)(columnItem->data);
      
      /* Put the header widget in an event box so that we can color its background. */
      GtkWidget *eventBox = gtk_event_box_new();
      gtk_container_add(GTK_CONTAINER(eventBox), columnInfo->headerWidget);
      
      /* Put the event box into the header bar. If it's the sequence column, set the expand property. */
      gboolean expand = (g_list_find(columnInfo->columnIds, GINT_TO_POINTER(BLXCOL_SEQUENCE)) != NULL);
      
      gtk_box_pack_start(GTK_BOX(headerBar), eventBox, expand, TRUE, 0);
    }
}


/* Sort comparison function for sorting by group */
static gint sortByGroupCompareFunc(const MSP *msp1, const MSP *msp2, GtkWidget *tree)
{
  gint result = 0;
  
  /* Get the order number out of the group and sort on that. If the sequence
   * is not in a group, its order number is UNSET_INT, and it gets sorted after
   * any sequences that are in groups. */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  const int msp1Order = sequenceGetGroupOrder(blxWindow, msp1->sSequence);
  const int msp2Order = sequenceGetGroupOrder(blxWindow, msp2->sSequence);
  
  if (msp1Order == UNSET_INT && msp2Order != UNSET_INT)
    {
      result = 1;
    }
  else if (msp1Order != UNSET_INT && msp2Order == UNSET_INT)
    {
      result = -1;
    }
  else
    {
      result = msp1Order - msp2Order;
    }

  return result;
}


/* Sort comparison function for sorting by string values. Allows NULL values and
 * sorts them AFTER non-null values. Comparison is case-insensitive. */
static int sortByStringCompareFunc(char *str1, char *str2)
{
  int result = 0;

  if (!str1 && !str2)
    {
      result = 0;
    }
  else if (!str1)
    {
      result = 1;
    }
  else if (!str2)
    {
      result = -1;
    }
  else
    {
      const int len = min(strlen(str1), strlen(str2));
      result = g_ascii_strncasecmp(str1, str2, len);
    }
  
  return result;
}


/* Sort comparison function.  Returns a negative value if the first row appears before the second,
 * positive if the second appears before the first, or 0 if they are equivalent according to
 * the current search criteria. */
static gint sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data)
{
  gint result = UNSET_INT;

  /* Extract the sort column and sort order */
  gint col;
  gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model), &col, NULL);

  BlxColumnId sortColumn = (BlxColumnId)col; 

  /* Extract the MSP lists from the tree rows */
  GList *mspGList1 = treeGetMsps(model, iter1);
  GList *mspGList2 = treeGetMsps(model, iter2);

  /* Get the first MSP in each list. */
  MSP *msp1 = (MSP*)(mspGList1->data);
  MSP *msp2 = (MSP*)(mspGList2->data);

  /* Check whether either row has more than one MSP. If so, it means some options aren't applicable. */
  const gboolean multipleMsps = g_list_length(mspGList1) > 1 || g_list_length(mspGList2) > 1;
  
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewColumnInfo *columnInfo = detailViewGetColumnInfo(detailView, sortColumn);
  
  switch (sortColumn)
    {
      case BLXCOL_SEQNAME:
        result = sortByStringCompareFunc(msp1->sname, msp2->sname);
        break;

      case BLXCOL_SOURCE:
        result = sortByStringCompareFunc(msp1->source, msp2->source);
        break;
        
      case BLXCOL_SCORE:
        result = multipleMsps ? 0 : (int)(msp1->score - msp2->score);
        break;
	
      case BLXCOL_ID:
        result = multipleMsps ? 0 : (int)(msp1->id - msp2->id);
        break;

      case BLXCOL_START:
        result = multipleMsps ? 0 : msp1->qRange.min - msp2->qRange.min;
        break;

      case BLXCOL_GROUP:
        result = sortByGroupCompareFunc(msp1, msp2, tree);
        break;

      case BLXCOL_ORGANISM:
        result = sortByStringCompareFunc(mspGetOrganism(msp1), mspGetOrganism(msp2));
        break;

      case BLXCOL_GENE_NAME:
        result = sortByStringCompareFunc(mspGetGeneName(msp1), mspGetGeneName(msp2));
        break;

      case BLXCOL_TISSUE_TYPE:
        result = sortByStringCompareFunc(mspGetTissueType(msp1), mspGetTissueType(msp2));
        break;

      case BLXCOL_STRAIN:
        result = sortByStringCompareFunc(mspGetStrain(msp1), mspGetStrain(msp2));
        break;
        
      default:
        g_warning("Sort function not implemented for column '%s'.\n", columnInfo->title);
        break;
    };
  
  return result;
}


/* Create the base data store for a detail view tree */
void treeCreateBaseDataModel(GtkWidget *tree, gpointer data)
{
  /* Create the data store for the tree view */
  GtkListStore *store = gtk_list_store_new(BLXCOL_NUM_COLUMNS, TREE_COLUMN_TYPE_LIST);
  
  /* Set the sort function for each column */
  int colNum = 0;
  for ( ; colNum < BLXCOL_NUM_COLUMNS; ++colNum)
    {
      gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), colNum,	sortColumnCompareFunc, tree, NULL);
    }
  
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
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_NONE);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_MULTIPLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, FALSE);
  
  /* Set the background color for the rows to be the same as the widget's background color */
  gtk_widget_modify_base(GTK_WIDGET(tree), GTK_STATE_NORMAL, GTK_WIDGET(tree)->style->bg);

  /* The default text color when rows are selected is white. This doesn't work 
   * well against our default background color of cyan, so use the same text color
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


/* Create the widget that will contain all the header widgets. */
static GtkWidget *createDetailViewTreeHeader(GtkWidget *detailView, 
					     const gboolean includeSnpTrack,
					     const BlxStrand strand,
					     GtkWidget *columnHeaderBar)
{
  GtkWidget *treeHeader = NULL; /* the outermost container for the header widgets */
  
  if (includeSnpTrack)
    {
      /* Pack the snp track and the column header bar into a vbox */
      treeHeader = gtk_vbox_new(FALSE, 0);
      gtk_widget_set_name(treeHeader, HEADER_CONTAINER_NAME);
      
      createSnpTrackHeader(GTK_BOX(treeHeader), detailView, strand);
      gtk_box_pack_start(GTK_BOX(treeHeader), columnHeaderBar, FALSE, FALSE, 0);
    }
  else
    {
      /* Only include the column headers */
      treeHeader = columnHeaderBar;
    }
  
  return treeHeader;
}


GtkWidget* createDetailViewTree(GtkWidget *grid, 
				GtkWidget *detailView, 
				GtkCellRenderer *renderer,
				GList **treeList,
				GList *columnList,
				BlxSeqType seqType,
				const char const *refSeqName,
				const int frame,
				const gboolean includeSnpTrack)
{
  /* Find the strand the tree belongs to from the grid */
  const BlxStrand strand = gridGetStrand(grid);
  
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

  GtkAdjustment *adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrollWin));
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollChangedTree), tree);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollChangedTree), tree);

  /* Create a header, and put the tree and header in a vbox. This vbox is the outermost
   * container for the tree */
  GtkWidget *columnHeaderBar = gtk_hbox_new(FALSE, 0);
  GtkWidget *treeHeader = createDetailViewTreeHeader(detailView, includeSnpTrack, strand, columnHeaderBar);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_widget_set_name(vbox, DETAIL_VIEW_TREE_CONTAINER_NAME);
  gtk_box_pack_start(GTK_BOX(vbox), treeHeader, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), scrollWin, TRUE, TRUE, 0);
  
  /* Create the tree columns */
  GList *treeColumnHeaderList = createTreeColumns(tree, detailView, renderer, seqType, columnList, columnHeaderBar, refSeqName, frame, strand);
  
  /* Add the columns to the tree header */
  addColumnsToTreeHeader(columnHeaderBar, treeColumnHeaderList);
  
  /* Set the essential tree properties */
  treeCreateProperties(tree, grid, detailView, frame, treeHeader, treeColumnHeaderList, includeSnpTrack);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);
  g_signal_connect(G_OBJECT(tree), "button-press-event",    G_CALLBACK(onButtonPressTree),	NULL);
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
  
  /* Add the tree's outermost container to the given tree list. Also increase its ref count so that
   * we can add/remove it from its parent (which we do to switch panes when we toggle strands)
   * without worrying about it being destroyed. */
  *treeList = g_list_append(*treeList, vbox);
  g_object_ref(vbox);
  
  return vbox;
}

