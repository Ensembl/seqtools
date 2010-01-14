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
#include <SeqTools/utilities.h>
#include <string.h>

#define DETAIL_VIEW_TREE_NAME		"DetailViewTreeName"
#define NO_SUBJECT_SELECTED_TEXT	"<no subject selected>"


enum {SORT_BY_NAME, SORT_BY_ID, SORT_BY_SCORE, SORT_BY_POS} SortType;


/* Local function declarations */
static void		onSelectionChangedTree(GObject *selection, gpointer data);
static GtkTreePath*	treeConvertBasePathToVisiblePath(GtkTreeView *tree, GtkTreePath *basePath);
//static void		sequenceColHeaderDataFunc(GtkCellLayout *cellLayout, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data);

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

int treeGetCharWidth(GtkWidget *tree)
{
  GtkCellRenderer *renderer = treeGetRenderer(tree);
  return SEQUENCE_CELL_RENDERER(renderer)->charWidth;
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
	  
	  /* Update the feedback box and the trees */
	  updateFeedbackBox(properties->detailView);
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

PangoFontDescription* treeGetFontDesc(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetFontDesc(detailView);
}


///* Calls the given function on the given widget if it is a detail-view-tree, or, if it is a
// * container, calls the function on all children/grandchildren/etc that are dettail-view-trees */
//static void callFuncRecursivelyOnChildTrees(GtkWidget *widget, gpointer data)
//{
//  GtkCallback func = (GtkCallback)data;
//  
//  const gchar *name = gtk_widget_get_name(widget);
//  if (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0)
//    {
//      func(widget, NULL);
//    }
//  else if (GTK_IS_CONTAINER(widget))
//    {
//      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllDetailViewTrees, func);
//    }
//}


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
      GtkWidget *fwdTree = detailViewGetFrameTree(detailView, TRUE, frame);
      GtkWidget *revTree = detailViewGetFrameTree(detailView, FALSE, frame);
      
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


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
void treeUpdateFontSize(GtkWidget *tree, gpointer data)
{
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);  
  gtk_widget_modify_font(tree, fontDesc);
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
  double leftEdge = (double)renderer->xpad + renderer->xalign;
  double rightEdge = colWidth - renderer->xpad;
  
  if (x >= leftEdge && x <= rightEdge)
    {
      /* Calculate the number of bases from the left edge */
      SequenceCellRenderer *seqRenderer = SEQUENCE_CELL_RENDERER(renderer);
      gint charWidth = seqRenderer->charWidth;
      double distFromLeft = (double)x - leftEdge;
      int numBasesFromLeft = trunc(distFromLeft / (double)charWidth);
      
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
      
      if (charIdx != UNSET_INT)
	{
	  GtkAdjustment *adjustment = treeGetAdjustment(tree);
	  baseIdx = charIdx + adjustment->value + 1; /* Adjustment is 0-based, display range is 1-based */
	}
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


/***********************************************************
 *			      Updates			   *
 ***********************************************************/

/* Queue a redraw for the given tree, and the grid that corresponds to it */
void refreshTreeAndGrid(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  if (properties->grid)
    gtk_widget_queue_draw(properties->grid);
 
  /* Re-draw the tree */
  //  gtk_widget_queue_draw(tree); /* re-rendering seems to happen before this... */
}


/* Updates the feedback box with info about any currently-selected row/base in the
 * given tree.  This currently assumes single-row selections, but could be extended
 * in the future to display, say, summary information about multiple rows. Returns true
 * if the given tree has rows selected. */
gboolean updateFeedbackBoxForTree(GtkWidget *tree)
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


/***********************************************************
 *			  Selections			   *
 ***********************************************************/

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


/***********************************************************
 *                   Sorting and filtering                 *
 ***********************************************************/

/* Refilter the data for the given tree */
void refilterTree(GtkWidget *tree, gpointer data)
{
  assertTree(tree);
  GtkTreeModelFilter *filter = GTK_TREE_MODEL_FILTER(gtk_tree_view_get_model(GTK_TREE_VIEW(tree)));
  gtk_tree_model_filter_refilter(filter);
}


/* Filter function. Returns true if the given row in the given tree model should be visible. */
gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  gboolean bDisplay = FALSE;
  
  /* Find the MSP in this row */
  const MSP *msp = NULL;
  gtk_tree_model_get(model, iter, MSP_COL, &msp, -1);
  
  /* Don't show introns */
  if (msp && !mspIsIntron(msp))
    {
      /* Only show this row if part of the MSP's range is inside the displayed range */
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
  gtk_tree_sortable_set_sort_column_id(model, MSP_COL, GTK_SORT_ASCENDING);
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
	int baseIdx = getBaseIndexAtTreeCoords(tree, event->x, event->y);
	
	if (baseIdx != UNSET_INT)
	  {
	    treeSetSelectedBaseIdx(tree, baseIdx);
	  }
	
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
  GtkWidget *detailView = treeGetDetailView(tree);

  /* Update the feedback box to tell the user which sequence is selected. */
  updateFeedbackBox(detailView);
  
  /* Redraw the corresponding grid */
  TreeProperties *properties = treeGetProperties(tree);
  if (properties->grid)
    gtk_widget_queue_draw(properties->grid);
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


/* Add the given msp as a row in the given model in the given tree. If the sequence
 * column header is supplied, it indicates that that widget should also point to the
 * row for this msp. */
void addMspToTreeModel(GtkTreeModel *model, MSP *msp, GtkWidget *tree)
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


///* Cell data function for the custom sequence-column header, if there is one. This is used
// * in protein matches to display the DNA triplets for each codon. */
//static void sequenceColHeaderDataFunc(GtkCellLayout *cellLayout,
//				      GtkCellRenderer *renderer, 
//				      GtkTreeModel *model, 
//				      GtkTreeIter *iter, 
//				      gpointer data)
//{
////  g_object_set(renderer, "msp", "test text", NULL);
//}


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

  if (mspIsFake(msp))
    {
      /* For the reference sequence, show the start and end of the display range. (Temp for debug.) */
      IntRange *displayRange = treeGetDisplayRange(tree);
      coord = rightToLeft ? displayRange->max : displayRange->min;
    }
  
  char displayText[numDigitsInInt(coord) + 1];
  sprintf(displayText, "%d", coord);
  g_object_set(renderer, "start", displayText, NULL);
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
  
  if (mspIsFake(msp))
    {
      /* For the reference sequence, show the start and end of the display range. (Temp for debug.) */
      IntRange *displayRange = treeGetDisplayRange(tree);
      coord = rightToLeft ? displayRange->min : displayRange->max;
    }

  char displayText[numDigitsInInt(coord) + 1];
  sprintf(displayText, "%d", coord);
  g_object_set(renderer, "end", displayText, NULL);
}


/* Create a single column in the tree. */
static void initColumn(GtkWidget *tree, 
		       GtkCellRenderer *renderer, 
		       const char *title,
		       char *rendererProperty, 
		       const int colNum, 
		       const int width,
		       const BlxSeqType seqType,
		       const gboolean hasHeaders)
{
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(title, 
								       renderer, 
								       rendererProperty, colNum, 
								       "data", MSP_COL, /* always set msp so each column has access to all the data */
								       NULL);
  
  gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_fixed_width(column, width);
  gtk_tree_view_column_set_resizable(column, TRUE);
  
  switch (colNum)
  {
    case MSP_COL:
      gtk_tree_view_column_set_expand(column, TRUE);
      break;
      
    case START_COL:
      /* Set custom data function for start col (so that coords flip when display is reversed) */
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionStartCol, tree, NULL);
      break;
      
    case END_COL:
      /* Set custom data function for end col (so that coords flip when display is reversed) */
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionEndCol, tree, NULL);
      break;
      
    default:
      break;
  }
}


/* Create the columns. Returns the header widget for the sequence column */
static void addTreeColumns(GtkWidget *tree, 
			   GtkCellRenderer *renderer, 
			   const BlxSeqType seqType,
			   const gboolean hasHeaders)
{
  /* Add the columns */
  initColumn(tree,  renderer,  "Name",     "name",  S_NAME_COL,	NAME_COLUMN_DEFAULT_WIDTH,  seqType, hasHeaders);
  initColumn(tree,  renderer,  "Score",    "score", SCORE_COL,	SCORE_COLUMN_DEFAULT_WIDTH, seqType, hasHeaders);
  initColumn(tree,  renderer,  "%ID",      "id",    ID_COL,	ID_COLUMN_DEFAULT_WIDTH,    seqType, hasHeaders);
  initColumn(tree,  renderer,  "Start",    "start", START_COL,	START_COLUMN_DEFAULT_WIDTH, seqType, hasHeaders);
  initColumn(tree,  renderer,  "Sequence", "msp",   MSP_COL,	SEQ_COLUMN_DEFAULT_WIDTH,   seqType, hasHeaders);
  initColumn(tree,  renderer,  "End",      "end",   END_COL,	END_COLUMN_DEFAULT_WIDTH,   seqType, hasHeaders);
}


static gint sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data)
{
  gint result = UNSET_INT;

  /* Extract the sort column and sort order */
  gint sortColumn;
  GtkSortType sortOrder;
  gtk_tree_sortable_get_sort_column_id(GTK_TREE_SORTABLE(model), &sortColumn, &sortOrder);
  gboolean ascending = (sortOrder == GTK_SORT_ASCENDING);

  /* Extract the MSPs from the tree rows */
  const MSP *msp1 = treeGetMsp(model, iter1);
  const MSP *msp2 = treeGetMsp(model, iter2);

  gint sortType = (gint)data;
  
  /* Always put "fake" sequences (i.e. the reference sequence) at the top */
  if (mspIsFake(msp1) && mspIsFake(msp2))
    {
      result = 0;
    }
  else if (mspIsFake(msp1))
    {
      result = ascending ? -1 : 1;
    }
  else if (mspIsFake(msp2))
    {
      result = ascending ? 1 : -1;
    }
  else
    {
      /* Otherwise, use standard string/int comparison */
      switch (sortType)
	{
	  case SORT_BY_NAME:
	    {
	      result = strcmp(msp1->sname, msp2->sname);
	      break;
	    }
	    
	  case SORT_BY_SCORE:
	    { 
	      result = msp1->score - msp2->score;
	      break;
	    }
	    
	  case SORT_BY_ID:
	    {
	      result = msp1->id - msp2->id;
	      break;
	    }
	    
	  case SORT_BY_POS:
	    {
	      /* Use the low end of the reference sequence range */
	      int qMin1, qMin2;
	      getMspRangeExtents(msp1, &qMin1, NULL, NULL, NULL);
	      getMspRangeExtents(msp2, &qMin2, NULL, NULL, NULL);
	      result = qMin1 - qMin2;
	      break;
	    }
	    
	  default:
	    break;
	}
      
      /* If values are the same, further sort by name, then position, then length */
      if (!result && sortType != SORT_BY_NAME)
	{
	  result = strcmp(msp1->sname, msp2->sname);
	}
      
      if (!result && sortType != SORT_BY_POS)
	{
	  int sMin1, sMin2;
	  getMspRangeExtents(msp1, NULL, NULL, &sMin1, NULL);
	  getMspRangeExtents(msp2, NULL, NULL, &sMin2, NULL);
	  result = sMin1 - sMin2;
	}
      
      if (!result)
	{
	  result = msp1->slength - msp2->slength;
	}
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
  
  /* Add a 'fake' msp for the ref sequence and add it to the model */
  char *refSeq = treeGetRefSeq(tree);
  
  MSP *refSeqMsp = createRefSeqMsp(tree,
				   treeGetStrand(tree) == FORWARD_STRAND, 
				   refSeq, 
				   treeGetDisplayRange(tree));
  
  addMspToTreeModel(GTK_TREE_MODEL(store), refSeqMsp, tree);

  /* Set the sort functions for each column */
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), S_NAME_COL, sortColumnCompareFunc, (gpointer)SORT_BY_NAME, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), ID_COL, sortColumnCompareFunc, (gpointer)SORT_BY_ID, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), SCORE_COL, sortColumnCompareFunc, (gpointer)SORT_BY_SCORE, NULL);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), MSP_COL, sortColumnCompareFunc, (gpointer)SORT_BY_POS, NULL);
  
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


static void setTreeStyle(GtkTreeView *tree, const gboolean hasHeaders)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_VERTICAL);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_SINGLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, hasHeaders);
  gtk_tree_view_set_headers_clickable(tree, TRUE);

  /* Set the expander size to 0 so that we can have tiny rows (otherwise the min is 12pt) */
  char parseString[500];
  sprintf(parseString, "style \"packedTree\"\n"
	  "{\n"
	  "GtkTreeView::expander-size	      = 0\n"
	  "GtkTreeView::vertical-separator    = %d\n"
	  "GtkTreeView::horizontal-separator  = 0\n"
	  "}"
	  "widget \"*%s*\" style \"packedTree\"", VERTICAL_SEPARATOR_HEIGHT, DETAIL_VIEW_TREE_NAME);
  gtk_rc_parse_string(parseString);
}


GtkWidget* createDetailViewTree(GtkWidget *grid, 
				GtkWidget *detailView, 
				GtkCellRenderer *renderer,
				GList **treeList,
				const gboolean hasHeaders,
				BlxSeqType seqType)
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
  g_object_ref(scrollWin); /* so we can remove/re-add without worrying about it getting destroyed */
  
  /* Add the tree to the given list. We add its overall container, so we can treat the thing as a whole. */
  *treeList = g_list_append(*treeList, scrollWin);

  /* Add the columns */
  addTreeColumns(tree, renderer, seqType, hasHeaders);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);

  /* The tree needs to know which grid and renderer it corresponds to, and vice versa */
  treeCreateProperties(tree, grid, detailView, renderer);
  gridGetProperties(grid)->tree = tree;
  
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
  
  return scrollWin;
}