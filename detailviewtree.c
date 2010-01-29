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
//static GtkTreePath*	treeConvertBasePathToVisiblePath(GtkTreeView *tree, GtkTreePath *basePath);
static int		calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree);

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

int treeGetSeqType(GtkWidget *tree)
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
	  detailViewSetSelectedBaseIdx(properties->detailView, selectedBaseIdx);
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
      GtkWidget *fwdTree = detailViewGetFrameTree(detailView, FORWARD_STRAND, frame);
      GtkWidget *revTree = detailViewGetFrameTree(detailView, REVERSE_STRAND, frame);
      
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
MSP* treeGetMsp(GtkTreeModel *model, GtkTreeIter *iter)
{
  MSP *msp = NULL;
  gtk_tree_model_get(model, iter, MSP_COL, &msp, -1);
  return msp;
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
      gint charWidth = treeGetCharWidth(tree);
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


/* Select the given row if its MSP is marked as selected. */
static gboolean selectRowIfMspSelected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *mainWindow = treeGetMainWindow(tree);
  MSP *msp = treeGetMsp(model, iter);
  
  if (mainWindowIsMspSelected(mainWindow, msp))
    {
      GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
      gtk_tree_selection_select_path(selection, path);
    }
  
  return FALSE;
}


/* Select all rows in this tree whose MSPs are marked as selected */
void selectRowsForSelectedMsps(GtkWidget *tree, gpointer data)
{
  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
  gtk_tree_model_foreach(model, selectRowIfMspSelected, tree);
}


/* Mark the given row's MSP as selected in the main window's list of selected MSPs.
 * Does not allow the main window to re-update the tree selection */
static void markRowMspSelected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *mainWindow = treeGetMainWindow(tree);
  
  MSP *msp = treeGetMsp(model, iter);
  
  mainWindowSelectMsp(mainWindow, msp, FALSE);
  
  gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(tree), path, NULL, FALSE, 0.0, 0.0);
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

  selectRowsForSelectedMsps(tree, NULL);
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

	  int qSeqMin = min(msp->displayStart, msp->displayEnd);
	  int qSeqMax = max(msp->displayStart, msp->displayEnd);
	  
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
	/* Left button: select row */
	
	mainWindowDeselectAllMsps(treeGetMainWindow(tree), FALSE);
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
		  mainWindowDeselectAllMsps(treeGetMainWindow(tree), FALSE);
		}
	    }
	  else
	    {
	      /* For some reason gtk_tree_path_next doesn't indicate if next is null, so we have to get an iter... */
	      GtkTreeIter iter;
	      gtk_tree_model_get_iter(model, &iter, path);
	      
	      if (gtk_tree_model_iter_next(model, &iter))
		{
		  mainWindowDeselectAllMsps(treeGetMainWindow(tree), FALSE);
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
  
  /* Middle button: scrolls */
  if (event->button == 2)
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
      /* Note that getSequenceSegment will reverse complement the ref seq if it is the rev
       * reverse strand. This means that where there is no gaps array the comparison is trivial 
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
					       !qForward);
      
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


/* Create a dummy MSP for the reference sequence so that we can display it in the
 * detail-view tree along with all the other sequences. Note that the q coords are
 * in the DNA sequence and the s coords are in the actual displayed sequence (which
 * may be a peptide sequence if we're doing protein matches). */
//static MSP* createRefSeqMsp(GtkWidget *tree, 
//			    gboolean fwd, 
//			    char *displaySeq, 
//			    const IntRange const *refSeqRange,
//			    const IntRange *displaySeqRange,
//			    const int frame)
//{
//  MSP *msp = g_malloc(sizeof(MSP));
//  
//  /* Convert the frame number to a char. It should only be one digit long. */
//  int frameCharLen = numDigitsInInt(frame);
//  char frameChar[frameCharLen];
//  sprintf(frameChar, "%d", frame);
//  
//  if (frameCharLen > 1)
//    messerror("Frame should only be one digit but frame = '%d'. Frame may be incorrect.", frame);
//
//  msp->next = NULL;
//  msp->type = BLX_MSP_INVALID;
//  msp->score = 0;
//  msp->id = -1;
//
//  msp->qname = REFERENCE_SEQUENCE_NAME;
//  msp->qframe[0] = '(';
//  msp->qframe[1] = fwd ? '+' : '-';
//  msp->qframe[2] = frameChar[0];
//  msp->qframe[3] = ')';
//  msp->qstart = refSeqRange->min;
//  msp->qend = refSeqRange->max;
//  
//  msp->displayStart = displaySeqRange->min;
//  msp->displayEnd = displaySeqRange->max;
//  
//  msp->sname = REFERENCE_SEQUENCE_NAME;
//  msp->sframe[0] = '(';
//  msp->sframe[1] = fwd ? '+' : '-';
//  msp->sframe[2] = frameChar[0];
//  msp->sframe[3] = ')';
//  msp->slength = displaySeqRange->max - displaySeqRange->min + 1;
//  msp->sstart = displaySeqRange->min;
//  msp->send = displaySeqRange->max;
//  msp->sseq = displaySeq;
//
//  return msp;
//}
//

/* Add the given msp as a row in the given model in the given tree. If the sequence
 * column header is supplied, it indicates that that widget should also point to the
 * row for this msp. */
void addMspToTreeModel(GtkTreeModel *model, MSP *msp, GtkWidget *tree)
{
  /* First, calculate the id and set the "real" reference sequence coords (the
   * ones relative to the ref seq that is actually displayed - different to the
   * DNA ref sequence if we're displaying protein matches). */
  calcID(msp, tree);
  
  const int numReadingFrames = treeGetNumReadingFrames(tree);
  msp->displayStart = convertDnaToPeptide(msp->qstart, numReadingFrames);
  msp->displayEnd = convertDnaToPeptide(msp->qend, numReadingFrames);
  
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

  char displayText[numDigitsInInt(coord) + 1];
  sprintf(displayText, "%d", coord);
  g_object_set(renderer, START_COLUMN_PROPERTY_NAME, displayText, NULL);
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
  g_object_set(renderer, END_COLUMN_PROPERTY_NAME, displayText, NULL);
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
								       "data", MSP_COL, /* always set msp so each column has access to all the data */
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


/* Refresh the sequence column header. This header shows the section of reference
 * sequence for the current display range, so it needs to be refreshed after
 * scrolling, zooming etc. */
static void refreshSequenceColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *detailView = treeGetDetailView(tree);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  
  /* Find the segment of the ref sequence to display (complemented if this tree is
   * displaying the reverse strand, and reversed if the display is toggled */
  IntRange *displayRange = treeGetDisplayRange(tree);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(mainWindow);

  gchar *segmentToDisplay = getSequenceSegment(mainWindow, 
					       mainWindowGetRefSeq(mainWindow),
					       mainWindowGetRefSeqRange(mainWindow),
					       displayRange->min, 
					       displayRange->max, 
					       treeGetStrand(tree), 
					       mainWindowGetSeqType(mainWindow),
					       treeGetFrame(tree), 
					       mainWindowGetNumReadingFrames(mainWindow),
					       rightToLeft);
  
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  if (selectedBaseIdx == UNSET_INT)
    {
      /* Just draw plain text */
      gtk_label_set_markup(GTK_LABEL(headerWidget), segmentToDisplay);
    }
  else
    {
      /* Markup the text so the selected base is highlighted in a different colour */
      int charIdx = rightToLeft ? displayRange->max - selectedBaseIdx : selectedBaseIdx - displayRange->min;
      const int segLen = strlen(segmentToDisplay);
      
      if (charIdx >= 0 && charIdx < segLen)
	{
	  char text1[charIdx];
	  char text2[segLen - charIdx + 200];
	  char text3[2];
	  
	  int i = 0;
	  for ( ; i < charIdx; ++i)
	    {
	      text1[i] = segmentToDisplay[i];
	    }
	  text1[i] = '\0';
	  
	  text3[0] = segmentToDisplay[i];
	  text3[1] = '\0';

	  int j = 0;
	  ++i;
	  for ( ; i < segLen; ++i, ++j)
	    {
	      text2[j] = segmentToDisplay[i];
	    }
	  text2[j] = '\0';

	  GdkColor *selectedColour = treeGetRefSeqColour(tree, TRUE);

	  char markupText[segLen + 200];
	  sprintf(markupText, "%s<span bgcolor='#%x'>%s</span>%s", 
		  text1, selectedColour->pixel, text3, text2);

	  gtk_label_set_markup(GTK_LABEL(headerWidget), markupText);
	}
      else
	{
	  gtk_label_set_markup(GTK_LABEL(headerWidget), segmentToDisplay);
	}
    }
  
  
  g_free(segmentToDisplay);
}


/* Refresh the start column header. This displays the start index of the current display range */
static void refreshStartColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  if (GTK_IS_LABEL(headerWidget))
    {
      int displayVal = getStartDnaCoord(treeGetDisplayRange(tree), treeGetSeqType(tree), treeGetStrandsToggled(tree), treeGetNumReadingFrames(tree));
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
      int displayVal = getEndDnaCoord(treeGetDisplayRange(tree), treeGetSeqType(tree), treeGetStrandsToggled(tree), treeGetNumReadingFrames(tree));
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
	
      case MSP_COL:
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
      gtk_box_pack_start(GTK_BOX(headerBar), parent, (columnInfo->columnId == MSP_COL), TRUE, 0);
      
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

  /* Extract the MSPs from the tree rows */
  const MSP *msp1 = treeGetMsp(model, iter1);
  const MSP *msp2 = treeGetMsp(model, iter2);

  gint sortType = (gint)data;
  
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
    };
  
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


static void setTreeStyle(GtkTreeView *tree)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_VERTICAL);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_SINGLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, FALSE);

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
  
  messout("Created detail-view tree [%x] in container [%x]", tree, vbox);
  
  return vbox;
}

