/*
 *  detailview.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/utilities.h>
#include <gtk/gtk.h>

#define DETAIL_VIEW_TOOLBAR_NAME	"DetailViewToolbarName"
#define SORT_BY_SCORE_STRING		"Score"
#define SORT_BY_ID_STRING		"Identity"
#define SORT_BY_NAME_STRING		"Name"
#define SORT_BY_POSITION_STRING		"Position"
#define FONT_INCREMENT_SIZE		1
#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20


/* Local function declarations */
static GtkWidget*	    detailViewGetFirstTree(GtkWidget *detailView);
static GtkWidget*	    detailViewGetBigPicture(GtkWidget *detailView);
static GtkWidget*	    detailViewGetHeader(GtkWidget *detailView);

static GtkToolItem*	    addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget);
static gboolean		    widgetIsTree(GtkWidget *widget);
static gboolean		    widgetIsTreeContainer(GtkWidget *widget);
static void		    updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc);
static void		    refreshDetailViewHeader(GtkWidget *detailView);

/***********************************************************
 *		       Utility functions                   *
 ***********************************************************/



/* Tries to return a fixed font from the list given in pref_families, returns
 * TRUE if it succeeded in finding a matching font, FALSE otherwise.
 * The list of preferred fonts is treated with most preferred first and least
 * preferred last.  The function will attempt to return the most preferred font
 * it finds.
 *
 * @param widget         Needed to get a context, ideally should be the widget you want to
 *                       either set a font in or find out about a font for.
 * @param pref_families  List of font families (as text strings).
 * @param points         Size of font in points.
 * @param weight         Weight of font (e.g. PANGO_WEIGHT_NORMAL)
 * @param font_out       If non-NULL, the font is returned.
 * @param desc_out       If non-NULL, the font description is returned.
 * @return               TRUE if font found, FALSE otherwise.
 */
static const char* findFixedWidthFontFamily(GtkWidget *widget, GList *pref_families)
{
  /* Find all the font families available */
  PangoContext *context = gtk_widget_get_pango_context(widget) ;
  PangoFontFamily **families;
  gint n_families;
  pango_context_list_families(context, &families, &n_families) ;
  
  /* Loop through the available font families looking for one in our preferred list */
  gboolean found_most_preferred = FALSE;
  gint most_preferred = g_list_length(pref_families);
  PangoFontFamily *match_family = NULL;

  gint displayTextPos;
  for (displayTextPos = 0 ; (displayTextPos < n_families && !found_most_preferred) ; displayTextPos++)
    {
      const gchar *name = pango_font_family_get_name(families[displayTextPos]) ;

      /* Look for this font family in our list of preferred families */
      GList *pref = g_list_first(pref_families) ;
      gint current = 1;
      
      while(pref)
	{
	  char *pref_font = (char *)pref->data ;
	  
	  if (g_ascii_strncasecmp(name, pref_font, strlen(pref_font)) == 0
#if GLIB_MAJOR_VERSION == 2 && GLIB_MINOR_VERSION >= 6
	      && pango_font_family_is_monospace(families[displayTextPos])
#endif
	      )
	    {
	      /* We prefer ones nearer the start of the list */
              if(current <= most_preferred)
                {
		  most_preferred = current;
		  match_family = families[displayTextPos];

                  if(most_preferred == 1)
		    {
		      found_most_preferred = TRUE;
		    }
                }

	      break;
	    }
	  
	  pref = g_list_next(pref);
	  ++current;
	}
    }

  const char *result = NULL;
  if (match_family)
    {
      result = pango_font_family_get_name(match_family);
    }
  else
    {
      messerror("Could not find a fixed-width font. Alignments may not be displayed correctly.");
    }
  
  return result;
}


/* Update the scroll position of the given adjustment to the given value. Does
 * bounds checking. */
void setDetailViewScrollPos(GtkAdjustment *adjustment, int value)
{  
  /* bounds checking */
  if (value < adjustment->lower)
    {
      value = adjustment->lower;
    }
  else
    {
      int maxValue = adjustment->upper - adjustment->page_size;
      if (value > maxValue)
	value = maxValue;
      
      if (value < 0)
	value = 0; /* Can happen if page size is larger than upper value */
    }
  
  adjustment->value = value;
  gtk_adjustment_value_changed(adjustment);
}


/* Scroll by one base */
void scrollDetailViewLeft1(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  
  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value + 1);
  else
    setDetailViewScrollPos(adjustment, adjustment->value - 1);
}

void scrollDetailViewRight1(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value - 1);
  else
    setDetailViewScrollPos(adjustment, adjustment->value + 1);
}

/* Scroll by one step increment */
void scrollDetailViewLeftStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value + adjustment->step_increment);
  else
    setDetailViewScrollPos(adjustment, adjustment->value - adjustment->step_increment);
}

void scrollDetailViewRightStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value - adjustment->step_increment);
  else
    setDetailViewScrollPos(adjustment, adjustment->value + adjustment->step_increment);
}

/* Scroll by one page size */
void scrollDetailViewLeftPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value + adjustment->page_increment);
  else
    setDetailViewScrollPos(adjustment, adjustment->value - adjustment->page_increment);
}

void scrollDetailViewRightPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (detailViewGetStrandsToggled(detailView))
    setDetailViewScrollPos(adjustment, adjustment->value - adjustment->page_increment);
  else
    setDetailViewScrollPos(adjustment, adjustment->value + adjustment->page_increment);
}


/* Calculate the number of bases that can be displayed in the sequence column */
static int calcNumBasesInSequenceColumn(GtkWidget *tree, int colWidth)
{
  /* Find the width of the sequence column (if we weren't already passed it) */
  if (colWidth == UNSET_INT)
    {
      GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), MSP_COL);
      colWidth = gtk_tree_view_column_get_width(sequenceCol);
    }
  
  /* Don't include the cell padding area */
  GtkCellRenderer *renderer = treeGetRenderer(tree);
  colWidth -= (2 * renderer->xpad) + (2 * renderer->xalign);
  
  /* Return the number of whole characters that fit in the column. */
  SequenceCellRenderer *seqRenderer = SEQUENCE_CELL_RENDERER(renderer);
  gint charWidth = seqRenderer->charWidth;
  int numChars = (int)((double)colWidth / (double)charWidth);
  
  return numChars;
}


/* This should be called when the width of the sequence column has changed. The
 * given tree can be any of the trees in the detail view - they all have the same
 * column sizes, so this function only needs to be called once. This function will
 * adjust the scroll range of our custom scroll adjustment so that it represents
 * the range that can be displayed in the new column width. */
static void updateSeqColumnSize(GtkWidget *tree, int colWidth)
{
  GtkAdjustment *adjustment = treeGetAdjustment(tree);
  
  if (adjustment)
    {
      int newPageSize = calcNumBasesInSequenceColumn(tree, colWidth);
      
      /* Only trigger updates if things have actually changed */
      if (newPageSize != adjustment->page_size)
	{
	  adjustment->page_size = newPageSize;
	  adjustment->page_increment = newPageSize;
	  
	  /* Reset the display range so that it is between the scrollbar min and max. Try to keep
	   * it centred on the same base. Note that display range is 1-based but adjustment is 0-based */
	  IntRange *displayRange = treeGetDisplayRange(tree);
	  int selectedBaseIdx = treeGetSelectedBaseIdx(tree);
	  
	  int centre = selectedBaseIdx == UNSET_INT ? getRangeCentre(displayRange) : selectedBaseIdx;
	  int offset = round((double)adjustment->page_size / 2.0);
	  adjustment->value = centre - offset - 1;
	  
	  /* If there is a gap at the end, shift the scrollbar back so that we show as much
	   * of the reference sequence as possible (unless the whole thing is already shown - 
	   * note that in this case we can end up with a gap at the start of the display if strands
	   * are toggled and the ref sequence is shorter than the page size, but this is unlikely 
	   * to happen in the real world so not worth worrying about. */
	  if (adjustment->value > adjustment->upper - adjustment->page_size)
	    {
	      adjustment->value = adjustment->upper - adjustment->page_size;
	    }
	  
	  if (adjustment->value < adjustment->lower)
	    {
	      adjustment->value = adjustment->lower;
	    }
	  
	  gtk_adjustment_changed(adjustment);
	}
    }
}


/* Find a child widget of this widget that is a paned widget. Returns null if it does not
 * have one. If there are two, it returns the first one it finds. */
static GtkWidget *findNestedPanedWidget(GtkContainer *parentWidget)
{
  GtkWidget *result = NULL;
  
  GList *children = gtk_container_get_children(parentWidget);
  GList *child = children;
  
  while (child)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);
      if (GTK_IS_PANED(childWidget))
	{
	  result = childWidget;
	  break;
	}
      
      child = child->next;
    }
  
  return result;
}


/* The given widget should be a detail-view-tree or a container (or series of nested 
 * containers) containing a single tree. This function extracts the tree and sets the headers
 * to be visible or hidden as indicated by the 'visible' flag. */
static void setTreeHeadersVisible(GtkWidget *widget, const gboolean visible)
{
  if (widgetIsTree(widget))
    {
      gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(widget), visible);
    }
  else if (widgetIsTreeContainer(widget))
    {
      /* Get the child and recurse. */
      GList *child = gtk_container_get_children(GTK_CONTAINER(widget));
      setTreeHeadersVisible(GTK_WIDGET(child->data), visible);
    }
  else 
    {
      messerror("Unexpected widget type [%d]. Expected a tree or container containing a single tree. Could not set tree header visibility.", widget);
    }
}


/* Add all trees in the given list to the detail view. Note that the treeList may not contain
 * actual trees, it may contain their containers (e.g. if they live in scrolled windows). 'first'
 * indicates that this is the first (set of) tree(s) added. This information is used to determine
 * which pane to put the tree in and/or whether the first tree in the list should have headers. */
static void addTreesToDetailView(GtkContainer *detailView, 
				 GList *treeList, 
				 const gboolean first)
{
  /* Currently we expect the detail view to be a paned widget */
  if (GTK_IS_PANED(detailView))
    {
      int numTrees = g_list_length(treeList);
      
      if (numTrees == 1)
	{
	  /* Two panes, one tree. Use 'first' flags to decide which pane to put it in and
	   * whether to show headers. */
	  GtkWidget *tree1 = GTK_WIDGET(treeList->data);
//	  setTreeHeadersVisible(tree1, first);
	  
	  if (first)
	    {
	      gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
	    }
	  else
	    {
	      gtk_paned_pack2(GTK_PANED(detailView), tree1, TRUE, TRUE);
	    }
	}
      else if (numTrees == 2)
	{
	  /* Two panes, two trees. Easy - put one in each. */
	  GtkWidget *tree1 = GTK_WIDGET(treeList->data);
	  GtkWidget *tree2 = GTK_WIDGET(treeList->next->data);
	  
	  gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
	  gtk_paned_pack2(GTK_PANED(detailView), tree2, TRUE, TRUE);
	  
//	  setTreeHeadersVisible(tree1, first);
//	  setTreeHeadersVisible(tree2, FALSE);
	}
      else if (numTrees > 2)
	{
	  /* Two panes, three trees. Put one tree in pane 1. There should be a
	   * nested widget in pane 2 that we can put the remaining trees in. */
	  GtkWidget *tree1 = GTK_WIDGET(treeList->data);
	  gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
//	  setTreeHeadersVisible(tree1, first);
	  
	  GtkWidget *nestedPanedWidget = findNestedPanedWidget(detailView);
	  
	  if (nestedPanedWidget)
	    {
	      /* Create a new list containing the remaining trees and call this 
	       * function again on the nested paned widget. */
	      GList *remainingTrees = NULL;
	      GList *nextTree = treeList->next;
	      while (nextTree)
		{
		  remainingTrees = g_list_append(remainingTrees, nextTree->data);
		  nextTree = nextTree->next;
		}
	      
	      /* Pass 'first' as false to make sure we don't add any more headers */
	      addTreesToDetailView(GTK_CONTAINER(nestedPanedWidget), remainingTrees, FALSE);
	    }
	}
    }
  else
    {
      messcrash("Unexpected detail view type: expected a paned widget. Could not add child trees.");
    }
}


/* Remove the given widget from the detail view. Does NOT remove nested paned
 * widgets, but does remove THEIR children. */
static void removeFromDetailView(GtkWidget *widget, gpointer data)
{
  GtkContainer *parent = GTK_CONTAINER(data);
  
  if (!GTK_IS_PANED(widget))
    {
      gtk_container_remove(parent, widget);
    }
  else
    {
      /* This widget is a nested paned widget. Leave it in the heirarchy, but 
       * recurse to remove all its child widgets. */
      gtk_container_foreach(GTK_CONTAINER(widget), removeFromDetailView, widget);
    }
}


/* Returns true if this widget is a detail-view-tree */
static gboolean widgetIsTree(GtkWidget *widget)
{
  gboolean result = FALSE;
  
  const gchar *name = gtk_widget_get_name(widget);
  if (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0)
    {
      result = TRUE;
    }
  
  return result;
}


/* Returns true if this widget is a detail-view-tree-container (that is, 
 * a container that has only one child which is a detail-view-tree). This is
 * intended so that we can identify the outermost parent of a single tree (so, for example,
 * we can identify the scrolled window that a tree lives in). It returns false if the widget 
 * is a container that contains multiple children and therefore might contain multiple trees. */
static gboolean widgetIsTreeContainer(GtkWidget *widget)
{
  gboolean result = FALSE;
  
  if (GTK_IS_CONTAINER(widget))
    {
      /* See if this container has an only child */
      GList *children = gtk_container_get_children(GTK_CONTAINER(widget));      
      if (children && g_list_length(children) == 1 && GTK_IS_WIDGET(children->data))
	{
	  GtkWidget *childWidget = GTK_WIDGET(children->data);
	  
	  /* See if this child is a tree, or if it is another container containing only a single tree. */
	  if (widgetIsTree(childWidget) || widgetIsTreeContainer(childWidget))
	    {
	      result = TRUE;
	    }
	}
    }
  
  return result;
}


/* This function removes all detail-view-trees from the given container widget. */
static void removeAllTreesFromContainer(GtkWidget *widget, gpointer data)
{
  /* See if this widget is a tree (or a container containing only a single tree - the
   * detail view contains the outermost single-tree container, e.g. typically the trees
   * will live in a scrolled window, and it is the scrolled windows that the detail view knows about) */
  if (widgetIsTree(widget) || widgetIsTreeContainer(widget))
    {
      GtkContainer *parent = GTK_CONTAINER(data);
      gtk_container_remove(parent, widget);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      /* It is a nested container, containing multiple children. Recurse. */
      gtk_container_foreach(GTK_CONTAINER(widget), removeAllTreesFromContainer, widget);
    }
}


/* This function removes the trees from the detail view and re-adds them in the
 * correct order according to the strandsToggled flag. It should be called every
 * time the strands are toggled. It assumes the trees are already in the 
 * detailView container, and that the properties have been set for all 3 widgets. */
static void refreshTreeOrder(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  gboolean toggled = detailViewGetStrandsToggled(detailView);

  gtk_container_foreach(GTK_CONTAINER(detailView), removeAllTreesFromContainer, detailView);

  if (properties->seqType == BLXSEQ_DNA)
    {
      /* Add both trees - the order they're added will depend on whether the display is toggled or not. */
      addTreesToDetailView(GTK_CONTAINER(detailView), properties->fwdStrandTrees, !toggled);
      addTreesToDetailView(GTK_CONTAINER(detailView), properties->revStrandTrees, toggled);
    }
  else if (properties->seqType == BLXSEQ_PEPTIDE)
    {
      /* Only add one set of trees - the reverse strand if toggled, the forward strand if not. */
      GList *treeList = toggled ? properties->revStrandTrees : properties->fwdStrandTrees;
      addTreesToDetailView(GTK_CONTAINER(detailView), treeList, TRUE);
    }
  
  gtk_widget_show_all(detailView);
}


/* Update the font for the custom sequence column header, if there is one */
static void updateDetailViewHeaderFont(GtkWidget *detailView, PangoFontDescription *fontDesc)
{
  GtkWidget *header = detailViewGetHeader(detailView);
  
  /* The sequence column header is a vbox containing 3 labels, one for each reading frame. */
  if (header != NULL && GTK_IS_CONTAINER(header))
    {
      GList *lines = gtk_container_get_children(GTK_CONTAINER(header));
      GList *line = lines;
      
      while (line)
	{
	  GtkWidget *lineWidget = GTK_WIDGET(line->data);
	  gtk_widget_modify_font(lineWidget, fontDesc);
	  line = line->next;
	}
    }
}


/* Refresh a specific line in the sequence column header */
static void refreshSequenceColHeaderLine(GtkWidget *detailView,
					 GtkWidget *lineWidget, 
					 const int frame, 
					 char *refSeq)
{
  if (!GTK_IS_LABEL(lineWidget))
    {
      messcrash("Unexpected widget type in detail view header: expected label");
    }

  
  GtkLabel *label = GTK_LABEL(lineWidget);
  
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  char displayText[displayRange->max - displayRange->min + 1 + 10];
  
  /* Loop forward/backward through the display range depending on which strand we're viewing.
   * The display range is the range for the peptides. We have to adjust by the number of reading
   * frames to get the index into the reference sequence (which is DNA bases). */
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  const int numFrames = detailViewGetNumReadingFrames(detailView);

  IntRange qRange;
  qRange.min = convertPeptideToDna(displayRange->min, frame, numFrames);
  qRange.max = convertPeptideToDna(displayRange->max, frame, numFrames);
  
  int qIdx = rightToLeft ? qRange.max : qRange.min;

  int incrementValue = numFrames;
  int displayTextPos = 0;
  gboolean done = FALSE;
  
  while (!done)
    {
      /* Get the base at this index */
      displayText[displayTextPos] = getRefSeqBase(refSeq,
						  qIdx, 
						  rightToLeft, /* we've got the reverse strand if display is toggled */
						  &qRange, 
						  BLXSEQ_DNA /* header always shows DNA bases */
						 );
      ++displayTextPos;

      /* Get the next index and check if we've finished */
      if (rightToLeft)
	qIdx -= incrementValue;
      else
	qIdx += incrementValue;
      
      if (qIdx < qRange.min || qIdx > qRange.max)
	done = TRUE;
    }

  /* to do: temp for testing to fill up seq col header while it's the wrong size (otherwise
   * it doesn't align left)... Need to fix the header size really. */
  displayText[displayTextPos++] = ' ';
  displayText[displayTextPos++] = ' ';
  displayText[displayTextPos++] = ' ';
  displayText[displayTextPos++] = ' ';
  displayText[displayTextPos++] = ' ';

  /* Make sure the string is terminated properly */
  displayText[displayTextPos] = '\0';
  
  gtk_label_set_text(label, displayText);
}


/* Refresh the detail view header */
static void refreshDetailViewHeader(GtkWidget *detailView)
{
  GtkWidget *header = detailViewGetHeader(detailView);
  
  /* The sequence column header is a vbox containing 3 labels, one for each reading frame. */
  if (header && GTK_IS_CONTAINER(header))
    {
      char* refSeq = detailViewGetRefSeq(detailView);

      GList *lines = gtk_container_get_children(GTK_CONTAINER(header));
      GList *line = lines;
      int frame = 1;
      
      while (line)
	{
	  GtkWidget *lineWidget = GTK_WIDGET(line->data);
	  refreshSequenceColHeaderLine(detailView, lineWidget, frame, refSeq);
	  
	  ++frame;
	  line = line->next;
	}
    }
}


/* Update the font description for all relevant components of the detail view */
static void updateDetailViewFontDesc(GtkWidget *detailView, PangoFontDescription *fontDesc)
{
  updateCellRendererFont(detailView, fontDesc);
  updateDetailViewHeaderFont(detailView, fontDesc);
  callFuncOnAllDetailViewTrees(detailView, treeUpdateFontSize);
}


static void incrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) + FONT_INCREMENT_SIZE;
  
  if (newSize <= MAX_FONT_SIZE)
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView, properties->fontDesc);
    }
}


static void decrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) - FONT_INCREMENT_SIZE;
  
  if (newSize >= MIN_FONT_SIZE)
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView, properties->fontDesc);
    }
}


/* Updates the feedback box with info about the currently-selected row/base. Takes
 * into account selections in any of the trees. This should be called every time 
 * the row selection or base selection is changed. 
 * TNote that there shouldn't be rows selected in different trees at the same time,
 * so this function doesn't really deal with that situation - what will happen is
 * that the contents of the feedback box will be overwritten by the last tree the
 * update function is called for. */
void updateFeedbackBox(GtkWidget *detailView)
{
  /* Loop through all of the trees. Stop if we find one that has selected rows. */
  int numFrames = detailViewGetNumReadingFrames(detailView);
  gboolean done = FALSE;
  
  int frame = 1;
  for ( ; frame <= numFrames && !done; ++frame)
    {
      /* Forward strand tree for this reading frame */
      GtkWidget *curTree = detailViewGetFrameTree(detailView, FORWARD_STRAND, frame);
      done = updateFeedbackBoxForTree(curTree);
      
      if (!done)
	{
	  /* Reverse strand tree for this reading frame */
	  GtkWidget *curTree = detailViewGetFrameTree(detailView, REVERSE_STRAND, frame);
	  done = updateFeedbackBoxForTree(curTree);
	}
    }
}

/***********************************************************
 *                    Detail view events                   *
 ***********************************************************/

static void onSizeAllocateDetailView(GtkWidget *detailView, GtkAllocation *allocation, gpointer data)
{
  GtkWidget *firstTree = detailViewGetFirstTree(detailView);
  if (firstTree)
    updateSeqColumnSize(firstTree, UNSET_INT);
}



/***********************************************************
 *              Detail view scrollbar events               *
 ***********************************************************/

static void onScrollRangeChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  int newStart = adjustment->value + 1;
  int newEnd = adjustment->value + adjustment->page_size;
  
  /* Only update if something has changed */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  if (displayRange->min != newStart || displayRange->max != newEnd)
    {
      /* Adjustment is zero-based but display range is 1-based */
      displayRange->min = newStart;
      displayRange->max = newEnd;
      
      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      callFuncOnAllDetailViewTrees(detailView, refilterTree);

      /* Scroll big picture if necessary to keep highlight box in view */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture);

      /* Recalculate the borders for all the grids and the header in the big picture */
      GtkWidget *header = bigPictureGetGridHeader(bigPicture);
      calculateGridHeaderBorders(header);
      callFuncOnAllBigPictureGrids(bigPicture, calculateGridBorders);
      
      /* Redraw all trees (and their corresponding grids) */
      callFuncOnAllDetailViewTrees(detailView, refreshTreeAndGrid);
      
      refreshDetailViewHeader(detailView);
    }
}


static void onScrollPosChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Set the display range so that it starts at the new scroll pos */
  int newStart = adjustment->value + 1;
  int newEnd = adjustment->value + adjustment->page_size;

  /* Only update if something has changed */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  if (displayRange->min != newStart || displayRange->max != newEnd)
    {
      displayRange->min = newStart;
      displayRange->max = newEnd;

      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      callFuncOnAllDetailViewTrees(detailView, refilterTree);

      /* Scroll big picture if necessary to keep highlight box in view */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture);
      
      /* Update the highlight box position for all grids in the big picture */
      callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders);
            
      /* Redraw all trees (and their corresponding grids) */
      callFuncOnAllDetailViewTrees(detailView, refreshTreeAndGrid);
      
      refreshDetailViewHeader(detailView);
    }
}


/***********************************************************
 *            Detail View toolbar events                   *
 ***********************************************************/

static void onZoomInDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Remember the size of the sequence column. This gets reset to 0 when we change the font. */
  GtkWidget *firstTree = detailViewGetFirstTree(detailView);
  int colWidth = UNSET_INT;
  if (firstTree)
    {
      GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(firstTree), MSP_COL);
      colWidth = gtk_tree_view_column_get_width(sequenceCol);
    }
  
  incrementFontSize(detailView);
  
  if (firstTree && colWidth != UNSET_INT)
    {
      updateSeqColumnSize(firstTree, colWidth);
    }
}

static void onZoomOutDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Remember the size of the sequence column. This gets reset to 0 when we change the font. */
  GtkWidget *firstTree = detailViewGetFirstTree(detailView);
  int colWidth = UNSET_INT;
  if (firstTree)
    {
      GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(firstTree), MSP_COL);
      colWidth = gtk_tree_view_column_get_width(sequenceCol);
    }
  
  decrementFontSize(detailView);
  
  if (firstTree && colWidth != UNSET_INT)
    {
      updateSeqColumnSize(firstTree, colWidth);
    }
}


/* Set the font for the detail view cell renderer. Should be called after
 * the font size is changed by zooming in/out of the detail view. */
static void updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc)
{
  /* Calculate the row height from the font metrics */
  PangoContext *context = gtk_widget_get_pango_context(detailView);
  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
  
  gint charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  gint charWidth = pango_font_metrics_get_approximate_char_width(metrics) / PANGO_SCALE;
  
  pango_font_metrics_unref(metrics);
  
  /* Cache these results in the sequence renderer for easier calculations when we render the sequence column */
  GtkCellRenderer *renderer = detailViewGetRenderer(detailView);
  SequenceCellRenderer *seqRenderer = SEQUENCE_CELL_RENDERER(renderer);
  seqRenderer->charHeight = charHeight;
  seqRenderer->charWidth = charWidth;
  
  /* Set the row height. Subtract the size of the vertical separators . This makes the
   * row smaller than it is, but we will render it at normal size, so we end up drawing over
   * the gaps, hence giving the impression that there are no gaps. */
  int verticalSeparator = detailViewGetVerticalSeparator(detailView);
  gint rowHeight = charHeight - verticalSeparator * 2;
  gtk_cell_renderer_set_fixed_size(renderer, 0, rowHeight);
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static void assertDetailView(GtkWidget *detailView)
{
  /* Check it's a valid detail-view tree type */
  if (!detailView)
    messcrash("Detail-view widget is null");
  
  if (!GTK_IS_WIDGET(detailView))
    messcrash("Detail-view is not a valid widget [%x]", detailView);
  
  if (!GTK_IS_CONTAINER(detailView))
    messcrash("Detail-view is not a valid container [%x]", detailView);
  
  if (!detailViewGetProperties(detailView))
    messcrash("Tree properties not set [widget=%x]", detailView);
}

GtkWidget* detailViewGetMainWindow(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->mainWindow : NULL;
}

GtkAdjustment* detailViewGetAdjustment(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->adjustment : NULL;
}

GtkCellRenderer* detailViewGetRenderer(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->renderer : NULL;
}

/* Get the reference sequence (always the forward strand and always the DNA seq) */
char* detailViewGetRefSeq(GtkWidget *detailView)
{
  assertDetailView(detailView);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetRefSeq(mainWindow);
}

/* Get the actual displayed reference sequence (may be the DNA sequence or peptide 
 * sequence - whichever we're displaying according to the BlxSeqType). Always the
 * forward strand.) */
char* detailViewGetDisplaySeq(GtkWidget *detailView)
{
  assertDetailView(detailView);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetDisplaySeq(mainWindow);
}

char** detailViewGetGeneticCode(GtkWidget *detailView)
{
  assertDetailView(detailView);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetGeneticCode(mainWindow);
}

static GtkWidget* detailViewGetHeader(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->header : NULL;
}

gboolean detailViewGetStrandsToggled(GtkWidget *detailView)
{
  return mainWindowGetStrandsToggled(detailViewGetMainWindow(detailView));
}

PangoFontDescription *detailViewGetFontDesc(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->fontDesc : NULL;
}

int detailViewGetVerticalSeparator(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->verticalSeparator : UNSET_INT;
}

GdkColor* detailViewGetRefSeqColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->refSeqColour;
}

GdkColor* detailViewGetRefSeqSelectedColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->refSeqSelectedColour;
}

GdkColor* detailViewGetMatchColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->matchColour;
}

GdkColor* detailViewGetMatchSelectedColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->matchSelectedColour;
}

GdkColor* detailViewGetMismatchColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->mismatchColour;
}

GdkColor* detailViewGetMismatchSelectedColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->mismatchSelectedColour;
}

GdkColor* detailViewGetExonColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->exonColour;
}

GdkColor* detailViewGetExonSelectedColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->exonSelectedColour;
}

GdkColor* detailViewGetGapColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->gapColour;
}

GdkColor* detailViewGetGapSelectedColour(GtkWidget *tree)
{
  DetailViewProperties *properties = detailViewGetProperties(tree);
  return &properties->gapSelectedColour;
}


GList *detailViewGetStrandTrees(GtkWidget *detailView, const Strand strand)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return (strand == FORWARD_STRAND) ? properties->fwdStrandTrees : properties->revStrandTrees;
}


/* Extract the tree view for the given frame on the given strand */
GtkWidget* detailViewGetFrameTree(GtkWidget *detailView, const Strand strand, const int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *result = NULL;
  
  /* Get the list of trees for the forward or reverse strand, as per the input argument */
  GList *list = (strand == FORWARD_STRAND) ? properties->fwdStrandTrees : properties->revStrandTrees;

  /* Extract the tree for this given frame number. The list should be sorted in order of frame
   * number, and we should not be requesting a frame number greater than the number of items in the list */
  if (frame <= g_list_length(list))
    {
      GList *listItem = list;
      int count = 1;
      
      for ( ; listItem && count < frame; ++count)
	{
	  listItem = listItem->next;
	}
      
      if (count == frame)
	{
	  /* The item in the list is probably a container (i.e. a scrolled window), so extract the actual tree */
	  if (GTK_IS_TREE_VIEW(listItem->data))
	    {
	      result = GTK_WIDGET(listItem->data);
	    }
	  else if (GTK_IS_CONTAINER(listItem->data))
	    {
	      GList *children = gtk_container_get_children(GTK_CONTAINER(listItem->data));
	      
	      if (children && g_list_length(children) == 1 && GTK_IS_WIDGET(children->data))
		{
		  result = GTK_WIDGET(children->data);
		}
	    }
	}
    }

  if (!result)
    {
      messerror("Tree not found for '%s' strand, frame '%d'. Number of trees = '%d'. Returning NULL.", ((strand == FORWARD_STRAND) ? "forward" : "reverse"), frame, g_list_length(list));
    }
  
  return result;
}

/* Get the first tree in the 'current' list of trees (i.e. the forward strand
 * list by default, or the reverse strand list if strands are toggled). */
static GtkWidget* detailViewGetFirstTree(GtkWidget *detailView)
{
  const Strand strand = detailViewGetStrandsToggled(detailView) ? REVERSE_STRAND : FORWARD_STRAND;
  return detailViewGetFrameTree(detailView, strand, 1);
}

GList* detailViewGetFwdStrandTrees(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties->fwdStrandTrees;
}

GList* detailViewGetRevStrandTrees(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties->revStrandTrees;
}

int detailViewGetNumReadingFrames(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->numReadingFrames : UNSET_INT;
}

IntRange* detailViewGetDisplayRange(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? &properties->displayRange : NULL;
}

IntRange* detailViewGetFullRange(GtkWidget *detailView)
{
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetFullRange(mainWindow);
}

IntRange* detailViewGetRefSeqRange(GtkWidget *detailView)
{
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetRefSeqRange(mainWindow);
}

BlxSeqType detailViewGetSeqType(GtkWidget *detailView)
{
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetSeqType(mainWindow);
}

int detailViewGetSelectedBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedBaseIdx : UNSET_INT;
}

BlxBlastMode detailViewGetBlastMode(GtkWidget *detailView)
{
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetBlastMode(mainWindow);
}

static GtkWidget *detailViewGetBigPicture(GtkWidget *detailView)
{
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetBigPicture(mainWindow);
}


DetailViewProperties* detailViewGetProperties(GtkWidget *widget)
{
  return widget ? (DetailViewProperties*)(g_object_get_data(G_OBJECT(widget), "DetailViewProperties")) : NULL;
}


static void onDestroyDetailView(GtkWidget *widget)
{
  DetailViewProperties *properties = detailViewGetProperties(widget);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "DetailViewProperties", NULL);
    }
}


static void detailViewCreateProperties(GtkWidget *detailView,
				       GtkWidget *mainWindow,
				       GtkCellRenderer *renderer,
				       GList *fwdStrandTrees,
				       GList *revStrandTrees,
				       GtkWidget *feedbackBox,
				       GtkWidget *header,
				       GtkAdjustment *adjustment, 
				       BlxSeqType seqType,
				       int numReadingFrames,
				       IntRange *displayRange,
				       PangoFontDescription *fontDesc)
{
  if (detailView)
    { 
      DetailViewProperties *properties = g_malloc(sizeof *properties);

      properties->mainWindow = mainWindow;
      properties->renderer = renderer;
      properties->fwdStrandTrees = fwdStrandTrees;
      properties->revStrandTrees = revStrandTrees;
      properties->feedbackBox = feedbackBox;
      properties->header = header;
      properties->adjustment = adjustment;
      properties->seqType = seqType;
      properties->numReadingFrames = numReadingFrames;
      properties->displayRange.min = displayRange->min;
      properties->displayRange.max = displayRange->max;
      properties->selectedBaseIdx = UNSET_INT;
      properties->fontDesc = fontDesc;
      properties->verticalSeparator = VERTICAL_SEPARATOR_HEIGHT;

      properties->refSeqColour		  = getGdkColor(GDK_YELLOW);
      properties->refSeqSelectedColour	  = getGdkColor(GDK_DARK_YELLOW);
      properties->matchColour		  = getGdkColor(GDK_CYAN);
      properties->matchSelectedColour	  = getGdkColor(GDK_DARK_CYAN);
      properties->mismatchColour	  = getGdkColor(GDK_GREY);
      properties->mismatchSelectedColour  = getGdkColor(GDK_DARK_GREY);
      properties->exonColour		  = getGdkColor(GDK_YELLOW);
      properties->exonSelectedColour	  = getGdkColor(GDK_DARK_YELLOW);
      properties->gapColour		  = getGdkColor(GDK_GREY);
      properties->gapSelectedColour	  = getGdkColor(GDK_DARK_GREY);
      
      g_object_set_data(G_OBJECT(detailView), "DetailViewProperties", properties);
      g_signal_connect(G_OBJECT(detailView), "destroy", G_CALLBACK(onDestroyDetailView), NULL); 
    }
}


/***********************************************************
 *                     Callbacks                           *
 ***********************************************************/

static void blixemSettings(void)
{
  //    if (!graphActivate(settingsGraph)) {
  //	settingsGraph = graphCreate(TEXT_FIT, "Blixem Settings", 0, 0, .6, .3);
  //	graphRegister(PICK, settingsPick);
  //	settingsRedraw();
  //    }
  //    else graphPop();
}

static void ToggleStrand(GtkWidget *detailView)
{
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(detailViewGetMainWindow(detailView));
  mainWindowProperties->strandsToggled = !mainWindowProperties->strandsToggled;
  
  /* Toggle the position of the scrollbar */
//  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
//  IntRange *displayRange = detailViewGetDisplayRange(detailView);
//  if (mainWindowProperties->strandsToggled)
//    {
//      adjustment->value = adjustment->upper - displayRange->max - 1;
//    }
//    else
//    {
//      adjustment->value = displayRange->min - 1;
//    }
//  
//  gtk_adjustment_changed(adjustment);
  
  /* Refresh the tree and grid order (i.e. switch them based on the new toggle status) */
  refreshTreeOrder(detailView);
  
  GtkWidget *bigPicture = mainWindowGetBigPicture(detailViewGetMainWindow(detailView));
  refreshGridOrder(bigPicture);
  
  /* Redraw all trees (and their corresponding grids) */
  callFuncOnAllDetailViewTrees(mainWindowProperties->detailView, refreshTreeAndGrid);
  
  /* Redraw the grid header and detail view header */
  gtk_widget_queue_draw(bigPictureGetGridHeader(mainWindowProperties->bigPicture));
  refreshDetailViewHeader(detailView);
}


static void Goto(void)
{
  //  static char dfault[32] = "";
  //  int i = 0;
  //  ACEIN pos_in;
  //  
  //  /*sprintf(posstr, "%d", dispstart + qoffset);*/
  //  
  //  if ((pos_in = messPrompt ("Goto which position: ", dfault, "t", 0)))
  //    {
  //      aceInInt(pos_in, &i);
  //      aceInDestroy (pos_in);
  //      
  //      dispstart = i - qoffset ;
  //      sprintf(dfault, "%d", i) ;
  //      
  //      blviewRedraw();
  //    }
  //  
  //  return ;
}

/* Sort the match entries by..... */

static void sortByName(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByName);
}

static void sortByScore(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByScore);
}

static void sortByPos(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByPos);
}

static void sortById(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortById);
}



static void prevMatch(void)
{
  //  gotoMatch(-1);
}

static void nextMatch(void)
{
  //  gotoMatch(1);
}


static void comboChange(GtkEditable *editBox, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Get the value to sort by from the combo box */
  gchar *val = gtk_editable_get_chars(editBox, 0, -1);
  
  if (GTK_WIDGET_REALIZED(detailView))
    {
      if (strcmp(val, SORT_BY_SCORE_STRING) == 0)
	sortByScore(detailView);
      else if (strcmp(val, SORT_BY_ID_STRING) == 0)
	sortById(detailView);
      else if (strcmp(val, SORT_BY_NAME_STRING) == 0)
	sortByName(detailView);
      else if (strcmp(val, SORT_BY_POSITION_STRING) == 0)
	sortByPos(detailView);
    }
  
  g_free(val);
}


static void GHelp(GtkButton *button, gpointer data)
{
  blxHelp();
}

static void GSettings(GtkButton *button, gpointer data)
{
  blixemSettings();
}

static void GGoto(GtkButton *button, gpointer data)
{
  Goto();
}

static void GprevMatch(GtkButton *button, gpointer data)
{
  prevMatch();
}

static void GnextMatch(GtkButton *button, gpointer data)
{
  nextMatch();
}

static void GscrollLeftBig(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  scrollDetailViewLeftPage(detailView);
}

static void GscrollRightBig(GtkButton *button, gpointer data)
{  
  GtkWidget *detailView = GTK_WIDGET(data);
  scrollDetailViewRightPage(detailView);
}

static void GscrollLeft1(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  scrollDetailViewLeft1(detailView);
}

static void GscrollRight1(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  scrollDetailViewRight1(detailView);
}

static void GToggleStrand(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  ToggleStrand(detailView);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

//static GtkWidget *createTripletCell(GtkCellRenderer *renderer)
//{
//  GtkWidget *cellView = gtk_cell_view_new();
//  
//  GtkCellLayout *layout = GTK_CELL_LAYOUT(cellView);
//  gtk_cell_layout_pack_start(layout, renderer, TRUE);
//  gtk_cell_layout_add_attribute(layout, renderer, "header", MSP_COL);
//  gtk_cell_layout_add_attribute(layout, renderer, "msp", MSP_COL);
//  gtk_cell_layout_add_attribute(layout, renderer, "data", MSP_COL);
//  //      gtk_cell_layout_set_cell_data_func(layout, renderer, sequenceColHeaderDataFunc, NULL, NULL);  
//  
//  return cellView;
//}


static GtkWidget* createTripletCell(const int frame)
{
  GtkWidget *label = gtk_label_new("");
  gtk_label_set_width_chars(GTK_LABEL(label), 1);
  return label;
}


/* Create the header bar for the detail view trees. Uses the custom header
 * widget for the sequence column if supplied, otherwise just creates a normal
 * label for this column header. */
static GtkWidget* createDetailViewHeader(GtkWidget *sequenceColHeader)
{
  GtkWidget *header = gtk_hbox_new(FALSE, 0);
  
  GtkWidget *nameLabel = gtk_label_new("Name");
  GtkWidget *scoreLabel = gtk_label_new("Score");
  GtkWidget *idLabel = gtk_label_new("%ID");
  GtkWidget *startLabel = gtk_label_new("Start");
  GtkWidget *sequenceLabel = sequenceColHeader ? sequenceColHeader : gtk_label_new("Sequence");
  GtkWidget *endLabel = gtk_label_new("End");  

  gtk_widget_set_size_request(nameLabel, NAME_COLUMN_DEFAULT_WIDTH, -1);
  gtk_widget_set_size_request(scoreLabel, SCORE_COLUMN_DEFAULT_WIDTH, -1);
  gtk_widget_set_size_request(idLabel, ID_COLUMN_DEFAULT_WIDTH, -1);
  gtk_widget_set_size_request(startLabel, START_COLUMN_DEFAULT_WIDTH, -1);
  gtk_widget_set_size_request(sequenceLabel, SEQ_COLUMN_DEFAULT_WIDTH, -1);
  gtk_widget_set_size_request(endLabel, END_COLUMN_DEFAULT_WIDTH, -1);
  
  gtk_box_pack_start(GTK_BOX(header), nameLabel, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(header), scoreLabel, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(header), idLabel, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(header), startLabel, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(header), sequenceLabel, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(header), endLabel, FALSE, TRUE, 0);
  
  return header;
}


/* Create a custom header widget for the sequence column for protein matches. This is where
 * we will display the triplets that make up the codons. Returns null for DNA matches. */
static GtkWidget* createSequenceColHeader(GtkWidget *detailView,
				 	 const BlxSeqType seqType,
					  const int numReadingFrames)
{
  GtkWidget *header = NULL;
  
  if (seqType == BLXSEQ_PEPTIDE)
    {
      header = gtk_vbox_new(FALSE, 0);
      
      /* Create a line for each reading frame */
      int frame = 0;
      for ( ; frame < numReadingFrames; ++frame)
	{
	  GtkWidget *line = createTripletCell(frame);
	  gtk_box_pack_start(GTK_BOX(header), line, FALSE, TRUE, 0);
	}
    }
  
  return header;
}


 GtkWidget* createDetailViewScrollBar(GtkAdjustment *adjustment, GtkWidget *detailView)
{
  GtkWidget *scrollBar = gtk_hscrollbar_new(adjustment);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollRangeChangedDetailView), detailView);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollPosChangedDetailView), detailView);
  
  return scrollBar;
}


static void buttonAttach(GtkHandleBox *handlebox, GtkWidget *toolbar, gpointer data)
{
  gtk_widget_set_usize(toolbar, 1, -2);
}


static void buttonDetach(GtkHandleBox *handlebox, GtkWidget *toolbar, gpointer data)
{
  gtk_widget_set_usize(toolbar, -1, -2);
}


/* Create an empty toolbar with our prefered settings. Sets the given
 * pointer to the actual toolbar and returns a pointer to the container of
 * the toolbar, if different (i.e. the widget that will be packed into the
 * parent container). */
static GtkWidget* createEmptyButtonBar(GtkWidget *parent, GtkToolbar **toolbar)
{
  /* Create a handle box for the toolbar and add it to the window */
  GtkWidget *handleBox = gtk_handle_box_new();
  
  /* Create the toolbar */
  *toolbar = GTK_TOOLBAR(gtk_toolbar_new());
  gtk_toolbar_set_tooltips(*toolbar, TRUE);
  gtk_toolbar_set_show_arrow(*toolbar, TRUE);
  gtk_toolbar_set_icon_size(*toolbar, GTK_ICON_SIZE_SMALL_TOOLBAR);

  /* Set the style property that controls the spacing */
  gtk_widget_set_name(GTK_WIDGET(*toolbar), DETAIL_VIEW_TOOLBAR_NAME);
  char parseString[500];
  sprintf(parseString, "style \"packedToolbar\"\n"
	  "{\n"
	  "GtkToolbar::space-size = 0\n"
	  "GtkToolbar::button-relief = GTK_RELIEF_NONE\n"
	  "}"
	  "widget \"*%s*\" style \"packedToolbar\"", DETAIL_VIEW_TOOLBAR_NAME);
  gtk_rc_parse_string(parseString);
  
  
  /* next three lines stop toolbar forcing the size of a blixem window */
  g_signal_connect(GTK_OBJECT(handleBox), "child-attached", G_CALLBACK(buttonAttach), NULL);
  g_signal_connect(GTK_OBJECT(handleBox), "child-detached", G_CALLBACK(buttonDetach), NULL);
  gtk_widget_set_usize(GTK_WIDGET(*toolbar), 1, -2);
  
  /* Add the toolbar to the handle box */
  gtk_container_add(GTK_CONTAINER(handleBox), GTK_WIDGET(*toolbar));
  
  return handleBox;
}


/* Create the combo box used for selecting sort criteria */
static void createSortBox(GtkToolbar *toolbar, GtkWidget *detailView)
{
  /* Add a label, to make it obvious what the combo box is for */
  GtkWidget *label = gtk_label_new(" <i>Sort HSPs by:</i>");
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  addToolbarWidget(toolbar, label);

  /* Create the combo box */
  GtkWidget *combo = gtk_combo_new();
  addToolbarWidget(toolbar, combo);

  gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(combo)->entry), FALSE);
  gtk_widget_set_usize(GTK_COMBO(combo)->entry, 80, -2);
  gtk_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", (GtkSignalFunc)comboChange, detailView);
  
  /* Create the list of strings the user can choose to sort by */
  GList *sortList = NULL;
  sortList = g_list_append(sortList, SORT_BY_SCORE_STRING);
  sortList = g_list_append(sortList, SORT_BY_ID_STRING);
  sortList = g_list_append(sortList, SORT_BY_NAME_STRING);
  sortList = g_list_append(sortList, SORT_BY_POSITION_STRING);
  gtk_combo_set_popdown_strings(GTK_COMBO(combo), sortList);
  
  /* Set the identity field as the default */
  gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), SORT_BY_ID_STRING);
}


/* Create the feedback box. (This feeds back info to the user about the currently-
 * selected base/sequence.) */
static GtkWidget* createFeedbackBox(GtkToolbar *toolbar)
{
  GtkWidget *feedbackBox = gtk_entry_new() ;

  /* User can copy text out but not edit contents */
  gtk_editable_set_editable(GTK_EDITABLE(feedbackBox), FALSE);

  /* Make it expandable so we use all available space. Set minimum size to be 0
   * because it's better to show it small than not at all. */
  gtk_widget_set_size_request(feedbackBox, 0, -1) ;
  GtkToolItem *item = addToolbarWidget(toolbar, feedbackBox) ;
  gtk_tool_item_set_expand(item, TRUE); /* make as big as possible */
  
  return feedbackBox;
}


/* Creates a single button for the detail-view toolbar */
static void makeToolbarButton(GtkToolbar *toolbar,
				 char *label,
				 char *tooltip,
				 GtkSignalFunc callback_func,
				 gpointer data)
{
  GtkToolItem *tool_button = gtk_tool_button_new(NULL, label);
  gtk_toolbar_insert(toolbar, tool_button, -1);	    /* -1 means "append" to the toolbar. */

  gtk_tool_item_set_homogeneous(tool_button, FALSE);
  gtk_tool_item_set_tooltip(tool_button, toolbar->tooltips, tooltip, NULL);
  
  gtk_signal_connect(GTK_OBJECT(tool_button), "clicked", GTK_SIGNAL_FUNC(callback_func), data);
}


/* Makes the given widget into a toolbar item on the given toolbar */
static GtkToolItem* addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget)
{
  GtkToolItem *toolItem = gtk_tool_item_new();
  gtk_container_add(GTK_CONTAINER(toolItem), widget);
  gtk_toolbar_insert(toolbar, toolItem, -1);	    /* -1 means "append" to the toolbar. */
  
  return toolItem;
}


/* Create the detail view toolbar */
static GtkWidget* createDetailViewButtonBar(GtkWidget *detailView, 
					    BlxBlastMode mode,
					    GtkWidget **feedbackBox)
{
  GtkToolbar *toolbar = NULL;
  GtkWidget *toolbarContainer = createEmptyButtonBar(detailView, &toolbar);
  
  /* Zoom buttons */
  makeToolbarButton(toolbar, "+", "Zoom in", (GtkSignalFunc)onZoomInDetailView, detailView);
  makeToolbarButton(toolbar, "-", "Zoom out", (GtkSignalFunc)onZoomOutDetailView, detailView);
  
  /* Help button */
  makeToolbarButton(toolbar, "Help",	 "Don't Panic",			 (GtkSignalFunc)GHelp,		  NULL);

  /* Combo box for sorting */
  createSortBox(toolbar, detailView);
  
  /* Settings button */
  makeToolbarButton(toolbar, "Settings", "Open the Preferences Window",  (GtkSignalFunc)GSettings,	  NULL);
  
  /* Navigation buttons */
  makeToolbarButton(toolbar, "Goto",	 "Go to specified co-ordinates", (GtkSignalFunc)GGoto,		  NULL);
  makeToolbarButton(toolbar, "< match",  "Next (leftward) match",	 (GtkSignalFunc)GprevMatch,	  NULL);
  makeToolbarButton(toolbar, "match >",  "Next (rightward) match",	 (GtkSignalFunc)GnextMatch,	  NULL);
  makeToolbarButton(toolbar, "<<",	 "Scroll leftward lots",	 (GtkSignalFunc)GscrollLeftBig,	  detailView);
  makeToolbarButton(toolbar, ">>",       "Scroll rightward lots",	 (GtkSignalFunc)GscrollRightBig,  detailView);
  makeToolbarButton(toolbar, "<",	 "Scroll leftward one base",	 (GtkSignalFunc)GscrollLeft1,	  detailView);
  makeToolbarButton(toolbar, ">",	 "Scroll rightward one base",	 (GtkSignalFunc)GscrollRight1,	  detailView);
  
  if (mode == BLXMODE_BLASTX || mode == BLXMODE_TBLASTX || mode == BLXMODE_BLASTN)
    {
      makeToolbarButton(toolbar, "Strand^v", "Toggle strand", (GtkSignalFunc)GToggleStrand, detailView);
    }
  
  *feedbackBox = createFeedbackBox(toolbar);
  
  return toolbarContainer;
}


/* Create two detail-view trees and place them in the 2 panes of the container (if one is
 * given). The first tree gets associated with grid1 and appended to list1, and the second
 * with grid2 and list2. If hasHeaders is true, the first tree will have visible headers. */
static void createTwoPanedTrees(GtkWidget *detailView,
				GtkWidget *container, 
				GtkCellRenderer *renderer,
				GtkWidget *grid1, 
				GtkWidget *grid2,
				GList **list1,
				GList **list2,
				BlxSeqType seqType,
				const gboolean hasHeaders,
				const int firstFrame)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, detailView, renderer, list1, hasHeaders, seqType, firstFrame);
  GtkWidget *tree2 = createDetailViewTree(grid2, detailView, renderer, list2, FALSE, seqType, firstFrame + 1);
  
  if (container)
    {
      gtk_paned_pack1(GTK_PANED(container), tree1, FALSE, TRUE);
      gtk_paned_pack2(GTK_PANED(container), tree2, TRUE, TRUE);
    }
}


/* Create three detail-view trees and place them in the paned widget. We only have
 * two panes in the widget, so two of the trees will be placed in a second, nested
 * paned widget. */
static void createThreePanedTrees(GtkWidget *detailView,
				  GtkCellRenderer *renderer,
				  GtkWidget *grid,
				  GList **list,
				  BlxSeqType seqType,
				  const gboolean addToDetailView,
				  const gboolean hasHeaders)
{
  /* Create a tree for pane1 (but only add it to the detailView if instructed to).
   * The first tree has headers. */
  GtkWidget *tree1 = createDetailViewTree(grid, detailView, renderer, list, hasHeaders, seqType, 1);
  
  GtkWidget *nestedPanedWidget = NULL;
  
  if (addToDetailView)
    {
      gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
      
      /* Create another paned widget and place it in pane 2 */
      nestedPanedWidget = gtk_vpaned_new();
      gtk_paned_pack2(GTK_PANED(detailView), nestedPanedWidget, TRUE, TRUE);
    }
  
  /* Create two more trees (and place them in the nested paned widget, if it is not null). 
   * Neither of these trees should have headers. */
  createTwoPanedTrees(detailView, nestedPanedWidget, renderer, grid, grid, list, list, seqType, FALSE, 2);
}


/* Create the trees for the detail view, creating sub-panes if necessary depending
 * on how many trees we need */
static void createDetailViewPanes(GtkWidget *detailView, 
				  GtkCellRenderer *renderer,
				  GtkWidget *fwdStrandGrid, 
				  GtkWidget *revStrandGrid, 
				  const int numReadingFrames,
				  GList **fwdStrandTrees,
				  GList **revStrandTrees,
				  BlxSeqType seqType,
				  const gboolean allowHeaders)
{
  if (numReadingFrames == 1)
    {
      /* DNA matches: we need 2 trees, one for the forward strand and one for the reverse */
      createTwoPanedTrees(detailView, 
			  detailView, 
			  renderer,
			  fwdStrandGrid, 
			  revStrandGrid, 
			  fwdStrandTrees, 
			  revStrandTrees, 
			  seqType,
			  allowHeaders,
			  1);
    }
  else if (numReadingFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames for EACH strand (although only
       * one set of trees will be displayed at a time). */
      createThreePanedTrees(detailView, renderer, fwdStrandGrid, fwdStrandTrees, seqType, TRUE, allowHeaders);
      createThreePanedTrees(detailView, renderer, revStrandGrid, revStrandTrees, seqType, FALSE, allowHeaders);
    }
}


/* Add the MSPs to the detail-view trees */
void detailViewAddMspData(GtkWidget *detailView, MSP *mspList)
{
  /* First, create a data store for each tree so we have something to add our data to. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateBaseDataModel);

  /* Loop through each MSP, and add it to the correct tree based on its strand and reading frame */
  MSP *msp = mspList;
  for ( ; msp; msp = msp->next)
    {
      Strand strand = (msp->qframe[1] == '+') ? FORWARD_STRAND : REVERSE_STRAND;
      int frame = atoi(&msp->qframe[2]);
      
      GtkWidget *tree = detailViewGetFrameTree(detailView, strand, frame);
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      
      addMspToTreeModel(model, msp, tree);
    }
  
  /* Finally, create a custom-filtered version of the data store for each tree. We do 
   * this AFTER adding the data so that it doesn't try to re-filter every time we add a row. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateFilteredDataModel);
}


/* We need a fixed-width font for displaying alignments. Find one from a 
 * list of possibilities. */
static const char* findDetailViewFont(GtkWidget *detailView)
{
  GList *fixed_font_list = NULL ;
  
  fixed_font_list = g_list_append(fixed_font_list, "Monaco");
  fixed_font_list = g_list_append(fixed_font_list, "Courier");
  fixed_font_list = g_list_append(fixed_font_list, "Courier New");
  fixed_font_list = g_list_append(fixed_font_list, "Monospace");
  fixed_font_list = g_list_append(fixed_font_list, "fixed");
  
  const char *fontFamily = findFixedWidthFontFamily(detailView, fixed_font_list);
  g_list_free(fixed_font_list);
  
  return fontFamily;
}


GtkWidget* createDetailView(GtkWidget *mainWindow,
			    GtkWidget *panedWidget,
			    GtkAdjustment *adjustment, 
			    GtkWidget *fwdStrandGrid, 
			    GtkWidget *revStrandGrid,
			    MSP *mspList,
			    BlxBlastMode mode,
			    BlxSeqType seqType,
			    int numReadingFrames,
			    IntRange *initDisplayRange)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so that we can identify them.) */
  GtkWidget *detailView = gtk_vpaned_new();
  
  /* Find a fixed-width font */
  const char *fontFamily = findDetailViewFont(detailView);
  PangoFontDescription *fontDesc = pango_font_description_copy(detailView->style->font_desc);
  pango_font_description_set_family(fontDesc, fontFamily);
  
  /* Create a custom cell renderer to render the sequences in the detail view */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();
  SequenceCellRenderer *seqRenderer = SEQUENCE_CELL_RENDERER(renderer);
  seqRenderer->detailView = detailView;

  /* Create the trees. */
  GList *fwdStrandTrees = NULL, *revStrandTrees = NULL;
  createDetailViewPanes(detailView, 
			renderer,
			fwdStrandGrid, 
			revStrandGrid, 
			numReadingFrames, 
			&fwdStrandTrees,
			&revStrandTrees,
			seqType,
			FALSE);
  
  /* Create the toolbar. We need to remember the feedback box. */
  GtkWidget *feedbackBox = NULL;
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, mode, &feedbackBox);

  GtkWidget *sequenceColHeader = createSequenceColHeader(detailView, seqType, numReadingFrames);
  GtkWidget *header = createDetailViewHeader(sequenceColHeader);
  
  /* Put everything in a vbox. */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), header, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), detailView, TRUE, TRUE, 0);
  
  /* Put the whole lot into the main window */
  gtk_paned_pack2(GTK_PANED(panedWidget), vbox, TRUE, TRUE);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(detailView), "size-allocate", G_CALLBACK(onSizeAllocateDetailView), NULL);
  
  /* Set the required properties */
  detailViewCreateProperties(detailView, 
			     mainWindow, 
			     renderer,
			     fwdStrandTrees,
			     revStrandTrees,
			     feedbackBox, 
			     sequenceColHeader,
			     adjustment, 
			     seqType,
			     numReadingFrames, 
			     initDisplayRange,
			     fontDesc);
  
  /* Set the initial font. (Requires detail view properties to be set) */
  updateDetailViewFontDesc(detailView, fontDesc);

  return detailView;
}

