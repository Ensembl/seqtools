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
#define SORT_BY_GROUP_ORDER_STRING	"Group"
#define FONT_INCREMENT_SIZE		1
#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20
#define NO_SUBJECT_SELECTED_TEXT	"<no subject selected>"
#define MULTIPLE_SUBJECTS_SELECTED_TEXT	"<multiple subjects selected>"


typedef struct 
  {
    const IntRange const *displayRange; 
    const gboolean searchRight;
    const int searchDirection;
    const gboolean rightToLeft;
    
    int smallestOffset;
  } MatchSearchData;


/* Local function declarations */
static GtkWidget*	    detailViewGetFirstTree(GtkWidget *detailView);
static GtkWidget*	    detailViewGetBigPicture(GtkWidget *detailView);
static GList*		    detailViewGetColumnList(GtkWidget *detailView);
static GtkWidget*	    detailViewGetFeedbackBox(GtkWidget *detailView);
static int		    detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView);
static int		    detailViewGetSelectedFrame(GtkWidget *detailView);
static GdkColor*	    detailViewGetTripletHighlightColour(GtkWidget *detailView, const gboolean isSelectedDnaBase);

static void		    detailViewCacheFontSize(GtkWidget *detailView, int charWidth, int charHeight);
static GtkToolItem*	    addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget);
static gboolean		    widgetIsTree(GtkWidget *widget);
static gboolean		    widgetIsTreeContainer(GtkWidget *widget);
static void		    updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc);
static void		    refreshDetailViewHeaders(GtkWidget *detailView);
static GtkWidget*	    createSequenceColHeader(GtkWidget *detailView, const BlxSeqType seqType, const int numReadingFrames);
static const char*	    findDetailViewFont(GtkWidget *detailView);
static void		    setDetailViewScrollPos(GtkAdjustment *adjustment, int value);

/***********************************************************
 *		       Utility functions                   *
 ***********************************************************/

/* Return the width of the column with the given column id */
int getDetailViewColumnWidth(GtkWidget *detailView, const ColumnId columnId)
{
  GList *listItem = detailViewGetColumnList(detailView);
  int result = UNSET_INT;
  
  for ( ; listItem; listItem = listItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)listItem->data;
      
      if (columnInfo && columnInfo->columnId == columnId)
	{
	  result = columnInfo->width;
	  break;
	}
    }
  
  return result;
}


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


/* Scroll the detail view so that the given coord is at the start of the display
 * range (within bounds). coordSeqType specifies the type of sequence that the coord
 * is an index into. */
void setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  /* If the given coord is on the DNA sequence but we're viewing peptide sequences, convert it */
  if (coordSeqType == BLXSEQ_DNA && detailViewGetSeqType(detailView) == BLXSEQ_PEPTIDE)
    {
      coord = convertDnaToPeptide(coord, 1, detailViewGetNumReadingFrames(detailView), NULL);
    }

  /* Adjust to a 0-based value */
  const IntRange const *fullRange = detailViewGetFullRange(detailView);
  coord = coord - fullRange->min;
  
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, coord);
}


/* Scroll the detail view so that the given coord is at the end of the display
 * range (within bounds). coordSeqType specifies the type of sequence that the coord
 * is an index into. */
static void setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  /* If the given coord is on the DNA sequence but we're viewing peptide sequences, convert it */
  if (coordSeqType == BLXSEQ_DNA && detailViewGetSeqType(detailView) == BLXSEQ_PEPTIDE)
    {
      coord = convertDnaToPeptide(coord, 1, detailViewGetNumReadingFrames(detailView), NULL);
    }

  /* Get the new start coord */
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  const int displayLen = getRangeLength(displayRange);
  int newStart = coord - displayLen + 1;

  /* Convert to 0-based value for the scrollbar */
  const IntRange const *fullRange = detailViewGetFullRange(detailView);
  newStart = newStart - fullRange->min;
  
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, newStart);
}


/* Update the scroll position of the adjustment to the given (0-based) value. Does
 * bounds checking. */
static void setDetailViewScrollPos(GtkAdjustment *adjustment, int value)
{  
  /* bounds checking */
  int maxValue = adjustment->upper - adjustment->page_size;
  
  if (value > maxValue)
    {
      value = maxValue;
    }
  
  if (value < adjustment->lower)
    {
      value = adjustment->lower;
    }
  
  adjustment->value = value;
  
  /* Emit notification that the scroll pos has changed */
   gtk_adjustment_value_changed(adjustment);
}


/* Scroll left/right by one base. */
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
      GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), SEQUENCE_COL);
      colWidth = gtk_tree_view_column_get_width(sequenceCol);
    }
  
  /* Don't include the cell padding area */
  GtkCellRenderer *renderer = treeGetRenderer(tree);
  colWidth -= (2 * renderer->xpad) + (2 * renderer->xalign);
  
  /* Return the number of whole characters that fit in the column. */
  gint charWidth = treeGetCharWidth(tree);
  int numChars = (int)((double)colWidth / (double)charWidth);
  
  return numChars;
}


/* This should be called when the width of the sequence column has changed (or the
 * size of the font has changed). The given tree can be any of the trees in the 
 * detail view - they all have the same column sizes, so this function only needs 
 * to be called once. This function will adjust the scroll range of our custom scroll
 * adjustment so that it represents the range that can be displayed in the new column 
 * width. */
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
	   * it centred on the same base. The base index is in terms of the display range coords, 
	   * so the sequence type of the coord is whatever the display sequence type is. */
	  GtkWidget *detailView = treeGetDetailView(tree);
	  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
	  const BlxSeqType seqType = detailViewGetSeqType(detailView);

	  int centre = getRangeCentre(displayRange);
	  int offset = roundNearest((double)adjustment->page_size / 2.0);
	  setDetailViewStartIdx(detailView, centre - offset, seqType);
	  
	  gtk_adjustment_changed(adjustment); /* signal that the scroll range has changed */
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
	  
	  if (first)
	    {
	      gtk_paned_pack1(GTK_PANED(detailView), tree1, TRUE, TRUE);
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
	  
	  gtk_paned_pack1(GTK_PANED(detailView), tree1, TRUE, TRUE);
	  gtk_paned_pack2(GTK_PANED(detailView), tree2, TRUE, TRUE);
	}
      else if (numTrees > 2)
	{
	  /* Two panes, three trees. Put one tree in pane 1. There should be a
	   * nested widget in pane 2 that we can put the remaining trees in. */
	  GtkWidget *tree1 = GTK_WIDGET(treeList->data);
	  gtk_paned_pack1(GTK_PANED(detailView), tree1, TRUE, TRUE);

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
  const char *name = gtk_widget_get_name(widget);
  return (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0);
}


/* Returns true if this widget is a detail-view-tree-container (that is, 
 * a container that has only one child that is a detail-view-tree.) */
static gboolean widgetIsTreeContainer(GtkWidget *widget)
{
  const char *name = gtk_widget_get_name(widget);
  return (strcmp(name, DETAIL_VIEW_TREE_CONTAINER_NAME) == 0);
}


/* Given a widget contains a single child that is a detail-view-tree, return
 * that tree. */
static GtkWidget *treeContainerGetTree(GtkContainer *container)
{
  GtkWidget *result = NULL;
  
  GList *child = gtk_container_get_children(container);
  
  for ( ; child && !result; child = child->next)
    {
      if (GTK_IS_WIDGET(child->data))
	{
	  GtkWidget *childWidget = GTK_WIDGET(child->data);
	  
	  if (widgetIsTree(childWidget))
	    {
	      result = childWidget;
	    }
	  else if (GTK_IS_CONTAINER(childWidget))
	   {
	     result = treeContainerGetTree(GTK_CONTAINER(childWidget));
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
   * will live in a vbox or scrolled window, and it is this container that the detail view
   * knows about) */
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
  
  /* Must show all child widgets because some of them may not have been in this parent before.
   * (Just calling gtk_widget_show on the individual trees doesn't seem to work.)
   * However, we then need to re-hide any that may have been previously hidden by the user. */
  gtk_widget_show_all(detailView);
  gtk_container_foreach(GTK_CONTAINER(detailView), hideUserHiddenWidget, NULL);
}


/* Set the markup in the given label, so that the label displays the given sequence
 * segment, with the given base index is highlighted in the given colour. If the given
 * base index is UNSET_INT, or is out of the current display range, then the label 
 * is just given the plain, non-marked-up sequence. */
static void createMarkupForLabel(GtkLabel *label, 
				 const char *segmentToDisplay, 
				 const int selectedCharIdx,
				 GdkColor *colour)
{
  const int segLen = strlen(segmentToDisplay);

  if (selectedCharIdx >= 0 && selectedCharIdx < segLen)
    {
      /* Cut the sequence into 4 sections: the text before the selected bases,
       * the text after them, the selected bases themselves. */
      char textBefore[selectedCharIdx];
      char textAfter[segLen - selectedCharIdx + 200];
      char selectedBaseText[2];
      
      int i = 0;
      for ( ; i < selectedCharIdx; ++i)
	{
	  textBefore[i] = segmentToDisplay[i];
	}
      textBefore[i] = '\0';
      
      selectedBaseText[0] = segmentToDisplay[i];
      selectedBaseText[1] = '\0';
      ++i;

      int j = 0;
      for ( ; i < segLen; ++i, ++j)
	{
	  textAfter[j] = segmentToDisplay[i];
	}
      textAfter[j] = '\0';
      
      /* Create the marked-up text */
      char markupText[segLen + 200];
      sprintf(markupText, "%s<span background='#%06x'>%s</span>%s", 
	      textBefore, colour->pixel, selectedBaseText, textAfter);
      
      gtk_label_set_markup(label, markupText);
    }
  else
    {
      /* Selected index is out of range, so no markup required */
      gtk_label_set_markup(label, segmentToDisplay);
    }
}


/* Refresh a specific line in the sequence column header */
static void refreshSequenceColHeaderLine(GtkWidget *detailView,
					 GtkWidget *lineWidget, 
					 int frame, /* local copy */
					 char *refSeq)
{
  if (!GTK_IS_LABEL(lineWidget))
    {
      messcrash("Unexpected widget type in detail view header: expected label");
    }

  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  const int selectedDnaBaseIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  const int selectedFrame = detailViewGetSelectedFrame(detailView);
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);

  /* We need to switch the order in which frames are displayed if the display is toggled.
   * Rather than bother moving widgets around, we'll just locally change the frame number 
   * here to pretend we're a different frame. */
  if (rightToLeft)
    {
      frame = numFrames - frame + 1;
    }
  
  /* Loop forward/backward through the display range depending on which strand we're viewing.
   * The display range is the range for the peptides. We have to adjust by the number of reading
   * frames to get the index into the reference sequence (which is DNA bases). */

  IntRange qRange;
  qRange.min = convertPeptideToDna(displayRange->min, frame, 1, numFrames);	    /* 1st base in frame */
  qRange.max = convertPeptideToDna(displayRange->max, frame, numFrames, numFrames); /* last base in frame */
  
  char displayText[displayRange->max - displayRange->min + 1 + 10];
  int incrementValue = numFrames;
  int displayTextPos = 0;
  int selectedCharIdx = UNSET_INT;
  GdkColor *highlightColour = NULL;

  int qIdx = rightToLeft ? qRange.max : qRange.min;

  while (qIdx >= qRange.min && qIdx <= qRange.max)
    {
      /* Get the base at this index */
      displayText[displayTextPos] = getRefSeqBase(refSeq,
						  qIdx, 
						  rightToLeft, /* we've got the reverse strand if display is toggled */
						  &qRange, 
						  BLXSEQ_DNA /* header always shows DNA bases */
						 );

      /* If the base (or its peptide) is selected, we need to highlight it */
      if (convertDnaToPeptide(qIdx, selectedFrame, numFrames, NULL) == selectedBaseIdx)
	{
	  /* Highlight colour depends on whether this actual DNA base is selected or just the peptide that it's in */
	  highlightColour = detailViewGetTripletHighlightColour(detailView, qIdx == selectedDnaBaseIdx);
	  selectedCharIdx = displayTextPos;
	}
      
      /* Increment indices */
      ++displayTextPos;

      if (rightToLeft)
	qIdx -= incrementValue;
      else
	qIdx += incrementValue;
    }

  /* Make sure the string is terminated properly */
  displayText[displayTextPos] = '\0';
  
  /* Mark up the text to highlight the selected base, if there is one */
  createMarkupForLabel(GTK_LABEL(lineWidget), displayText, selectedCharIdx, highlightColour);
}


/* Refresh the header label for the start-coord column. It shows 'Start' for normal
 * orientation or 'End' if the display is reversed. */
static void refreshDetailViewStartColHeader(GtkWidget *header, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  if (detailViewGetStrandsToggled(detailView))
    {
      gtk_label_set_text(GTK_LABEL(header), END_COLUMN_HEADER_TEXT);
    }
  else
    {
      gtk_label_set_text(GTK_LABEL(header), START_COLUMN_HEADER_TEXT);
    }
}


/* Refresh the header label for the end-coord column. It shows 'End' for normal
 * orientation or 'Start' if the display is reversed. */
static void refreshDetailViewEndColHeader(GtkWidget *header, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  if (detailViewGetStrandsToggled(detailView))
    {
      gtk_label_set_text(GTK_LABEL(header), START_COLUMN_HEADER_TEXT);
    }
  else
    {
      gtk_label_set_text(GTK_LABEL(header), END_COLUMN_HEADER_TEXT);
    }
}


/* Refresh the sequence column header in the detail view. This column header
 * contains the DNA sequence for the currently displayed section of peptide 
 * sequence for protein matches. It currently displays nothing for DNA matches, 
 * so this callback is not required in that case. */
static void refreshDetailViewSequenceColHeader(GtkWidget *header, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  PangoFontDescription *fontDesc = detailViewGetFontDesc(detailView);
  
  /* If the header is a normal label, there is nothing to do. If it is a container
   * containing multiple labels, then refresh each of the labels  */
  if (GTK_IS_CONTAINER(header))
    {
      char* refSeq = detailViewGetRefSeq(detailView);

      GList *lines = gtk_container_get_children(GTK_CONTAINER(header));
      GList *line = lines;
      int frame = 1;
      
      while (line)
	{
	  GtkWidget *lineWidget = GTK_WIDGET(line->data);
	  
	  /* update the font, in case it has changed */
	  gtk_widget_modify_font(lineWidget, fontDesc);

	  /* refresh the displayed text */
	  refreshSequenceColHeaderLine(detailView, lineWidget, frame, refSeq);
	  
	  ++frame;
	  line = line->next;
	}
    }
}


/* Refresh the headers for the detail view. */
static void refreshDetailViewHeaders(GtkWidget *detailView)
{
  /* Loop through all columns and call the refresh callback function on each column
   * header widget, if one exists. */
  
  GList	*columnList = detailViewGetColumnList(detailView);
  GList *column = columnList;
  
  for ( ; column; column = column->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)column->data;
      if (columnInfo)
	{
	  if (columnInfo->headerWidget && columnInfo->refreshFunc)
	    {
	      columnInfo->refreshFunc(columnInfo->headerWidget, detailView);
	    }
	}
      else
	{
	  messerror("refreshDetailViewHeaders: Invalid column data for detail view header. Header may not be refreshed correctly.");
	}
    }
}


/* Update the font description for all relevant components of the detail view */
void updateDetailViewFontDesc(GtkWidget *detailView)
{
  PangoFontDescription *fontDesc = detailViewGetFontDesc(detailView);
  
  updateCellRendererFont(detailView, fontDesc);
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
  callFuncOnAllDetailViewTrees(detailView, treeUpdateFontSize);
  gtk_widget_queue_draw(detailView);
}


static void incrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) + FONT_INCREMENT_SIZE;
  
  if (newSize <= MAX_FONT_SIZE)
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView);
    }
}


static void decrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) - FONT_INCREMENT_SIZE;
  
  /* Note that the font needs to be big enough to cover the vertical separators and padding around the
   * cell, otherwise we will end up with visible gaps between rows. */
  if (newSize >= MIN_FONT_SIZE && newSize >= (detailViewGetCellYPadding(detailView) * 2))
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView);
    }
}


/* Zoom the detail view in/out */
void zoomDetailView(GtkWidget *detailView, const gboolean zoomIn)
{
  /* Remember the size of the sequence column. This gets reset to 0 when we change the font. */
  GtkWidget *firstTree = detailViewGetFirstTree(detailView);
  int colWidth = UNSET_INT;  
  if (firstTree)
    {
      GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(firstTree), SEQUENCE_COL);
      colWidth = gtk_tree_view_column_get_width(sequenceCol);
    }
  
  if (zoomIn)
    {
      incrementFontSize(detailView);
    }
  else
    {
      decrementFontSize(detailView);
    }
  
  if (firstTree && colWidth != UNSET_INT)
    {
      updateSeqColumnSize(firstTree, colWidth);
    }  
}


/* Get the text displayed in the user feedback box based on the given MSPs sequence name
 * (if an MSP is given), and also the currently-selected base index (if there is one). 
 * The string returned by this function must be free'd with g_free. */
static char* getFeedbackText(GtkWidget *detailView, const char *seqName, const int numSeqsSelected)
{
  /* The info we need to find... */
  int qIdx = UNSET_INT; /* index into the ref sequence. Ref seq is always a DNA seq */
  int sIdx = UNSET_INT; /* index into the match sequence. Will be coords into the peptide sequence if showing peptide matches */
  
  const BlxSeqType seqType = detailViewGetSeqType(detailView);
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const int selectedDnaBaseIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  
  /* See if a base is selected. */
  qIdx = selectedDnaBaseIdx;
  
  /* Make sure we have enough space for all the bits to go in the result text */
  int msgLen = numDigitsInInt(qIdx);
  
  /* Find the sequence name text (or some default text to indicate that a sequence is not selected) */
  const char *noSeqText = numSeqsSelected > 0 ? MULTIPLE_SUBJECTS_SELECTED_TEXT : NO_SUBJECT_SELECTED_TEXT;
  
  if (!seqName)
    {
      msgLen += strlen(noSeqText);
    }
  else
    {
      msgLen += strlen(seqName);
      
      /* If a q index is selected, see if there is a valid base at that index 
       * for any of the MSPs for the selected sequence. */
      if (qIdx != UNSET_INT)
	{
	  GList *mspListItem = detailViewGetSequenceMsps(detailView, seqName);
	  
	  for ( ; mspListItem; mspListItem = mspListItem->next)
	    {
	      MSP *msp = (MSP*)(mspListItem->data);
	      sIdx = getMatchIdxFromDisplayIdx(msp, 
					       qIdx, 
					       mspGetRefFrame(msp, seqType), 
					       mspGetRefStrand(msp), 
					       rightToLeft, 
					       BLXSEQ_DNA, /* q index has already been converted to DNA coords */ 
					       numFrames);
	      
	      if (sIdx != UNSET_INT)
		{
		  msgLen += numDigitsInInt(sIdx);
		  break;
		}
	    }
	}
    }
  
  msgLen += 10; /* for format text */
  char *messageText = g_malloc(sizeof(char) * msgLen);
  
  /* Create the message text. */
  if (seqName != NULL && qIdx != UNSET_INT)
    {
      sprintf(messageText, "%d   %s: %d", qIdx, seqName, sIdx);		/* display all */
    }
  else if (seqName != NULL)
    {
      sprintf(messageText, "%s", seqName);				/* just the name */
    }
  else if (qIdx != UNSET_INT)
    {
      sprintf(messageText, "%d   %s", qIdx, noSeqText);  /* just the index */
    }
  else
    {
      sprintf(messageText, " ");
    }
  
  return messageText;
}


/* Updates the feedback box with info about any currently-selected MSPs. This 
 * currently assumes single-MSP selections, but could be extended in the future 
 * to display, say, summary information about multiple MSPs. */
void updateFeedbackBox(GtkWidget *detailView)
{
  char *messageText = NULL;

  GList *selectedSeqs = mainWindowGetSelectedSeqs(detailViewGetMainWindow(detailView));
  const int numSeqsSelected = g_list_length(selectedSeqs);
  
  if (numSeqsSelected == 1) /* currently we only properly handle single sequence selection */
    {
      const char *seqName = (const char*)(selectedSeqs->data);
      messageText = getFeedbackText(detailView, seqName, numSeqsSelected);
    }
  else
    {
      /* 0 or multiple MSPs selected. Just see if a base index is selected. */
      messageText = getFeedbackText(detailView, NULL, numSeqsSelected);
    }
  
  GtkWidget *feedbackBox = detailViewGetFeedbackBox(detailView);
  gtk_entry_set_text(GTK_ENTRY(feedbackBox), messageText);
  gtk_widget_queue_draw(feedbackBox);
  
  g_free(messageText);
}


/* If the selected base index is outside the current display range, scroll to keep it in range. */
static void scrollToKeepSelectionInRange(GtkWidget *detailView)
{
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  
  if (selectedBaseIdx != UNSET_INT)
    {
      if (selectedBaseIdx < displayRange->min)
	{
	  setDetailViewStartIdx(detailView, selectedBaseIdx, detailViewGetSeqType(detailView));
	}
      else if (selectedBaseIdx > displayRange->max)
	{
	  setDetailViewEndIdx(detailView, selectedBaseIdx, detailViewGetSeqType(detailView));
	}
    }
}


/* Squash matches switches to the condensed view of the trees (where multiple
* MSPs in the same sequence appear on the same row) if squash is true, or reverts
* to the expanded version (where each MSP has its own row) if squash is false. */
void detailViewSquashMatches(GtkWidget *detailView, const gboolean squash)
{
  if (squash)
    {
      callFuncOnAllDetailViewTrees(detailView, treeSquashMatches);
    }
  else
    {
      callFuncOnAllDetailViewTrees(detailView, treeUnsquashMatches);
    }

  gtk_widget_queue_draw(detailView);
}


/* Returns true if the matches are squashed */
gboolean detailViewGetMatchesSquashed(GtkWidget *detailView)
{
  /* Just check the state of any one of the trees */
  return treeGetMatchesSquashed(detailViewGetFirstTree(detailView));
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

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
      gtk_widget_queue_draw(detailView);

      /* Update the big picture because the highlight box has moved (and changed size) */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture, TRUE);
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

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
      gtk_widget_queue_draw(detailView);

      /* Update the big picture because the highlight box has moved */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture, TRUE);
    }
}


/***********************************************************
 *            Detail View toolbar events                   *
 ***********************************************************/

static void onZoomInDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  zoomDetailView(detailView, TRUE);
}

static void onZoomOutDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  zoomDetailView(detailView, FALSE);
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
  
  /* Cache these results, because we use them often for calculations */
  detailViewCacheFontSize(detailView, charWidth, charHeight);
  
  /* Set the row height. Subtract the padding between the cell's actual area and
   * its background area. We will render at the background area's height, so that
   * we draw over the "gaps" between the cells, giving the impression of no gaps. */
  gint rowHeight = charHeight - (detailViewGetCellYPadding(detailView) * 2);

  GtkCellRenderer *renderer = detailViewGetRenderer(detailView);
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

GHashTable *detailViewGetSeqTable(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->seqTable : NULL;
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

static GList* detailViewGetColumnList(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->columnList : NULL;
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

GdkColor* detailViewGetRefSeqColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->refSeqColourSelected : &properties->refSeqColour;
}

GdkColor* detailViewGetMatchColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->matchColourSelected : &properties->matchColour;
}

GdkColor* detailViewGetMismatchColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->mismatchColourSelected : &properties->mismatchColour;
}

GdkColor* detailViewGetConsColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->consColourSelected : &properties->consColour;
}

GdkColor* detailViewGetExonColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->exonColourSelected : &properties->exonColour;
}

GdkColor* detailViewGetGapColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->gapColourSelected : &properties->gapColour;
}

static GdkColor* detailViewGetTripletHighlightColour(GtkWidget *detailView, const gboolean isSelectedDnaBase)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return isSelectedDnaBase ? &properties->highlightDnaBaseColour : &properties->highlightTripletColour;
}

GdkColor* detailViewGetExonBoundaryColour(GtkWidget *detailView, const gboolean isStart)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return isStart ? &properties->exonBoundaryColourStart : &properties->exonBoundaryColourEnd;
}

int detailViewGetExonBoundaryWidth(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties->exonBoundaryLineWidth;
}

GdkLineStyle detailViewGetExonBoundaryStyle(GtkWidget *detailView, const gboolean isStart)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return isStart ? properties->exonBoundaryLineStyleStart : properties->exonBoundaryLineStyleEnd;
}


GList *detailViewGetStrandTrees(GtkWidget *detailView, const Strand activeStrand)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return (activeStrand == FORWARD_STRAND) ? properties->fwdStrandTrees : properties->revStrandTrees;
}


int detailViewGetCellXPadding(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->cellXPadding : 0;
}

int detailViewGetCellYPadding(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->cellYPadding : 0;
}

/* Get the character width. */
int detailViewGetCharWidth(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->charWidth : 0;
}

int detailViewGetCharHeight(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->charHeight : 0;
}

static void detailViewCacheFontSize(GtkWidget *detailView, int charWidth, int charHeight)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->charWidth = charWidth;
  properties->charHeight = charHeight;
}


/* Get the outermost container for the GtkTreeView for the given frame on the given strand */
GtkWidget* detailViewGetTreeContainer(GtkWidget *detailView, const Strand activeStrand, const int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *result = NULL;
  
  /* Get the list of trees for the relevant strand */
  GList *list = (activeStrand == FORWARD_STRAND) ? properties->fwdStrandTrees : properties->revStrandTrees;

  /* The list should be sorted in order of frame number, and we should not be requesting a frame number
   * greater than the number of items in the list */
  if (frame <= g_list_length(list))
    {
      GList *listItem = list;
      int count = 1;
      
      for ( ; listItem && count < frame; ++count)
	{
	  listItem = listItem->next;
	}

      if (count == frame && 
	  GTK_IS_WIDGET(listItem->data) && 
	  (widgetIsTree(GTK_WIDGET(listItem->data)) || widgetIsTreeContainer(GTK_WIDGET(listItem->data))))
	{
	  result = GTK_WIDGET(listItem->data);
	}
    }
  
  return result;
}


/* Extract the actual GtkTreeView for the given frame on the given strand */
GtkWidget* detailViewGetTree(GtkWidget *detailView, const Strand activeStrand, const int frame)
{
  GtkWidget *result = NULL;
  GtkWidget *treeContainer = detailViewGetTreeContainer(detailView, activeStrand, frame);
  
  if (treeContainer && widgetIsTree(treeContainer))
    {
      result = treeContainer;
    }
  else if (treeContainer)
    {
      result = treeContainerGetTree(GTK_CONTAINER(treeContainer));
    }
  
  if (!result)
    {
      printf("Tree not found for '%s' strand, frame '%d'. Returning NULL.", ((activeStrand == FORWARD_STRAND) ? "forward" : "reverse"), frame);
    }
  
  return result;
}

/* Get the first visible tree in the 'current' list of trees (i.e. the forward 
 * strand list by default, or the reverse strand list if strands are toggled). */
static GtkWidget* detailViewGetFirstTree(GtkWidget *detailView)
{
  const gboolean toggled = detailViewGetStrandsToggled(detailView);
  const Strand activeStrand = toggled ? REVERSE_STRAND : FORWARD_STRAND;
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  
  GtkWidget *result = NULL;
  
  if (detailViewGetSeqType(detailView) == BLXSEQ_PEPTIDE)
    {
      /* First tree might be hidden, so loop until we find a visible one */
      int frame = 1;
      
      for ( ; frame <= numFrames; ++frame)
	{
	  GtkWidget *tree = detailViewGetTree(detailView, activeStrand, frame);
	  
	  if (tree && GTK_WIDGET_VISIBLE(tree))
	    {
	      result = tree;
	      break;
	    }
	}
    }
  else
    {
      /* Try the active strand, and if that's hidden try the other strand */
      result = detailViewGetTree(detailView, activeStrand, 1);
      
      if (!result || !GTK_WIDGET_VISIBLE(result))
	{
	  result = detailViewGetTree(detailView, toggled ? FORWARD_STRAND : REVERSE_STRAND, 1);
	}
    }
  
  return result;
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

static GtkWidget* detailViewGetFeedbackBox(GtkWidget *detailView)
{
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  return detailViewProperties->feedbackBox;
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

static int detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedDnaBaseIdx : UNSET_INT;
}

static int detailViewGetSelectedFrame(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedFrame : UNSET_INT;
}

GtkSortType detailViewGetSortType(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->sortType : GTK_SORT_ASCENDING;
}

/* Return a list of all MSPs that have the given match sequence name */
GList *detailViewGetSequenceMsps(GtkWidget *detailView, const char *seqName)
{
  SubjectSequence *subjectSeq = detailViewGetSequenceFromName(detailView, seqName);
  GList *result = subjectSeq ? subjectSeq->mspList : NULL;
  return result;
}


/* Return the SubjectSequence struct for the sequence with the given name */
SubjectSequence* detailViewGetSequenceFromName(GtkWidget *detailView, const char *seqName)
{
  SubjectSequence *result = NULL;
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  if (properties)
    {
      result = (SubjectSequence*)(g_hash_table_lookup(properties->seqTable, seqName));
    }
  
  return result;
}


/* Set the selected base index. Performs any required refreshes */
void detailViewSetSelectedBaseIdx(GtkWidget *detailView, const int selectedBaseIdx, const int frame, const int baseNum, const gboolean allowScroll)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  properties->selectedBaseIdx = selectedBaseIdx;
  properties->selectedFrame = frame;
  properties->selectedBaseNum = baseNum;
  
  /* For protein matches, calculate the base index in terms of the DNA sequence and cache it */
  properties->selectedDnaBaseIdx = convertPeptideToDna(selectedBaseIdx, frame, baseNum, detailViewGetNumReadingFrames(detailView));

  if (allowScroll)
    {
      scrollToKeepSelectionInRange(detailView);
    }

  /* Update the feedback box */
  updateFeedbackBox(detailView);

  /* Update the headers so that the newly-selected index is highlighted */
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
  gtk_widget_queue_draw(detailView);
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
				       GList *columnList,
				       GtkAdjustment *adjustment, 
				       BlxSeqType seqType,
				       int numReadingFrames)
{
  if (detailView)
    { 
      DetailViewProperties *properties = g_malloc(sizeof *properties);

      /* Find a fixed-width font */
      const char *fontFamily = findDetailViewFont(detailView);
      PangoFontDescription *fontDesc = pango_font_description_copy(detailView->style->font_desc);
      pango_font_description_set_family(fontDesc, fontFamily);
      
      properties->mainWindow = mainWindow;
      properties->renderer = renderer;
      properties->fwdStrandTrees = fwdStrandTrees;
      properties->revStrandTrees = revStrandTrees;
      properties->feedbackBox = feedbackBox;
      properties->columnList = columnList;
      properties->adjustment = adjustment;
      properties->seqType = seqType;
      properties->numReadingFrames = numReadingFrames;
      properties->displayRange.min = 1;
      properties->displayRange.max = 1;
      properties->selectedBaseIdx = UNSET_INT;
      properties->selectedBaseNum = UNSET_INT;
      properties->selectedFrame = UNSET_INT;
      properties->selectedDnaBaseIdx = UNSET_INT;
      properties->fontDesc = fontDesc;
      properties->charWidth = 0;
      properties->charHeight = 0;
      properties->seqTable = g_hash_table_new(g_str_hash, g_str_equal);
      properties->sortType = GTK_SORT_ASCENDING;

      /* Find the padding between the background area of the tree cells and the actual
       * drawing area. This is used to render the full height of the background area, so
       * that we don't have gaps between rows. */
      if (fwdStrandTrees)
	{
	  GtkWidget *widget = GTK_WIDGET(fwdStrandTrees->data);

	  GtkWidget *tree = widgetIsTree(widget) ? widget : treeContainerGetTree(GTK_CONTAINER(widget));
	  gtk_widget_realize(tree); /* must realize tree to pick up any overriden style properties */
	    
	  /* I can't get this to work properly using gtk_tree_view_get_cell_area etc. The cell
	   * area and background area come out wrong - perhaps because I'm not using a real
	   * row path. After a bit of experimentation it seems that the y padding is always
	   * related to the vertical-separator as follows: ypad = trunc(vseparator / 2) + 1 */
	  gint vertical_separator, horizontal_separator;
	  gtk_widget_style_get (tree, "vertical-separator", &vertical_separator, NULL);
	  gtk_widget_style_get (tree, "horizontal-separator", &horizontal_separator, NULL);
	  
	  properties->cellXPadding = (horizontal_separator / 2) + 1;
	  properties->cellYPadding = (vertical_separator / 2) + 1;
	}
      
      properties->refSeqColour		  = getGdkColor(GDK_YELLOW);
      properties->refSeqColourSelected	  = getGdkColor(GDK_DARK_YELLOW);
      properties->matchColour		  = getGdkColor(GDK_TURQUOISE);
      properties->matchColourSelected	  = getGdkColor(GDK_DARK_CYAN);
      properties->mismatchColour	  = getGdkColor(GDK_GREY);
      properties->mismatchColourSelected  = getGdkColor(GDK_DARK_GREY);
      properties->consColour		  = getGdkColor(GDK_LIGHT_STEEL_BLUE);
      properties->consColourSelected	  = getGdkColor(GDK_ROYAL_BLUE);
      properties->exonColour		  = getGdkColor(GDK_YELLOW);
      properties->exonColourSelected	  = getGdkColor(GDK_DARK_YELLOW);
      properties->gapColour		  = getGdkColor(GDK_GREY);
      properties->gapColourSelected	  = getGdkColor(GDK_DARK_GREY);
      properties->exonBoundaryColourStart = getGdkColor(GDK_BLUE);
      properties->exonBoundaryColourEnd	  = getGdkColor(GDK_DARK_BLUE);
      properties->highlightTripletColour  = getGdkColor(GDK_GREEN);
      properties->highlightDnaBaseColour  = getGdkColor(GDK_RED);

      properties->exonBoundaryLineWidth	  = 1;
      properties->exonBoundaryLineStyleStart = GDK_LINE_SOLID;
      properties->exonBoundaryLineStyleEnd   = GDK_LINE_SOLID;
      
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


/* If there are two trees and one is visible and the other hidden, toggle their hidden states */
static void swapTreeVisibility(GtkWidget *detailView)
{
  /* Only do anything for DNA matches, where we have 1 frame and both strands are visible. */
  if (detailViewGetSeqType(detailView) == BLXSEQ_DNA)
    {
      GtkWidget *tree1 = detailViewGetTreeContainer(detailView, FORWARD_STRAND, 1);
      GtkWidget *tree2 = detailViewGetTreeContainer(detailView, REVERSE_STRAND, 1);
      
      if (widgetGetHidden(tree1) != widgetGetHidden(tree2))
	{
	  widgetSetHidden(tree1, !widgetGetHidden(tree1));
	  widgetSetHidden(tree2, !widgetGetHidden(tree2));
	}
    }
}


/* If one grid is visible and the other hidden, toggle their hidden states */
static void swapGridVisibility(GtkWidget *bigPicture)
{
  GtkWidget *grid1 = bigPictureGetFwdGrid(bigPicture);
  GtkWidget *grid2 = bigPictureGetRevGrid(bigPicture);
  
  if (widgetGetHidden(grid1) != widgetGetHidden(grid2))
    {
      widgetSetHidden(grid1, !widgetGetHidden(grid1));
      widgetSetHidden(grid2, !widgetGetHidden(grid2));
    }
}


void toggleStrand(GtkWidget *detailView)
{
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(detailViewGetMainWindow(detailView));
  GtkWidget *bigPicture = mainWindowProperties->bigPicture;

  /* Update the flag */
  mainWindowProperties->strandsToggled = !mainWindowProperties->strandsToggled;

  /* If one grid/tree is hidden and the other visible, toggle which is hidden */
  swapTreeVisibility(detailView);
  swapGridVisibility(bigPicture);
  
  /* Toggle the order of the trees and grids. */
  refreshTreeOrder(detailView);
  refreshGridOrder(bigPicture);

  /* Redraw the tree and detail-view headers */
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
  gtk_widget_queue_draw(detailView);
  
  /* Redraw the grids and grid headers */
  callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders); /* must recalculate the highlight box first */
  callFuncOnAllBigPictureGrids(bigPicture, widgetClearCachedDrawable);
  widgetClearCachedDrawable(bigPictureGetGridHeader(bigPicture));
  gtk_widget_queue_draw(bigPicture);
}


/* Go to a user-specified coord (in terms of the DNA or peptide sequence as specified
 * by coordSeqType). */
void goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType)
{
  static gchar defaultInput[32] = "";
  
  /* Pop up a dialog to request a coord from the user */
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Goto which position: ", 
						  NULL, 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_NO_SEPARATOR,
						  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
						  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
						  NULL);

  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  GtkWidget *entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(contentArea), entry, TRUE, TRUE, 0);

  gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  gtk_widget_show(entry);
  

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const gchar *inputText = gtk_entry_get_text(GTK_ENTRY(entry));
      
      /* Convert the input string to an int */
      int requestedCoord = atoi(inputText);
      
      if (requestedCoord)
	{
	  /* Remember this input for next time */
	  sprintf(defaultInput, "%d", requestedCoord);

	  setDetailViewStartIdx(detailView, requestedCoord, coordSeqType);
	}
    }

  gtk_widget_destroy(dialog);
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

static void sortByGroupOrder(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByGroupOrder);
}

/* Find the next MSP (out the MSPs in this tree) whose start is the next closest to the
 * current start position of the display, searching only in the direction specified by 
 * the search criteria passed in the user data. */
static gboolean findNextMatchInTree(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  /* Loop through all MSPs in this tree row. */
  GList* mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      if (mspIsBlastMatch(msp) || mspIsExon(msp))
	{
	  MatchSearchData *searchData = (MatchSearchData*)data;

          int qSeqMin = min(msp->displayStart, msp->displayEnd);
          int qSeqMax = max(msp->displayStart, msp->displayEnd);
	  
	  /* Get the offset of this MSP from the current display start position */
	  int curOffset = UNSET_INT;
	  if (searchData->rightToLeft) /* display reversed */
	    {
	      curOffset = (qSeqMax - searchData->displayRange->max) * searchData->searchDirection;
	    }
	  else
	    {
	      curOffset = (qSeqMin - searchData->displayRange->min) * searchData->searchDirection;
	    }
	  
	  if (curOffset > 0 && (curOffset < searchData->smallestOffset || searchData->smallestOffset == UNSET_INT))
	    {
	      searchData->smallestOffset = curOffset;
	    }
	}
    }
  
  return FALSE;
}


/* Find and go to the next match (either left or right depending on the search flag) */
static void goToNextMatch(GtkWidget *detailView, const gboolean searchRight)
{
  /* If the display is toggled (i.e. showing right-to-left), then the active
   * strand is the reverse strand. This function currently only searchs for 
   * MSPs in the active strand. */
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  const Strand activeStrand = rightToLeft ? REVERSE_STRAND : FORWARD_STRAND;
  const int searchDirection = (searchRight != rightToLeft) ? 1 : -1;
  
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  
  /* Loop through each tree associated with this strand */
  MatchSearchData searchData = {displayRange, searchRight, searchDirection, rightToLeft, UNSET_INT};

  int frame = 1;
  for (  ; frame <= numFrames; ++frame)
    {
      /* Loop through all the MSPs in this tree and find the smallest offset
       * in the direction we're searching. */
      GtkWidget *tree = detailViewGetTree(detailView, activeStrand, frame);
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
    }
  
  if (searchData.smallestOffset != UNSET_INT)
    {
      /* Scroll the detail view by this offset amount. The offset coords are of
       * the same type as the display sequence type. */
      const int newStart = displayRange->min + (searchData.smallestOffset * searchDirection);
      const BlxSeqType seqType = detailViewGetSeqType(detailView);

      setDetailViewStartIdx(detailView, newStart, seqType);
    }
}

static void prevMatch(GtkWidget *detailView)
{
  goToNextMatch(detailView, FALSE);
}

static void nextMatch(GtkWidget *detailView)
{
  goToNextMatch(detailView, TRUE);
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
      else if (strcmp(val, SORT_BY_GROUP_ORDER_STRING) == 0)
	sortByGroupOrder(detailView);
    }
  
  g_free(val);
}


static void GHelp(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  displayHelp(detailViewGetMainWindow(detailView));
}

static void GSettings(GtkButton *button, gpointer data)
{
  blixemSettings();
}

static void GGoto(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  goToDetailViewCoord(detailView, BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
}

static void GprevMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  prevMatch(detailView);
}

static void GnextMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  nextMatch(detailView);
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
  toggleStrand(detailView);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* Adds the given column-header widget to the container, and adds its column info
 * to the columnList. Also sets the default width of the widget and the alignment. */
static void addHeaderColumn(GtkBox *container, 
			    ColumnId columnId, 
			    GtkWidget *specialWidget,
			    GtkCallback callbackFn, 
			    char *title,
			    char *propertyName,
			    const int defaultWidth,
			    const gboolean expand,
			    GList **columnList)
{
  /* Create a simple label for the header (unless already passed a special header widget) */
  GtkWidget *headerWidget = specialWidget ? specialWidget : gtk_label_new(title);
  gtk_box_pack_start(container, headerWidget, expand, TRUE, 0);
  gtk_widget_set_size_request(headerWidget, defaultWidth, -1);
  
  if (GTK_IS_MISC(headerWidget))
    {
      /* Align the text bottom-left */
      gtk_misc_set_alignment(GTK_MISC(headerWidget), 0.0, 1.0);
    }

  /* Create the column info and place it in the list */
  DetailViewColumnInfo *columnInfo = g_malloc(sizeof(DetailViewColumnInfo));
  columnInfo->columnId = columnId;
  columnInfo->headerWidget = headerWidget;
  columnInfo->refreshFunc = callbackFn,
  columnInfo->title = title;
  columnInfo->propertyName = propertyName;
  columnInfo->width = defaultWidth;
  *columnList = g_list_append(*columnList, columnInfo);
}


/* Create the header bar for the detail view. This contains the labels for the
 * detail-view trees (since we only want one label at the top, rather than 
 * individual labels for each tree). For protein sequence matches, the header
 * for the sequence column will also show the DNA sequence (separated into reading
 * frames). 
 * Column data is compiled into the detailViewColumns return argument. */
static GtkWidget* createDetailViewHeader(GtkWidget *detailView, 
					 const BlxSeqType seqType, 
					 const int numReadingFrames,
					 GList **columnList)
{
  GtkBox *header = GTK_BOX(gtk_hbox_new(FALSE, 0));

  /* The sequence column has a special header and callback when we're dealing with peptide sequences */
  GtkWidget *seqHeader = createSequenceColHeader(detailView, seqType, numReadingFrames);
  GtkCallback seqCallback = (seqType == BLXSEQ_PEPTIDE) ? refreshDetailViewSequenceColHeader : NULL;
  
  /* The start and end columns have callbacks to switch the start/end text when display is toggled */
  GtkCallback startCallback = refreshDetailViewStartColHeader;
  GtkCallback endCallback = refreshDetailViewEndColHeader;
  
  const int endColWidth = END_COLUMN_DEFAULT_WIDTH;
  
  addHeaderColumn(header, S_NAME_COL, NULL,	  NULL,		 NAME_COLUMN_HEADER_TEXT,  NAME_COLUMN_PROPERTY_NAME,  NAME_COLUMN_DEFAULT_WIDTH,  FALSE, columnList);
  addHeaderColumn(header, SCORE_COL,  NULL,	  NULL,		 SCORE_COLUMN_HEADER_TEXT, SCORE_COLUMN_PROPERTY_NAME, SCORE_COLUMN_DEFAULT_WIDTH, FALSE, columnList);
  addHeaderColumn(header, ID_COL,     NULL,	  NULL,		 ID_COLUMN_HEADER_TEXT,    ID_COLUMN_PROPERTY_NAME,    ID_COLUMN_DEFAULT_WIDTH,    FALSE, columnList);
  addHeaderColumn(header, START_COL,  NULL,	  startCallback, START_COLUMN_HEADER_TEXT, START_COLUMN_PROPERTY_NAME, START_COLUMN_DEFAULT_WIDTH, FALSE, columnList);
  addHeaderColumn(header, SEQUENCE_COL, seqHeader, seqCallback,	 SEQ_COLUMN_HEADER_TEXT,   SEQ_COLUMN_PROPERTY_NAME,   SEQ_COLUMN_DEFAULT_WIDTH,   TRUE,  columnList);
  addHeaderColumn(header, END_COL,    NULL,	  endCallback,	 END_COLUMN_HEADER_TEXT,   END_COLUMN_PROPERTY_NAME,   endColWidth,		  FALSE, columnList);

  return GTK_WIDGET(header);
}


/* Create a custom header widget for the sequence column for protein matches (this is where
 * we will display the triplets that make up the codons.) Returns NULL for DNA matches. */
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
	  GtkWidget *line = createLabel(NULL, 0.0, 0.0, TRUE, TRUE);
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
  sortList = g_list_append(sortList, SORT_BY_GROUP_ORDER_STRING);
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
  makeToolbarButton(toolbar, "Help",	 "Don't Panic",			 (GtkSignalFunc)GHelp,		  detailView);

  /* Combo box for sorting */
  createSortBox(toolbar, detailView);
  
  /* Settings button */
  makeToolbarButton(toolbar, "Settings", "Open the Preferences Window",  (GtkSignalFunc)GSettings,	  NULL);
  
  /* Navigation buttons */
  makeToolbarButton(toolbar, "Goto",	 "Go to specified co-ordinates", (GtkSignalFunc)GGoto,		  detailView);
  makeToolbarButton(toolbar, "< match",  "Next (leftward) match",	 (GtkSignalFunc)GprevMatch,	  detailView);
  makeToolbarButton(toolbar, "match >",  "Next (rightward) match",	 (GtkSignalFunc)GnextMatch,	  detailView);
  makeToolbarButton(toolbar, "<<",	 "Scroll leftward lots",	 (GtkSignalFunc)GscrollLeftBig,	  detailView);
  makeToolbarButton(toolbar, ">>",       "Scroll rightward lots",	 (GtkSignalFunc)GscrollRightBig,  detailView);
  makeToolbarButton(toolbar, "<",	 "Scroll leftward one base",	 (GtkSignalFunc)GscrollLeft1,	  detailView);
  makeToolbarButton(toolbar, ">",	 "Scroll rightward one base",	 (GtkSignalFunc)GscrollRight1,	  detailView);
  
  /* Strand toggle button */
  if (mode == BLXMODE_BLASTX || mode == BLXMODE_TBLASTX || mode == BLXMODE_BLASTN)
    {
      makeToolbarButton(toolbar, "Strand^v", "Toggle strand", (GtkSignalFunc)GToggleStrand, detailView);
    }
  
  /* Feedback box */
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
				GList *columnList,
				const char const *refSeqName,
				const int frame1,
				const int frame2)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, detailView, renderer, list1, columnList, seqType, refSeqName, frame1);
  GtkWidget *tree2 = createDetailViewTree(grid2, detailView, renderer, list2, columnList, seqType, refSeqName, frame2);
  
  if (container)
    {
      gtk_paned_pack1(GTK_PANED(container), tree1, TRUE, TRUE);
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
				  GList *columnList,
				  const char const *refSeqName)
{
  /* Add the frames in increasing order for the forward strand, decreasing order for
   * the reverse strand. This is the order they'll be added to the list, and they
   * are currently placed in the display in the same order. */
  const gboolean forward = gridGetStrand(grid) == FORWARD_STRAND;
  const int frame1 = forward ? 1 : 3;
  const int frame2 = 2;
  const int frame3 = forward ? 3 : 1;
  
  /* Create a tree for pane1 (but only add it to the detailView if instructed to).
   * The first tree has headers. */
  GtkWidget *tree1 = createDetailViewTree(grid, detailView, renderer, list, columnList, seqType, refSeqName, frame1);

  GtkWidget *nestedPanedWidget = NULL;
  if (addToDetailView)
    {
      gtk_paned_pack1(GTK_PANED(detailView), tree1, TRUE, TRUE);
      
      /* Create another paned widget and place it in pane 2 */
      nestedPanedWidget = gtk_vpaned_new();
      gtk_paned_pack2(GTK_PANED(detailView), nestedPanedWidget, TRUE, TRUE);
    }
  
  /* Create two more trees (and place them in the nested paned widget, if it is not null). 
   * Neither of these trees should have headers. */
  createTwoPanedTrees(detailView, nestedPanedWidget, renderer, grid, grid, list, list, seqType, columnList, refSeqName, frame2, frame3);
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
				  GList *columnList,
				  const char const *refSeqName)
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
			  columnList,
			  refSeqName,
			  1,
			  1);
    }
  else if (numReadingFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames for EACH strand (although only
       * one set of trees will be displayed at a time). */
      createThreePanedTrees(detailView, renderer, fwdStrandGrid, fwdStrandTrees, seqType, TRUE, columnList, refSeqName);
      createThreePanedTrees(detailView, renderer, revStrandGrid, revStrandTrees, seqType, FALSE, columnList, refSeqName);
    }
}


/* Add the MSPs to the detail-view trees */
void detailViewAddMspData(GtkWidget *detailView, MSP *mspList)
{
  /* First, create a data store for each tree so we have something to add our data to. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateBaseDataModel);

  /* Loop through each MSP, and add it to the correct tree based on its strand and reading frame */
  const BlxSeqType seqType = detailViewGetSeqType(detailView);
  MSP *msp = mspList;
  
  for ( ; msp; msp = msp->next)
    {
      /* Find the tree that this MSP should belong to based on its reading frame and strand */
      Strand activeStrand = (msp->qframe[1] == '+') ? FORWARD_STRAND : REVERSE_STRAND;
      const int frame = mspGetRefFrame(msp, seqType);
      GtkWidget *tree = detailViewGetTree(detailView, activeStrand, frame);

      addMspToTree(tree, msp);
      
      /* Add the MSP to the hash table that will group MSPs by sequence name. */
      GHashTable *seqTable = detailViewGetSeqTable(detailView);
      addMspToHashTable(seqTable, msp, msp->sname);
    }
  
  /* Finally, create a custom-filtered version of the data store for each tree. We do 
   * this AFTER adding the data so that it doesn't try to re-filter every time we add a row. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateFilteredDataModel);
  
  /* Also create a second data store that will store one sequence per row (as opposed to one
   * MSP per row). This data store will be switched in when the user selects 'squash matches'. */
  callFuncOnAllDetailViewTrees(detailView, addSequencesToTree);
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
			    GtkWidget *container,
			    GtkAdjustment *adjustment, 
			    GtkWidget *fwdStrandGrid, 
			    GtkWidget *revStrandGrid,
			    MSP *mspList,
			    BlxBlastMode mode,
			    BlxSeqType seqType,
			    int numReadingFrames,
			    const char const *refSeqName)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so that we can identify them.) */
  GtkWidget *detailView = gtk_vpaned_new();

  /* Create a custom cell renderer to render the sequences in the detail view */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();

  /* Create the toolbar. We need to remember the feedback box. */
  GtkWidget *feedbackBox = NULL;
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, mode, &feedbackBox);

  /* Create the header, and compile a list of columns */
  GList *columnList = NULL;
  GtkWidget *header = createDetailViewHeader(detailView, seqType, numReadingFrames, &columnList);
  
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
			columnList,
			refSeqName);
  
  /* Put everything in a vbox, and pack it into the main window. */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), header, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), detailView, TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(container), vbox, TRUE, TRUE, 0);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(detailView), "size-allocate", G_CALLBACK(onSizeAllocateDetailView), NULL);
  
  detailViewCreateProperties(detailView, 
			     mainWindow, 
			     renderer,
			     fwdStrandTrees,
			     revStrandTrees,
			     feedbackBox, 
			     columnList,
			     adjustment, 
			     seqType,
			     numReadingFrames);
  
  return detailView;
}

