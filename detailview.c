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
#define SORT_BY_POS_STRING		"Position"
#define SORT_BY_GROUP_ORDER_STRING	"Group"
#define FONT_INCREMENT_SIZE		1
#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20
#define NO_SUBJECT_SELECTED_TEXT	"<no subject selected>"
#define MULTIPLE_SUBJECTS_SELECTED_TEXT	"<multiple subjects selected>"

typedef enum {SORT_TYPE_COL, SORT_TEXT_COL, N_SORT_COLUMNS} SortColumns;


typedef struct 
  {
    const int startCoord;		/* the display coord to start searching from */
    const gboolean searchRight;		/* search towards the right or left */
    const int searchDirection;		/* multiplier to add/subtract depending on whether searching right/left */
    const gboolean rightToLeft;		/* true if the display is reversed */
    const BlxSeqType seqType;		/* whether viewing DNA or peptide seqs */
    const int numFrames;		/* number of reading frames */
    const IntRange const *refSeqRange;	/* the full range of the reference sequence */
    GList *seqNameList;			/* only search matches in these sequences */

    int frame;				/* the frame of the current tree we're looking at MSPs in */
    int offset;				/* the offset of the found MSP */
    int isStart;			/* whether the offset is to the start or end coord of the MSP */
    int foundFrame;			/* the frame of the MSP we chose */
  } MatchSearchData;


/* Local function declarations */
static GtkWidget*	      detailViewGetFirstTree(GtkWidget *detailView);
static GtkWidget*	      detailViewGetBigPicture(GtkWidget *detailView);
static GtkWidget*	      detailViewGetFeedbackBox(GtkWidget *detailView);
static int		      detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView);
static int		      detailViewGetSelectedFrame(GtkWidget *detailView);
static GdkColor*	      detailViewGetTripletHighlightColour(GtkWidget *detailView, const gboolean isSelectedDnaBase);

static int		      seqColHeaderGetFrameNumber(GtkWidget *header);

static void		      detailViewCacheFontSize(GtkWidget *detailView, int charWidth, int charHeight);
static GtkToolItem*	      addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget);
static gboolean		      widgetIsTree(GtkWidget *widget);
static gboolean		      widgetIsTreeContainer(GtkWidget *widget);
static void		      updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc);
static void		      refreshDetailViewHeaders(GtkWidget *detailView);
static GtkWidget*	      createSeqColHeader(GtkWidget *detailView, const BlxSeqType seqType, const int numReadingFrames);
static const char*	      findDetailViewFont(GtkWidget *detailView);
static void		      setDetailViewScrollPos(GtkAdjustment *adjustment, int value);

/***********************************************************
 *		       Utility functions                   *
 ***********************************************************/

/* Return the width of the column with the given column id */
int detailViewGetColumnWidth(GtkWidget *detailView, const ColumnId columnId)
{
  int result = 0;

  DetailViewColumnInfo *columnInfo = detailViewGetColumnInfo(detailView, columnId);
  if (columnInfo)
    {
      result = columnInfo->width;
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

  gint family;
  for (family = 0 ; (family < n_families && !found_most_preferred) ; family++)
    {
      const gchar *name = pango_font_family_get_name(families[family]) ;
      
      /* Look for this font family in our list of preferred families */
      GList *pref = g_list_first(pref_families) ;
      gint current = 1;
      
      while(pref)
	{
	  char *pref_font = (char *)pref->data ;
	  
	  if (g_ascii_strncasecmp(name, pref_font, strlen(pref_font)) == 0
#if GLIB_MAJOR_VERSION >= 2 && GLIB_MINOR_VERSION >= 6
	      && pango_font_family_is_monospace(families[family])
#endif
	      )
	    {
	      /* We prefer ones nearer the start of the list */
              if(current <= most_preferred)
                {
		  most_preferred = current;
		  match_family = families[family];

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
 * range (within bounds). the coord should be in terms of display coords */
void setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, coord);
}


/* Scroll the detail view so that the given coord is at the end of the display
 * range (within bounds). the coord should be in terms of display coords */
void setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  /* Get the new start coord */
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  const int displayLen = getRangeLength(displayRange);
  int newStart = coord - displayLen + 1;

  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, newStart);
}


/* Update the scroll position of the adjustment to the given value. Does bounds checking. */
static void setDetailViewScrollPos(GtkAdjustment *adjustment, int value)
{  
  /* bounds checking */
  int maxValue = adjustment->upper - adjustment->page_size + 1;
  
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
  setDetailViewScrollPos(adjustment, adjustment->value - 1);
}

void scrollDetailViewRight1(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + 1);
}

/* Scroll by one step increment */
void scrollDetailViewLeftStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value - adjustment->step_increment);
}

void scrollDetailViewRightStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + adjustment->step_increment);
}

/* Scroll by one page size */
void scrollDetailViewLeftPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value - adjustment->page_increment);
}

void scrollDetailViewRightPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + adjustment->page_increment);
}


/* Calculate the number of bases that can be displayed in the sequence column */
static int calcNumBasesInSequenceColumn(GtkWidget *detailView)
{
  /* Find the width of the sequence column */
  int colWidth = detailViewGetColumnWidth(detailView, SEQUENCE_COL);
  
  /* Don't include the cell padding area */
  GtkCellRenderer *renderer = detailViewGetRenderer(detailView);
  colWidth -= (2 * renderer->xpad) + (2 * renderer->xalign);
  
  /* Return the number of whole characters that fit in the column. */
  gint charWidth = detailViewGetCharWidth(detailView);
  int numChars = (int)((double)colWidth / (double)charWidth);
  
  return numChars;
}


/* This should be called when the width of the sequence column has changed (or the
 * size of the font has changed). This function will adjust the scroll range of our
 * custom scroll adjustment so that it represents the range that can be displayed 
 * in the new column width. */
void updateSeqColumnSize(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  
  if (adjustment)
    {
      int newPageSize = calcNumBasesInSequenceColumn(detailView);
      
      /* Only trigger updates if things have actually changed */
      if (newPageSize != adjustment->page_size)
	{
	  adjustment->page_size = newPageSize;
	  adjustment->page_increment = newPageSize;
	  
	  /* Reset the display range so that it is between the scrollbar min and max. Try to keep
	   * it centred on the same base. The base index is in terms of the display range coords, 
	   * so the sequence type of the coord is whatever the display sequence type is. */
	  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);

	  /* First time through, both coords are set to the initial start coord. So, if
	   * they are the same, just use that coord as the start coord */
	  int newStart = displayRange->min;
	  
	  if (displayRange->min != displayRange->max)
	    {
	      int centre = getRangeCentre(displayRange);
	      int offset = roundNearest((double)adjustment->page_size / 2.0);
	      newStart = centre - offset;
	    }
	      
	  const BlxSeqType seqType = detailViewGetSeqType(detailView);
	  setDetailViewStartIdx(detailView, newStart, seqType);
	  
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
	  /* Two panes, one tree. Use 'first' flags to decide which pane to put it in */
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
static void refreshSeqColHeaderLine(GtkWidget *detailView,
					 GtkWidget *lineWidget, 
					 int frame, /* local copy */
					 char *refSeq)
{
  if (!GTK_IS_LABEL(lineWidget))
    {
      messcrash("Unexpected widget type in detail view header: expected label");
    }

  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  const IntRange const *refSeqRange = detailViewGetRefSeqRange(detailView);
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  const int selectedDnaBaseIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  const int selectedFrame = detailViewGetSelectedFrame(detailView);
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  const BlxSeqType seqType = detailViewGetSeqType(detailView);

  /* Loop forward/backward through the display range depending on which strand we're viewing.
   * We need to convert display coords to actual coords on the ref sequence */
  const int qIdx1 = convertDisplayIdxToDnaIdx(displayRange->min, seqType, frame, 1, numFrames, rightToLeft, refSeqRange);	      /* 1st base in frame */
  const int qIdx2 = convertDisplayIdxToDnaIdx(displayRange->max, seqType, frame, numFrames, numFrames, rightToLeft, refSeqRange); /* last base in frame */
  
  IntRange qRange;
  qRange.min = min(qIdx1, qIdx2);
  qRange.max = max(qIdx1, qIdx2);
  
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
						  refSeqRange, 
						  BLXSEQ_DNA /* header always shows DNA bases */
						 );

      /* If the base (or its peptide) is selected, we need to highlight it */
      const int peptideIdx = convertDnaIdxToDisplayIdx(qIdx, seqType, selectedFrame, numFrames, rightToLeft, refSeqRange, NULL);
      if (peptideIdx == selectedBaseIdx)
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
static void refreshStartColHeader(GtkWidget *header, gpointer data)
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
static void refreshEndColHeader(GtkWidget *header, gpointer data)
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
static void refreshSeqColHeader(GtkWidget *header, gpointer data)
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
	  refreshSeqColHeaderLine(detailView, lineWidget, frame, refSeq);
	  
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


/* Resize the detail view header widgets. Should be called whenever a column is resized. */
void resizeDetailViewHeaders(GtkWidget *detailView)
{
  GList *listItem = detailViewGetColumnList(detailView);
  
  for ( ; listItem; listItem = listItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)listItem->data;

      if (columnInfo->columnId == SEQUENCE_COL)
	{
	  /* For the sequence column, don't set the size request, or we won't be
	   * able to shrink the window. The sequence col header will be resized
	   * dynamically to fit its text, which is the width that we want anyway. */
	  gtk_widget_set_size_request(columnInfo->headerWidget, SEQ_COLUMN_DEFAULT_WIDTH, -1);
	}
      else
	{
	  /* For other columns, we can set the size request: they're small enough
	   * that we can live without the need to shrink below their sum widths. */
	  gtk_widget_set_size_request(columnInfo->headerWidget, columnInfo->width, -1);
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
  if (zoomIn)
    {
      incrementFontSize(detailView);
    }
  else
    {
      decrementFontSize(detailView);
    }
  
  updateSeqColumnSize(detailView);
}


/* Get the text displayed in the user feedback box based on the given MSPs sequence name
 * (if an MSP is given), and also the currently-selected base index (if there is one). 
 * The string returned by this function must be free'd with g_free. */
static char* getFeedbackText(GtkWidget *detailView, const char *seqName, const int numSeqsSelected)
{
  /* The info we need to find... */
  int qIdx = UNSET_INT; /* index into the ref sequence. Ref seq is always a DNA seq */
  int sIdx = UNSET_INT; /* index into the match sequence. Will be coords into the peptide sequence if showing peptide matches */
  int sLen = UNSET_INT; /* the length of the match sequence */
  
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  const BlxSeqType seqType = mainWindowGetSeqType(mainWindow);
  const int numFrames = mainWindowGetNumReadingFrames(mainWindow);
  const int selectedDnaBaseIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(mainWindow);
  
  /* See if a base is selected. */
  qIdx = selectedDnaBaseIdx;
  
  /* Find the sequence name text (or some default text to indicate that a sequence is not selected) */
  const char *noSeqText = numSeqsSelected > 0 ? MULTIPLE_SUBJECTS_SELECTED_TEXT : NO_SUBJECT_SELECTED_TEXT;
  
  if (seqName)
    {
      GList *mspList = detailViewGetSequenceMsps(detailView, seqName);
      if (g_list_length(mspList) > 0)
	{
	  MSP *firstMsp = (MSP*)(mspList->data);
	  
	  if (firstMsp->sseq && firstMsp->sseq != mainWindowGetPaddingSeq(mainWindow))
	    {
	      sLen = strlen(firstMsp->sseq);
	    }

	  /* If a q index is selected, see if there is a valid base at that index 
	   * for any of the MSPs for the selected sequence. */
	  if (qIdx != UNSET_INT)
	    {
	      GList *mspListItem = mspList;
	      for ( ; mspListItem; mspListItem = mspListItem->next)
		{
		  MSP *msp = (MSP*)(mspListItem->data);
		  
		  sIdx = gapCoord(msp, qIdx, numFrames, mspGetRefFrame(msp, seqType), rightToLeft, NULL);

		  if (sIdx != UNSET_INT)
		    {
		      break;
		    }
		}
	    }
	}
    }

  /* Add all the bits into a text string */
  GString *resultString = g_string_sized_new(200); /* will be extended if we need more space */
  
  if (qIdx != UNSET_INT)
    {
      /* Adjust by the offset to give the real value for displaying */
      qIdx += mainWindowGetOffset(mainWindow);
      g_string_printf(resultString, "%d   ", qIdx);
    }
  
  if (seqName != NULL)
    { 
      g_string_append_printf(resultString, "%s", seqName);
    }
  else if (qIdx != UNSET_INT)
    {
      g_string_append_printf(resultString, "%s", noSeqText); 
    }
    
  if (sLen != UNSET_INT)
    {
      g_string_append_printf(resultString, "(%d)", sLen);
    }

  if (sIdx != UNSET_INT)
    {
      g_string_append_printf(resultString, " : %d", sIdx);
    }
  
  char *messageText = resultString->str;
  g_string_free(resultString, FALSE);
  
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


/* If the selected base index is outside the current display range, scroll to
 * keep it in range. We scroll by the minimum number of bases possible if 
 * scrollMinimum is true; otherwise we re-centre on the selection. */
static void scrollToKeepSelectionInRange(GtkWidget *detailView, const gboolean scrollMinimum)
{
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  
  if (selectedBaseIdx != UNSET_INT && !valueWithinRange(selectedBaseIdx, displayRange))
    {
      const BlxSeqType seqType = detailViewGetSeqType(detailView);
      
      if (scrollMinimum)
	{
	  if (selectedBaseIdx < displayRange->min)
	    {
	      setDetailViewStartIdx(detailView, selectedBaseIdx, seqType);
	    }
	  else if (selectedBaseIdx > displayRange->max)
	    {
	      setDetailViewEndIdx(detailView, selectedBaseIdx, seqType);
	    }
	}
      else
	{
	  const int newStart = selectedBaseIdx - (getRangeLength(displayRange) / 2);
	  setDetailViewStartIdx(detailView, newStart, seqType);
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


/* Set the value of the 'invert sort order' flag */
void detailViewSetSortInverted(GtkWidget *detailView, const gboolean invert)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->sortInverted = invert;

  callFuncOnAllDetailViewTrees(detailView, resortTree);
}


/* Set the value of the 'highlight differences' flag */
void detailViewSetHighlightDiffs(GtkWidget *detailView, const gboolean highlightDiffs)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->highlightDiffs = highlightDiffs;

  /* No data to recalculate, but we need to make all the trees redraw themselves */
  gtk_widget_queue_draw(detailView);
}


/***********************************************************
 *                    Detail view events                   *
 ***********************************************************/

static void onSizeAllocateDetailView(GtkWidget *detailView, GtkAllocation *allocation, gpointer data)
{
  /* The sequence column should be the only one that dynamically resizes as the window
   * window size changes. Find its new width and cache it. All trees should resize
   * in the same manner, so their columns should be the same size. */
  GtkWidget *tree = detailViewGetFirstTree(detailView);
  GtkTreeViewColumn *column = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), SEQUENCE_COL);

  DetailViewColumnInfo *columnInfo = detailViewGetColumnInfo(detailView, SEQUENCE_COL);
  columnInfo->width = gtk_tree_view_column_get_width(column);
      
  /* Perform updates required on the sequence col after its size has changed */
  updateSeqColumnSize(detailView);
}



/***********************************************************
 *              Detail view scrollbar events               *
 ***********************************************************/

static void onScrollRangeChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  int newStart = adjustment->value;
  int newEnd = adjustment->value + adjustment->page_size - 1;
  
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
  int newStart = adjustment->value;
  int newEnd = adjustment->value + adjustment->page_size - 1;

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

char** detailViewGetGeneticCode(GtkWidget *detailView)
{
  assertDetailView(detailView);
  GtkWidget *mainWindow = detailViewGetMainWindow(detailView);
  return mainWindowGetGeneticCode(mainWindow);
}

/* Get the list of columns */
GList* detailViewGetColumnList(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->columnList : NULL;
}

/* Get the column info for a particular column */
DetailViewColumnInfo *detailViewGetColumnInfo(GtkWidget *detailView, const ColumnId columnId)
{
  DetailViewColumnInfo *result = NULL;
  
  GList *listItem = detailViewGetColumnList(detailView);
  for ( ; listItem; listItem = listItem->next)
  {
    DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(listItem->data);
    if (columnInfo && columnInfo->columnId == columnId)
      {
	result = columnInfo;
	break;
      }
  }
  
  return result;
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

GdkColor* detailViewGetInsertionColour(GtkWidget *detailView, const gboolean selected)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return selected ? &properties->insertionColourSelected : &properties->insertionColour;
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

gboolean detailViewGetSortInverted(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->sortInverted : FALSE;
}

gboolean detailViewGetHighlightDiffs(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->highlightDiffs : FALSE;
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


/* Set the selected base index. Performs any required refreshes. Scrolls the view to
 * keep the selected base in view if allowScroll is true. (Such scrolling is by the minimum
 * number of bases necessary if scrollMinimum is true.) */
void detailViewSetSelectedBaseIdx(GtkWidget *detailView, 
				  const int selectedBaseIdx, 
				  const int frame, 
				  const int baseNum, 
				  const gboolean allowScroll,
				  const gboolean scrollMinimum)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  properties->selectedBaseIdx = selectedBaseIdx;
  properties->selectedFrame = frame;
  properties->selectedBaseNum = baseNum;
  
  /* For protein matches, calculate the base index in terms of the DNA sequence and cache it */
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  const IntRange const *refSeqRange = detailViewGetRefSeqRange(detailView);
  const BlxSeqType seqType = detailViewGetSeqType(detailView);
  
  properties->selectedDnaBaseIdx = convertDisplayIdxToDnaIdx(selectedBaseIdx, seqType, frame, baseNum, numFrames, rightToLeft, refSeqRange);

  if (allowScroll)
    {
      scrollToKeepSelectionInRange(detailView, scrollMinimum);
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

  /* N.B. Don't free the cell renderer, or it causes memory corruption. I'm not
   * sure what owns it - the columns it is added to? */
  
  if (properties->columnList > 0)
    {
      g_list_free(properties->columnList);
      properties->columnList = NULL;
    }
    
  if (properties->seqTable)
    {
      g_hash_table_unref(properties->seqTable);
      properties->seqTable = NULL;
    }
  
  if (properties->fwdStrandTrees)
    {
      g_list_free(properties->fwdStrandTrees);
      properties->fwdStrandTrees = NULL;
    }
  
  if (properties->revStrandTrees)
    {
      g_list_free(properties->revStrandTrees);
      properties->revStrandTrees = NULL;
    }

  if (properties->fontDesc)
    {
      pango_font_description_free(properties->fontDesc);
      properties->fontDesc = NULL;
    }
  
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
				       int numReadingFrames,
				       const int startCoord,
				       const gboolean sortInverted)
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
      properties->selectedBaseIdx = UNSET_INT;
      properties->selectedBaseNum = UNSET_INT;
      properties->selectedFrame = UNSET_INT;
      properties->selectedDnaBaseIdx = UNSET_INT;
      properties->fontDesc = fontDesc;
      properties->charWidth = 0;
      properties->charHeight = 0;
      properties->seqTable = g_hash_table_new(g_str_hash, g_str_equal);
      properties->sortInverted = sortInverted;
      properties->highlightDiffs = FALSE;

      /* Set initial display range to something valid but only 1 base wide. Then if we try to do any 
       * calculations on the range before it gets set properly, it won't have much work to do. */
      properties->displayRange.min = startCoord;
      properties->displayRange.max = startCoord;

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
      
      properties->refSeqColour		  = getGdkColor(DEFAULT_REF_SEQ_BG_COLOUR);
      properties->refSeqColourSelected	  = getGdkColor(GDK_DARK_YELLOW);
      properties->matchColour		  = getGdkColor(GDK_TURQUOISE);
      properties->matchColourSelected	  = getGdkColor(GDK_DARK_TURQUOISE);
      properties->mismatchColour	  = getGdkColor(GDK_GREY);
      properties->mismatchColourSelected  = getGdkColor(GDK_DARK_GREY);
      properties->consColour		  = getGdkColor(GDK_LIGHT_STEEL_BLUE);
      properties->consColourSelected	  = getGdkColor(GDK_STEEL_BLUE);
      properties->exonColour		  = getGdkColor(GDK_YELLOW);
      properties->exonColourSelected	  = getGdkColor(GDK_DARK_YELLOW);
      properties->insertionColour	  = getGdkColor(GDK_YELLOW);
      properties->insertionColourSelected = getGdkColor(GDK_DARK_YELLOW);
      properties->exonBoundaryColourStart = getGdkColor(GDK_BLUE);
      properties->exonBoundaryColourEnd	  = getGdkColor(GDK_DARK_BLUE);
      properties->highlightTripletColour  = getGdkColor(GDK_GREEN);
      properties->highlightDnaBaseColour  = getGdkColor(GDK_DARK_GREEN);

      properties->exonBoundaryLineWidth	  = 1;
      properties->exonBoundaryLineStyleStart = GDK_LINE_SOLID;
      properties->exonBoundaryLineStyleEnd   = GDK_LINE_SOLID;
      
      g_object_set_data(G_OBJECT(detailView), "DetailViewProperties", properties);
      g_signal_connect(G_OBJECT(detailView), "destroy", G_CALLBACK(onDestroyDetailView), NULL); 
    }
}


/* Get/set functions for the sequence column header frame number. (There is only
 * one property to set, so it's not worth having a struct for the properties.) */
static void seqColHeaderSetFrameNumber(GtkWidget *header, const int frame)
{
  g_object_set_data(G_OBJECT(header), "seqColHeaderFrameNumber", GINT_TO_POINTER(frame));
}

static int seqColHeaderGetFrameNumber(GtkWidget *header)
{
  return header ? GPOINTER_TO_INT(g_object_get_data(G_OBJECT(header), "seqColHeaderFrameNumber")) : UNSET_INT;
}


/***********************************************************
 *                     Events                              *
 ***********************************************************/

static int getBaseIndexAtColCoords(const int x, const int y, const int charWidth, const IntRange const *displayRange)
{
  int result = UNSET_INT;
  
  const int leftEdge = 0; /* to do: allow for padding? check upper bound? */
  
  if (x > leftEdge)
    {
      result = (int)(((double)x - leftEdge) / (double)charWidth);
    }
  
  result += displayRange->min;
  
  return result;
}


/* Select the nucleotide at the given x/y position in the given header (protein
 * matches only). */
static void selectClickedNucleotide(GtkWidget *header, GtkWidget *detailView, const int x, const int y)
{
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  
  int coord = getBaseIndexAtColCoords(x, y, detailViewGetCharWidth(detailView), detailViewGetDisplayRange(detailView));
  
  if (coord != UNSET_INT)
    {	
      /* Get the active (i.e. currently-selected) frame. Use frame 1 if none selected. */
      int frame = detailViewGetSelectedFrame(detailView);
      if (frame == UNSET_INT)
	{
	  frame = 1;
	}
      
      /* Get the base number of the clicked base within the active frame. The
       * base number is determined by which row in the header the mouse pointer 
       * is over: for frame 1, row 1 will give base 1; for frame 2, row 2 will
       * give base 1, etc. Start by getting the frame number for the clicked row: */
      int row = seqColHeaderGetFrameNumber(header);
      
      /* The header widget passed to this function is the originally-clicked widget.
       * If the pointer has dragged onto another row in the header, we can work out
       * which one based on the y coord. If it is outside altogether, take the 
       * topmost/bottommost row. (NB this assumes all rows in the header have the 
       * same height as this one - should be true because they're just simple labels). */
      if (y < 0 || y > header->allocation.height)
	{
	  const int offsetRows = floor((double)y / (double)header->allocation.height); /* can be negative */
	  row += offsetRows;
	  
	  if (row < 1)
	    {
	      row = 1;
	    }
	  else if (row > numFrames)
	    {
	      row = numFrames;
	    }
	}
      
      /* Calculate the base number based on the row and the currently-active frame */
      int baseNum = row - frame + 1;
      
      if (baseNum < 1)
	{
	  /* Cycle round if gone below the min base number */
	  baseNum += numFrames;
	  --coord;
	}
      
      detailViewSetSelectedBaseIdx(detailView, coord, frame, baseNum, FALSE, TRUE);
    }
}


/* Handler for when a mouse button is pressed on a sequence column header widget.
 * This signal only gets connected for protein matches, where there are 3 header
 * widgets, each showing the DNA nucleotides for a specific reading frame. */
static gboolean onButtonPressSeqColHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
  {
    case 2:
      {
	/* Middle button: select the nucleotide that was clicked on. */
	GtkWidget *detailView = GTK_WIDGET(data);
	selectClickedNucleotide(header, detailView, event->x, event->y);
	handled = TRUE;
      }
  }
  
  return handled;
}


static gboolean onButtonReleaseSeqColHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
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


static gboolean onMouseMoveSeqColHeader(GtkWidget *header, GdkEventMotion *event, gpointer data)
{
  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base */
      GtkWidget *detailView = GTK_WIDGET(data);
      selectClickedNucleotide(header, detailView, event->x, event->y);
    }
  
  return TRUE;
}


/***********************************************************
 *                     Callbacks                           *
 ***********************************************************/

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


/* If one exon view is visible and the other hidden, toggle their hidden states */
static void swapExonViewVisibility(GtkWidget *bigPicture)
{
  GtkWidget *exonView1 = bigPictureGetFwdExonView(bigPicture);
  GtkWidget *exonView2 = bigPictureGetRevExonView(bigPicture);
  
  if (widgetGetHidden(exonView1) != widgetGetHidden(exonView2))
    {
      widgetSetHidden(exonView1, !widgetGetHidden(exonView1));
      widgetSetHidden(exonView2, !widgetGetHidden(exonView2));
    }
}


void toggleStrand(GtkWidget *detailView)
{
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(detailViewGetMainWindow(detailView));
  GtkWidget *bigPicture = mainWindowProperties->bigPicture;

  /* Update the flag */
  mainWindowProperties->strandsToggled = !mainWindowProperties->strandsToggled;
  
  /* Invert the display range */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const IntRange const *fullRange = &mainWindowProperties->fullDisplayRange;
  const int newStart = fullRange->max - displayRange->max + fullRange->min;
  setDetailViewStartIdx(detailView, newStart, mainWindowProperties->seqType);

  /* Invert the currently-selected index, if there is one. We want to select
   * the same index but counting from the other end. */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  if (properties->selectedBaseIdx != UNSET_INT)
    {
      const int newIdx = fullRange->max - properties->selectedBaseIdx + fullRange->min;
      const int newBaseNum = mainWindowProperties->numReadingFrames - properties->selectedBaseNum + 1;
      
      detailViewSetSelectedBaseIdx(detailView, newIdx, properties->selectedFrame, newBaseNum, FALSE, TRUE);
    }
  
  /* If one grid/tree is hidden and the other visible, toggle which is hidden */
  swapTreeVisibility(detailView);
  swapGridVisibility(bigPicture);
  swapExonViewVisibility(bigPicture);
  
  /* Toggle the order of the trees and grids. */
  refreshTreeOrder(detailView);
  refreshGridOrder(bigPicture);

  /* Redraw the tree and detail-view headers */
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders);
  gtk_widget_queue_draw(detailView);
  
  /* Redraw the grids and grid headers */
  refreshBigPictureDisplayRange(bigPicture, TRUE);
}


/* Go to a user-specified coord on the reference sequence (in terms of the DNA
 * or peptide sequence as specified by coordSeqType). Selects the base index
 * and centres on it. */
void goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType)
{
  static gchar defaultInput[32] = "";
  
  /* Pop up a dialog to request a coord from the user */
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Go to position: ", 
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

	  /* Convert the input coord to display coords. */
	  const int activeFrame = detailViewGetSelectedFrame(detailView);
	  int baseNum;
	  
	  const int displayIdx = convertDnaIdxToDisplayIdx(requestedCoord, 
							   detailViewGetSeqType(detailView), 
							   activeFrame,
							   detailViewGetNumReadingFrames(detailView), 
							   detailViewGetStrandsToggled(detailView), 
							   detailViewGetRefSeqRange(detailView),
							   &baseNum);
	  
	  /* Select the base index. */
	  detailViewSetSelectedBaseIdx(detailView, displayIdx, activeFrame, baseNum, TRUE, FALSE);
	}
    }

  gtk_widget_destroy(dialog);
}


/* Sort the match entries by..... */

static void detailViewSortByName(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByName);
}

static void detailViewSortByScore(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByScore);
}

static void detailViewSortByPos(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByPos);
}

static void detailViewSortById(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortById);
}

static void detailViewSortByGroupOrder(GtkWidget *detailView)
{
  callFuncOnAllDetailViewTrees(detailView, treeSortByGroupOrder);
}

void detailViewSortByType(GtkWidget *detailView, const SortByType sortByType)
{
  switch (sortByType)
    {
      case SORTBYNAME:
	detailViewSortByName(detailView);
	break;
    
      case SORTBYSCORE:
	detailViewSortByScore(detailView);
	break;

      case SORTBYID:
	detailViewSortById(detailView);
	break;

      case SORTBYPOS:
	detailViewSortByPos(detailView);
	break;

      case SORTBYGROUPORDER:
	detailViewSortByGroupOrder(detailView);
	break;
	
      default:
	break;
    };
}


/* Find the next MSP (out the MSPs in this tree) whose start/end is the next closest to the
 * current start position of the display, searching only in the direction specified by the
 * search criteria passed in the user data. Updates the searchData with the offset found. */
static gboolean findNextMatchInTree(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  MatchSearchData *searchData = (MatchSearchData*)data;

  /* Loop through all MSPs in this tree row. */
  GList* mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      /* Check its in the sequence list (if given), and is not an intron */
      if (!mspIsIntron(msp) && (!searchData->seqNameList || findStringInList(searchData->seqNameList, msp->sname)))
	{
	  /* Get the msp start/end in terms of the display coords */
	  const int coord1 = convertDnaIdxToDisplayIdx(msp->qstart, searchData->seqType, searchData->frame, searchData->numFrames, searchData->rightToLeft, searchData->refSeqRange, NULL);
	  const int coord2 = convertDnaIdxToDisplayIdx(msp->qend, searchData->seqType, searchData->frame, searchData->numFrames, searchData->rightToLeft, searchData->refSeqRange, NULL);
	  
	  /* Get the offset of the coords from the given display range and find the smallest
	   * ignorning zero and negative offsets: negative means its the wrong direction) */
	  const int offset1 = (coord1 - searchData->startCoord) * searchData->searchDirection;
	  const int offset2 = (coord2 - searchData->startCoord) * searchData->searchDirection;
	  
	  int currentBest = UNSET_INT;
	  gboolean isStart = TRUE;
	  
	  if (offset1 > 0 && (offset2 <= 0 || offset2 >= offset1))
	    {
	      currentBest = offset1;
	      isStart = TRUE;
	    }
	  else
	    {
	      currentBest = offset2;
	      isStart = FALSE;
	    }

	  if (currentBest > 0)
	    {
	      /* Use the new offset if it is smaller than the previous. If it is the same, give
	       * preference to the one with the lowest frame number (i.e. so that we select the
	       * topmost frame by preference, in the hope that this will be least confusing). */
	      
	      if (searchData->offset == UNSET_INT  ||
		  currentBest < searchData->offset ||
		  (currentBest == searchData->offset && searchData->frame < searchData->foundFrame))
		{
		  searchData->offset = currentBest;
		  searchData->isStart = isStart;
		  searchData->foundFrame = searchData->frame;
		}
	    }
	}
    }
  
  return FALSE;
}


/* Find and go to the next match (either left or right depending on the search flag).
 * If a list of sequence names is given, only look at matches in those sequences.
 * startCoord determines where to start searching from. Sets the found match coord as
 * the currently-selected base index */
static void goToNextMatch(GtkWidget *detailView, const int startCoord, const gboolean searchRight, GList *seqNameList)
{
  const gboolean rightToLeft = detailViewGetStrandsToggled(detailView);
  const int numFrames = detailViewGetNumReadingFrames(detailView);
  const int searchDirection = searchRight ? 1 : -1;
  
  MatchSearchData searchData = {startCoord, 
				searchRight, 
				searchDirection, 
				rightToLeft, 
				detailViewGetSeqType(detailView),
				numFrames,
				detailViewGetRefSeqRange(detailView),
				seqNameList,
				UNSET_INT,
				UNSET_INT,
				TRUE,
				UNSET_INT};

  /* Loop through the MSPs in all visible trees */
  searchData.frame = 1;

  for (  ; searchData.frame <= searchData.numFrames; ++searchData.frame)
    {
      GtkWidget *treeContainer = detailViewGetTreeContainer(detailView, FORWARD_STRAND, searchData.frame);
      GtkWidget *tree = treeContainerGetTree(GTK_CONTAINER(treeContainer));
      
      if (GTK_WIDGET_VISIBLE(tree) && gtk_widget_get_parent(treeContainer)) /* ignore if not currently included in view */
	{
	  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
	  gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
	}

      treeContainer = detailViewGetTreeContainer(detailView, REVERSE_STRAND, searchData.frame);
      tree = treeContainerGetTree(GTK_CONTAINER(treeContainer));
      
      if (GTK_WIDGET_VISIBLE(tree) && gtk_widget_get_parent(treeContainer)) /* ignore if not currently included in view */
	{
	  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
	  gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
	}
    }
  
  if (searchData.offset != UNSET_INT)
    {
      /* Offset the start coord by the found amount. */
      const int newCoord = searchData.startCoord + (searchData.offset * searchDirection);
      
      /* Select the first base in the frame if moving to the start of the MSP 
       * and the last base if moving to the end */
      const int selectedFrame = searchData.foundFrame;
      const int selectedBaseNum = searchData.isStart ? 1 : numFrames;
      
      detailViewSetSelectedBaseIdx(detailView, newCoord, selectedFrame, selectedBaseNum, TRUE, FALSE);
    }
}


/* Go to the previous match (optionally limited to matches in the given list)  */
void prevMatch(GtkWidget *detailView, GList *seqNameList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  
  if (startCoord == UNSET_INT || !valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
    }
  
  goToNextMatch(detailView, startCoord, FALSE, seqNameList);
}


/* Go to the next match (optionally limited to matches in the given list)  */
void nextMatch(GtkWidget *detailView, GList *seqNameList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  
  if (startCoord == UNSET_INT || !valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
    }
  
  goToNextMatch(detailView, startCoord, TRUE, seqNameList);
}


/* Go to the first match (optionally limited to matches in the given list)  */
void firstMatch(GtkWidget *detailView, GList *seqNameList)
{
  /* Jump to the nearest match to the whole display start */
  const IntRange const *fullRange = detailViewGetFullRange(detailView);
  goToNextMatch(detailView, fullRange->min, TRUE, seqNameList);
}


/* Go to the last match (optionally limited to matches in the given list)  */
void lastMatch(GtkWidget *detailView, GList *seqNameList)
{
  /* Jump to the nearest match to the whole display end */
  const IntRange const *fullRange = detailViewGetFullRange(detailView);
  goToNextMatch(detailView, fullRange->max, FALSE, seqNameList);
}


/* Callback called when the sort order has been changed in the drop-down box */
static void onSortOrderChanged(GtkComboBox *combo, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  GtkTreeIter iter;
  
  if (GTK_WIDGET_REALIZED(detailView) && gtk_combo_box_get_active_iter(combo, &iter))
    {
      GtkTreeModel *model = gtk_combo_box_get_model(combo);
      
      GValue val = {0};
      gtk_tree_model_get_value(model, &iter, SORT_TYPE_COL, &val);
      
      SortByType sortType = g_value_get_int(&val);
      
      switch (sortType)
	{
	  case SORTBYSCORE:
	    detailViewSortByScore(detailView);
	    break;
	    
	  case SORTBYID:
	    detailViewSortById(detailView);
	    break;
	    
	  case SORTBYNAME:
	    detailViewSortByName(detailView);
	    break;
	    
	  case SORTBYPOS:
	    detailViewSortByPos(detailView);
	    break;
	    
	  case SORTBYGROUPORDER:
	    detailViewSortByGroupOrder(detailView);
	    break;
	    
	  default:
	    break;
	};
    }
}


static void GHelp(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  displayHelp(detailViewGetMainWindow(detailView));
}

static void GSettings(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  showSettingsDialog(detailViewGetMainWindow(detailView));
}

static void GFind(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  showFindDialog(detailViewGetMainWindow(detailView));
}

//static void GGroup(GtkButton *button, gpointer data)
//{
//  GtkWidget *detailView = GTK_WIDGET(data);
//  showGroupsDialog(detailViewGetMainWindow(detailView), TRUE);
//}

//static void GCopy(GtkButton *button, gpointer data)
//{
//  GtkWidget *detailView = GTK_WIDGET(data);
//  copySelectionToClipboard(detailViewGetMainWindow(detailView));
//}
//
//static void GPaste(GtkButton *button, gpointer data)
//{
//  GtkWidget *detailView = GTK_WIDGET(data);
//  requestDefaultClipboardText(findSeqsFromClipboard, detailViewGetMainWindow(detailView));
//}

static void GGoto(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* We currently only accept input in terms of DNA coords on the ref seq */
  const BlxSeqType seqType = BLXSEQ_DNA;
  
  goToDetailViewCoord(detailView, seqType);
}

static void GprevMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqNameList = mainWindowGetSelectedSeqs(detailViewGetMainWindow(detailView));
  prevMatch(detailView, seqNameList);
}

static void GnextMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqNameList = mainWindowGetSelectedSeqs(detailViewGetMainWindow(detailView));
  nextMatch(detailView, seqNameList);
}

static void GfirstMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqNameList = mainWindowGetSelectedSeqs(detailViewGetMainWindow(detailView));
  firstMatch(detailView, seqNameList);
}

static void GlastMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqNameList = mainWindowGetSelectedSeqs(detailViewGetMainWindow(detailView));
  lastMatch(detailView, seqNameList);
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
  GtkWidget *seqHeader = createSeqColHeader(detailView, seqType, numReadingFrames);
  GtkCallback seqCallback = (seqType == BLXSEQ_PEPTIDE) ? refreshSeqColHeader : NULL;
  
  /* The start and end columns have callbacks to switch the start/end text when display is toggled */
  GtkCallback startCallback = refreshStartColHeader;
  GtkCallback endCallback = refreshEndColHeader;
  
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
static GtkWidget* createSeqColHeader(GtkWidget *detailView,
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
	  
	  seqColHeaderSetFrameNumber(line, frame + 1);
	  
	  gtk_drag_dest_set(line, GTK_DEST_DEFAULT_MOTION, NULL, 0, GDK_ACTION_PRIVATE);
	  
	  gtk_widget_add_events(line, GDK_ENTER_NOTIFY_MASK);
	  gtk_widget_add_events(line, GDK_LEAVE_NOTIFY_MASK);
	  g_signal_connect(G_OBJECT(line), "button-press-event", G_CALLBACK(onButtonPressSeqColHeader), detailView);
	  g_signal_connect(G_OBJECT(line), "button-release-event", G_CALLBACK(onButtonReleaseSeqColHeader), detailView);
	  g_signal_connect(G_OBJECT(line), "motion-notify-event", G_CALLBACK(onMouseMoveSeqColHeader), detailView);
	}
    }
  
  return header;
}


/* Create the horizontal scrollbar for custom scrolling. This takes ownership of
 * the GtkAdjustment. */
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
  gtk_toolbar_set_icon_size(*toolbar, GTK_ICON_SIZE_MENU);
  gtk_toolbar_set_style(*toolbar, GTK_TOOLBAR_ICONS);

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


/* Add an option for the sorting drop-down box */
static GtkTreeIter* addSortBoxItem(GtkTreeStore *store, 
				  GtkTreeIter *parent, 
				  SortByType sortType, 
				  const char *sortName,
				  SortByType initSortType,
				  GtkComboBox *combo)
{
  GtkTreeIter iter;
  gtk_tree_store_append(store, &iter, parent);

  gtk_tree_store_set(store, &iter, SORT_TYPE_COL, sortType, SORT_TEXT_COL, sortName, -1);

  if (sortType == initSortType)
    {
      gtk_combo_box_set_active_iter(combo, &iter);
    }
  
  return NULL;
}


/* Create the combo box used for selecting sort criteria */
static void createSortBox(GtkToolbar *toolbar, GtkWidget *detailView, const SortByType initSortType)
{
  /* Add a label, to make it obvious what the combo box is for */
  GtkWidget *label = gtk_label_new(" <i>Sort by:</i>");
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  addToolbarWidget(toolbar, label);

  /* Create the data for the drop-down box. Use a tree so that we can sort by
   * multiple criteria. */
  GtkTreeStore *store = gtk_tree_store_new(N_SORT_COLUMNS, G_TYPE_INT, G_TYPE_STRING);
  GtkComboBox *combo = GTK_COMBO_BOX(gtk_combo_box_new_with_model(GTK_TREE_MODEL(store)));
  g_object_unref(store);
  addToolbarWidget(toolbar, GTK_WIDGET(combo));

  /* Create a cell renderer to display the sort text. */
  GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
  gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, FALSE);
  gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", SORT_TEXT_COL, NULL);

  /* Add the options */
  GtkTreeIter *iter = NULL;
  iter = addSortBoxItem(store, iter, SORTBYNAME, SORT_BY_NAME_STRING, initSortType, combo);
  iter = addSortBoxItem(store, iter, SORTBYSCORE, SORT_BY_SCORE_STRING, initSortType, combo);
  iter = addSortBoxItem(store, iter, SORTBYID, SORT_BY_ID_STRING, initSortType, combo);
  iter = addSortBoxItem(store, iter, SORTBYPOS, SORT_BY_POS_STRING, initSortType, combo);
  iter = addSortBoxItem(store, iter, SORTBYGROUPORDER, SORT_BY_GROUP_ORDER_STRING, initSortType, combo);

  g_signal_connect(G_OBJECT(combo), "changed", G_CALLBACK(onSortOrderChanged), detailView);
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


/* Creates a single button for the detail-view toolbar. */
static void makeToolbarButton(GtkToolbar *toolbar,
			      char *label,
			      char *stockId,
			      char *tooltip,
			      GtkSignalFunc callback_func,
			      gpointer data)
{
  GtkStockItem stockItem;
  GtkToolItem *tool_button = NULL;
  
  if (stockId && gtk_stock_lookup(stockId, &stockItem))
    {
      tool_button = gtk_tool_button_new_from_stock(stockId);
      gtk_tool_button_set_label(GTK_TOOL_BUTTON(tool_button), label);
    }
  else
    {
      tool_button = gtk_tool_button_new(NULL, label);
    }
  
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


static void insertToolbarSeparator(GtkToolbar *toolbar)
{
  GtkToolItem *separator = gtk_separator_tool_item_new();
  gtk_separator_tool_item_set_draw(GTK_SEPARATOR_TOOL_ITEM(separator), TRUE);
  gtk_toolbar_insert(toolbar, separator, -1);
}


/* Create the detail view toolbar */
static GtkWidget* createDetailViewButtonBar(GtkWidget *detailView, 
					    BlxBlastMode mode,
					    const SortByType sortByType,
					    GtkWidget **feedbackBox)
{
  GtkToolbar *toolbar = NULL;
  GtkWidget *toolbarContainer = createEmptyButtonBar(detailView, &toolbar);
  
  /* Help */
  makeToolbarButton(toolbar, "Help", GTK_STOCK_HELP,	    "Help (Ctrl-H)",	 (GtkSignalFunc)GHelp,		  detailView);

  /* Combo box for sorting */
  createSortBox(toolbar, detailView, sortByType);

  /* Settings button */
  makeToolbarButton(toolbar, "Settings", GTK_STOCK_PREFERENCES,  "Settings (Ctrl-S)",		 (GtkSignalFunc)GSettings,	  detailView);

  /* Zoom buttons */
  makeToolbarButton(toolbar, "Zoom in",		GTK_STOCK_ZOOM_IN,  "Zoom in (=)",  (GtkSignalFunc)onZoomInDetailView, detailView);
  makeToolbarButton(toolbar, "Zoom out",	GTK_STOCK_ZOOM_OUT, "Zoom out (-)", (GtkSignalFunc)onZoomOutDetailView, detailView);
  
  /* Navigation buttons */
  makeToolbarButton(toolbar, "Go to",		GTK_STOCK_JUMP_TO,    "Go to position (p)",		  (GtkSignalFunc)GGoto,		  detailView);

  makeToolbarButton(toolbar, "First match",	GTK_STOCK_GOTO_FIRST, "First match (Ctrl-Home)",	  (GtkSignalFunc)GfirstMatch,	  detailView);
  makeToolbarButton(toolbar, "Previous match",	GTK_STOCK_GO_BACK,    "Previous match (Ctrl-left)",	  (GtkSignalFunc)GprevMatch,	  detailView);
  makeToolbarButton(toolbar, "Next match",	GTK_STOCK_GO_FORWARD, "Next match (Ctrl-right)",	  (GtkSignalFunc)GnextMatch,	  detailView);
  makeToolbarButton(toolbar, "Last match",	GTK_STOCK_GOTO_LAST,  "Last match (Ctrl-End)",		  (GtkSignalFunc)GlastMatch,	  detailView);
  insertToolbarSeparator(toolbar);
  
  makeToolbarButton(toolbar, "<<", NULL,	"Scroll back lots",	  (GtkSignalFunc)GscrollLeftBig,	  detailView);
  makeToolbarButton(toolbar, "<",  NULL,	"Scroll back one (,)",	  (GtkSignalFunc)GscrollLeft1,	  detailView);
  makeToolbarButton(toolbar, ">",  NULL,	"Scroll forward one (.)", (GtkSignalFunc)GscrollRight1,  detailView);
  makeToolbarButton(toolbar, ">>", NULL,	"Scroll forward lots",	  (GtkSignalFunc)GscrollRightBig,  detailView);
  
  /* Find/copy/paste */
  makeToolbarButton(toolbar, "Find", GTK_STOCK_FIND,    "Find sequences (f, Ctrl-F)",		  (GtkSignalFunc)GFind,  detailView);
//  makeToolbarButton(toolbar, "Groups", NULL,  "Group sequences (g, Ctrl-G)",		  (GtkSignalFunc)GGroup,  detailView);

  /* Strand toggle button */
  if (mode == BLXMODE_BLASTX || mode == BLXMODE_TBLASTX || mode == BLXMODE_BLASTN)
    {
      makeToolbarButton(toolbar, "Toggle strand", GTK_STOCK_REFRESH, "Toggle strand (t)", (GtkSignalFunc)GToggleStrand, detailView);
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
  const int frame1 = 1, frame2 = 2, frame3 = 3;
  
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

  fixed_font_list = g_list_append(fixed_font_list, "andale");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida sans typewriter");
  fixed_font_list = g_list_append(fixed_font_list, "Bitstream vera");
  fixed_font_list = g_list_append(fixed_font_list, "monaco");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida console");
  fixed_font_list = g_list_append(fixed_font_list, "Courier");
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
			    const char const *refSeqName,
			    const int startCoord,
			    const gboolean sortInverted,
			    const SortByType sortByType)
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
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, mode, sortByType, &feedbackBox);

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
			     numReadingFrames,
			     startCoord,
			     sortInverted);
  
  return detailView;
}

