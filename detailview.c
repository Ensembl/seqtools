/*
 *  detailview.c
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/exonview.h>
#include <SeqTools/utilities.h>
#include <gtk/gtk.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define DETAIL_VIEW_TOOLBAR_NAME	"DetailViewToolbarName"
#define DETAIL_VIEW_WIDGET_NAME		"DetailViewWidget"
#define SORT_BY_NAME_STRING		"Name"
#define SORT_BY_SCORE_STRING		"Score"
#define SORT_BY_ID_STRING		"Identity"
#define SORT_BY_POS_STRING		"Position"
#define SORT_BY_GROUP_ORDER_STRING	"Group"
#define FONT_INCREMENT_SIZE		1
#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20
#define NO_SUBJECT_SELECTED_TEXT	"<no subject selected>"
#define MULTIPLE_SUBJECTS_SELECTED_TEXT	"<multiple subjects selected>"
#define DEFAULT_SNP_CONNECTOR_HEIGHT	0
#define DEFAULT_NUM_UNALIGNED_BASES     5     /* the default number of additional bases to show if displaying unaligned parts of the match sequence */
#define POLYA_SIG_BASES_UPSTREAM        50    /* the number of bases upstream from a polyA tail to search for polyA signals */

/* Define the columns' default widths and titles. */
#define BLXCOL_INT_COLUMN_WIDTH		40    /* default width for ordinary integer columns */
#define BLXCOL_HIDDEN_COLUMN_WIDTH      0     /* default width for columns that are initially hidden */
#define BLXCOL_SEQNAME_WIDTH            120   /* default width for the name column */
#define BLXCOL_START_WIDTH              50    /* default width for the start coord column */
#define BLXCOL_END_WIDTH                80    /* default width for end coord column (bigger because it also spans the scrollbar) */
#define BLXCOL_SEQUENCE_WIDTH           40    /* default width for sequence column */


typedef enum {SORT_TYPE_COL, SORT_TEXT_COL, N_SORT_COLUMNS} SortColumns;


typedef struct 
  {
    const int startDnaIdx;		/* the DNA coord to start searching from */
    const gboolean searchRight;		/* search towards the right or left */
    const int searchDirection;		/* multiplier to add/subtract depending on whether searching right/left */
    const gboolean displayRev;		/* true if the display is reversed */
    const BlxSeqType seqType;		/* whether viewing DNA or peptide seqs */
    const int numFrames;		/* number of reading frames */
    const IntRange const *refSeqRange;	/* the full range of the reference sequence */
    GList *seqList;			/* only search matches in these sequences */

    int offset;				/* the offset of the found MSP */
    int foundFrame;			/* which ref seq frame the MSP we chose is in */
    int foundBase;			/* the base number of the DNA coord we chose within the foundFrame */
  } MatchSearchData;


/* Local function declarations */
static BlxViewContext*	      detailViewGetContext(GtkWidget *detailView);
static GtkWidget*	      detailViewGetFirstTree(GtkWidget *detailView);
static GtkWidget*	      detailViewGetBigPicture(GtkWidget *detailView);
static GtkWidget*	      detailViewGetHeader(GtkWidget *detailView);
static GtkWidget*	      detailViewGetFeedbackBox(GtkWidget *detailView);
static int		      detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView);
static int		      detailViewGetSnpConnectorHeight(GtkWidget *detailView);

static void		      snpTrackSetStrand(GtkWidget *snpTrack, const int strand);
static int		      snpTrackGetStrand(GtkWidget *snpTrack);
static void		      getVariationDisplayRange(const MSP *msp, const gboolean expand, const BlxSeqType seqType, const int numFrames, const gboolean displayRev, const int activeFrame, const IntRange const *refSeqRange, IntRange *displayRange, IntRange *expandedRange);

static void		      detailViewCacheFontSize(GtkWidget *detailView, int charWidth, int charHeight);
static GtkToolItem*	      addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget);
static gboolean		      widgetIsTree(GtkWidget *widget);
static gboolean		      widgetIsTreeContainer(GtkWidget *widget);
static void		      updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc);
static GtkWidget*	      createSeqColHeader(GtkWidget *detailView, const BlxSeqType seqType, const int numFrames);
static const char*	      findDetailViewFont(GtkWidget *detailView);
static void		      setDetailViewScrollPos(GtkAdjustment *adjustment, int value);
static const char*            spliceSiteGetBases(const BlxSpliceSite *spliceSite, const gboolean donor, const gboolean reverse);

/***********************************************************
 *		       Utility functions                   *
 ***********************************************************/

void detailViewRedrawAll(GtkWidget *detailView)
{
  /* Redraw all the trees */
  callFuncOnAllDetailViewTrees(detailView, widgetClearCachedDrawable, NULL);
  gtk_widget_queue_draw(detailView);
}


/* Return the width of the column with the given column id */
int detailViewGetColumnWidth(GtkWidget *detailView, const BlxColumnId columnId)
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
#if GLIB_MAJOR_VERSION >= 1 && GLIB_MINOR_VERSION >= 4
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
      g_debug("Using fixed-width font '%s'\n", result);
    }
  else
    {
      g_critical("Could not find a fixed-width font. Alignments may not be displayed correctly.\n");
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
  int colWidth = detailViewGetColumnWidth(detailView, BLXCOL_SEQUENCE);
  
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
static void addTreesToDetailViewPane(GtkPaned *panedWin, 
				     GList *treeList, 
				     const gboolean first)
{
  int numTrees = g_list_length(treeList);
  
  if (numTrees == 1)
    {
      /* Two panes, one tree. Use 'first' flags to decide which pane to put it in */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      
      if (first)
	{
	  gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);
	}
      else
	{
	  gtk_paned_pack2(panedWin, tree1, TRUE, TRUE);
	}
    }
  else if (numTrees == 2)
    {
      /* Two panes, two trees. Easy - put one in each. */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      GtkWidget *tree2 = GTK_WIDGET(treeList->next->data);
      
      gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);
      gtk_paned_pack2(panedWin, tree2, TRUE, TRUE);
    }
  else if (numTrees > 2)
    {
      /* Two panes, three trees. Put one tree in pane 1. There should be a
       * nested widget in pane 2 that we can put the remaining trees in. */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);

      GtkWidget *nestedPanedWidget = findNestedPanedWidget(GTK_CONTAINER(panedWin));
      
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
	  addTreesToDetailViewPane(GTK_PANED(nestedPanedWidget), remainingTrees, FALSE);
	}
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
 * correct order according to the displayRev flag. It should be called every
 * time the strands are toggled. It assumes the trees are already in the 
 * detailView container, and that the properties have been set for all 3 widgets. */
static void refreshTreeOrder(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  BlxViewContext *bc = detailViewGetContext(detailView);
  gboolean toggled = detailViewGetDisplayRev(detailView);

  /* Get the direct parent of the trees */
  GtkWidget *detailViewContainer = getNamedChildWidget(detailView, DETAIL_VIEW_WIDGET_NAME);

  gtk_container_foreach(GTK_CONTAINER(detailView), removeAllTreesFromContainer, detailViewContainer);

  /* Extract the named detail-view widget, which is the direct parent of the trees. This 
   * should be a paned window. */
  GtkWidget *panedWin = getNamedChildWidget(detailView, DETAIL_VIEW_WIDGET_NAME);

  if (!GTK_IS_PANED(panedWin))
    {
      g_error("Unexpected detail view type: expected a paned widget. Could not add child trees.\n");
    }

  if (bc->seqType == BLXSEQ_DNA)
    {
      /* Add both trees - the order they're added will depend on whether the display is toggled or not. */
      addTreesToDetailViewPane(GTK_PANED(panedWin), properties->fwdStrandTrees, !toggled);
      addTreesToDetailViewPane(GTK_PANED(panedWin), properties->revStrandTrees, toggled);
    }
  else if (bc->seqType == BLXSEQ_PEPTIDE)
    {
      /* Only add one set of trees - the reverse strand if toggled, the forward strand if not. */
      GList *treeList = toggled ? properties->revStrandTrees : properties->fwdStrandTrees;
      addTreesToDetailViewPane(GTK_PANED(panedWin), treeList, TRUE);
    }
  
  /* Must show all child widgets because some of them may not have been in this parent before.
   * (Just calling gtk_widget_show on the individual trees doesn't seem to work.)
   * However, we then need to re-hide any that may have been previously hidden by the user. */
  gtk_widget_show_all(detailView);
  gtk_container_foreach(GTK_CONTAINER(detailView), hideUserHiddenWidget, NULL);
}


/* Refresh the header label for the start-coord column. It shows 'Start' for normal
 * orientation or 'End' if the display is reversed. */
static void refreshStartColHeader(GtkWidget *header, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  if (detailViewGetDisplayRev(detailView))
    {
      gtk_label_set_text(GTK_LABEL(header), "End");
    }
  else
    {
      gtk_label_set_text(GTK_LABEL(header), "Start");
    }
}


/* Refresh the header label for the end-coord column. It shows 'End' for normal
 * orientation or 'Start' if the display is reversed. */
static void refreshEndColHeader(GtkWidget *header, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  if (detailViewGetDisplayRev(detailView))
    {
      gtk_label_set_text(GTK_LABEL(header), "Start");
    }
  else
    {
      gtk_label_set_text(GTK_LABEL(header), "End");
    }
}


/* Refresh a header that contains fixed-width text that needs to be resized whenever a zoom
 * happens; updates the height of the widget and the font description, and clears its cached
 * drawable, if it has one. This just sets things up ready for the next expose event, where 
 * the actual drawing will take place. */
void refreshTextHeader(GtkWidget *header, gpointer data)
{
  const char *widgetName = gtk_widget_get_name(header);
  
  if (GTK_IS_LABEL(header))
    {
      /* For straightforward labels, just update the font size. */
      GtkWidget *detailView = GTK_WIDGET(data);
      
      PangoFontDescription *dvFont = detailViewGetFontDesc(detailView);
      int newSize = pango_font_description_get_size(dvFont);

      pango_font_description_set_size(header->style->font_desc, newSize);
      gtk_widget_modify_font(header, header->style->font_desc);
    }
  else if (!strcmp(widgetName, SNP_TRACK_HEADER_NAME) || !strcmp(widgetName, DNA_TRACK_HEADER_NAME))
    {
      /* Clear the cached drawable so that it gets recreated on the next expose */
      widgetClearCachedDrawable(header, NULL);

      /* Update the font and the widget height, in case the font-size has changed. */
      GtkWidget *detailView = GTK_WIDGET(data);
      gtk_widget_modify_font(header, detailViewGetFontDesc(detailView));

      const int charHeight = detailViewGetCharHeight(detailView);
      
      if (!strcmp(widgetName, SNP_TRACK_HEADER_NAME))
	{
	  /* SNP track */
	  BlxViewContext *bc = detailViewGetContext(detailView);
	
	  if (!bc->flags[BLXFLAG_SHOW_SNP_TRACK])
	    {
	      /* SNP track is hidden, so set the height to 0 */
	      gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, 0);
	      gtk_widget_set_size_request(header, -1, 0);
	    }
	  else
	    {
	      /* Set the height to the character height plus the connector line height */
	      const int height = charHeight + detailViewGetSnpConnectorHeight(detailView);
	      gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, height);
	      gtk_widget_set_size_request(header, -1, height);
	    }
	}
      else if (GTK_IS_LAYOUT(header))
	{
	  /* Normal text header. Set the height to the character height */
	  gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, charHeight);
	  gtk_widget_set_size_request(header, -1, charHeight);
	}

      gtk_widget_queue_draw(header);
    }
  
  /* If this is a container of header widgets, recurse over each child */
  if (!strcmp(widgetName, HEADER_CONTAINER_NAME) && GTK_IS_CONTAINER(header))
    {
      gtk_container_foreach(GTK_CONTAINER(header), refreshTextHeader, data);
    }
}


/* Refresh the headers for the detail view. */
void refreshDetailViewHeaders(GtkWidget *detailView)
{
  /* Loop through all widgets in the header and call refreshTextHeader. This
   * updates the font etc. if it is a type of widget that requires that. */
  GtkWidget *header = detailViewGetHeader(detailView);
  gtk_container_foreach(GTK_CONTAINER(header), refreshTextHeader, detailView);
  
  /* Loop through all columns and call the individual refresh callbacks for
   * each of their headers. This updates the specific information within the header. */
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
	  g_critical("Invalid column data for detail view header; header may not be refreshed correctly.\n");
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

      /* For the sequence column, don't set the size request, or we won't be able
       * to shrink the window. (The sequence col header will be resized dynamically.) */
      if (columnInfo->columnId != BLXCOL_SEQUENCE)
	{
	  /* For other columns, we can set the size request: they're small enough
	   * that we can live without the need to shrink below their sum widths. */
          if (columnInfo->width > 0)
            {
              widgetSetHidden(columnInfo->headerWidget, FALSE);
              gtk_widget_set_size_request(columnInfo->headerWidget, columnInfo->width, -1);
            }
          else
            {
              /* Zero width: we must hide the widget because it always seems to have
               * some width even if we set the size request to 0 */
              widgetSetHidden(columnInfo->headerWidget, TRUE);
            }
	}
    }
}


/* Update the font description for all relevant components of the detail view */
void updateDetailViewFontDesc(GtkWidget *detailView)
{
  PangoFontDescription *fontDesc = detailViewGetFontDesc(detailView);
  
  updateCellRendererFont(detailView, fontDesc);
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  callFuncOnAllDetailViewTrees(detailView, treeUpdateFontSize, NULL);
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
static char* getFeedbackText(GtkWidget *detailView, const BlxSequence *seq, const int numSeqsSelected)
{
  /* The info we need to find... */
  int qIdx = UNSET_INT; /* index into the ref sequence. Ref seq is always a DNA seq */
  int sIdx = UNSET_INT; /* index into the match sequence. Will be coords into the peptide sequence if showing peptide matches */
  int sLen = UNSET_INT; /* the length of the match sequence */

  BlxViewContext *bc = detailViewGetContext(detailView);
  const int selectedDnaBaseIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  
  /* See if a base is selected. */
  qIdx = selectedDnaBaseIdx;
  
  /* Find the sequence name text (or some default text to indicate that a sequence is not selected) */
  const char *noSeqText = numSeqsSelected > 0 ? MULTIPLE_SUBJECTS_SELECTED_TEXT : NO_SUBJECT_SELECTED_TEXT;
  
  if (seq)
    {
      if (g_list_length(seq->mspList) > 0)
	{
	  MSP *firstMsp = (MSP*)(seq->mspList->data);
	  
	  if (mspGetMatchSeq(firstMsp))
	    {
	      sLen = strlen(mspGetMatchSeq(firstMsp));
	    }

	  /* If a q index is selected, see if there is a valid base at that index 
	   * for any of the MSPs for the selected sequence. */
	  if (qIdx != UNSET_INT)
	    {
	      GList *mspListItem = seq->mspList;
              const int numUnalignedBases = detailViewGetNumUnalignedBases(detailView);
	      
	      for ( ; mspListItem; mspListItem = mspListItem->next)
		{
		  MSP *msp = (MSP*)(mspListItem->data);
		  
		  sIdx = gapCoord(msp, qIdx, bc->numFrames, mspGetRefFrame(msp, bc->seqType), bc->displayRev, TRUE, numUnalignedBases, bc->flags, bc->mspList);

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
      g_string_printf(resultString, "%d   ", qIdx);
    }
  
  if (seq)
    {
      const char *seqName = blxSequenceGetDisplayName(seq);
      
      if (seqName)
	{
	  g_string_append_printf(resultString, "%s", seqName);
	}
        
      if (seq->type == BLXSEQUENCE_VARIATION && seq->sequence && seq->sequence->str)
        {
          g_string_append_printf(resultString, " : %s", seq->sequence->str);
        }
    }
  else if (qIdx != UNSET_INT)
    {
      g_string_append_printf(resultString, "%s", noSeqText); 
    }
    
  if (sLen != UNSET_INT && (!seq || seq->type != BLXSEQUENCE_VARIATION))
    {
      g_string_append_printf(resultString, "(%d)", sLen);
    }

  if (sIdx != UNSET_INT && (!seq || seq->type != BLXSEQUENCE_VARIATION))
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

  BlxViewContext *bc = detailViewGetContext(detailView);
  const int numSeqsSelected = g_list_length(bc->selectedSeqs);
  
  if (numSeqsSelected == 1) /* currently we only properly handle single sequence selection */
    {
      const BlxSequence *seq = (const BlxSequence*)(bc->selectedSeqs->data);
      messageText = getFeedbackText(detailView, seq, numSeqsSelected);
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


/* This is called when the Squash matches option has been toggled. It switches to the condensed 
 * view of the trees (where multiple MSPs in the same sequence appear on the same row) if squash
 * is true, or reverts to the expanded version (where each MSP has its own row) if squash is false. */
void detailViewUpdateSquashMatches(GtkWidget *detailView, const gboolean squash)
{
  if (squash)
    {
      callFuncOnAllDetailViewTrees(detailView, treeSquashMatches, NULL);
    }
  else
    {
      callFuncOnAllDetailViewTrees(detailView, treeUnsquashMatches, NULL);
    }

  gtk_widget_queue_draw(detailView);
}


/* Set the value of the 'invert sort order' flag */
void detailViewUpdateSortInverted(GtkWidget *detailView, const gboolean invert)
{
  callFuncOnAllDetailViewTrees(detailView, resortTree, NULL);
}


/* Set the value of the 'Show SNP track' flag */
void detailViewUpdateShowSnpTrack(GtkWidget *detailView, const gboolean showSnpTrack)
{
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  detailViewRedrawAll(detailView);
  
  refreshDialog(BLXDIALOG_SETTINGS, detailViewGetBlxWindow(detailView));
}


/* Set the value of the 'Limit Unaligned Bases' flag */
void detailViewUpdateLimitUnalignedBases(GtkWidget *detailView, const gboolean limitUnalignedBases)
{
  /* Refilter and re-draw */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  gtk_widget_queue_draw(detailView);
}


/* Set the number of additional bases to show for the show-unaligned-sequence option */
void detailViewSetNumUnalignedBases(GtkWidget *detailView, const int numBases)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->numUnalignedBases = numBases;
  
  /* Refilter and re-draw */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  gtk_widget_queue_draw(detailView);
}


int getBaseIndexAtColCoords(const int x, const int y, const int charWidth, const IntRange const *displayRange)
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
  int coord = getBaseIndexAtColCoords(x, y, detailViewGetCharWidth(detailView), detailViewGetDisplayRange(detailView));
  
  if (coord != UNSET_INT)
    {	
      /* Get the base number of the clicked base within the active frame. The
       * base number is determined by which row in the header the mouse pointer 
       * is over: for frame 1, row 1 will give base 1; for frame 2, row 2 will
       * give base 1, etc. Start by getting the frame number for the clicked row: */
      int frame = detailViewGetActiveFrame(detailView);
      int row = seqColHeaderGetRow(header);
      const int numFrames = detailViewGetNumFrames(detailView);
      
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


/* If the user clicked on a SNP in the SNP header, select it. This updates the 
 * feedback box with the SNP name and coord, selects the SNP coord and recentres 
 * on it (if allowScroll is true), clearing any previous selections. If expandSnps
 * is true it means that the SNPs are drawn full-width (i.e. showing all the sequence
 * in the MSP) and we should take into account when clicking; otherwise, the SNP
 * is just shown as a marker at its actual coord, and we only need to check the
 * position of that coord to see if the user clicked there. snpTrack is the widget
 * that was clicked on and the optional 'colHeader' argument can be passed if coordinates
 * should be converted to the coordinate system of that column header (i.e. if the SNP track
 * widget is not the same width of the sequence column, the coords will need converting). */
void selectClickedSnp(GtkWidget *snpTrack,
		      GtkWidget *colHeader,
		      GtkWidget *detailView, 
		      const int xIn, 
		      const int yIn, 
		      const gboolean reCentre,
		      const gboolean expandSnps,
		      const int clickedBase)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  /* Convert x coord to sequence-column coords */
  int x = xIn, y = yIn;
  
  if (colHeader)
    {
      gtk_widget_translate_coordinates(snpTrack, colHeader, xIn, 0, &x, NULL);
    }
  
  int clickedDisplayIdx = getBaseIndexAtColCoords(x, y, properties->charWidth, &properties->displayRange);
  
  if (clickedDisplayIdx != UNSET_INT)
    {
      GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
      BlxViewContext *bc = blxWindowGetContext(blxWindow);
      const int activeFrame = detailViewGetActiveFrame(detailView);
      
      /* See if there are any SNPs at this displayIdx */
      GList *snpList = NULL;
      const MSP *msp = bc->mspList;
      
      for ( ; msp; msp = msp->next)
	{
	  if (mspIsVariation(msp))
	    {
	      /* Get the variation coords in terms of display coords, and also the 'expanded' range
               * of coords where the variation is displayed */
	      IntRange mspDisplayRange, mspExpandedRange;
	      getVariationDisplayRange(msp, expandSnps, bc->seqType, bc->numFrames, bc->displayRev, activeFrame, &bc->refSeqRange, &mspDisplayRange, &mspExpandedRange);

              const int clickedDnaIdx = convertDisplayIdxToDnaIdx(clickedDisplayIdx, bc->seqType, 1, clickedBase, bc->numFrames, bc->displayRev, &bc->refSeqRange);
              gboolean found = FALSE;
              int idxToSelect = UNSET_INT, baseToSelect = UNSET_INT;
              
              if (expandSnps && valueWithinRange(clickedDisplayIdx, &mspExpandedRange))
                {
                  /* We clicked inside this MSP on the variation track. Select the first coord in
                   * the MSP. */
                  found = TRUE;
                  idxToSelect = mspDisplayRange.min;
                  baseToSelect = 1;
                }
              else if (!expandSnps && valueWithinRange(clickedDnaIdx, &msp->qRange))
                {
                  /* We clicked on a coord in the ref seq header that is affected by this variation.
                   * Select the clicked coord. */
		  found = TRUE;
                  idxToSelect = clickedDisplayIdx;
                  baseToSelect = clickedBase;
                }
	      
	      if (found)
		{
		  snpList = g_list_prepend(snpList, msp->sSequence);
		  detailViewSetSelectedBaseIdx(detailView, idxToSelect, activeFrame, baseToSelect, TRUE, FALSE);

		  if (reCentre)
		    {
		      const int newStart = idxToSelect - (getRangeLength(&properties->displayRange) / 2);
		      setDetailViewStartIdx(detailView, newStart, bc->seqType);
		    }
		}
	    }
	}
      
      /* Clear any existing selections and select the new SNP(s) */
      blxWindowSetSelectedSeqList(blxWindow, snpList);
    }
}


/***********************************************************
 *			    Drawing			   *
 ***********************************************************/

/* Draw a vertical separator line at the edges of the given widget */
void drawColumnSeparatorLine(GtkWidget *widget, 
                             GdkDrawable *drawable, 
                             GdkGC *gc, 
                             const BlxViewContext *bc)
{
  const int lineWidth = 1; /* width of the separator lines */
  
  if (GTK_WIDGET_VISIBLE(widget) && widget->allocation.width > lineWidth)
    {
      GdkColor *color = getGdkColor(BLXCOLOR_TREE_GRID_LINES, bc->defaultColors, FALSE, bc->usePrintColors);
      gdk_gc_set_foreground(gc, color);
      
      /* We want dashed lines */
      gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      /* Make the dashes very short and closely packed (i.e. dash length of 1 pixel and gaps of 1 pixel) */
      int listLen = 1;
      gint8 dashList[listLen];
      dashList[0] = 1;
      gdk_gc_set_dashes(gc, 0, dashList, listLen);
      
      /* Draw vertical lines. Draw the rightmost edge of each widget (because we don't really want a
       * line at the leftmost edge of the first cell of the tree view.) */
      const int x = widget->allocation.x + widget->allocation.width - lineWidth;
      const int y1 = widget->allocation.y;
      const int y2 = y1 + widget->allocation.height;
      
      gdk_draw_line(drawable, gc, x, y1, x, y2);
    }
}


/* Populates the 'bases' string with the two bases from the ref seq adjacent to the start/end
 * of the given MSP (if it is an exon/match). 'Start' here means at the lower value coord end
 * of the MSP. 'revStrand' is true if this MSP is on the reverse strand of the ref seq */
static void mspGetAdjacentBases(const MSP const *msp, char *bases, const gboolean start, const gboolean revStrand, const BlxViewContext *bc)
{
  if (start)
    {
      bases[0] = getRefSeqBase(bc->refSeq, msp->qRange.min - 2, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[1] = getRefSeqBase(bc->refSeq, msp->qRange.min - 1, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[2] = '\0';
    }  
  else
    {
      bases[0] = getRefSeqBase(bc->refSeq, msp->qRange.max + 1, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[1] = getRefSeqBase(bc->refSeq, msp->qRange.max + 2, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[2] = '\0';
    }
}


/* This function returns the prev/next MSP in the given sequence in order of start position on
 * the ref seq. Returns the next MSP (with the next-highest start coord) if 'getNext' is true; 
 * otherwise returns the previous one. The return value will be NULL if there is no next/prev MSP.
 * The error will be set if there was any problem (e.g. if the given MSP does not exist the the sequence).
 * Only considers exons and matches. */
static const MSP* sequenceGetNextMsp(const MSP const *msp, 
                                     const BlxSequence *blxSeq, 
                                     const gboolean getNext, 
                                     GError **error)
{
  const MSP *result = NULL;

  GList *mspItem = blxSeq->mspList;
  const MSP *prevMsp = NULL;
  gboolean found = FALSE;
  
  /* Loop through all MSPs and look for the given one. */
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *curMsp = (const MSP*)(mspItem->data);
      
      if (curMsp == msp)
        {
          found = TRUE;
          
          /* If looking for the previous MSP, we already know it. */
          if (!getNext)
            {
              result = prevMsp;
              break;
            }
          else
            {
              /* Continue looping to find the next MSP. */
              continue;
            }
        }
      
      if (found && getNext)
        {
          /* We've found the input MSP and we're on the next one. See if it is an exon/match,
           * otherwise continue looping. */
          if (mspIsExon(curMsp) || mspIsBlastMatch(msp))
            {
              result = curMsp;
              break;
            }
        }
      
      /* Set the previous MSP (but only if it's an exon/match) */
      if (mspIsExon(curMsp) || mspIsBlastMatch(curMsp))
        {
          prevMsp = curMsp;
        }
    }
  
  if (!found && error)
    {
      g_set_error(error, BLX_ERROR, 1, "The given MSP '%s' was not found in the given sequence '%s'.\n", mspGetSName(msp), blxSeq->fullName);
    }
  
  return result;
}


/* Determine whether the bases at the other end of the intron for the start/end of 
 * the given MSP are canonical. 'canonicalStart' and 'canonicalEnd' contain the two bases to look
 * for at the start/end of the next/previous MSP to determine if the result is canonical. */
static gboolean mspIsSpliceSiteCanonicalOtherEnd(const MSP const *msp,
                                                 const BlxSequence *blxSeq,
                                                 const gboolean isMinCoord, 
                                                 const gboolean revStrand,
                                                 const gboolean donor,
                                                 const BlxViewContext *bc,
                                                 const BlxSpliceSite *spliceSite)
{
  gboolean result = FALSE;
 
  /* Get the previous MSP if we're looking at the min coord of the current MSP, or 
   * the next MSP if we're at the end */
  const MSP *nextMsp = sequenceGetNextMsp(msp, blxSeq, !isMinCoord, NULL);
  
  if (nextMsp)
    {
      /* Get the two bases at the end of the previous MSP / start of the next MSP */
      char bases[3];
      mspGetAdjacentBases(nextMsp, bases, !isMinCoord, revStrand, bc);
      
      /* If the original site was a donor site, look for an acceptor site (or vice versa) */
      const char *canonicalBases = spliceSiteGetBases(spliceSite, !donor, revStrand);

      if (stringsEqual(bases, canonicalBases, FALSE))
        {
          result = TRUE;
        }
    }
  
  return result;
}


/* Create a BlxSpliceSite and add it to the given list */
static void addBlxSpliceSite(GSList **spliceSites, char *donorSite, char *acceptorSite, const gboolean bothReqd)
{
  BlxSpliceSite *spliceSite = g_malloc(sizeof(BlxSpliceSite));

  if (strlen(donorSite) < 2 || strlen(acceptorSite) < 2)
    {
      g_critical("Error adding splice site info ['%s', '%s'].\n", donorSite, acceptorSite);
      return;
    }
  
  spliceSite->donorSite[0] = donorSite[0];
  spliceSite->donorSite[1] = donorSite[1];
  spliceSite->donorSite[2] = '\0';

  spliceSite->donorSiteRev[0] = donorSite[1];
  spliceSite->donorSiteRev[1] = donorSite[0];
  spliceSite->donorSiteRev[2] = '\0';

  spliceSite->acceptorSite[0] = acceptorSite[0];
  spliceSite->acceptorSite[1] = acceptorSite[1];
  spliceSite->acceptorSite[2] = '\0';
  
  spliceSite->acceptorSiteRev[0] = acceptorSite[1];
  spliceSite->acceptorSiteRev[1] = acceptorSite[0];
  spliceSite->acceptorSiteRev[2] = '\0';
  
  spliceSite->bothReqd = bothReqd;
  
  /* Add it to the list */
  *spliceSites = g_slist_append(*spliceSites, spliceSite);
}


/* Frees the memory used by a BlxSpliceSite */
static void destroyBlxSpliceSite(gpointer listItemData, gpointer data)
{
  BlxSpliceSite *spliceSite = (BlxSpliceSite*)listItemData;
  
  g_free(spliceSite);
}


/* Return the canonical bases at the donor/acceptor end of the given BlxSpliceSite. Returns them in
 * reverse order if 'reverse' is true, e.g. for a GC-AG intron, the dono is GC and acceptor is AG.
 * If reverse is true the donor returns CG and acceptor returns GA. The result is a pointer to the
 * string in the splice site, which is owned by the Detail View properties. */
static const char* spliceSiteGetBases(const BlxSpliceSite *spliceSite, const gboolean donor, const gboolean reverse)
{
  const char *result = NULL;
  
  if (donor)
    {
      result = reverse ? spliceSite->donorSiteRev : spliceSite->donorSite;
    }
  else
    {
      result = reverse ? spliceSite->acceptorSiteRev : spliceSite->acceptorSite;
    }
  
  return result;
}


/* Determines whether the intron splice sites for the given MSP are canonical or not.
 * We look for GC-AG or AT-AC introns. The latter must have both ends of the intron matching
 * to be canonical. GC is the donor site and AG is the acceptor site. */
static gboolean isMspSpliceSiteCanonical(const MSP const *msp, 
                                         const BlxSequence *blxSeq, 
                                         const gboolean isMinCoord, 
                                         const gboolean revStrand,
                                         const BlxViewContext *bc,
                                         GSList *spliceSites)
{
  gboolean result = FALSE;

  GSList *item = spliceSites;
  for ( ; item; item = item->next)
    {
      const BlxSpliceSite *spliceSite = (const BlxSpliceSite*)(item->data);

      /* Get the two adjacent bases at the start/end of the given MSP */
      char bases[3];
      mspGetAdjacentBases(msp, bases, isMinCoord, revStrand, bc);
      
      /* The min coord end of the MSP is the acceptor site of the intron and the max coord the 
       * donor site (or vice versa if the strand is reversed). */
      const gboolean donor = (isMinCoord == revStrand);
      const char *canonicalBases = spliceSiteGetBases(spliceSite, donor, revStrand);

      if (stringsEqual(bases, canonicalBases, FALSE))
        {
          if (spliceSite->bothReqd)
            {
              /* Must also check if the other end of the intron matches */
              result = mspIsSpliceSiteCanonicalOtherEnd(msp, blxSeq, isMinCoord, revStrand, donor, bc, spliceSite);
            }
          else
            {
              result = TRUE;
            }
        }
    }

  return result;
}


/* If the given MSP is an exon/match, get the 2 nucleotides from the reference sequence at the
 * splice sites of the adjacent introns. Inserts the results into the given hash table with
 * an enum indicating whether the nucleotides are canonical/non-canonical. Only considers
 * MSPs that are within the given ref seq range. */
static void mspGetSpliceSiteCoords(const MSP const *msp, 
                                   const BlxSequence *blxSeq, 
                                   const IntRange const *qRange, 
                                   const BlxViewContext *bc, 
                                   GSList *spliceSites,
                                   GHashTable *result)
{
  const gboolean revStrand = (mspGetRefStrand(msp) != BLXSTRAND_FORWARD);
  
  /* Ignore the termini (i.e. 5' end of first exon and 3' end of last exon) */
  const gboolean getMin = (msp->sSequence && msp->qRange.min != msp->sSequence->qRange.min);
  const gboolean getMax = (msp->sSequence && msp->qRange.max != msp->sSequence->qRange.max);
  
  /* See if the min coord is within the given range */
  if (getMin && valueWithinRange(msp->qRange.min, qRange) && msp->qRange.min >= bc->refSeqRange.min + 2)
    {
      /* Find out if the adjacent bases are canonical/non-canonical. */
      const gboolean canonical = isMspSpliceSiteCanonical(msp, blxSeq, TRUE, revStrand, bc, spliceSites);
      BlxColorId colorId = canonical ? BLXCOLOR_CANONICAL : BLXCOLOR_NON_CANONICAL;

      /* Insert the two coords into the hash table */
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.min - 2), GINT_TO_POINTER(colorId));
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.min - 1), GINT_TO_POINTER(colorId));
    }
  
  /* See if the max coord is within the given range */
  if (getMax && valueWithinRange(msp->qRange.max, qRange) && msp->qRange.max <= bc->refSeqRange.max - 2)
    {
      /* Find out if they are canonical/non-canonical */
      const gboolean canonical = isMspSpliceSiteCanonical(msp, blxSeq, FALSE, revStrand, bc, spliceSites);
      BlxColorId colorId = canonical ? BLXCOLOR_CANONICAL : BLXCOLOR_NON_CANONICAL;
      
      /* Insert the two coords into the hash table */
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.max + 1), GINT_TO_POINTER(colorId));
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.max + 2), GINT_TO_POINTER(colorId));
    }
}


/* This looks through all of the polyA signals and sees if any of them are in the given display
 * range (in nucleotide coords) and whether we want to display them. If the 'display polyA signals for
 * selected sequences only' option is enabled, then any selected sequence MSPs that have polyA tails
 * should be passed in the given GSList. Any coords within a relevant polyA signal range are added to
 * the given hash table with the BlxColorId that they should be drawn with. */
static void getPolyASignalBasesToHighlight(const BlxViewContext *bc, GSList *polyATailMsps, const BlxStrand qStrand, const IntRange const *qRange, GHashTable *result)
{
  /* We've more work to do if showing polyA signals (unless we're only showing them for selected
   * sequences and there were no relevant selected sequences) */
  if (bc->flags[BLXFLAG_SHOW_POLYA_SIG] && (!bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED] || g_slist_length(polyATailMsps) > 0))
    {
      /* Loop through all polyA signals */
      const MSP const *sigMsp = bc->mspList;
      BlxColorId colorId = BLXCOLOR_POLYA_TAIL;
      const int direction = (qStrand == BLXSTRAND_REVERSE ? -1 : 1);
      
      for ( ; sigMsp; sigMsp = sigMsp->next)
        {
          /* Only interested if it's a polyA signal with the correct strand and is within the display range. */
          if (sigMsp->type == BLXMSP_POLYA_SIGNAL && sigMsp->qStrand == qStrand && rangesOverlap(&sigMsp->qRange, qRange))
            {
              gboolean addSignal = FALSE;
              
              if (bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED])
                {
                  /* We're only interested in polyA signals that are 50 bases upstream of one of our polyA-tail MSPs */
                  GSList *item = polyATailMsps;
                  
                  for ( ; item && !addSignal; item = item->next)
                    {
                      const MSP const *tailMsp = (const MSP const*)(item->data);
                      const int qEnd = mspGetQEnd(tailMsp);
                      const int qStart = qEnd - direction * POLYA_SIG_BASES_UPSTREAM;
                      
                      IntRange upstreamRange;
                      intrangeSetValues(&upstreamRange, qStart, qEnd); /* sorts out which is min and which is max */
                                        
                      addSignal = rangesOverlap(&sigMsp->qRange, &upstreamRange);
                    }
                }
              else
                {
                  /* Add all signals that are in the display range */
                  addSignal = TRUE;
                }
            
              if (addSignal)
                {
                  /* Add each base in the polyA signal range to the hash table. This may overwrite
                   * splice site bases that we previously found, which is fine (something has to take priority). */
                  int i = sigMsp->qRange.min;
                  for ( ; i <= sigMsp->qRange.max; ++i)
                    {
                      g_hash_table_insert(result, GINT_TO_POINTER(i), GINT_TO_POINTER(colorId));
                    }
                }
            }
        }
    }
}


/* This function looks for special bases in the reference sequence header to highlight and stores their
 * coords in the returned hash table with the BlxColorId (converted to a gpointer with GINT_TO_POINTER)
 * of the fill color they should be drawn with. These special coords include splice sites and polyA signals. */
GHashTable* getRefSeqBasesToHighlight(GtkWidget *detailView, 
                                      const IntRange const *qRange, 
                                      const BlxSeqType seqType,
                                      const BlxStrand qStrand)
{
  GHashTable *result = g_hash_table_new(g_direct_hash, g_direct_equal);

  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  const BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* We only highlight nucleotides, so there is nothing to do if we're showing a peptide sequence,
   * or if the show-splice-sites or show-polyA-signals options are disabled. */
  if (seqType == BLXSEQ_PEPTIDE || 
      (!bc->flags[BLXFLAG_SHOW_SPLICE_SITES] && !bc->flags[BLXFLAG_SHOW_POLYA_SIG]))
    {
      return result;
    }
  
  /* Loop through the selected sequences. If showing polyA signals for selected sequences, we'll 
   * compile a list of MSPs we're interested in seeing polyA signals for (i.e. those with polyA tails) */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GList *seqItem = blxWindowGetSelectedSeqs(blxWindow);
  GSList *polyATailMsps = NULL;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      const BlxSequence *blxSeq = (const BlxSequence*)(seqItem->data);
      GList *mspItem = blxSeq->mspList;
     
      /* Loop through all MSPs for this sequence */ 
      for ( ; mspItem; mspItem = mspItem->next)
        {
          MSP *msp = (MSP*)(mspItem->data);
          
          /* Only look at matches/exons on the correct strand */
          if ((mspIsBlastMatch(msp) || msp->type == BLXMSP_EXON) && mspGetRefStrand(msp) == qStrand)
            {
              if (bc->flags[BLXFLAG_SHOW_SPLICE_SITES])
                {
                  mspGetSpliceSiteCoords(msp, blxSeq, qRange, bc, properties->spliceSites, result);
                }
              
              if (bc->flags[BLXFLAG_SHOW_POLYA_SIG] && bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED])
                {
                   if (mspHasPolyATail(msp, bc->mspList))
                     {
                       polyATailMsps = g_slist_append(polyATailMsps, msp);
                     }
                }
            }
        }
    }
  
  /* Now check the polyA signals and see if any of them are in range. */
  getPolyASignalBasesToHighlight(bc, polyATailMsps, qStrand, qRange, result);
  
  return result;
}


static void drawDnaTrack(GtkWidget *dnaTrack, GtkWidget *detailView, const BlxStrand strand, const int frame)
{
  GdkDrawable *drawable = createBlankPixmap(dnaTrack);
  
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Find the segment of the ref sequence to display (complemented if this tree is
   * displaying the reverse strand, and reversed if the display is toggled) */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  GError *error = NULL;

  gchar *segmentToDisplay = getSequenceSegment(bc,
					       bc->refSeq,
					       displayRange->min, 
					       displayRange->max, 
					       strand, 
					       bc->seqType,
					       frame, 
					       bc->displayRev,
					       bc->displayRev,
					       bc->displayRev,
					       FALSE,
					       &error);
  
  if (!segmentToDisplay)
    {
      g_assert(error);
      prefixError(error, "Could not draw DNA header. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
    
    
  GdkGC *gc = gdk_gc_new(drawable);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  const BlxStrand activeStrand = blxWindowGetActiveStrand(blxWindow);
  const gboolean showSnpTrack = TRUE; /* always highlight SNP positions in the DNA header */

  gtk_layout_set_size(GTK_LAYOUT(dnaTrack), dnaTrack->allocation.width, properties->charHeight);
  
  /* Loop forward/backward through the display range depending on which strand we're viewing.
   * We need to convert display coords to actual coords on the ref sequence */
  const int qIdx1 = convertDisplayIdxToDnaIdx(displayRange->min, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);	  /* 1st base in frame */
  const int qIdx2 = convertDisplayIdxToDnaIdx(displayRange->max, bc->seqType, frame, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange); /* last base in frame */
  
  IntRange qRange = {min(qIdx1, qIdx2), max(qIdx1, qIdx2)};
  
  /* Find out if there are any bases in the introns that need highlighting. */
  GHashTable *basesToHighlight = getRefSeqBasesToHighlight(detailView, &qRange, BLXSEQ_DNA, activeStrand);
  
  int incrementValue = bc->displayRev ? -bc->numFrames : bc->numFrames;
  int displayLen = qRange.max - qRange.min + 1;
  char displayText[displayLen + 1];
  int displayTextPos = 0;
  
  int qIdx = bc->displayRev ? qRange.max : qRange.min;
  int displayIdx = convertDnaIdxToDisplayIdx(qIdx, bc->seqType, activeFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int y = 0;
  
  while (qIdx >= qRange.min && qIdx <= qRange.max)
    {
      /* Get the character to display at this index */
      displayText[displayTextPos] = getRefSeqBase(bc->refSeq, qIdx, bc->displayRev, &bc->refSeqRange, BLXSEQ_DNA);
      
      /* Color the base depending on whether it is selected or affected by a SNP */
      const gboolean displayIdxSelected = (displayIdx == properties->selectedBaseIdx);
      const gboolean dnaIdxSelected = (qIdx == properties->selectedDnaBaseIdx);
      const int x = displayTextPos * properties->charWidth;
      const char base = displayText[displayTextPos];

      drawHeaderChar(bc, properties, qIdx, base, strand, UNSET_INT, BLXSEQ_DNA, displayIdxSelected, dnaIdxSelected, FALSE, showSnpTrack, TRUE, BLXCOLOR_BACKGROUND, drawable, gc, x, y, basesToHighlight);
      
      /* Increment indices */
      ++displayTextPos;
      ++displayIdx;
      qIdx += incrementValue;
    }
  
  /* Make sure the string is terminated properly */
  displayText[displayTextPos] = '\0';
  
  /* Draw the text */
  PangoLayout *layout = gtk_widget_create_pango_layout(detailView, displayText);
  pango_layout_set_font_description(layout, detailViewGetFontDesc(detailView));
  
  if (layout)
    {
      gtk_paint_layout(dnaTrack->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, 0, 0, layout);
      g_object_unref(layout);
    }
  
  drawColumnSeparatorLine(dnaTrack, drawable, gc, bc);
  
  g_free(segmentToDisplay);
  g_object_unref(gc);
}


/* Utility to find the display coord(s) that a variation lies on. If expand is true, 
 * it means that if the variation contains more than one alternative, and these will be 
 * displayed horizontally across the display, taking up more width than where 
 * the actual variation coords lie. expandedRange will take this into account. displayRange 
 * always returns the actual msp start/end in display coords. */
static void getVariationDisplayRange(const MSP *msp, 
                                     const gboolean expand,
                                     const BlxSeqType seqType, 
                                     const int numFrames,
                                     const gboolean displayRev, 
                                     const int activeFrame,
                                     const IntRange const *refSeqRange,
                                     IntRange *displayRange,
                                     IntRange *expandedRange)
{
  /* Convert the MSP coords to display coords */
  int base1, base2;
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, seqType, activeFrame, numFrames, displayRev, refSeqRange, &base1);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, seqType, activeFrame, numFrames, displayRev, refSeqRange, &base2);

  intrangeSetValues(displayRange, coord1, coord2);
  
  /* The calculated display index will be different depending on which reading frame is 
   * active. However, we always want to display the variation in the same position, so adjust 
   * the display range so that we have coords in terms of frame 1. Return the real variation display
   * coords too, though - we need this so that we can highlight the selected display index and
   * still have the correct frame and base highlighted. */
  int adjustedIdx1 = coord1;
  if (base1 > (numFrames - activeFrame + 1))
    {
      ++adjustedIdx1;
    }
  
  int adjustedIdx2 = coord2;
  if (base2 > (numFrames - activeFrame + 1))
    {
      ++adjustedIdx2;
    }
  
  intrangeSetValues(expandedRange, adjustedIdx1, adjustedIdx2);
  
  if (expand && mspGetMatchSeq(msp))
    {
      /* Expand the variation range so that we display its entire sequence. We'll position 
       * the variation so that the middle of its sequence lies at the centre coord of its ref seq range */
      const int numChars = strlen(mspGetMatchSeq(msp));
      
      int offset = (int)((double)numChars / 2.0);
      
      if (numChars % 2 == 0)
        {
          /* Shift to the right by one if we've got an even number of chars */
          --offset;
        }
      
      expandedRange->min = getRangeCentre(expandedRange) - offset;
      expandedRange->max = expandedRange->min + numChars - 1;
    }
}


/* Determine whether the given coord in the given frame/strand is affected by
 * a variation. Sets whether to draw the start/end/top/bottom boundaries for the coord
 * when outlining it by the extent of the found variation. Also sets drawBackground to true
 * if the background of the base should be highlighted in the variation color (i.e. if the base
 * is part of a selected variation). */
gboolean coordAffectedByVariation(const int dnaIdx,
                                  const BlxStrand strand, 
                                  const MSP *mspList, 
                                  BlxViewContext *bc,
                                  const MSP **mspOut, /* the variation we found */
                                  gboolean *drawStartBoundary, 
                                  gboolean *drawEndBoundary, 
                                  gboolean *drawTopBoundary, 
                                  gboolean *drawBottomBoundary,
                                  gboolean *drawBackground)
{
  gboolean result = FALSE;
  
  /* Loop through the MSPs and check any that are SNPs */
  const MSP *msp = mspList;
  
  for ( ; msp; msp = msp->next)
    {
      if (mspIsVariation(msp) && mspGetRefStrand(msp) == strand && valueWithinRange(dnaIdx, &msp->qRange))
        {
          result = TRUE;
          
          if (mspOut)
            *mspOut = msp;
          
          /* We draw the start boundary of this base if it's the first base in the variation's range 
           * and the end boundary if it's the last, unless it's a zero-length feature, i.e. an insertion
           * site; in this case we draw a line to the right of the coord in the direction of the
           * reference sequence (i.e. to give just a vertical line at the site, not an outline around the bases) */
          const gboolean isFirst = (dnaIdx == msp->qRange.min);
          const gboolean isLast = (dnaIdx == msp->qRange.max);
          const gboolean isZeroLen = mspIsZeroLenVariation(msp);

          if (drawStartBoundary)
            *drawStartBoundary = (isZeroLen != bc->displayRev ? isLast : isFirst);
            
          if (drawEndBoundary)
            *drawEndBoundary = (isZeroLen != bc->displayRev ? isFirst : isLast);

          if (drawTopBoundary)
            *drawTopBoundary = !isZeroLen;
            
          if (drawBottomBoundary)
            *drawBottomBoundary = !isZeroLen;

          /* Highlight the background if the variation is selected (unless it's a zero-length variation) */
          if (drawBackground)
            *drawBackground = !isZeroLen && blxContextIsSeqSelected(bc, msp->sSequence);
          
          break;
	}
    }
  
  return result;
}


/* Draw a rectangle with an outline and a fill color */
static void drawRectangle(GdkDrawable *drawable, 
			  GdkGC *gc, 
			  GdkColor *fillColor,
			  GdkColor *outlineColor,
			  const int x, 
			  const int y,
			  const int width, 
			  const int height,
                          const gboolean drawLeft,
                          const gboolean drawRight,
                          const gboolean drawTop,
                          const gboolean drawBottom)
{
  const int lineWidth = 1;
  
  if (fillColor)
    {
      gdk_gc_set_foreground(gc, fillColor);
      gdk_draw_rectangle(drawable, gc, TRUE, x, y, width, height);
    }

  if (outlineColor)
    {
      gdk_gc_set_foreground(gc, outlineColor);
      gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      int w = width - lineWidth;
      int h = height - lineWidth;
      
      if (drawLeft)
        gdk_draw_line(drawable, gc, x, y, x, y + h);
        
      if (drawRight)
        gdk_draw_line(drawable, gc, x + w, y, x + w, y + h);

      if (drawTop)
        gdk_draw_line(drawable, gc, x, y, x + w, y);

      if (drawBottom)
        gdk_draw_line(drawable, gc, x, y + h, x + w, y + h);

      //gdk_draw_rectangle(drawable, gc, FALSE, x, y, width - lineWidth, height - lineWidth);
    }
}


/* Returns true if the given coord (which must be in nucleotide coords) is within the
 * range of any MSP that is selected. Only considers MSPs that are on the given ref seq
 * frame and strand (if given; otherwise considers all MSPs). */
static gboolean isCoordInSelectedMspRange(const BlxViewContext *bc,
                                          const int dnaIdx, 
                                          const BlxStrand refSeqStrand, 
                                          const int refSeqFrame,
                                          const BlxSeqType seqType)
{
  gboolean inSelectedMspRange = FALSE;
  
  /* Loop through all the selected sequences */
  GList *seqItem = bc->selectedSeqs;
  
  for ( ; seqItem && !inSelectedMspRange; seqItem = seqItem->next)
    {
      /* Loop through all the MSPs in this sequence */
      const BlxSequence const *blxSeq = (const BlxSequence const *)(seqItem->data);
      GList *mspItem = blxSeq->mspList;
      
      for ( ; mspItem && !inSelectedMspRange; mspItem = mspItem->next)
        {
          const MSP const *msp = (const MSP const *)(mspItem->data);
          const int mspFrame = mspGetRefFrame(msp, seqType);
          
          if ((mspIsBlastMatch(msp) || mspIsExon(msp)) && 
              (refSeqStrand == BLXSTRAND_NONE || mspGetRefStrand(msp) == refSeqStrand) &&
              (refSeqFrame == UNSET_INT || mspFrame == refSeqFrame))
            {
              /* Passed-in coord is the nucleotide for base 1 in frame 1 of the current display index.
               * We need to compare against base1 for the frame for the individual MSP, though, so 
               * we need to convert it: just add 1 base if we're in frame 2 or 2 bases if we're in frame3. */
              int idx = dnaIdx;
              
              if (bc->displayRev)
                {
                  idx -= (mspFrame - 1);
                }
              else
                {
                  idx += (mspFrame - 1);
                }
              
              inSelectedMspRange = valueWithinRange(idx, &msp->qRange);
            }
        }
    }

  return inSelectedMspRange;
}

/* Draw a given nucleotide or peptide. Determines the color depending on various
 * parameters */
void drawHeaderChar(BlxViewContext *bc,
		    DetailViewProperties *properties,
		    const int dnaIdx,
		    const char baseChar,
		    const BlxStrand strand,
                    const int frame,
		    const BlxSeqType seqType,
		    const gboolean displayIdxSelected,
		    const gboolean dnaIdxSelected,
		    const gboolean showBackground,    /* whether to use default background color or leave blank */
		    const gboolean showSnps,	      /* whether to show SNPs */
		    const gboolean showCodons,	      /* whether to highlight DNA bases within the selected codon, for protein matches */
                    const BlxColorId defaultBgColor,  /* the default background color for the header */
		    GdkDrawable *drawable,
		    GdkGC *gc,
		    const int x,
		    const int y,
                    GHashTable *basesToHighlight)
{
  GdkColor *fillColor = NULL;
  GdkColor *outlineColor = NULL;

  /* If drawing an outline, these can be set to false to omit certain parts of the outline */
  gboolean drawLeft = TRUE, drawRight = TRUE, drawTop = TRUE, drawBottom = TRUE;

  /* Shade the background if the base is selected XOR if the base is within the range of a 
   * selected sequence. (If both conditions are true we don't shade, to give the effect of an
   * inverted selection color.) */
  gboolean inSelectedMspRange = isCoordInSelectedMspRange(bc, dnaIdx, strand, frame, seqType);
  const gboolean shadeBackground = (displayIdxSelected != inSelectedMspRange);
  
  /* Check if this coord already has a special color stored for it */
  gpointer hashValue = g_hash_table_lookup(basesToHighlight, GINT_TO_POINTER(dnaIdx));
  
  if (hashValue)
    {
      BlxColorId colorId = GPOINTER_TO_INT(hashValue);
      fillColor = getGdkColor(colorId, bc->defaultColors, shadeBackground, bc->usePrintColors);
    }
  
  if (seqType == BLXSEQ_DNA && showCodons && (dnaIdxSelected || displayIdxSelected || shadeBackground))
    {
      if (dnaIdxSelected || displayIdxSelected)
        {
          /* The coord is a nucleotide in the currently-selected codon. The color depends
           * on whether the actual nucleotide itself is selected, or just the codon that it 
           * belongs to. */
          fillColor = getGdkColor(BLXCOLOR_CODON, bc->defaultColors, dnaIdxSelected, bc->usePrintColors);
        }
      else
        {
          /* The coord is not selected but this coord is within the range of a selected MSP, so 
           * shade the background. */
          fillColor = getGdkColor(defaultBgColor, bc->defaultColors, shadeBackground, bc->usePrintColors);
        }
    }

  if (showSnps)
    {
      gboolean drawBackground = TRUE;
    
      if (coordAffectedByVariation(dnaIdx, strand, bc->mspList, bc, NULL, &drawLeft, &drawRight, &drawTop, &drawBottom, &drawBackground))
        {
          /* The coord is affected by a SNP. Outline it in the "selected" SNP color
           * (which is darker than the normal color) */
          outlineColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, TRUE, bc->usePrintColors);
          
          /* If the SNP is selected, also fill it with the SNP color (using the
           * "unselected" SNP color, which is lighter than the outline). */
          if (drawBackground)
            {
              fillColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, shadeBackground, bc->usePrintColors);
            }
        }
    }
  
  if (!fillColor)
    {
      if (seqType == BLXSEQ_PEPTIDE && baseChar == SEQUENCE_CHAR_MET)
	{
	  /* The coord is a MET codon */
	  fillColor = getGdkColor(BLXCOLOR_MET, bc->defaultColors, shadeBackground, bc->usePrintColors);
	}
      else if (seqType == BLXSEQ_PEPTIDE && baseChar == SEQUENCE_CHAR_STOP)
	{
	  /* The coord is a STOP codon */
	  fillColor = getGdkColor(BLXCOLOR_STOP, bc->defaultColors, shadeBackground, bc->usePrintColors);
	}
      else if (showBackground)
	{
	  /* Use the default background color for the reference sequence */
	  fillColor = getGdkColor(defaultBgColor, bc->defaultColors, shadeBackground, bc->usePrintColors);
	}
    }
  
  
  drawRectangle(drawable, gc, fillColor, outlineColor, x, y, properties->charWidth, properties->charHeight, drawLeft, drawRight, drawTop, drawBottom);
}


/* Function that does the drawing for the SNP track */
static void drawSnpTrack(GtkWidget *snpTrack, GtkWidget *detailView)
{
  /* Create the drawable for the widget (whether we're actually going to do any drawing or not) */
  GdkDrawable *drawable = createBlankPixmap(snpTrack);

  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (!bc->flags[BLXFLAG_SHOW_SNP_TRACK])
    {
      return;
    }
  
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  
  BlxStrand strand = (snpTrackGetStrand(snpTrack) == UNSET_INT)
    ? blxWindowGetActiveStrand(blxWindow)
    : (BlxStrand)snpTrackGetStrand(snpTrack);
  
  GdkGC *gc = gdk_gc_new(drawable);
  
  /* Find the left margin. It will be at the same x coord as the left edge of
   * the sequence column header. */
  int leftMargin = UNSET_INT;
  DetailViewColumnInfo *seqColInfo = detailViewGetColumnInfo(detailView, BLXCOL_SEQUENCE);
  gtk_widget_translate_coordinates(seqColInfo->headerWidget, snpTrack, 0, 0, &leftMargin, NULL);
  
  /* Loop through all the MSPs looking for SNPs in the current display range */
  const MSP* msp = blxWindowGetMspList(blxWindow);
  const int y = 0;
  
  for ( ; msp; msp = msp->next)
    {
      if (mspIsVariation(msp) && mspGetRefStrand(msp) == strand && mspGetMatchSeq(msp))
	{
	  /* Get the range of display coords where this SNP will appear */
          IntRange mspDisplayRange, mspExpandedRange;
	  getVariationDisplayRange(msp, TRUE, bc->seqType, bc->numFrames, bc->displayRev, activeFrame, &bc->refSeqRange, &mspDisplayRange, &mspExpandedRange);

	  /* See if the variation is in the current display range */
	  if (rangesOverlap(&mspExpandedRange, &properties->displayRange))
	    {
	      int x = leftMargin + ((mspExpandedRange.min - properties->displayRange.min) * properties->charWidth);
	      const int width = strlen(mspGetMatchSeq(msp)) * properties->charWidth;
	      const gboolean isSelected = blxWindowIsSeqSelected(blxWindow, msp->sSequence);
	      
	      /* Draw the outline in the default SNP color. If the SNP is selected, also
	       * fill in the rectangle in the SNP color (use the selected color for the
	       * outline and the unselected color for the fill, so that the outline is darker). */
	      GdkColor *outlineColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, TRUE, bc->usePrintColors);
	      GdkColor *fillColor = isSelected ? getGdkColor(BLXCOLOR_SNP, bc->defaultColors, FALSE, bc->usePrintColors) : NULL;
	      
	      /* Draw the background rectangle for the char */
	      drawRectangle(drawable, gc, fillColor, outlineColor, x, y, width, properties->charHeight, TRUE, TRUE, TRUE, TRUE);
	      
	      /* Draw the text */
	      PangoLayout *layout = gtk_widget_create_pango_layout(detailView, mspGetMatchSeq(msp));
	      pango_layout_set_font_description(layout, properties->fontDesc);

	      if (layout)
		{
		  gtk_paint_layout(snpTrack->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, x, y, layout);
		  g_object_unref(layout);
		}	      
	    }
	}
    }
  
  drawColumnSeparatorLine(snpTrack, drawable, gc, bc);
  
  g_object_unref(gc);
}


/* This function centres the detail view display on the currently-selected base index. Does 
 * nothing if there isn't a selected base. */
static void detailViewCentreOnSelection(GtkWidget *detailView)
{
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


/* Gets the x coords at the start/end of the given column and populate them into the range
 * return argument. */
void detailViewGetColumnXCoords(GtkWidget *detailView, const BlxColumnId columnId, IntRange *xRange)
{
  xRange->min = 0;
  xRange->max = 0;
  
  /* Loop through all the columns up to the sequence column, summing their widths. */
  GList *columnItem = detailViewGetColumnList(detailView);
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(columnItem->data);
      
      if (columnInfo->columnId != columnId)
        {
          xRange->min += columnInfo->width;
        }
      else
        {
          xRange->max = xRange->min + columnInfo->width;
          break;
        }
    }
}


/* In the detail view, get the base index at the given coords, if those coords lie within the
 * sequence column; returns UNSET_INT otherwise. */
static int getBaseIndexAtDetailViewCoords(GtkWidget *detailView, const int x, const int y)
{
  int baseIdx = UNSET_INT;

  /* Get the x coords at the start/end of the sequence column */
  IntRange xRange;
  detailViewGetColumnXCoords(detailView, BLXCOL_SEQUENCE, &xRange);
  
  /* See if our x coord lies inside the sequence column */
  if (x >= xRange.min && x <= xRange.max)
    {
      /* Get the 0-based char index at x */
      gint charWidth = detailViewGetCharWidth(detailView);
      int charIdx = (int)(((double)x - xRange.min) / (double)charWidth);

      /* Add the start of the scroll range to convert this to the display index */
      GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
      baseIdx = charIdx + adjustment->value;
    }
  
  return baseIdx;
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
  
  if (tree && GTK_WIDGET_VISIBLE(tree))
    {
      GtkTreeViewColumn *column = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), BLXCOL_SEQUENCE);

      DetailViewColumnInfo *columnInfo = detailViewGetColumnInfo(detailView, BLXCOL_SEQUENCE);
      columnInfo->width = gtk_tree_view_column_get_width(column);
          
      /* Perform updates required on the sequence col after its size has changed */
      updateSeqColumnSize(detailView);
    }
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
      callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
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
      callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
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
  gint charWidth = pango_font_metrics_get_approximate_digit_width(metrics) / PANGO_SCALE;
  
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
    g_error("Detail-view widget is null\n");
  
  if (!GTK_IS_WIDGET(detailView))
    g_error("Detail-view is not a valid widget [%p]\n", detailView);
  
  if (!GTK_IS_CONTAINER(detailView))
    g_error("Detail-view is not a valid container [%p]\n", detailView);
  
  if (!detailViewGetProperties(detailView))
    g_error("Tree properties not set [widget=%p]\n", detailView);
}

GtkWidget* detailViewGetBlxWindow(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->blxWindow : NULL;
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

/* Get the list of columns */
GList* detailViewGetColumnList(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->columnList : NULL;
}

/* Get the column info for a particular column */
DetailViewColumnInfo *detailViewGetColumnInfo(GtkWidget *detailView, const BlxColumnId columnId)
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


gboolean detailViewGetDisplayRev(GtkWidget *detailView)
{
  return blxWindowGetDisplayRev(detailViewGetBlxWindow(detailView));
}

PangoFontDescription *detailViewGetFontDesc(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->fontDesc : NULL;
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
GtkWidget* detailViewGetTreeContainer(GtkWidget *detailView, const BlxStrand activeStrand, const int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *result = NULL;
  
  /* Get the list of trees for the relevant strand */
  GList *list = (activeStrand == BLXSTRAND_FORWARD) ? properties->fwdStrandTrees : properties->revStrandTrees;

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
GtkWidget* detailViewGetTree(GtkWidget *detailView, const BlxStrand activeStrand, const int frame)
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
      printf("Tree not found for '%s' strand, frame '%d'. Returning NULL.\n", ((activeStrand == BLXSTRAND_FORWARD) ? "forward" : "reverse"), frame);
    }
  
  return result;
}

/* Get the first visible tree in the 'current' list of trees (i.e. the forward 
 * strand list by default, or the reverse strand list if strands are toggled). */
static GtkWidget* detailViewGetFirstTree(GtkWidget *detailView)
{
  const gboolean toggled = detailViewGetDisplayRev(detailView);
  const BlxStrand activeStrand = toggled ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
  const int numFrames = detailViewGetNumFrames(detailView);
  
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
	  result = detailViewGetTree(detailView, toggled ? BLXSTRAND_FORWARD : BLXSTRAND_REVERSE, 1);
	}
    }
  
  return result;
}


int detailViewGetNumFrames(GtkWidget *detailView)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  return bc->numFrames;
}

static GtkWidget* detailViewGetHeader(GtkWidget *detailView)
{
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  return detailViewProperties->header;
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

IntRange* detailViewGetRefSeqRange(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetRefSeqRange(blxWindow);
}

BlxSeqType detailViewGetSeqType(GtkWidget *detailView)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  return bc->seqType;
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

/* Get the active frame. Returns the last-selected frame, or 1 if no frame is selected. */
int detailViewGetActiveFrame(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return (!properties || properties->selectedFrame == UNSET_INT) ? 1 : properties->selectedFrame;
}

/* Get the strand of the tree that was last selected (which defaults to the active strand if none is selected). */
BlxStrand detailViewGetSelectedStrand(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedStrand : BLXSTRAND_FORWARD;
}

/* Set the strand of the tree that was last selected. */
void detailViewSetSelectedStrand(GtkWidget *detailView, BlxStrand strand)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->selectedStrand = strand;
}


int detailViewGetNumUnalignedBases(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->numUnalignedBases : FALSE;
}


static int detailViewGetSnpConnectorHeight(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->snpConnectorHeight : UNSET_INT;
}


/* Perform required updates after the selected base has changed. */
static void updateFollowingBaseSelection(GtkWidget *detailView,
					 const gboolean allowScroll,
					 const gboolean scrollMinimum)
{
  if (allowScroll)
    {
      scrollToKeepSelectionInRange(detailView, scrollMinimum);
    }
  
  /* Update the feedback box */
  updateFeedbackBox(detailView);
  
  /* Update the headers so that the newly-selected index is highlighted */
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  gtk_widget_queue_draw(detailView);
}


/* Set the selected base index to a specific DNA index. Performs any required 
 * refreshes. Scrolls the view to keep the selected base in view if allowScroll 
 * is true. (Such scrolling is by the minimum number of bases necessary if 
 * scrollMinimum is true.) */
static void detailViewSetSelectedDnaBaseIdx(GtkWidget *detailView,
					    const int selectedDnaBaseIdx,
					    const int frame,
					    const gboolean allowScroll,
					    const gboolean scrollMinimum)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  properties->selectedDnaBaseIdx = selectedDnaBaseIdx;
  properties->selectedFrame = frame;
  
  /* For protein matches, calculate the display index and base number of this dna idx */
  BlxViewContext *bc = detailViewGetContext(detailView);
  properties->selectedBaseIdx = convertDnaIdxToDisplayIdx(selectedDnaBaseIdx, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &properties->selectedBaseNum);
  
  updateFollowingBaseSelection(detailView, allowScroll, scrollMinimum);
}


/* Set the selected base index to the given display index and base/frame number.
 *  Performs any required refreshes. Scrolls the view to keep the selected base 
 * in view if allowScroll is true. (Such scrolling is by the minimum
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
  BlxViewContext *bc = detailViewGetContext(detailView);
  properties->selectedDnaBaseIdx = convertDisplayIdxToDnaIdx(selectedBaseIdx, bc->seqType, frame, baseNum, bc->numFrames, bc->displayRev, &bc->refSeqRange);

  updateFollowingBaseSelection(detailView, allowScroll, scrollMinimum);
}

static GtkWidget *detailViewGetBigPicture(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetBigPicture(blxWindow);
}


DetailViewProperties* detailViewGetProperties(GtkWidget *widget)
{
  return widget ? (DetailViewProperties*)(g_object_get_data(G_OBJECT(widget), "DetailViewProperties")) : NULL;
}

static BlxViewContext* detailViewGetContext(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetContext(blxWindow);
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
  
  if (properties->spliceSites)
    {
      g_slist_foreach(properties->spliceSites, destroyBlxSpliceSite, NULL);
      g_slist_free(properties->spliceSites);
      properties->spliceSites = NULL;
    }
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "DetailViewProperties", NULL);
    }
}


static void detailViewCreateProperties(GtkWidget *detailView,
				       GtkWidget *blxWindow,
				       GtkCellRenderer *renderer,
				       GList *fwdStrandTrees,
				       GList *revStrandTrees,
				       GtkWidget *header,
				       GtkWidget *feedbackBox,
				       GtkWidget *statusBar,
				       GList *columnList,
				       GtkAdjustment *adjustment, 
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
      
      properties->blxWindow = blxWindow;
      properties->renderer = renderer;
      properties->fwdStrandTrees = fwdStrandTrees;
      properties->revStrandTrees = revStrandTrees;
      properties->header = header;
      properties->feedbackBox = feedbackBox;
      properties->statusBar = statusBar;
      properties->columnList = columnList;
      properties->adjustment = adjustment;
      properties->selectedBaseIdx = UNSET_INT;
      properties->selectedBaseNum = UNSET_INT;
      properties->selectedFrame = UNSET_INT;
      properties->selectedDnaBaseIdx = UNSET_INT;
      properties->fontDesc = fontDesc;
      properties->charWidth = 0;
      properties->charHeight = 0;
      properties->snpConnectorHeight = DEFAULT_SNP_CONNECTOR_HEIGHT;
      properties->numUnalignedBases = DEFAULT_NUM_UNALIGNED_BASES;
      
      /* Add the splice sites that we want Blixem to identify as canonical */
      properties->spliceSites = NULL;
      addBlxSpliceSite(&properties->spliceSites, "GT", "AG", FALSE);
      addBlxSpliceSite(&properties->spliceSites, "GC", "AG", FALSE);
      addBlxSpliceSite(&properties->spliceSites, "AT", "AC", TRUE);

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
      properties->exonBoundaryLineWidth	  = 1;
      properties->exonBoundaryLineStyleStart = GDK_LINE_SOLID;
      properties->exonBoundaryLineStyleEnd   = GDK_LINE_SOLID;
      
      g_object_set_data(G_OBJECT(detailView), "DetailViewProperties", properties);
      g_signal_connect(G_OBJECT(detailView), "destroy", G_CALLBACK(onDestroyDetailView), NULL); 
    }
}


/* Get/set functions for the sequence column header's row number. (There is only
 * one property to set, so it's not worth having a struct for the properties.) */
void seqColHeaderSetRow(GtkWidget *header, const int frame)
{
  g_object_set_data(G_OBJECT(header), "seqColHeaderFrameNumber", GINT_TO_POINTER(frame));
}

int seqColHeaderGetRow(GtkWidget *header)
{
  return header ? GPOINTER_TO_INT(g_object_get_data(G_OBJECT(header), "seqColHeaderFrameNumber")) : UNSET_INT;
}

/* Gets the base number within the given reading frame that this header row corresponds to*/
int seqColHeaderGetBase(GtkWidget *header, const int frame, const int numFrames)
{
  const int row = seqColHeaderGetRow(header);
  
  int baseNum = UNSET_INT;
  
  if (row != UNSET_INT)
    {
      baseNum = row - (frame - 1);

      /* Cyclic decrement */
      if (baseNum < 1)
	{
	  baseNum += numFrames;
	}
    }
  
  return baseNum;
}


/* Get/set functions for the SNP track header strand. (There is only one property
 * to set, so it's not worth having a struct for the properties.) */
static void snpTrackSetStrand(GtkWidget *snpTrack, const int strand)
{
  g_object_set_data(G_OBJECT(snpTrack), "snpTrackStrand", GINT_TO_POINTER(strand));
}

static int snpTrackGetStrand(GtkWidget *snpTrack)
{
  return snpTrack ? GPOINTER_TO_INT(g_object_get_data(G_OBJECT(snpTrack), "snpTrackStrand")) : UNSET_INT;
}


/***********************************************************
 *                     Events                              *
 ***********************************************************/

/* Expose function for generic text headers. Draws vertical separator lines and then lets the
 * default handler continue.  */
gboolean onExposeGenericHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkGC *gc = gdk_gc_new(headerWidget->window);
  GtkWidget *detailView = GTK_WIDGET(data);
  const BlxViewContext *bc = detailViewGetContext(detailView);
  
  drawColumnSeparatorLine(headerWidget, headerWidget->window, gc, bc);
  
  return FALSE;
}


/* expose function to push a cached bitmap to screen */
gboolean onExposeDnaTrack(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkWindow *window = GTK_IS_LAYOUT(headerWidget) ? GTK_LAYOUT(headerWidget)->bin_window : headerWidget->window;
  
  if (window)
    {
      /* See if there's a cached drawable and, if not, create it */
      GdkDrawable *bitmap = widgetGetDrawable(headerWidget);

      if (!bitmap)
	{
	  /* There isn't a bitmap yet. Create it now. */
	  GtkWidget *detailView = GTK_WIDGET(data);
	  const BlxStrand activeStrand = detailViewGetDisplayRev(detailView) ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
	  
	  drawDnaTrack(headerWidget, detailView, activeStrand, seqColHeaderGetRow(headerWidget));
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


/* Expose handler for the SNP track header */
static gboolean onExposeSnpTrack(GtkWidget *snpTrack, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(snpTrack)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(snpTrack);
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  GtkWidget *detailView = GTK_WIDGET(data);
          drawSnpTrack(snpTrack, detailView);
	  
	  bitmap = widgetGetDrawable(snpTrack);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc2 = gdk_gc_new(window);
          gdk_draw_drawable(window, gc2, bitmap, 0, 0, 0, 0, -1, -1);
        }
      else
	{
	  g_critical("Failed to draw SNP track [%p] - could not create bitmap.\n", snpTrack);
	}
    }
  
  return TRUE;
}


/* Handler for when the mouse button is pressed in the SNP track header */
static gboolean onButtonPressSnpTrack(GtkWidget *snpTrack, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  switch (event->button)
  {
    case 1:
    {
      GtkWidget *detailView = GTK_WIDGET(data);
      GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
      
      if (event->type == GDK_BUTTON_PRESS) /* first click */
        {
          /* Select the variation that was clicked on.  */
          blxWindowDeselectAllSeqs(blxWindow);

          /* The SNP track is not the same width as the sequence column, so pass the
           * sequence column header so that we can convert to the correct coords */
          DetailViewColumnInfo *seqColInfo = detailViewGetColumnInfo(detailView, BLXCOL_SEQUENCE);

          selectClickedSnp(snpTrack, seqColInfo->headerWidget, detailView, event->x, event->y, FALSE, TRUE, UNSET_INT); /* SNPs are always expanded in the SNP track */
          
          refreshDetailViewHeaders(detailView);
          callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
        }      
      else if (event->type == GDK_2BUTTON_PRESS) /* double-click */
        {
          /* If a variation was double-clicked, open its URL in a browser. If multiple are selected,
           * use the last-selected one. */
          GList *seqItem = g_list_last(blxWindowGetSelectedSeqs(blxWindow));
          
          if (seqItem)
            {
              BlxSequence *seq = (BlxSequence*)(seqItem->data);
              
              /* (The url lives in the MSP, and there is only one MSP in a BlxSequence for a variation,
               * so perhaps the url should be moved to the BlxSequence...) */
              if (seq->type == BLXSEQUENCE_VARIATION && g_list_length(seq->mspList) > 0)
                {
                  const MSP const *msp = (const MSP const *)(seq->mspList->data);
                  
                  if (msp->url)
                    {
                      GError *error = NULL;
                      seqtoolsLaunchWebBrowser(msp->url, &error);
                      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
                    }
                  else
                    {
                      g_warning("Variation '%s' does not have a URL.\n", mspGetSName(msp));
                    }
                }
            }
        }

      handled = TRUE;
      break;
    }
  }
  
  return handled;
}


/* Parent handler for button presses anywhere on the detail view. */
static gboolean onButtonPressDetailView(GtkWidget *detailView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
  {
    case 2:
    {
      /* Middle button: select the base index at the clicked coords */
      int baseIdx = getBaseIndexAtDetailViewCoords(detailView, event->x, event->y);
      
      if (baseIdx != UNSET_INT)
        {
          /* For protein matches, select the first base in the triplet */
          const int baseNum = 1;
          const int frame = detailViewGetActiveFrame(detailView);
          detailViewSetSelectedBaseIdx(detailView, baseIdx, frame, baseNum, FALSE, TRUE);
        }
      
      handled = TRUE;
      break;
    }
      
    case 3:
    {
      /* Right button: show context menu. */
      GtkWidget *mainmenu = blxWindowGetMainMenu(detailViewGetBlxWindow(detailView));
      gtk_menu_popup (GTK_MENU(mainmenu), NULL, NULL, NULL, NULL, event->button, event->time);
      handled = TRUE;
      break;
    }
      
    default:
      break;
  };
  
  return handled;
}


static gboolean onMouseMoveDetailView(GtkWidget *detailView, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base
       * (but don't re-centre on the selected base until the mouse button is released). */
      const int baseIdx = getBaseIndexAtDetailViewCoords(detailView, event->x, event->y);
      
      if (baseIdx != UNSET_INT)
        {
          /* For protein matches, get the 1st base in the peptide */
          const int baseNum = 1;
          const int frame = detailViewGetActiveFrame(detailView);
          detailViewSetSelectedBaseIdx(detailView, baseIdx, frame, baseNum, FALSE, TRUE);
        }
      
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseDetailView(GtkWidget *detailView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Right button: show context menu (handled in button-press event) */
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


/* Handler for when a mouse button is pressed on a sequence column header widget.
 * This signal only gets connected for protein matches, where there are 3 header
 * widgets, each showing the DNA nucleotides for a specific reading frame. */
static gboolean onButtonPressSeqColHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
  {
    case 1:
      {
	GtkWidget *detailView = GTK_WIDGET(data);

	if (event->type == GDK_BUTTON_PRESS)
	  {
	    /* Select the SNP that was clicked on.  */
	    blxWindowDeselectAllSeqs(detailViewGetBlxWindow(detailView));
	    const int clickedBase = seqColHeaderGetBase(header, detailViewGetActiveFrame(detailView), detailViewGetNumFrames(detailView));

	    selectClickedSnp(header, NULL, detailView, event->x, event->y, FALSE, FALSE, clickedBase); /* SNPs are always un-expanded in the DNA track */
	    
	    refreshDetailViewHeaders(detailView);
	    callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
	    /* Double-click: toggle SNP track visibility */
	    BlxViewContext *bc = detailViewGetContext(detailView);
            gboolean showSnpTrack = !bc->flags[BLXFLAG_SHOW_SNP_TRACK];
            
	    bc->flags[BLXFLAG_SHOW_SNP_TRACK] = showSnpTrack;
            detailViewUpdateShowSnpTrack(detailView, showSnpTrack);
	  }

	handled = TRUE;
	break;
      }
      
    case 2:
      {
	/* Middle button: select the nucleotide that was clicked on. */
	GtkWidget *detailView = GTK_WIDGET(data);
	selectClickedNucleotide(header, detailView, event->x, event->y);
	handled = TRUE;
	break;
      }
      
    default:
      break;
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
          detailViewCentreOnSelection(detailView);
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
      GtkWidget *tree1 = detailViewGetTreeContainer(detailView, BLXSTRAND_FORWARD, 1);
      GtkWidget *tree2 = detailViewGetTreeContainer(detailView, BLXSTRAND_REVERSE, 1);
      
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
  BlxViewContext *blxContext = detailViewGetContext(detailView);
  GtkWidget *bigPicture = detailViewGetBigPicture(detailView);

  /* Update the flag */
  blxContext->displayRev = !blxContext->displayRev;
  
  /* Invert the display range */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const IntRange const *fullRange = &blxContext->fullDisplayRange;
  const int newStart = fullRange->max - displayRange->max + fullRange->min;
  setDetailViewStartIdx(detailView, newStart, blxContext->seqType);

  /* Invert the currently-selected index, if there is one. We want to select
   * the same index but counting from the other end. */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  if (properties->selectedBaseIdx != UNSET_INT)
    {
      const int newIdx = fullRange->max - properties->selectedBaseIdx + fullRange->min;
      const int newBaseNum = blxContext->numFrames - properties->selectedBaseNum + 1;
      
      detailViewSetSelectedBaseIdx(detailView, newIdx, properties->selectedFrame, newBaseNum, FALSE, TRUE);
    }
  
  /* If one grid/tree is hidden and the other visible, toggle which is hidden */
  swapTreeVisibility(detailView);
  swapGridVisibility(bigPicture);
  swapExonViewVisibility(bigPicture);
  
  /* Toggle the order of the trees and grids. */
  refreshTreeOrder(detailView);
  refreshGridOrder(bigPicture);

  /* Redraw the tree and detail-view headers. Also need to resort because if the display is
   * reversed it affects sorting by position. */
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, resortTree, NULL);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
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
	  const int activeFrame = detailViewGetActiveFrame(detailView);
	  const BlxViewContext *bc = detailViewGetContext(detailView);
	  int baseNum;
	  
	  const int displayIdx = convertDnaIdxToDisplayIdx(requestedCoord, 
							   bc->seqType, 
							   activeFrame,
							   bc->numFrames, 
							   bc->displayRev, 
							   &bc->refSeqRange,
							   &baseNum);
	  
	  /* Select the base index. */
	  detailViewSetSelectedBaseIdx(detailView, displayIdx, activeFrame, baseNum, TRUE, FALSE);
	}
    }

  gtk_widget_destroy(dialog);
}


/* Sort the detail view trees by the given column */
void detailViewSetSortColumn(GtkWidget *detailView, const BlxColumnId sortColumn)
{
  callFuncOnAllDetailViewTrees(detailView, treeSetSortColumn, GINT_TO_POINTER(sortColumn));
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

      /* Check its in the sequence list (if given), is a valid match or exon, and is in a visible layer */
      if (mspLayerIsVisible(msp) &&
          (mspIsExon(msp) || mspIsBlastMatch(msp)) &&
	  (!searchData->seqList || g_list_find(searchData->seqList, msp->sSequence)))
	{
	  /* Get the offset of the msp coords from the given start coord and find the smallest,
	   * ignorning zero and negative offsets (negative means its the wrong direction) */
	  const int offset1 = (msp->qRange.min - searchData->startDnaIdx) * searchData->searchDirection;
	  const int offset2 = (msp->qRange.max - searchData->startDnaIdx) * searchData->searchDirection;
	  
	  int currentFrame = msp->qFrame;
	  int currentBest = UNSET_INT;
	  
	  if (offset1 > 0 && (offset2 <= 0 || offset2 >= offset1))
	    {
	      currentBest = offset1;
	    }
	  else
	    {
	      currentBest = offset2;
	    }
	  if (currentBest > 0)
	    {
	      gboolean useNew = FALSE;

	      /* Check if this is the first one we've looked at */
	      useNew |= (searchData->offset == UNSET_INT);
	      
	      /* Check if this offset is smaller than the previous best */
	      useNew |= (currentBest < searchData->offset);
	      
	      /* If it's the same as the previous best, use the one with the smallest 
	       * frame number, i.e. give preference to the topmost frame (in the hope
	       * that this will be the least confusing!) */
	      useNew |= (currentBest == searchData->offset && currentFrame < searchData->foundFrame);
	      
	      if (useNew)
		{
		  searchData->offset = currentBest;
		  searchData->foundFrame = currentFrame;
		}
	    }
	}
    }
  
  return FALSE;
}


/* Find and go to the next match (either left or right depending on the search flag).
 * If a list of sequence names is given, only look at matches in those sequences.
 * startDnaIdx determines where to start searching from. Sets the found match idx as
 * the currently-selected base index */
static void goToNextMatch(GtkWidget *detailView, const int startDnaIdx, const gboolean searchRight, GList *seqList)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  const int searchDirection = (searchRight != bc->displayRev) ? 1 : -1;

  MatchSearchData searchData = {startDnaIdx, 
				searchRight, 
				searchDirection, 
				bc->displayRev, 
				bc->seqType,
				bc->numFrames,
				&bc->refSeqRange,
				seqList,
				UNSET_INT,
				UNSET_INT,
				UNSET_INT};

  /* Loop through the MSPs in all visible trees */
  int frame = 1;

  for (  ; frame <= searchData.numFrames; ++frame)
    {
      GtkWidget *treeContainer = detailViewGetTreeContainer(detailView, BLXSTRAND_FORWARD, frame);
      GtkWidget *tree = treeContainerGetTree(GTK_CONTAINER(treeContainer));
      
      if (GTK_WIDGET_VISIBLE(tree) && gtk_widget_get_parent(treeContainer)) /* ignore if not currently included in view */
	{
	  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
	  gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
	}

      treeContainer = detailViewGetTreeContainer(detailView, BLXSTRAND_REVERSE, frame);
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
      int newDnaIdx = searchData.startDnaIdx + (searchData.offset * searchDirection);
      detailViewSetSelectedDnaBaseIdx(detailView, newDnaIdx, searchData.foundFrame, TRUE, FALSE);
    }
}


/* Go to the previous match (optionally limited to matches in the given list)  */
void prevMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startDnaIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  
  if (!valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
      
      /* Use base 1 within the currently selected frame for this display coord */
      int frame = detailViewGetActiveFrame(detailView);
      BlxViewContext *bc = detailViewGetContext(detailView);
      startDnaIdx = convertDisplayIdxToDnaIdx(startCoord, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  goToNextMatch(detailView, startDnaIdx, FALSE, seqList);
}


/* Go to the next match (optionally limited to matches in the given list)  */
void nextMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startDnaIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
  
  if (!valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
      
      /* Use base 1 within the currently selected frame for this display coord */
      int frame = detailViewGetActiveFrame(detailView);
      BlxViewContext *bc = detailViewGetContext(detailView);
      startDnaIdx = convertDisplayIdxToDnaIdx(startCoord, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  goToNextMatch(detailView, startDnaIdx, TRUE, seqList);
}


/* Go to the first match (optionally limited to matches in the given list)  */
void firstMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the start of the ref seq */
  const IntRange const *refSeqRange = detailViewGetRefSeqRange(detailView);
  const int startIdx = detailViewGetDisplayRev(detailView) ? refSeqRange->max : refSeqRange->min;
  
  goToNextMatch(detailView, startIdx, TRUE, seqList);
}


/* Go to the last match (optionally limited to matches in the given list)  */
void lastMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the end of the reference sequence */
  const IntRange const *refSeqRange = detailViewGetRefSeqRange(detailView);
  const int startIdx = detailViewGetDisplayRev(detailView) ? refSeqRange->min : refSeqRange->max;

  goToNextMatch(detailView, startIdx, FALSE, seqList);
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
      
      BlxColumnId sortColumn = g_value_get_int(&val);
      detailViewSetSortColumn(detailView, sortColumn);
    }
}


static void GHelp(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  showHelpDialog(detailViewGetBlxWindow(detailView), TRUE);
}

static void GSettings(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  showSettingsDialog(detailViewGetBlxWindow(detailView), TRUE);
}

static void GFind(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  showFindDialog(detailViewGetBlxWindow(detailView), TRUE);
}

//static void GInfo(GtkButton *button, gpointer data)
//{
//  GtkWidget *detailView = GTK_WIDGET(data);
//  showInfoDialog(detailViewGetBlxWindow(detailView));
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
  GList *seqList = blxWindowGetSelectedSeqs(detailViewGetBlxWindow(detailView));
  prevMatch(detailView, seqList);
}

static void GnextMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(detailViewGetBlxWindow(detailView));
  nextMatch(detailView, seqList);
}

static void GfirstMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(detailViewGetBlxWindow(detailView));
  firstMatch(detailView, seqList);
}

static void GlastMatch(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(detailViewGetBlxWindow(detailView));
  lastMatch(detailView, seqList);
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

/* Comparison function for two DetailViewColumnInfo structs
 * Returns : negative value if a < b; zero if a = b; positive value if a > b.  */
static gint columnCompareFunc(gconstpointer a, gconstpointer b)
{
  DetailViewColumnInfo *col1 = (DetailViewColumnInfo*)a;
  DetailViewColumnInfo *col2 = (DetailViewColumnInfo*)b;
  
  return col1->columnId - col2->columnId;
}

/* Creates a detail-view column from the given info and adds it to the columnList. */
static void createColumn(BlxColumnId columnId, 
                         GtkWidget *specialWidget,
                         GtkCallback callbackFn, 
                         char *title,
                         char *propertyName,
                         const int defaultWidth,
                         const gboolean optionalDataLoaded,
                         char *sortName,
                         GList **columnList,
                         GtkWidget *detailView)
{
  /* Create a simple label for the header (unless already passed a special header widget) */
  GtkWidget *headerWidget = specialWidget;
  
  if (headerWidget == NULL)
    {
      headerWidget = gtk_label_new(title);
      g_signal_connect(G_OBJECT(headerWidget), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
    }
  
  gtk_widget_set_size_request(headerWidget, defaultWidth, -1);
  
  if (GTK_IS_MISC(headerWidget))
    {
      /* Align the text bottom-left */
      gtk_misc_set_alignment(GTK_MISC(headerWidget), 0.0, 1.0);
      gtk_misc_set_padding(GTK_MISC(headerWidget), DEFAULT_LABEL_X_PAD, 0);
    }

  /* Create the column info */
  DetailViewColumnInfo *columnInfo = g_malloc(sizeof(DetailViewColumnInfo));
  
  columnInfo->columnId = columnId;
  columnInfo->headerWidget = headerWidget;
  columnInfo->refreshFunc = callbackFn,
  columnInfo->title = title;
  columnInfo->propertyName = propertyName;
  columnInfo->width = defaultWidth;
  columnInfo->sortName = sortName;
  columnInfo->dataLoaded = optionalDataLoaded;
  
  /* Place it in the list. Sort the list by BlxColumnId because the list must be sorted in the same
   * order as the variable types in the TREE_COLUMN_TYPE_LIST definition */
  *columnList = g_list_insert_sorted(*columnList, columnInfo, columnCompareFunc);
}


/* This creates DetailViewColumnInfo entries for each column required in the detail view. It
 * returns a list of the columns created. */
static GList* createColumns(GtkWidget *detailView, const BlxSeqType seqType, const int numFrames, const gboolean optionalDataLoaded)
{
  GList *columnList = NULL;
  
  /* The sequence column has a special header widget and callback when we're dealing 
   * with peptide sequences. This returns NULL for DNA sequences, in which case createColumn
   * will create a simple label header for us instead. */
  GtkWidget *seqHeader = createSeqColHeader(detailView, seqType, numFrames);
  GtkCallback seqCallback = (seqType == BLXSEQ_PEPTIDE) ? refreshTextHeader : NULL;
  
  /* The start and end columns have special callbacks to switch the start/end text when display is toggled */
  GtkCallback startCallback = refreshStartColHeader;
  GtkCallback endCallback = refreshEndColHeader;
  
  /* Create the column headers and pack them into the column header bar */
  createColumn(BLXCOL_SEQNAME,   NULL,       NULL,          "Name",      RENDERER_TEXT_PROPERTY,     BLXCOL_SEQNAME_WIDTH,        TRUE,  "Name",     &columnList, detailView);
  createColumn(BLXCOL_SCORE,     NULL,       NULL,          "Score",     RENDERER_TEXT_PROPERTY,     BLXCOL_INT_COLUMN_WIDTH,     TRUE,  "Score",    &columnList, detailView);
  createColumn(BLXCOL_ID,        NULL,       NULL,          "%Id",       RENDERER_TEXT_PROPERTY,     BLXCOL_INT_COLUMN_WIDTH,     TRUE,  "Identity", &columnList, detailView);
  createColumn(BLXCOL_START,     NULL,       startCallback, "Start",     RENDERER_TEXT_PROPERTY,     BLXCOL_START_WIDTH,          TRUE,  "Position", &columnList, detailView);
  createColumn(BLXCOL_SEQUENCE,  seqHeader,  seqCallback,   "Sequence",  RENDERER_SEQUENCE_PROPERTY, BLXCOL_SEQUENCE_WIDTH,       TRUE,  NULL,       &columnList, detailView);
  createColumn(BLXCOL_END,       NULL,       endCallback,   "End",       RENDERER_TEXT_PROPERTY,     BLXCOL_END_WIDTH,            TRUE,  NULL,       &columnList, detailView);
  createColumn(BLXCOL_SOURCE,    NULL,       NULL,          "Source",    RENDERER_TEXT_PROPERTY,     BLXCOL_HIDDEN_COLUMN_WIDTH,  TRUE,  NULL,       &columnList, detailView);
  createColumn(BLXCOL_GROUP,     NULL,       NULL,          "Group",     RENDERER_TEXT_PROPERTY,     BLXCOL_HIDDEN_COLUMN_WIDTH,  TRUE,  "Group",    &columnList, detailView);
  createColumn(BLXCOL_ORGANISM,  NULL,       NULL,          "Organism",  RENDERER_TEXT_PROPERTY,     BLXCOL_HIDDEN_COLUMN_WIDTH,  optionalDataLoaded,  "Organism", &columnList, detailView);
  createColumn(BLXCOL_GENE_NAME, NULL,       NULL,          "Gene Name", RENDERER_TEXT_PROPERTY,     BLXCOL_HIDDEN_COLUMN_WIDTH,  optionalDataLoaded, "Gene name", &columnList, detailView);
  createColumn(BLXCOL_TISSUE_TYPE,NULL,      NULL,          "Tissue Type",RENDERER_TEXT_PROPERTY,    BLXCOL_HIDDEN_COLUMN_WIDTH,  optionalDataLoaded, "Tissue type",&columnList, detailView);
  createColumn(BLXCOL_STRAIN,    NULL,       NULL,          "Strain",    RENDERER_TEXT_PROPERTY,     BLXCOL_HIDDEN_COLUMN_WIDTH,  optionalDataLoaded, "Strain",   &columnList, detailView);

  return columnList;
}


/* This loops through all the columns and adds each column's header widget to the header bar. */
static void addColumnsToHeaderBar(GtkBox *headerBar, GList *columnList)
{
  /* Loop through each column and add its header widget to the header bar */
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(columnItem->data);

      /* The sequence column is a special one that wants to fill any additional space, so
       * set the 'expand' property to true for that column only. */
      const gboolean expand = (columnInfo->columnId == BLXCOL_SEQUENCE );
      gtk_box_pack_start(headerBar, columnInfo->headerWidget, expand, TRUE, 0);
    }
}


/* Create the header bar for the detail view. This contains the labels for the
 * detail-view trees (since we only want one label at the top, rather than 
 * individual labels for each tree). For protein sequence matches, the header
 * for the sequence column will also show the DNA sequence (separated into reading
 * frames). 
 * Column data is compiled into the detailViewColumns return argument. */
static GtkWidget* createDetailViewHeader(GtkWidget *detailView, 
					 const BlxSeqType seqType, 
					 const int numFrames,
					 GList *columnList,
					 const gboolean includeSnpTrack)
{
  GtkBox *header = NULL; /* outermost container for the header */
  GtkBox *headerBar = GTK_BOX(gtk_hbox_new(FALSE, 0)); /* hbox for the column headers */
  
  /* Create a SNP track, if requested */
  if (includeSnpTrack)
    {
      header = GTK_BOX(gtk_vbox_new(FALSE, 0));
      createSnpTrackHeader(GTK_BOX(header), detailView, UNSET_INT);
      gtk_box_pack_start(header, GTK_WIDGET(headerBar), FALSE, FALSE, 0);
      gtk_widget_set_name(GTK_WIDGET(headerBar), HEADER_CONTAINER_NAME);
    }
  else
    {
      /* Only need the column header bar */
      header = headerBar;
    }

  /* Add all the column headers to the header bar */
  addColumnsToHeaderBar(headerBar, columnList);

  return GTK_WIDGET(header);
}


/* Create the SNP track header widget. The strand can be passed as UNSET_INT,
 * in which case the active strand will be used. */
GtkWidget* createSnpTrackHeader(GtkBox *parent, GtkWidget *detailView, const int strand)
{
  GtkWidget *snpTrack = gtk_layout_new(NULL, NULL);

  gtk_widget_set_name(snpTrack, SNP_TRACK_HEADER_NAME);  
  gtk_box_pack_start(parent, snpTrack, FALSE, TRUE, 0);
  snpTrackSetStrand(snpTrack, strand);

  gtk_widget_add_events(snpTrack, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(snpTrack), "expose-event", G_CALLBACK(onExposeSnpTrack), detailView);
  g_signal_connect(G_OBJECT(snpTrack), "button-press-event", G_CALLBACK(onButtonPressSnpTrack), detailView);
  
  return snpTrack;
}


/* Create a custom header widget for the sequence column for protein matches (this is where
 * we will display the triplets that make up the codons.) Returns NULL for DNA matches. */
static GtkWidget* createSeqColHeader(GtkWidget *detailView,
				     const BlxSeqType seqType,
				     const int numFrames)
{
  GtkWidget *header = NULL;
  
  if (seqType == BLXSEQ_PEPTIDE)
    {
      header = gtk_vbox_new(FALSE, 0);
      gtk_widget_set_name(header, HEADER_CONTAINER_NAME);
      
      int frame = 0;
      for ( ; frame < numFrames; ++frame)
	{
	  GtkWidget *line = gtk_layout_new(NULL, NULL);
	  gtk_box_pack_start(GTK_BOX(header), line, FALSE, TRUE, 0);

	  gtk_widget_set_name(line, DNA_TRACK_HEADER_NAME);
	  seqColHeaderSetRow(line, frame + 1);
	  
	  gtk_widget_add_events(line, GDK_BUTTON_PRESS_MASK);
	  gtk_widget_add_events(line, GDK_BUTTON_RELEASE_MASK);
	  gtk_widget_add_events(line, GDK_POINTER_MOTION_MASK);
	  g_signal_connect(G_OBJECT(line), "expose-event", G_CALLBACK(onExposeDnaTrack), detailView);
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
static GtkWidget* createEmptyButtonBar(GtkToolbar **toolbar)
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
				  BlxColumnId sortColumn, 
				  const char *sortName,
				  BlxColumnId initSortColumn,
				  GtkComboBox *combo)
{
  GtkTreeIter iter;
  gtk_tree_store_append(store, &iter, parent);

  gtk_tree_store_set(store, &iter, SORT_TYPE_COL, sortColumn, SORT_TEXT_COL, sortName, -1);

  if (sortColumn == initSortColumn)
    {
      gtk_combo_box_set_active_iter(combo, &iter);
    }
  
  return NULL;
}


/* Create the combo box used for selecting sort criteria */
static void createSortBox(GtkToolbar *toolbar, GtkWidget *detailView, const BlxColumnId initSortColumn, GList *columnList)
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

  /* Add an option to sort by each column that has the 'sortName' property set. */
  GtkTreeIter *iter = NULL;
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(columnItem->data);
      
      if (columnInfo->sortName)
        {
          iter = addSortBoxItem(store, iter, columnInfo->columnId, columnInfo->sortName, initSortColumn, combo);
        }
    }
  
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


/* Create the status bar for the detail-view toolbar. (This feeds back info to the user 
 * about the currently-moused-over sequence.) */
static GtkWidget* createStatusBar(GtkToolbar *toolbar)
{
  GtkWidget *statusBar = gtk_statusbar_new() ;
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), FALSE);

  /* Make it expandable so we use all available space. Set minimum size to be 0
   * because it's better to show it small than not at all. */
  gtk_widget_set_size_request(statusBar, 0, -1) ;
  GtkToolItem *item = addToolbarWidget(toolbar, statusBar) ;
  gtk_tool_item_set_expand(item, TRUE); /* make as big as possible */

  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");

  return statusBar;
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
					    const BlxColumnId sortColumn,
                                            GList *columnList,
					    GtkWidget **feedbackBox,
					    GtkWidget **statusBar)
{
  GtkToolbar *toolbar = NULL;
  GtkWidget *toolbarContainer = createEmptyButtonBar(&toolbar);
  
  /* Help */
  makeToolbarButton(toolbar, "Help", GTK_STOCK_HELP,	    "Help (Ctrl-H)",			(GtkSignalFunc)GHelp,		  detailView);

  /* Combo box for sorting */
  createSortBox(toolbar, detailView, sortColumn, columnList);

  /* Settings button */
  makeToolbarButton(toolbar, "Settings", GTK_STOCK_PREFERENCES,  "Settings (Ctrl-S)",		 (GtkSignalFunc)GSettings,	  detailView);

  /* Zoom buttons */
  makeToolbarButton(toolbar, "Zoom in",		GTK_STOCK_ZOOM_IN,  "Zoom in (=)",		 (GtkSignalFunc)onZoomInDetailView, detailView);
  makeToolbarButton(toolbar, "Zoom out",	GTK_STOCK_ZOOM_OUT, "Zoom out (-)",		 (GtkSignalFunc)onZoomOutDetailView, detailView);
  
  /* Navigation buttons */
  makeToolbarButton(toolbar, "Go to",		GTK_STOCK_JUMP_TO,    "Go to position (p)",		  (GtkSignalFunc)GGoto,		  detailView);

  makeToolbarButton(toolbar, "First match",	GTK_STOCK_GOTO_FIRST, "First match (Ctrl-Home)",	  (GtkSignalFunc)GfirstMatch,	  detailView);
  makeToolbarButton(toolbar, "Previous match",	GTK_STOCK_GO_BACK,    "Previous match (Ctrl-left)",	  (GtkSignalFunc)GprevMatch,	  detailView);
  makeToolbarButton(toolbar, "Next match",	GTK_STOCK_GO_FORWARD, "Next match (Ctrl-right)",	  (GtkSignalFunc)GnextMatch,	  detailView);
  makeToolbarButton(toolbar, "Last match",	GTK_STOCK_GOTO_LAST,  "Last match (Ctrl-End)",		  (GtkSignalFunc)GlastMatch,	  detailView);
  insertToolbarSeparator(toolbar);
  
  makeToolbarButton(toolbar, "<<", NULL,	"Scroll back one page (Ctrl-,)",    (GtkSignalFunc)GscrollLeftBig,  detailView);
  makeToolbarButton(toolbar, "<",  NULL,	"Scroll back one index (,)",	    (GtkSignalFunc)GscrollLeft1,    detailView);
  makeToolbarButton(toolbar, ">",  NULL,	"Scroll forward one index (.)",	    (GtkSignalFunc)GscrollRight1,   detailView);
  makeToolbarButton(toolbar, ">>", NULL,	"Scroll forward one page (Ctrl-.)", (GtkSignalFunc)GscrollRightBig, detailView);
  
  /* Find/Msp-info */
  makeToolbarButton(toolbar, "Find",          GTK_STOCK_FIND,    "Find sequences (f, Ctrl-F)",                      (GtkSignalFunc)GFind,  detailView);
//  makeToolbarButton(toolbar, "Sequence info", GTK_STOCK_INFO,    "Display info about the selected sequence(s) (i)", (GtkSignalFunc)GInfo,  detailView);

  /* Strand toggle button */
  if (mode == BLXMODE_BLASTX || mode == BLXMODE_TBLASTX || mode == BLXMODE_BLASTN)
    {
      makeToolbarButton(toolbar, "Toggle strand", GTK_STOCK_REFRESH, "Toggle strand (t)", (GtkSignalFunc)GToggleStrand, detailView);
    }

  *feedbackBox = createFeedbackBox(toolbar);
  *statusBar = createStatusBar(toolbar);

  return toolbarContainer;
}


/* Create two detail-view trees and place them in the 2 panes of the container (if one is
 * given). The first tree gets associated with grid1 and appended to list1, and the second
 * with grid2 and list2. If hasHeaders is true, the first tree will have visible headers. */
static void createTwoPanedTrees(GtkWidget *detailView,
				GtkPaned *panedWin, 
				GtkCellRenderer *renderer,
				GtkWidget *grid1, 
				GtkWidget *grid2,
				GList **list1,
				GList **list2,
				BlxSeqType seqType,
				GList *columnList,
				const char const *refSeqName,
				const int frame1,
				const int frame2,
				const gboolean includeSnpTrack)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, detailView, renderer, list1, columnList, seqType, refSeqName, frame1, includeSnpTrack);
  GtkWidget *tree2 = createDetailViewTree(grid2, detailView, renderer, list2, columnList, seqType, refSeqName, frame2, includeSnpTrack);
  
  if (panedWin)
    {
      gtk_paned_pack1(GTK_PANED(panedWin), tree1, TRUE, TRUE);
      gtk_paned_pack2(GTK_PANED(panedWin), tree2, TRUE, TRUE);
    }
}


/* Create three detail-view trees and place them in the paned widget. We only have
 * two panes in the widget, so two of the trees will be placed in a second, nested
 * paned widget. */
static void createThreePanedTrees(GtkWidget *detailView,
				  GtkPaned *panedWin,
				  GtkCellRenderer *renderer,
				  GtkWidget *grid,
				  GList **list,
				  BlxSeqType seqType,
				  const gboolean addToDetailView,
				  GList *columnList,
				  const char const *refSeqName,
				  const gboolean includeSnpTrack)
{
  const int frame1 = 1, frame2 = 2, frame3 = 3;
  
  /* Create a tree for pane1 (but only add it to the detailView if instructed to).
   * The first tree has headers. */
  GtkWidget *tree1 = createDetailViewTree(grid, detailView, renderer, list, columnList, seqType, refSeqName, frame1, includeSnpTrack);

  GtkPaned *nestedPanedWidget = NULL;
  if (addToDetailView)
    {
      gtk_paned_pack1(GTK_PANED(panedWin), tree1, TRUE, TRUE);
      
      /* Create another paned widget and place it in pane 2 */
      nestedPanedWidget = GTK_PANED(gtk_vpaned_new());
      gtk_paned_pack2(panedWin, GTK_WIDGET(nestedPanedWidget), TRUE, TRUE);
    }
  
  /* Create two more trees (and place them in the nested paned widget, if it is not null). 
   * Neither of these trees should have headers. */
  createTwoPanedTrees(detailView, nestedPanedWidget, renderer, grid, grid, list, list, seqType, columnList, refSeqName, frame2, frame3, includeSnpTrack);
}


/* Create the trees for the detail view, creating sub-panes if necessary depending
 * on how many trees we need */
static void createDetailViewPanes(GtkWidget *detailView, 
				  GtkPaned *panedWin,
				  GtkCellRenderer *renderer,
				  GtkWidget *fwdStrandGrid, 
				  GtkWidget *revStrandGrid, 
				  const int numFrames,
				  GList **fwdStrandTrees,
				  GList **revStrandTrees,
				  BlxSeqType seqType,
				  GList *columnList,
				  const char const *refSeqName,
				  const gboolean includeSnpTrack)
{
  if (numFrames == 1)
    {
      /* DNA matches: we need 2 trees, one for the forward strand and one for the reverse. */
      createTwoPanedTrees(detailView, 
			  panedWin, 
			  renderer,
			  fwdStrandGrid, 
			  revStrandGrid, 
			  fwdStrandTrees, 
			  revStrandTrees, 
			  seqType,
			  columnList,
			  refSeqName,
			  1,
			  1,
			  includeSnpTrack);
    }
  else if (numFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames for EACH strand (although only
       * one set of trees will be displayed at a time). */
      createThreePanedTrees(detailView, panedWin, renderer, fwdStrandGrid, fwdStrandTrees, seqType, TRUE, columnList, refSeqName, includeSnpTrack);
      createThreePanedTrees(detailView, panedWin, renderer, revStrandGrid, revStrandTrees, seqType, FALSE, columnList, refSeqName, includeSnpTrack);
    }
}



/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* Add the MSPs to the detail-view trees. Calculates and returns the lowest ID 
 * out of all the blast matches */
void detailViewAddMspData(GtkWidget *detailView, MSP *mspList)
{
  BlxViewContext *bc = detailViewGetContext(detailView);

  /* First, create a data store for each tree so we have something to add our data to. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateBaseDataModel, NULL);

  /* Loop through each MSP, and add it to the correct tree based on its strand and
   * reading frame. Also find the lowest ID value out of all the matches. */
  MSP *msp = mspList;
  
  for ( ; msp; msp = msp->next)
    {
      /* Only add matches/exons to trees */
      if (msp->type == BLXMSP_MATCH || msp->type == BLXMSP_EXON)
	{
	  /* Find the tree that this MSP should belong to based on its reading frame and strand */
	  BlxStrand strand = mspGetRefStrand(msp);
	  const int frame = mspGetRefFrame(msp, bc->seqType);
	  GtkWidget *tree = detailViewGetTree(detailView, strand, frame);

	  if (tree)
	    {
	      addMspToTree(tree, msp);
	    }
	  else
	    {
	      printf("Error: could not determine alignment list. Sequence may not be shown. (sequence = '%s', q range [%d-%d], s range [%d-%d], q frame=%s)\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max, msp->qframe);
	    }
	}
    }
    
  /* Finally, create a custom-filtered version of the data store for each tree. We do 
   * this AFTER adding the data so that it doesn't try to re-filter every time we add a row. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateFilteredDataModel, NULL);
  
  /* Also create a second data store that will store one sequence per row (as opposed to one
   * MSP per row). This data store will be switched in when the user selects 'squash matches'. */
  callFuncOnAllDetailViewTrees(detailView, addSequencesToTree, NULL);
}


/* We need a fixed-width font for displaying alignments. Find one from a 
 * list of possibilities. */
static const char* findDetailViewFont(GtkWidget *detailView)
{
  GList *fixed_font_list = NULL ;

  fixed_font_list = g_list_append(fixed_font_list, "andale mono");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida sans typewriter");
  fixed_font_list = g_list_append(fixed_font_list, "deja vu sans mono");
  fixed_font_list = g_list_append(fixed_font_list, "Bitstream vera sans mono");
  fixed_font_list = g_list_append(fixed_font_list, "monaco");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida console");
  fixed_font_list = g_list_append(fixed_font_list, "Courier 10 pitch");
  fixed_font_list = g_list_append(fixed_font_list, "Courier new");
  fixed_font_list = g_list_append(fixed_font_list, "Courier");
  fixed_font_list = g_list_append(fixed_font_list, "Monospace");
  fixed_font_list = g_list_append(fixed_font_list, "fixed");
  
  const char *fontFamily = findFixedWidthFontFamily(detailView, fixed_font_list);
  g_list_free(fixed_font_list);
  
  return fontFamily;
}


GtkWidget* createDetailView(GtkWidget *blxWindow,
			    GtkContainer *parent,
			    GtkAdjustment *adjustment, 
			    GtkWidget *fwdStrandGrid, 
			    GtkWidget *revStrandGrid,
			    MSP *mspList,
			    BlxBlastMode mode,
			    BlxSeqType seqType,
			    int numFrames,
			    const char const *refSeqName,
			    const int startCoord,
			    const gboolean sortInverted,
			    const BlxColumnId sortColumn,
                            const gboolean optionalDataLoaded)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so that we can identify them.) */
  GtkWidget *detailView = gtk_vbox_new(FALSE, 0);
  gtk_container_add(parent, detailView);
  
  GtkWidget *panedWin = gtk_vpaned_new();
  gtk_widget_set_name(panedWin, DETAIL_VIEW_WIDGET_NAME);

  /* Create a custom cell renderer to render the sequences in the detail view */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();

  /* Create the columns */
  GList *columnList = createColumns(detailView, seqType, numFrames, optionalDataLoaded);

  /* Create the header bar. If viewing protein matches include one SNP track in the detail 
   * view header; otherwise create SNP tracks in each tree header. */
  const gboolean singleSnpTrack = (seqType == BLXSEQ_PEPTIDE);
  GtkWidget *header = createDetailViewHeader(detailView, seqType, numFrames, columnList, singleSnpTrack);

  /* Create the toolbar. We need to remember the feedback box and status bar so we can set them in the properties. */
  GtkWidget *feedbackBox = NULL;
  GtkWidget *statusBar = NULL;
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, mode, sortColumn, columnList, &feedbackBox, &statusBar);
  
  /* Create the trees. */
  GList *fwdStrandTrees = NULL, *revStrandTrees = NULL;
  createDetailViewPanes(detailView, 
			GTK_PANED(panedWin),
			renderer,
			fwdStrandGrid, 
			revStrandGrid, 
			numFrames, 
			&fwdStrandTrees,
			&revStrandTrees,
			seqType,
			columnList,
			refSeqName,
			!singleSnpTrack);
    
  gtk_box_pack_start(GTK_BOX(detailView), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(detailView), header, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(detailView), panedWin, TRUE, TRUE, 0);

  /* Connect signals */
  gtk_widget_add_events(detailView, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(detailView), "size-allocate", G_CALLBACK(onSizeAllocateDetailView), NULL);
  g_signal_connect(G_OBJECT(detailView), "button-press-event", G_CALLBACK(onButtonPressDetailView), NULL);
  g_signal_connect(G_OBJECT(detailView), "button-release-event", G_CALLBACK(onButtonReleaseDetailView), NULL);
  g_signal_connect(G_OBJECT(detailView), "motion-notify-event", G_CALLBACK(onMouseMoveDetailView), NULL);

  detailViewCreateProperties(detailView, 
			     blxWindow, 
			     renderer,
			     fwdStrandTrees,
			     revStrandTrees,
			     header,
			     feedbackBox, 
			     statusBar,
			     columnList,
			     adjustment, 
			     startCoord,
			     sortInverted);
  
  return detailView;
}

