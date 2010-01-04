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
#include <gtk/gtk.h>

#define DETAIL_VIEW_TOOLBAR_NAME	"DetailViewToolbarName"
#define SORT_BY_SCORE_STRING		"Score"
#define SORT_BY_ID_STRING		"Identity"
#define SORT_BY_NAME_STRING		"Name"
#define SORT_BY_POSITION_STRING		"Position"


/* Local function declarations */
static GtkToolItem*	    addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget);
static GtkWidget*	    detailViewGetFirstTree(GtkWidget *detailView);
static GtkWidget*	    detailViewGetBigPicture(GtkWidget *detailView);


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

  gint i;
  for (i = 0 ; (i < n_families && !found_most_preferred) ; i++)
    {
      const gchar *name = pango_font_family_get_name(families[i]) ;

      /* Look for this font family in our list of preferred families */
      GList *pref = g_list_first(pref_families) ;
      gint current = 1;
      
      while(pref)
	{
	  char *pref_font = (char *)pref->data ;
	  
	  if (g_ascii_strncasecmp(name, pref_font, strlen(pref_font)) == 0
#if GLIB_MAJOR_VERSION == 2 && GLIB_MINOR_VERSION >= 6
	      && pango_font_family_is_monospace(families[i])
#endif
	      )
	    {
	      /* We prefer ones nearer the start of the list */
              if(current <= most_preferred)
                {
		  most_preferred = current;
		  match_family = families[i];

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
  /* Find the width of the sequence column */
  GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), MSP_COL);
  
  if (colWidth == UNSET_INT)
    colWidth = gtk_tree_view_column_get_width(sequenceCol);
  
  /* Find the character width of the font */
  PangoFontDescription *font_desc = pango_font_description_copy (tree->style->font_desc);
  PangoContext *context = gtk_widget_get_pango_context (tree);
  PangoFontMetrics *metrics = pango_context_get_metrics (context,
							 font_desc,
							 pango_context_get_language (context));
  pango_font_description_free (font_desc);
  
  int charWidth = pango_font_metrics_get_approximate_char_width(metrics) / PANGO_SCALE;
  pango_font_metrics_unref(metrics);
  
  /* Return the number of whole characters that fit in the column. */
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
	  gtk_adjustment_changed(adjustment);
	}
    }
}


/* Add all trees in the given list to the detail view */
static void addTreesToDetailView(GtkContainer *detailView, GList *treeList, gboolean first, gboolean last)
{
  if (GTK_IS_PANED(detailView))
    {
      int numTrees = g_list_length(treeList);
      if (numTrees > 0)
	{
	  GtkWidget *tree1 = GTK_WIDGET(treeList->data);
	  
	  if (first)
	    {
	      gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
	    }
	  else if (last)
	    {
	      gtk_paned_pack2(GTK_PANED(detailView), tree1, TRUE, TRUE);
	    }
	}
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


/* This function removes the trees from the detail view and re-adds them in the
 * correct order according to the strandsToggled flag. It should be called every
 * time the strands are toggled. It assumes the trees are already in the 
 * detailView container, and that the properties have been set for all 3 widgets. */
static void refreshTreeOrder(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  gboolean toggled = detailViewGetStrandsToggled(detailView);

  if (properties->seqType == BLXSEQ_DNA)
    {
      gtk_container_foreach(GTK_CONTAINER(detailView), removeFromDetailView, detailView);
      
      if (toggled)
	{
	  addTreesToDetailView(GTK_CONTAINER(detailView), properties->revStrandTrees, TRUE, FALSE);
	  addTreesToDetailView(GTK_CONTAINER(detailView), properties->fwdStrandTrees, FALSE, TRUE);
	}
      else
	{
	  addTreesToDetailView(GTK_CONTAINER(detailView), properties->fwdStrandTrees, TRUE, FALSE);
	  addTreesToDetailView(GTK_CONTAINER(detailView), properties->revStrandTrees, FALSE, TRUE);
	}
    }
  else if (properties->seqType == BLXSEQ_PEPTIDE)
    {
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
  
  /* Reset the display range so that it is between the scrollbar min and max. */
  int newStart = adjustment->value;
  int newEnd = adjustment->value + adjustment->page_size;
  
  /* If there is a gap at the end, shift the scrollbar back so that we show more
   * of the reference sequence (unless the whole thing is already shown). */
  int len = strlen(detailViewGetRefSeq(detailView));
  if (newEnd >= len)
    {
      newStart = len - adjustment->page_size;
      if (newStart < adjustment->lower)
	newStart = adjustment->lower;
      
      newEnd = newStart + adjustment->page_size;
      adjustment->value = newStart;
    }
  
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  if (displayRange->min != newStart + 1 || displayRange->max != newEnd + 1)
    {
      /* Adjustment is zero-based but display range is 1-based */
      displayRange->min = newStart + 1;
      displayRange->max = newEnd + 1;
      
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
    }
}


static void onScrollPosChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Set the display range so that it starts at the new scroll pos */
  int newEnd = adjustment->value + adjustment->page_size;

  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  if (displayRange->min != adjustment->value + 1 || displayRange->max != newEnd + 1)
    {
      displayRange->min = adjustment->value + 1;
      displayRange->max = newEnd + 1;

      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      callFuncOnAllDetailViewTrees(detailView, refilterTree);

      /* Scroll big picture if necessary to keep highlight box in view */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture);
      
      /* Update the highlight box position for all grids in the big picture */
      callFuncOnAllBigPictureGrids(bigPicture, calculateHighlightBoxBorders);
            
      /* Redraw all trees (and their corresponding grids) */
      callFuncOnAllDetailViewTrees(detailView, refreshTreeAndGrid);
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
  
  /* Call the increment-font-size function on all trees in the detail view */
  callFuncOnAllDetailViewTrees(detailView, treeIncrementFontSize);
  
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
  
  /* Call the decrement-font-size function on all trees in the detail view */
  callFuncOnAllDetailViewTrees(detailView, treeDecrementFontSize);
  
  if (firstTree && colWidth != UNSET_INT)
    {
      updateSeqColumnSize(firstTree, colWidth);
    }
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

GtkAdjustment* detailViewGetAdjustment(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->adjustment : NULL;
}

/* Get the forward strand of the reference sequence */
char* detailViewGetRefSeq(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties->refSeq;
}

static GtkWidget* detailViewGetMainWindow(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->mainWindow : NULL;
}

gboolean detailViewGetStrandsToggled(GtkWidget *detailView)
{
  return mainWindowGetStrandsToggled(detailViewGetMainWindow(detailView));
}

const char* detailViewGetFontFamily(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->fontFamily : NULL;
}


/* Extract the tree view for the given frame on the given strand */
GtkWidget* detailViewGetFrameTree(GtkWidget *detailView, gboolean forward, int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *result = NULL;
  
  /* Get the list of trees for the forward or reverse strand, as per the input argument */
  GList *list = forward ? properties->fwdStrandTrees : properties->revStrandTrees;

  /* Extract the tree for this given frame number. The list should be sorted in order of frame
   * number, and we should not be requesting a frame number greater than the number of items in the list */
  if (frame > g_list_length(list))
    {
      if (g_list_length(list) > 0)
	{
	  frame = 1;
	}
      else
	{
	  messcrash("Invalid frame number. Requested frame=%d on %s strand but number of detail view trees=%d.",
		    frame, (forward ? "forward" : "reverse"), g_list_length(list));
	}
    }

  GtkWidget *item = NULL;
  switch (frame)
  {
    case 1:
      item = list->data;
      break;
    case 2:
      item = list->next->data;
      break;
    case 3:
      item = list->next->next->data;
      break;
    default:
      messcrash("Invalid frame number. Requested frame=%d but max frame number=3", frame);
  }
  
  /* The item in the list is probably a container (i.e. a scrolled window), so extract the actual tree */
  if (GTK_IS_TREE_VIEW(item))
    {
      result = GTK_WIDGET(item);
    }
  else if (GTK_IS_CONTAINER(item))
    {
      GList *children = gtk_container_get_children(GTK_CONTAINER(item));
      if (children && g_list_length(children) == 1)
	{
	  result = GTK_WIDGET(children->data);
	}
    }

  if (!result)
    messcrash("Tree list for %s strand contains an invalid widget type for frame %d.", (forward ? "forward" : "reverse"), frame);
  
  return result;
}

/* Get the first tree in the 'current' list of trees (i.e. the forward strand
 * list by default, or the reverse strand list if strands are toggled). */
static GtkWidget* detailViewGetFirstTree(GtkWidget *detailView)
{
  gboolean forward = !detailViewGetStrandsToggled(detailView);
  return detailViewGetFrameTree(detailView, forward, 1);
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

//static BlxSeqType detailViewGetSeqType(GtkWidget *detailView)
//{
//  DetailViewProperties *properties = detailViewGetProperties(detailView);
//  return properties ? properties->seqType : UNSET_INT;
//}


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
				       GList *fwdStrandTrees,
				       GList *revStrandTrees,
				       GtkWidget *feedbackBox,
				       GtkAdjustment *adjustment, 
				       char *refSeq,
				       BlxSeqType seqType,
				       int numReadingFrames,
				       IntRange *displayRange,
				       const char *fontFamily)
{
  if (detailView)
    { 
      DetailViewProperties *properties = g_malloc(sizeof *properties);

      properties->mainWindow = mainWindow;
      properties->fwdStrandTrees = fwdStrandTrees;
      properties->revStrandTrees = revStrandTrees;
      properties->feedbackBox = feedbackBox;
      properties->adjustment = adjustment;
      properties->refSeq = refSeq;
      properties->seqType = seqType;
      properties->numReadingFrames = numReadingFrames;
      properties->displayRange.min = displayRange->min;
      properties->displayRange.max = displayRange->max;
      properties->selectedBaseIdx = UNSET_INT;
      properties->fontFamily = fontFamily;

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
  
  /* Redraw the grid header */
  gtk_widget_queue_draw(bigPictureGetGridHeader(mainWindowProperties->bigPicture));
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
static GtkWidget* createDetailViewButtonBar(GtkWidget *detailView, GtkWidget **feedbackBox)
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
  
  if (1) //blastx ||tblastx || blastn)
    {
      makeToolbarButton(toolbar, "Strand^v", "Toggle strand", (GtkSignalFunc)GToggleStrand, detailView);
    }
  
  *feedbackBox = createFeedbackBox(toolbar);
  
  return toolbarContainer;
}


/* Create two detail-view trees and place them in the 2 panes of the paned widget. The first
 * tree gets associated with grid1 and appended to list1, and the second with grid2 and list2. */
static void createTwoPanedTrees(GtkWidget *panedWidget, 
				GtkWidget *grid1, 
				GtkWidget *grid2,
				GList **list1,
				GList **list2,
				const char *fontFamily)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, panedWidget, list1, fontFamily);
  GtkWidget *tree2 = createDetailViewTree(grid2, panedWidget, list2, fontFamily);
  gtk_paned_pack1(GTK_PANED(panedWidget), tree1, FALSE, TRUE);
  gtk_paned_pack2(GTK_PANED(panedWidget), tree2, TRUE, TRUE);
}


/* Create the trees for the detail view, creating sub-panes if necessary depending
 * on how many trees we need */
static void createDetailViewPanes(GtkWidget *detailView, 
				  GtkWidget *fwdStrandGrid, 
				  GtkWidget *revStrandGrid, 
				  const int numReadingFrames,
				  GList **fwdStrandTrees,
				  GList **revStrandTrees,
				  const char *fontFamily)
{
  if (numReadingFrames == 1)
    {
      /* DNA matches: we need 2 trees, one for the forward strand and one for the reverse */
      createTwoPanedTrees(detailView, fwdStrandGrid, revStrandGrid, fwdStrandTrees, revStrandTrees, fontFamily);
    }
  else if (numReadingFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames. We only have
       * 2 panes in our parent container, so one pane will have to contain
       * a second, nested paned widget. */
      GtkWidget *tree1 = createDetailViewTree(fwdStrandGrid, detailView, fwdStrandTrees, fontFamily);
      GtkWidget *nestedPanedWidget = gtk_vpaned_new();
      
      gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
      gtk_paned_pack1(GTK_PANED(detailView), nestedPanedWidget, TRUE, TRUE);
      
      createTwoPanedTrees(nestedPanedWidget, fwdStrandGrid, fwdStrandGrid, fwdStrandTrees, fwdStrandTrees, fontFamily);
    }
}


/* Add the MSPs to the detail-view trees */
static void detailViewAddMspData(GtkWidget *detailView, MSP *mspList)
{
  /* First, create a data store for each tree so we have something to add our data to. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateBaseDataModel);

  /* Loop through each MSP, and add it to the correct tree based on its strand and reading frame */
  MSP *msp = mspList;
  for ( ; msp; msp = msp->next)
    {
      gboolean qForward = (msp->qframe[1] == '+');
      int frame = atoi(&msp->qframe[2]);
      
      GtkWidget *tree = detailViewGetFrameTree(detailView, qForward, frame);
      GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
      
      addMspToTreeModel(model, msp, tree);
    }
  
  /* Finally, create a custom-filtered version of the data store for each tree. We do 
   * this AFTER adding the data so that it doesn't try to re-filter every time we add a row. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateFilteredDataModel);
}


GtkWidget* createDetailView(GtkWidget *mainWindow,
			    GtkWidget *panedWidget,
			    GtkAdjustment *adjustment, 
			    GtkWidget *fwdStrandGrid, 
			    GtkWidget *revStrandGrid,
			    MSP *mspList,
			    char *refSeq,
			    BlxBlastMode mode,
			    BlxSeqType seqType,
			    int numReadingFrames,
			    IntRange *refSeqRange)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so we can identify them.) */
  GtkWidget *detailView = gtk_vpaned_new();
  
  /* We need a fixed-width font for displaying alignments. Find one from a 
   * list of possibilities. */
  GList *fixed_font_list = NULL ;

  fixed_font_list = g_list_append(fixed_font_list, "Monaco");
  fixed_font_list = g_list_append(fixed_font_list, "Courier");
  fixed_font_list = g_list_append(fixed_font_list, "Courier New");
  fixed_font_list = g_list_append(fixed_font_list, "Monospace");
  fixed_font_list = g_list_append(fixed_font_list, "fixed");

  const char *fontFamily = findFixedWidthFontFamily(detailView, fixed_font_list);
  g_list_free(fixed_font_list);
    
  /* Create the trees. */
  GList *fwdStrandTrees = NULL, *revStrandTrees = NULL;
  createDetailViewPanes(detailView, 
			fwdStrandGrid, 
			revStrandGrid, 
			numReadingFrames, 
			&fwdStrandTrees,
			&revStrandTrees,
			fontFamily);
  
  /* Create the toolbar. We need to remember the feedback box. */
  GtkWidget *feedbackBox = NULL;
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, &feedbackBox);
  
  /* Put everything in a vbox, because we need to have a toolbar at the top
   * of the detail view. */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), detailView, TRUE, TRUE, 0);
  
  /* Put the whole lot into the main window */
  gtk_paned_pack2(GTK_PANED(panedWidget), vbox, TRUE, TRUE);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(detailView), "size-allocate", G_CALLBACK(onSizeAllocateDetailView), NULL);
  
  /* Set the required properties */
  detailViewCreateProperties(detailView, 
			     mainWindow, 
			     fwdStrandTrees,
			     revStrandTrees,
			     feedbackBox, 
			     adjustment, 
			     refSeq, 
			     seqType,
			     numReadingFrames, 
			     refSeqRange,
			     fontFamily);

  /* Add the MSP's to the trees */
  detailViewAddMspData(detailView, mspList);
  
  return detailView;
}

