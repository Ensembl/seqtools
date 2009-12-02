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

#define MIN_FONT_SIZE			2
#define MAX_FONT_SIZE			20
#define FONT_INCREMENT_SIZE		2
#define DETAIL_VIEW_TOOLBAR_NAME	"DetailViewToolbarName"


/***********************************************************
 *		       Utility functions                   *
 ***********************************************************/


/* Increase the font size in the detail view trees (i.e. effectively zoom in) */
static void incrementFontSize(GtkWidget *tree, gpointer data)
{
  gint newSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE) + FONT_INCREMENT_SIZE;
  
  if (newSize <= MAX_FONT_SIZE)
    {
      GtkCellRenderer *renderer = treeGetRenderer(tree);
      setCellRendererFont(tree, renderer, FIXED_WIDTH_FONT, newSize);
    }
}


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
static void decrementFontSize(GtkWidget *tree, gpointer data)
{
  gint newSize = (pango_font_description_get_size(tree->style->font_desc) / PANGO_SCALE) - FONT_INCREMENT_SIZE;
  
  if (newSize >= MIN_FONT_SIZE)
    {
      GtkCellRenderer *renderer = treeGetRenderer(tree);
      setCellRendererFont(tree, renderer, FIXED_WIDTH_FONT, newSize);
    }
}

/* Update the scroll position of the given adjustment to the given value. Does
 * bounds checking. */
void setDetailViewScrollPos(GtkAdjustment *adjustment, int value)
{
  if (value < adjustment->lower)
    {
      value = adjustment->lower;
    }
  else
    {
      int maxValue = adjustment->upper - adjustment->page_size + 1;
      if (value > maxValue)
	value = maxValue;
      
      if (value < 0)
	value = 0; /* Can happen if page size is larger than upper value */
    }
  
  adjustment->value = value;
  gtk_adjustment_value_changed(adjustment);
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
  int newPageSize = calcNumBasesInSequenceColumn(tree, colWidth);
  
  /* Only perform an update if it's actually changed */
  if (newPageSize != adjustment->page_size)
    {
      adjustment->page_size = newPageSize;
      adjustment->page_increment = newPageSize;
      gtk_adjustment_changed(adjustment);
    }
}


/***********************************************************
 *                    Detail view events                   *
 ***********************************************************/

static void onSizeAllocateDetailView(GtkWidget *detailView, GtkAllocation *allocation, gpointer data)
{
  //  static int count = 0; ++count; printf("onSizeAllocateDetailView%d\n", count);
  
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  if (detailViewProperties->refTree)
    updateSeqColumnSize(detailViewProperties->refTree, UNSET_INT);
}



/***********************************************************
 *              Detail view scrollbar events               *
 ***********************************************************/

static void onScrollRangeChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *mainWindow = GTK_WIDGET(data);
  
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(mainWindow);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(mainWindowProperties->detailView);
  
  /* Reset the display range so that it is between the scrollbar min and max. */
  int newStart = adjustment->value;
  int newEnd = adjustment->value + adjustment->page_size;
  
  /* If there is a gap at the end, shift the scrollbar back so that we show more
   * of the reference sequence (unless the whole thing is already shown). */
  int len = strlen(detailViewProperties->refSeq);
  if (newEnd >= len)
    {
      newStart = len - adjustment->page_size;
      if (newStart < adjustment->lower)
	newStart = adjustment->lower;
      
      newEnd = newStart + adjustment->page_size;
      adjustment->value = newStart;
    }
  
  IntRange *displayRange = &detailViewProperties->displayRange;
  if (displayRange->min != newStart + 1 || displayRange->max != newEnd + 1)
    {
      /* Adjustment is zero-based but display range is 1-based */
      displayRange->min = newStart + 1;
      displayRange->max = newEnd + 1;
      
      /* Recalculate the borders for all the grids and the header in the big picture */
      BigPictureProperties *bigPictureProperties = bigPictureGetProperties(mainWindowProperties->bigPicture);
      calculateGridHeaderBorders(bigPictureProperties->header);
      callFuncOnAllBigPictureGrids(mainWindowProperties->bigPicture, calculateGridBorders);
      
      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      callFuncOnAllDetailViewTrees(mainWindowProperties->detailView, refilterTree);
      
      /* Redraw all trees (and their corresponding grids) */
      callFuncOnAllDetailViewTrees(mainWindowProperties->detailView, refreshTreeAndGrid);
    }
}


static void onScrollPosChangedDetailView(GtkObject *object, gpointer data)
{
  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *mainWindow = GTK_WIDGET(data);
  
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(mainWindow);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(mainWindowProperties->detailView);
  
  /* Set the display range so that it starts at the new scroll pos */
  int newEnd = adjustment->value + adjustment->page_size;

  IntRange *displayRange = &detailViewProperties->displayRange;
  if (displayRange->min != adjustment->value + 1 || displayRange->max != newEnd + 1)
    {
      displayRange->min = adjustment->value + 1;
      displayRange->max = newEnd + 1;
      
      /* Update the highlight box position for all grids in the big picture */
      callFuncOnAllBigPictureGrids(mainWindowProperties->bigPicture, calculateHighlightBoxBorders);
      
      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      callFuncOnAllDetailViewTrees(mainWindowProperties->detailView, refilterTree);
      
      /* Redraw all trees (and their corresponding grids) */
      callFuncOnAllDetailViewTrees(mainWindowProperties->detailView, refreshTreeAndGrid);
    }
}


/***********************************************************
 *            Detail View toolbar events                   *
 ***********************************************************/

static void onZoomInDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Remember the size of the sequence column. This gets reset to 0 when we change the font. */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(properties->refTree), MSP_COL);
  int colWidth = gtk_tree_view_column_get_width(sequenceCol);
  
  /* Call the increment-font-size function on all trees in the detail view */
  callFuncOnAllDetailViewTrees(detailView, incrementFontSize);
  
  if (properties->refTree)
    {
      updateSeqColumnSize(properties->refTree, colWidth);
    }
}

static void onZoomOutDetailView(GtkButton *button, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Remember the size of the sequence column. This gets reset to 0 when we change the font. */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkTreeViewColumn *sequenceCol = gtk_tree_view_get_column(GTK_TREE_VIEW(properties->refTree), MSP_COL);
  int colWidth = gtk_tree_view_column_get_width(sequenceCol);
  
  /* Call the decrement-font-size function on all trees in the detail view */
  callFuncOnAllDetailViewTrees(detailView, decrementFontSize);
  
  if (properties->refTree)
    {
      updateSeqColumnSize(properties->refTree, colWidth);
    }
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

char* detailViewGetRefSeq(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->refSeq : NULL;
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


static void detailViewCreateProperties(GtkWidget *widget,
				       GtkWidget *refTree, 
				       GtkAdjustment *adjustment, 
				       char *refSeq,
				       int numReadingFrames,
				       IntRange *displayRange)
{
  if (widget)
    { 
      DetailViewProperties *properties = g_malloc(sizeof *properties);
      
      properties->refTree = refTree;
      properties->adjustment = adjustment;
      properties->refSeq = refSeq;
      properties->numReadingFrames = numReadingFrames;
      properties->displayRange.min = displayRange->min;
      properties->displayRange.max = displayRange->max;
      properties->selectedBaseIdx = UNSET_INT;
      
      g_object_set_data(G_OBJECT(widget), "DetailViewProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyDetailView), NULL); 
    }
}


/***********************************************************
 *                     Callbacks                           *
 ***********************************************************/


static void blxHelp(void)
{
  //  graphMessage (messprintf("\
  //			   \
  //			   BLIXEM - BLast matches\n\
  //			   In an\n\
  //			   X-windows\n\
  //			   Embedded\n\
  //			   Multiple alignment\n\
  //			   \n\
  //			   LEFT MOUSE BUTTON:\n\
  //			   Pick on boxes and sequences.\n\
  //			   Fetch annotation by double clicking on sequence (Requires 'efetch' to be installed.)\n\
  //			   \n\
  //			   MIDDLE MOUSE BUTTON:\n\
  //			   Scroll horizontally.\n\
  //			   \n\
  //			   RIGHT MOUSE BUTTON:\n\
  //			   Menu.  Note that the buttons Settings and Goto have their own menus.\n\
  //			   \n\
  //			   \n\
  //			   Keyboard shortcuts:\n\
  //			   \n\
  //			   Cntl-Q          quit application\n\
  //			   Cntl-P          print\n\
  //			   Cntl-H          help\n\
  //			   \n\
  //			   \n\
  //			   m/M             for mark/unmark a set of matches from the cut buffer\n\
  //			   \n\
  //			   g/G             go to the match in the cut buffer\n\
  //			   \n\
  //			   \n\
  //			   \n\
  //			   RESIDUE COLOURS:\n\
  //			   Yellow = Query.\n\
  //			   See Settings Panel for matching residues (click on Settings button).\n\
  //			   \n\
  //			   version %s\n\
  //			   (c) Erik Sonnhammer", blixemVersion));
  //  
  //  return ;
}

static void blixemSettings(void)
{
  //    if (!graphActivate(settingsGraph)) {
  //	settingsGraph = graphCreate(TEXT_FIT, "Blixem Settings", 0, 0, .6, .3);
  //	graphRegister(PICK, settingsPick);
  //	settingsRedraw();
  //    }
  //    else graphPop();
}

static void ToggleStrand(void)
{
  //  dispstart += plusmin*(displen-1);
  //  
  //  plusmin = -plusmin;
  //  
  //  sprintf(actframe, "(%+d)", plusmin);
  //  blviewRedraw();
  //  
  //  return ;
}

static void scrollRightBig(void)
{
  //  dispstart = dispstart + plusmin*displen*.5;
  //  blviewRedraw();
}
static void scrollLeftBig(void)
{
  //  dispstart = dispstart - plusmin*displen*.5;
  //  blviewRedraw();
}
static void scrollRight1(void)
{
  //  dispstart = dispstart + plusmin*symbfact;
  //  blviewRedraw();
}
static void scrollLeft1(void)
{
  //  dispstart = dispstart - plusmin*symbfact;
  //  blviewRedraw();
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

static void prevMatch(void)
{
  //  gotoMatch(-1);
}

static void nextMatch(void)
{
  //  gotoMatch(1);
}


static void comboChange(GtkWidget *window, GtkEditable *edit, gpointer args)
{
  
  //  gchar *val = gtk_editable_get_chars(edit, 0, -1);
  //  
  //  if (GTK_WIDGET_REALIZED(window))
  //    {
  //      if (strcmp(val, "score") == 0)
  //	sortByScore();
  //      else if (strcmp(val, "identity") == 0)
  //	sortById();
  //      else if (strcmp(val, "name") == 0)
  //	sortByName();
  //      else if (strcmp(val, "position") == 0)
  //	sortByPos();
  //    }
  //  g_free(val);
}


static void GHelp(GtkButton *button, gpointer args)
{
  blxHelp();
}

static void GSettings(GtkButton *button, gpointer args)
{
  blixemSettings();
}

static void GGoto(GtkButton *button, gpointer args)
{
  Goto();
}

static void GprevMatch(GtkButton *button, gpointer args)
{
  prevMatch();
}

static void GnextMatch(GtkButton *button, gpointer args)
{
  nextMatch();
}

static void GscrollLeftBig(GtkButton *button, gpointer args)
{
  scrollLeftBig();
}

static void GscrollRightBig(GtkButton *button, gpointer args)
{
  scrollRightBig();
}

static void GscrollLeft1(GtkButton *button, gpointer args)
{
  scrollLeft1();
}

static void GscrollRight1(GtkButton *button, gpointer args)
{
  scrollRight1();
}

static void GToggleStrand(GtkButton *button, gpointer args)
{
  ToggleStrand();
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

GtkWidget* createDetailViewScrollBar(GtkAdjustment *adjustment, GtkWidget *mainWindow)
{
  GtkWidget *scrollBar = gtk_hscrollbar_new(adjustment);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollRangeChangedDetailView), mainWindow);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollPosChangedDetailView), mainWindow);
  
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
static GtkWidget* createEmptyButtonBar(GtkWidget *parent, GtkWidget **toolbar)
{
  /* Create a handle box for the toolbar and add it to the window */
  GtkWidget *handleBox = gtk_handle_box_new();
  
  /* Create the toolbar */
  *toolbar = gtk_toolbar_new();
  gtk_toolbar_set_tooltips(GTK_TOOLBAR(*toolbar), TRUE);
  gtk_toolbar_set_show_arrow(GTK_TOOLBAR(*toolbar), TRUE);
  gtk_toolbar_set_icon_size(GTK_TOOLBAR(*toolbar), GTK_ICON_SIZE_SMALL_TOOLBAR);

  /* Set the style property that controls the spacing */
  gtk_widget_set_name(*toolbar, DETAIL_VIEW_TOOLBAR_NAME);
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
  gtk_widget_set_usize(*toolbar, 1, -2);
  
  /* Add the toolbar to the handle box */
  gtk_container_add(GTK_CONTAINER(handleBox), *toolbar);
  
  return handleBox;
}


/* Create the combo box used for selecting sort criteria */
static GtkWidget* createSortBox()
{
  GtkWidget *combo = gtk_combo_new();

  gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(combo)->entry), FALSE);
  gtk_widget_set_usize(GTK_COMBO(combo)->entry, 80, -2);
  gtk_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", (GtkSignalFunc)comboChange, NULL);
  
  /* Create the list of strings the user can choose to sort by */
  GList *sortList = NULL;
  sortList = g_list_append(sortList, "score");
  sortList = g_list_append(sortList, "identity");
  sortList = g_list_append(sortList, "name");
  sortList = g_list_append(sortList, "position");
  gtk_combo_set_popdown_strings(GTK_COMBO(combo), sortList);
  
  /* Set the identity field as the default */
  gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), "identity");
  
  return combo;
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
static GtkWidget* createDetailViewButtonBar(GtkWidget *detailView)
{
  GtkWidget *toolbar = NULL;
  GtkWidget *toolbarContainer = createEmptyButtonBar(detailView, &toolbar);
  
  /* Zoom buttons */
  makeToolbarButton(GTK_TOOLBAR(toolbar), "+", "Zoom in", (GtkSignalFunc)onZoomInDetailView, detailView);
  makeToolbarButton(GTK_TOOLBAR(toolbar), "-", "Zoom out", (GtkSignalFunc)onZoomOutDetailView, detailView);
  
  /* Help button */
  makeToolbarButton(GTK_TOOLBAR(toolbar), "Help",	 "Don't Panic",			 (GtkSignalFunc)GHelp,		  NULL);

  /* Combo box for sorting (and a label to say what it is) */
  GtkWidget *label = gtk_label_new(" <i>Sort HSPs by:</i>");
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  GtkWidget *combo = createSortBox();
  addToolbarWidget(GTK_TOOLBAR(toolbar), label);
  addToolbarWidget(GTK_TOOLBAR(toolbar), combo);
  
  /* Settings button */
  makeToolbarButton(GTK_TOOLBAR(toolbar), "Settings", "Open the Preferences Window",  (GtkSignalFunc)GSettings,	  NULL);
  
  /* Navigation buttons */
  makeToolbarButton(GTK_TOOLBAR(toolbar), "Goto",	 "Go to specified co-ordinates", (GtkSignalFunc)GGoto,		  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), "< match",  "Next (leftward) match",	 (GtkSignalFunc)GprevMatch,	  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), "match >",  "Next (rightward) match",	 (GtkSignalFunc)GnextMatch,	  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), "<<",	 "Scroll leftward lots",	 (GtkSignalFunc)GscrollLeftBig,	  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), ">>",       "Scroll rightward lots",	 (GtkSignalFunc)GscrollRightBig,  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), "<",	 "Scroll leftward one base",	 (GtkSignalFunc)GscrollLeft1,	  NULL);
  makeToolbarButton(GTK_TOOLBAR(toolbar), ">",	 "Scroll rightward one base",	 (GtkSignalFunc)GscrollRight1,	  NULL);
  
  if (1) //blastx ||tblastx || blastn)
    {
      makeToolbarButton(GTK_TOOLBAR(toolbar), "Strand^v", "Toggle strand", (GtkSignalFunc)GToggleStrand, NULL);
    }
  
  /* Feedback box. (Feeds back info to the user about the currently-selected base/sequence. The 
   * user can copy text out of this box, but cannot edit the box's contents.) */
  GtkWidget *feedbackBox = gtk_entry_new() ;
  gtk_editable_set_editable(GTK_EDITABLE(feedbackBox), FALSE);
  gtk_widget_set_size_request(feedbackBox, 0, -1) ;
  GtkToolItem *item = addToolbarWidget(GTK_TOOLBAR(toolbar), feedbackBox) ;
  gtk_tool_item_set_expand(item, TRUE); /* make as big as possible */
  
  return toolbarContainer;
}


/* Create two detail-view trees and place them in the 2 panes of the paned widget. The first
 * tree gets associated with grid1 and the second with grid2. */
static void createTwoPanedTrees(GtkWidget *panedWidget, 
				GtkWidget *grid1, 
				GtkWidget *grid2,
				gboolean hasHeaders, 
				const MSP const *mspList)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, panedWidget, hasHeaders, mspList);
  GtkWidget *tree2 = createDetailViewTree(grid2, panedWidget, FALSE, mspList);
  gtk_paned_pack1(GTK_PANED(panedWidget), tree1, FALSE, TRUE);
  gtk_paned_pack2(GTK_PANED(panedWidget), tree2, TRUE, TRUE);
}


/* Create the trees for the detail view, creating sub-panes if necessary depending
 * on how many trees we need */
static void createDetailViewPanes(GtkWidget *detailView, 
				  GtkWidget *fwdStrandGrid, 
				  GtkWidget *revStrandGrid, 
				  const MSP const *mspList,
				  const int numReadingFrames)
{
  if (numReadingFrames == 1)
    {
      /* DNA matches: we need 2 trees, one for the forward strand and one for the reverse */
      createTwoPanedTrees(detailView, fwdStrandGrid, revStrandGrid, TRUE, mspList);
    }
  else if (numReadingFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames. We only have
       * 2 panes in our parent container, so one pane will have to contain
       * a second, nested paned widget. */
      GtkWidget *tree1 = createDetailViewTree(fwdStrandGrid, detailView, TRUE, mspList);
      GtkWidget *nestedPanedWidget = gtk_vpaned_new();
      
      gtk_paned_pack1(GTK_PANED(detailView), tree1, FALSE, TRUE);
      gtk_paned_pack1(GTK_PANED(detailView), nestedPanedWidget, TRUE, TRUE);
      
      createTwoPanedTrees(nestedPanedWidget, fwdStrandGrid, fwdStrandGrid, FALSE, mspList);
    }
}


GtkWidget* createDetailView(GtkWidget *container,
			    GtkAdjustment *adjustment, 
			    const MSP const *mspList,
			    GtkWidget *fwdStrandGrid, 
			    GtkWidget *revStrandGrid,
			    char *refSeq,
			    int numReadingFrames,
			    IntRange *refSeqRange)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so we can identify them.) */
  GtkWidget *detailView = gtk_vpaned_new();
  detailViewCreateProperties(detailView, NULL, adjustment, refSeq, numReadingFrames, refSeqRange);
  
  /* Create the trees */
  createDetailViewPanes(detailView, fwdStrandGrid, revStrandGrid, mspList, numReadingFrames);
  
  /* Create the toolbar */
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView);
  
  /* Put everything in a vbox, because we need to have a toolbar at the top
   * of the detail view. */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), detailView, TRUE, TRUE, 0);
  
  /* Put the vbox in the container we were passed */
  gtk_paned_pack2(GTK_PANED(container), vbox, TRUE, TRUE);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(detailView), "size-allocate", G_CALLBACK(onSizeAllocateDetailView), NULL);
  
  return detailView;
}

