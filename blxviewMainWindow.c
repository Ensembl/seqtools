/*
 *  blxviewMainWindow.c
 *  acedb
 *
 *  Created by Gemma Barson on 24/11/2009.
 *
 */

#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/blxdotter.h>
#include <SeqTools/utilities.h>
#include <gdk/gdkkeysyms.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      10   /* used to change the default border width around the main window */
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2   /* used to start with a smaller font than the default widget font */
#define DEFAULT_SCROLL_STEP_INCREMENT	 5    /* how many bases the scrollbar scrolls by for each increment */
#define DEFAULT_WINDOW_WIDTH_FRACTION	 0.9  /* what fraction of the screen size the blixem window width defaults to */
#define DEFAULT_WINDOW_HEIGHT_FRACTION	 0.8  /* what fraction of the screen size the blixem window height defaults to */

typedef struct _GroupDialogData
  {
    GtkWidget *mainWindow;
    GtkWidget *searchTextBox;
    GtkWidget *searchToggleButton;
  } GroupDialogData;

typedef struct _CompareSeqData
  {
    const char *searchStr;
    GList *matchList;
  } CompareSeqData;


/* Local function declarations */
static void			  onHelpMenu(GtkAction *action, gpointer data);
static void			  onPrintMenu(GtkAction *action, gpointer data);
static void			  onSettingsMenu(GtkAction *action, gpointer data);
static void			  onViewMenu(GtkAction *action, gpointer data);
static void			  onCreateGroupMenu(GtkAction *action, gpointer data);
static void			  onEditGroupsMenu(GtkAction *action, gpointer data);
static void			  onDotterMenu(GtkAction *action, gpointer data);
static void			  onSelectFeaturesMenu(GtkAction *action, gpointer data);
static void			  onDeselectAllRows(GtkAction *action, gpointer data);
static void			  onStatisticsMenu(GtkAction *action, gpointer data);

//static void			  onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
//static void			  onDrawPage(GtkPrintOperation *operation, GtkPrintContext *context, gint pageNum, gpointer data);

static Strand			  mainWindowGetActiveStrand(GtkWidget *mainWindow);
static Strand			  mainWindowGetInactiveStrand(GtkWidget *mainWindow);

static GList*			  findSelectedSeqInList(GList *list, const char *seqName);
static gint			  runConfirmationBox(GtkWidget *mainWindow, char *title, char *messageText);
static void			  onButtonClickedDeleteGroup(GtkWidget *button, gpointer data);
static void			  copySelectionToClipboard(GtkWidget *mainWindow);


/* Menu builders */
static const GtkActionEntry mainMenuEntries[] = {
  { "Quit",		NULL, "_Quit\t\t\t\tCtrl-Q",	      "<control>Q",	    "Quit the program",		  gtk_main_quit},
  { "Help",		NULL, "_Help\t\t\t\tCtrl-H",	      "<control>H",	    "Display Blixem help",	  G_CALLBACK(onHelpMenu)},
  { "Print",		NULL, "_Print\t\t\t\tCtrl-P",	      "<control>P",	    "Print",			  G_CALLBACK(onPrintMenu)},
  { "Settings",		NULL, "_Settings\t\t\t\tCtrl-S",      "<control>S",	    "Edit Blixem settings",	  G_CALLBACK(onSettingsMenu)},

  { "View",		NULL, "_View\t\t\t\tCtrl-V",		  "<control>V",	    "Edit view settings",	  G_CALLBACK(onViewMenu)},
  { "CreateGroup",	NULL, "Create Group\t\tShift-Ctrl-G",	  "<shift><control>G",  "Group sequences together",	  G_CALLBACK(onCreateGroupMenu)},
  { "EditGroups",	NULL, "Edit _Groups\t\t\tCtrl-G",	  "<control>G",	 "Groups",			  G_CALLBACK(onEditGroupsMenu)},
  { "DeselectAllRows",	NULL, "Deselect _all\t\t\tShift-Ctrl-A",  "<shift><control>A", "Deselect all",		  G_CALLBACK(onDeselectAllRows)},

  { "Dotter",		NULL, "_Dotter\t\t\t\tCtrl-D",		      "<control>D",	 "Start Dotter",	  G_CALLBACK(onDotterMenu)},
  { "SelectFeatures",	NULL, "Feature series selection tool",	NULL, "Feature series selection tool",		  G_CALLBACK(onSelectFeaturesMenu)},

  { "Statistics",	NULL, "Statistics",   NULL,		      "Show memory statistics",			  G_CALLBACK(onStatisticsMenu)}
};


/* This defines the layout of the menu for a standard user */
static const char *standardMenuDescription =
"<ui>"
"  <popup name='MainMenu'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
"      <menuitem action='Print'/>"
"      <menuitem action='Settings'/>"
"      <separator/>"
"      <menuitem action='View'/>"
"      <menuitem action='CreateGroup'/>"
"      <menuitem action='EditGroups'/>"
"      <menuitem action='DeselectAllRows'/>"
"      <separator/>"
"      <menuitem action='Dotter'/>"
"      <menuitem action='SelectFeatures'/>"
"  </popup>"
"</ui>";


/* This defines the layout of the menu for a developer user */
static const char *developerMenuDescription =
"<ui>"
"  <popup name='MainMenu'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
"      <menuitem action='Print'/>"
"      <menuitem action='Settings'/>"
"      <separator/>"
"      <menuitem action='View'/>"
"      <menuitem action='CreateGroup'/>"
"      <menuitem action='EditGroups'/>"
"      <menuitem action='DeselectAllRows'/>"
"      <separator/>"
"      <menuitem action='Dotter'/>"
"      <menuitem action='SelectFeatures'/>"
"      <separator/>"
"      <menuitem action='Statistics'/>"
"  </popup>"
"</ui>";



/***********************************************************
 *			   Utilities			   *
 ***********************************************************/

/* Return true if the current user is in our list of developers. */
static gboolean userIsDeveloper()
{
  gchar* developers[] = {"edgrif", "gb10"};

  gboolean result = FALSE;
  const gchar *user = g_get_user_name();
  int numDevelopers = sizeof(developers) / sizeof(gchar*);
  
  int i = 0;
  for (i = 0; i < numDevelopers; ++i)
    {
      if (strcmp(user, developers[i]) == 0)
        {
          result = TRUE;
          break;
        }
    }

  return result;
}


/* Copy a segment of the given sequence into a new string. The result must be
 * free'd with g_free by the caller. The given indices must be 0-based. */
static gchar *copySeqSegment(const char const *inputSeq, const int idx1, const int idx2)
{
  const int minIdx = min(idx1, idx2);
  const int maxIdx = max(idx1, idx2);
  
  const int segmentLen = maxIdx - minIdx + 1;
  gchar *segment = g_malloc(sizeof(gchar) * segmentLen + 1);
  
  strncpy(segment, inputSeq + minIdx, segmentLen);
  segment[segmentLen] = '\0';
  
  return segment;
}


/* Get a segment of the given sequence (which is always the DNA sequence and
 * always the forward strand) and reverse/complement it as required for the given
 * strand and reading frame (note that the sequence will only actually be reversed
 * if the 'reverse' argument is true). This function also translates it to a peptide
 * sequence if relevant (if the 'translate' flag allows it). */
gchar *getSequenceSegment(GtkWidget *mainWindow,
			  const char const *dnaSequence,
			  const IntRange const *dnaSequenceRange,
			  const int coord1, 
			  const int coord2,
			  const Strand strand,
			  const BlxSeqType inputCoordType,
			  const int frame,
			  const int numFrames,
			  const gboolean rightToLeft,
			  const gboolean reverseResult,
			  const gboolean translateResult)
{
  gchar *result = NULL;

  /* Convert input coord to ref seq coords and find the min/max */
  const int qIdx1 = convertDisplayIdxToDnaIdx(coord1, inputCoordType, frame, 1, numFrames, rightToLeft, dnaSequenceRange);	 /* 1st base in frame */
  const int qIdx2 = convertDisplayIdxToDnaIdx(coord2, inputCoordType, frame, numFrames, numFrames, rightToLeft, dnaSequenceRange); /* last base in frame */
  int qMin = min(qIdx1, qIdx2);
  int qMax = max(qIdx1, qIdx2);
  
  /* Check that the requested segment is within the sequence's range */
  if (qMin < dnaSequenceRange->min || qMax > dnaSequenceRange->max)
    {
      if (inputCoordType == BLXSEQ_PEPTIDE)
	printf ( "Requested query sequence %d - %d out of available range: %d - %d. Input coords on peptide sequence were %d - %d\n", qMin, qMax, dnaSequenceRange->min, dnaSequenceRange->max, coord1, coord2);
      else
	printf ( "Requested query sequence %d - %d out of available range: %d - %d\n", qMin, qMax, dnaSequenceRange->min, dnaSequenceRange->max);
	
      if (qMax > dnaSequenceRange->max)
	qMax = dnaSequenceRange->max;
	
      if (qMin < dnaSequenceRange->min)
	qMin = dnaSequenceRange->min;
    }
  
  /* Get 0-based indices into the sequence */
  const int idx1 = qMin - dnaSequenceRange->min;
  const int idx2 = qMax - dnaSequenceRange->min;
  
  /* Copy the portion of interest from the reference sequence and translate as necessary */
  const BlxBlastMode mode = mainWindowGetBlastMode(mainWindow);
  
  if (mode == BLXMODE_BLASTP || mode == BLXMODE_TBLASTN)
    {
      /* Just get a straight copy of this segment from the ref seq. Must pass 0-based
       * indices into the sequence */
      result = copySeqSegment(dnaSequence, idx1, idx2);
      
      if (reverseResult)
	{
	  g_strreverse(result);
	}
    }
  else
    {
      /* Get the segment of the ref seq, adjusted as necessary for this strand */
      gchar *segment = NULL;
      
      if (strand == FORWARD_STRAND)
	{
	  /* Straight copy of the ref seq segment */
	  segment = copySeqSegment(dnaSequence, idx1, idx2);
	}
      else
	{
	  /* Get the segment of the ref seq and then complement it */
	  segment = copySeqSegment(dnaSequence, idx1, idx2);
	  blxComplement(segment);
	  
	  if (!segment)
	    {
	      messcrash ("Error getting the reference sequence segment: Failed to complement the reference sequence for the range %d - %d.", qMin, qMax);
	    }
	}
      
      if (reverseResult)
	{
	  g_strreverse(segment);
	}
      
      if (mode == BLXMODE_BLASTN || !translateResult)
	{
	  /* Just return the segment of DNA sequence */
	  result = segment;
	}
      else
	{
	  /* Translate the DNA sequence to a peptide sequence */
	  result = blxTranslate(segment, mainWindowGetGeneticCode(mainWindow));
	  
	  g_free(segment); /* delete this because we're not returning it */
	  segment = NULL;
	  
	  if (!result)
	    {
	      messcrash ("Error getting the reference sequence segment: Failed to translate the DNA sequence for the reference sequence range %d - %d/\n", qMin, qMax) ;
	    }
	}
    }
  
  return result;
}


/* Returns the order number of the group that this sequence belongs to, or
 * UNSET_INT if it does not belong to a group. */
int sequenceGetGroupOrder(GtkWidget *mainWindow, const char *seqName)
{
  SequenceGroup *group = mainWindowGetSequenceGroup(mainWindow, seqName);
  return group ? group->order : UNSET_INT;
}


/* Scroll the detail view left/right by 1 base */
static void scrollDetailViewBy1(GtkWidget *window, const gboolean moveLeft)
{
  GtkWidget *detailView = mainWindowGetDetailView(window);
  
  if (moveLeft)
    {
      scrollDetailViewLeft1(detailView);
    }
  else
    {
      scrollDetailViewRight1(detailView);
    }
}


/* Move the selected base index 1 base to the left/right. Moves by individual
 * DNA bases (i.e. you have to move 3 bases in order to scroll a full peptide
 * if viewing protein matches). Scrolls the detail view if necessary to keep 
 * the new base in view. */
static void moveSelectedBaseIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  MainWindowProperties *properties = mainWindowGetProperties(window);
  GtkWidget *detailView = properties->detailView;
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  if (detailViewProperties->selectedBaseIdx != UNSET_INT)
    {
      /* Decrement the index if moving left decrease and increment it if moving right, 
       * unless the display is toggled, in which case do the opposite */
      const int numFrames = properties->numReadingFrames;
      const int minBaseNum = 1;
      const int maxBaseNum = numFrames;

      int newBaseNum = detailViewProperties->selectedBaseNum;
      int newSelectedBaseIdx = detailViewProperties->selectedBaseIdx;
      
      if (moveLeft)
	{
	  --newBaseNum;
	  
	  if (newBaseNum < minBaseNum)
	    {
	      newBaseNum = maxBaseNum;
	      --newSelectedBaseIdx;
	    }
	}
      else
	{
	  ++newBaseNum;
	  
	  if (newBaseNum > maxBaseNum)
	    {
	      newBaseNum = minBaseNum;
	      ++newSelectedBaseIdx;
	    }
	}

      IntRange *fullRange = mainWindowGetFullRange(window);
      boundsLimitValue(&newSelectedBaseIdx, fullRange);

      detailViewSetSelectedBaseIdx(detailView, newSelectedBaseIdx, detailViewProperties->selectedFrame, newBaseNum, TRUE);
    }
}


/* Move the selected display index 1 value to the left/right. Moves by full peptides
 * if viewing protein matches. Scrolls the detail view if necessary to keep the new 
 * index in view. */
static void moveSelectedDisplayIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  MainWindowProperties *properties = mainWindowGetProperties(window);
  GtkWidget *detailView = properties->detailView;
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  if (detailViewProperties->selectedBaseIdx != UNSET_INT)
    {
      /* Decrement the index if moving left decrease and increment it if moving right, 
       * unless the display is toggled, in which case do the opposite */
      int newSelectedBaseIdx = detailViewProperties->selectedBaseIdx;
      
      if (moveLeft)
	{
	  --newSelectedBaseIdx;
	}
      else
	{
	  ++newSelectedBaseIdx;
	}
      
      IntRange *fullRange = mainWindowGetFullRange(window);
      boundsLimitValue(&newSelectedBaseIdx, fullRange);
      
      detailViewSetSelectedBaseIdx(detailView, newSelectedBaseIdx, detailViewProperties->selectedFrame, detailViewProperties->selectedBaseNum, TRUE);
    }
}


/* Zooms the display in/out. if zoomOverview is true it zooms the big-picture section of
 * the display; otherwise it zooms the detail-view section. */
static void zoomMainWindow(GtkWidget *window, const gboolean zoomIn, const gboolean zoomOverview)
{
  if (zoomOverview)
    {
      zoomBigPicture(mainWindowGetBigPicture(window), zoomIn);
    }
  else
    {
      zoomDetailView(mainWindowGetDetailView(window), zoomIn);
    }
}


/* Force a redraw of all widgets. Clears cached bitmaps etc. first */
void mainWindowRedrawAll(GtkWidget *mainWindow)
{
  GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
  widgetClearCachedDrawable(bigPictureGetFwdGrid(bigPicture));
  widgetClearCachedDrawable(bigPictureGetRevGrid(bigPicture));
  
  gtk_widget_queue_draw(mainWindow);
}


/* Utility to create a vbox with the given border and pack it into the given container */
static GtkWidget* createVBoxWithBorder(GtkWidget *parent, const int borderWidth)
{
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), borderWidth);
  gtk_container_add(GTK_CONTAINER(parent), vbox);
  return vbox;
}


/***********************************************************
 *			   View panes menu                 *
 ***********************************************************/

/* Toggle visibility the n'th tree. This is the active strand's frame n if displaying
 * protein matches (where we only display one strand), or the forward or reverse
 * strand tree if displaying DNA matches (where both strands are displayed). */
static void toggleTreeVisibility(GtkWidget *mainWindow, const int number)
{
  const gboolean toggled = mainWindowGetStrandsToggled(mainWindow);
  const Strand activeStrand = toggled ? REVERSE_STRAND : FORWARD_STRAND;
  
  /* For protein matches, trees are always displayed in frame order (i.e. 1, 2, 3), 
   * so just use the number pressed for the frame, and the active strand for the
   * strand. */
  int frame = number;
  Strand strand = activeStrand;
  
  /* For DNA matches, the frame is always 1, but the strand depends on which number
   * was pressed: use 1 to toggle active strand, 2 for other strand */
  if (mainWindowGetSeqType(mainWindow) == BLXSEQ_DNA)
    {
      frame = 1;
      
      if (number == 1)
	{
	  strand = activeStrand;
	}
      else if (number == 2)
	{
	  strand = toggled ? FORWARD_STRAND : REVERSE_STRAND;
	}
    }
  
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (tree && gtk_widget_get_parent(tree))
    {
      widgetSetHidden(tree, !widgetGetHidden(tree));
    }
}


/* Toggle visibility of the active (1) or other (2) strand grid depending on the number pressed */
static void toggleGridVisibility(GtkWidget *mainWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
      const gboolean useFwdGrid = (number == 1) != mainWindowGetStrandsToggled(mainWindow);
      
      GtkWidget *grid = useFwdGrid ? bigPictureGetFwdGrid(bigPicture) : bigPictureGetRevGrid(bigPicture);
      widgetSetHidden(grid, !widgetGetHidden(grid));
    }
}


/* Toggle visibility of the active (1) or other (2) strand exon view depending on the number pressed */
static void toggleExonViewVisibility(GtkWidget *mainWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
      const gboolean useFwdExonView = (number == 1) != mainWindowGetStrandsToggled(mainWindow);
      
      GtkWidget *exonView = useFwdExonView ? bigPictureGetFwdExonView(bigPicture) : bigPictureGetRevExonView(bigPicture);
      widgetSetHidden(exonView, !widgetGetHidden(exonView));
    }
}


/* Toggle the visibility of tree/grid panes following a number key press */
static void togglePaneVisibility(GtkWidget *mainWindow, const int number, const gboolean modifier1, const gboolean modifier2)
{
  /* Affects big picture if modifier1 was pressed, the detail view otherwise */
  if (modifier1)
    {
      /* If modifier 2 was also pressed, affects the exon views; otherwise the grids */
      if (modifier2)
	{
	  toggleExonViewVisibility(mainWindow, number);
	}
      else
	{
	  toggleGridVisibility(mainWindow, number);
	}
    }
  else
    {
      toggleTreeVisibility(mainWindow, number);
    }
}


/* Called when the state of a check button is toggled */
static void onVisibilityButtonToggled(GtkWidget *button, gpointer data)
{
  GtkWidget *widgetToToggle = GTK_WIDGET(data);
  gboolean visible = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  widgetSetHidden(widgetToToggle, !visible);
}


/* Create a check button to control visibility of the given widget */
static void createVisibilityButton(GtkWidget *widgetToToggle, const char *mnemonic, GtkWidget *container)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_container_add(GTK_CONTAINER(container), button);

  /* Set the state depending on the widget's current visibility */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), GTK_WIDGET_VISIBLE(widgetToToggle));
  
  g_signal_connect(G_OBJECT(button), "toggled", G_CALLBACK(onVisibilityButtonToggled), widgetToToggle);
}


/* Create a check button to control visibility of the given tree. */
static void createTreeVisibilityButton(GtkWidget *detailView, const Strand strand, const int frame, GtkWidget *container)
{
  /* Some trees may have been removed from the main window if they are not on the active 
   * strand, so only show check boxes for those that are in the window (i.e. have a parent). 
   * Note that we get the tree container here, which consists of the tree itself plus any headers etc. */
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (gtk_widget_get_parent(tree))
    {
      if (detailViewGetSeqType(detailView) == BLXSEQ_DNA)
	{
	  /* We only have 1 frame, but trees are from both strands, so distinguish between strands */
	  const gboolean toggled = detailViewGetStrandsToggled(detailView);
	  gboolean isActiveStrand = ((strand == FORWARD_STRAND) != toggled);
	  
	  char text1[] = "Act_ive strand alignments";
	  char text2[] = "Othe_r strand alignments";
	  createVisibilityButton(tree, isActiveStrand ? text1 : text2, container);
	}
      else
	{
	  /* All the visible trees should be in the same strand, so just distinguish by frame number */
	  char formatStr[] = "Alignment list _%d";
	  char displayText[strlen(formatStr) + numDigitsInInt(frame) + 1];
	  sprintf(displayText, formatStr, frame);

	  createVisibilityButton(tree, displayText, container);
	}
    }
}


/* Shows the "View panes" dialog. This dialog allows the user to show/hide certain portions of the window. */
void showViewPanesDialog(GtkWidget *mainWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("View panes", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  
  /* Big picture */
  GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
  createVisibilityButton(bigPicture, "_Big picture", contentArea);
  
  GtkWidget *bigPictureSubBox = createVBoxWithBorder(contentArea, borderWidth);
  createVisibilityButton(bigPictureGetActiveGrid(bigPicture), "_Active strand grid", bigPictureSubBox);
  createVisibilityButton(bigPictureGetActiveExonView(bigPicture), "Active strand _exons", bigPictureSubBox);
  createVisibilityButton(bigPictureGetInactiveExonView(bigPicture), "Other strand e_xons", bigPictureSubBox);
  createVisibilityButton(bigPictureGetInactiveGrid(bigPicture), "O_ther strand grid", bigPictureSubBox);
  
  /* Detail view */
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  createVisibilityButton(detailView, "_Detail view", contentArea);
  
  GtkWidget *detailViewSubBox = createVBoxWithBorder(contentArea, borderWidth);
  const int numFrames = mainWindowGetNumReadingFrames(mainWindow);
  int frame = 1;
  for ( ; frame <= numFrames; ++frame)
    {
      createTreeVisibilityButton(detailView, mainWindowGetActiveStrand(mainWindow), frame, detailViewSubBox);
      createTreeVisibilityButton(detailView, mainWindowGetInactiveStrand(mainWindow), frame, detailViewSubBox);
    }
  
  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
}


/***********************************************************
 *		      Group sequences menu                 *
 ***********************************************************/

/* Free the memory used by the given sequence group and its members. */
static void destroySequenceGroup(SequenceGroup *seqGroup)
{
  if (seqGroup->groupName)
    g_free(seqGroup->groupName);
  
  if (seqGroup->seqList)
    g_list_free(seqGroup->seqList);
  
  g_free(seqGroup);
}


/* Delete a single group */
static void mainWindowDeleteSequenceGroup(GtkWidget *mainWindow, SequenceGroup *group)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  
  if (properties->sequenceGroups)
    {
      properties->sequenceGroups = g_list_remove(properties->sequenceGroups, group);
    }	 
  
  /* Refilter the trees (because hidden rows may now be visible again), and redraw */
  callFuncOnAllDetailViewTrees(mainWindowGetDetailView(mainWindow), refilterTree);
  mainWindowRedrawAll(mainWindow);
}


/* Delete all groups */
static void mainWindowDeleteAllSequenceGroups(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  GList *groupItem = properties->sequenceGroups;
  
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      destroySequenceGroup(group);
    }
  
  g_list_free(properties->sequenceGroups);
  properties->sequenceGroups = NULL;
  
  /* Refilter the trees (because hidden rows may now be visible again), and redraw */
  callFuncOnAllDetailViewTrees(mainWindowGetDetailView(mainWindow), refilterTree);
  mainWindowRedrawAll(mainWindow);
}


/* Create a new, empty sequence group with a unique ID and name (unique from all
 * others in the given list - does not actually add it to the list though). The result 
 * should be destroyed with destroySequenceGroup. */
static SequenceGroup* createSequenceGroup(GList *groupList)
{
  /* Create the new group */
  SequenceGroup *seqGroup = g_malloc(sizeof(SequenceGroup));
  
  seqGroup->seqList = NULL;
  seqGroup->hidden = FALSE;
  
  /* Find a unique ID */
  GList *lastItem = g_list_last(groupList);
  if (lastItem)
    {
      SequenceGroup *lastGroup = (SequenceGroup*)(lastItem->data);
      seqGroup->groupId = lastGroup->groupId + 1;
    }
  else
    {
      seqGroup->groupId = 1;
    }

  /* Create a default name based on the unique ID */
  char formatStr[] = "Group%d";
  const int nameLen = strlen(formatStr) + numDigitsInInt(seqGroup->groupId);
  seqGroup->groupName = g_malloc(nameLen * sizeof(*seqGroup->groupName));
  sprintf(seqGroup->groupName, formatStr, seqGroup->groupId);

  /* Set the order number. For simplicity, set the default order to be the same
   * as the ID number, so groups are sorted in the order they were added */
  seqGroup->order = seqGroup->groupId;

  /* Set the default highlight colour. */
  seqGroup->highlighted = TRUE;
  seqGroup->highlightColour = getGdkColor(GDK_RED);

  return seqGroup;
}


/* Make a group from the currently selection. Does nothing except give a warning
 * if no sequences are currently selected. */
static void makeGroupFromSelection(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  
  if (g_list_length(properties->selectedSeqs) > 0)
    {
      /* Create the group */
      SequenceGroup *seqGroup = createSequenceGroup(properties->sequenceGroups);
      seqGroup->seqList = g_list_copy(properties->selectedSeqs);
      
      /* Add it to the list of groups */
      properties->sequenceGroups = g_list_append(properties->sequenceGroups, seqGroup);
    }
  else
    {
      messout("Warning: cannot create group; no sequences are currently selected");
    }
}


/* This function is called when the sequence-group-name text entry widget's
 * value has changed. It sets the new group name in the group. */
static gboolean onGroupNameChanged(GtkWidget *widget, GdkEventFocus *event, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newName = gtk_entry_get_text(entry);
  
  if (!newName || strlen(newName) < 1)
    {
      messout("Invalid group name '%s' entered; reverting to previous group name '%s'.", newName, group->groupName);
      gtk_entry_set_text(entry, group->groupName);
    }
  else
    {
      if (group->groupName) 
	g_free(group->groupName);
      
      group->groupName = g_strdup(newName);
    }
  
  return FALSE;
}


/* This function is called when the sequence-group-order text entry widget's
 * value has changed. It sets the new order number in the group. */
static gboolean onGroupOrderChanged(GtkWidget *widget, GdkEventFocus *event, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newOrder = gtk_entry_get_text(entry);
  
  if (!newOrder || strlen(newOrder) < 1)
    {
      messout("Invalid order number '%s' entered; reverting to previous order number '%d'.", newOrder, group->order);
      char *orderStr = convertIntToString(group->order);
      gtk_entry_set_text(entry, orderStr);
      g_free(orderStr);
    }
  else
    {
      group->order = convertStringToInt(newOrder);
      
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(widget));
      GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
      
      /* Re-sort the trees and re-draw */
      callFuncOnAllDetailViewTrees(mainWindowGetDetailView(mainWindow), resortTree);
      mainWindowRedrawAll(mainWindow);
    }
  
  return FALSE;
}


/* This callback is called when the toggle button for a group's "hidden" flag is toggled.
 * It updates the group's hidden flag according to the button's new status. */
static void onGroupHiddenToggled(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  group->hidden = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

  /* Refilter trees and redraw all immediately show the new status */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  callFuncOnAllDetailViewTrees(mainWindowGetDetailView(mainWindow), refilterTree);
  mainWindowRedrawAll(mainWindow);
}


/* This callback is called when the toggle button for a group's "highlighted" flag is toggled.
 * It updates the group's highlighted flag according to the button's new status. */
static void onGroupHighlightedToggled(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  group->highlighted = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  /* Redraw the main window to immediately show the new status. The toplevel
   * parent of the button is the dialog, and the main window is the transient
   * parent of the dialog. */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  mainWindowRedrawAll(mainWindow);
}


/* Called when the user has changed the colour of a group in the 'edit groups' dialog */
static void onGroupColourChanged(GtkColorButton *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  gtk_color_button_get_color(button, &group->highlightColour);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), &group->highlightColour, TRUE, TRUE);
  
  /* Redraw everything in the new colours */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(button)));
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  mainWindowRedrawAll(mainWindow);
}


/* This function creates a widget that allows the user to edit the group
 * pointed to by the given list item, and adds it to the table container
 * widget at the given row. */
static void createEditGroupWidget(SequenceGroup *group, GtkTable *table, const int row, const int xpad, const int ypad)
{
  /* Show the group's name in a text box that the user can edit */
  GtkWidget *nameWidget = gtk_entry_new();
  gtk_entry_set_text(GTK_ENTRY(nameWidget), group->groupName);
  g_signal_connect(G_OBJECT(nameWidget), "focus-out-event", G_CALLBACK(onGroupNameChanged), group);
  
  /* Add a check box for the 'hidden' flag */
  GtkWidget *isHiddenWidget = gtk_check_button_new();
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHiddenWidget), group->hidden);
  g_signal_connect(G_OBJECT(isHiddenWidget), "toggled", G_CALLBACK(onGroupHiddenToggled), group);

  /* Add a check box for the 'highlighted' flag */
  GtkWidget *isHighlightedWidget = gtk_check_button_new();
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHighlightedWidget), group->highlighted);
  g_signal_connect(G_OBJECT(isHighlightedWidget), "toggled", G_CALLBACK(onGroupHighlightedToggled), group);
  
  /* Show the group's order number in an editable text box */
  GtkWidget *orderWidget = gtk_entry_new();
  char *orderStr = convertIntToString(group->order);
  gtk_entry_set_text(GTK_ENTRY(orderWidget), orderStr);
  g_free(orderStr);
  gtk_widget_set_size_request(orderWidget, 20, -1);
  g_signal_connect(G_OBJECT(orderWidget), "focus-out-event", G_CALLBACK(onGroupOrderChanged), group);

  /* Show the group's highlight colour in a button that will also launch a colour-picker */
  GtkWidget *colourButton = gtk_color_button_new_with_color(&group->highlightColour);
  g_signal_connect(G_OBJECT(colourButton), "color-set", G_CALLBACK(onGroupColourChanged), group);
  
  /* Create a button that will delete this group */
  GtkWidget *deleteButton = gtk_button_new_with_label("Delete");
  g_signal_connect(G_OBJECT(deleteButton), "clicked", G_CALLBACK(onButtonClickedDeleteGroup), group);
  
  /* Put everything in the table */
  gtk_table_attach(table, nameWidget,		1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, isHiddenWidget,	2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, isHighlightedWidget,	3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, orderWidget,		4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, colourButton,		5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, deleteButton,		6, 7, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
}


/* Called when the "from-name" radio button in the "create group" section of the
 * group sequences dialog box is toggled. It enables/disables the "from-name"
 * text entry box according to the state of the toggle button. */
static void onGroupSourceButtonToggled(GtkWidget *button, gpointer data)
{
  GtkWidget *widget = GTK_WIDGET(data);
  const gboolean enable = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  gtk_widget_set_sensitive(widget, enable);
  
  if (enable)
    {
      GtkWindow *window = GTK_WINDOW(gtk_widget_get_toplevel(button));
      gtk_window_set_focus(window, widget);
    }
}


/* Called for each entry in a hash table. Compares the key of the table, which is
 * a sequence name, to the search string given in the user data. If it matches, it
 * appends the sequence name to the result list in the user data. */
static void getSequencesThatMatch(gpointer key, gpointer value, gpointer data)
{
  CompareSeqData *searchData = (CompareSeqData*)data;
  char *seqName = (char *)(key);
  
  if (wildcardSearch(seqName, searchData->searchStr))
    {
      searchData->matchList = g_list_append(searchData->matchList, seqName);
    }
}


/* Called when the user has clicked "add group" on the "group sequences" dialog. */
void onButtonClickedAddGroup(GtkWidget *button, gpointer data)
{
  /* The text entry box was passed as the user data */
  GtkEntry *entry = GTK_ENTRY(data);
  
  /* Extract the main window from our parent window. */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);

  /* If the entry box is enabled, that means the user has selected the option to
   * use a search string. */
  if (GTK_WIDGET_SENSITIVE(entry))
    {
      const char *searchStr = gtk_entry_get_text(entry);

      /* Loop through all the sequences and see if the name matches the search string */
      GHashTable *seqTable = detailViewGetSeqTable(detailView);
      CompareSeqData searchData = {searchStr, NULL};

      g_hash_table_foreach(seqTable, getSequencesThatMatch, &searchData);
      
      /* If we found anything, create a group */
      if (g_list_length(searchData.matchList) > 0)
	{
	  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
	  SequenceGroup *group = createSequenceGroup(properties->sequenceGroups);
	  group->seqList = searchData.matchList;
	  properties->sequenceGroups = g_list_append(properties->sequenceGroups, group);
	}
      else
	{
	  messout("No sequences found");
	}
    }
  else
    {
      /* Currently the only other option is to use the current selection */
      makeGroupFromSelection(mainWindow);
    }
  
  /* Re-sort the trees, because the group ordering has changed, which may affect the sort order */
  callFuncOnAllDetailViewTrees(detailView, resortTree);
  
  /* Refresh dialog by closing and re-opening. Open it on the 'edit groups' page,
   * to display the group that we've just created. */
  gtk_widget_destroy(GTK_WIDGET(dialogWindow));
  showGroupSequencesDialog(mainWindow, TRUE);
}


/* Called when the user has clicked the "delete all groups" button on the "group sequences" dialog. */
static void onButtonClickedDeleteAllGroups(GtkWidget *button, gpointer data)
{
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

  /* Ask the user if they're sure */
  gint response = runConfirmationBox(mainWindow, "Delete groups", "This will delete ALL groups. Are you sure?");
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      mainWindowDeleteAllSequenceGroups(mainWindow);
      
      /* Close the dialog, because there are no groups left to display. */
      gtk_widget_destroy(GTK_WIDGET(dialogWindow));
    }
}


/* Called when the user has clicked the "delete group" button on the "group sequences" dialog. */
static void onButtonClickedDeleteGroup(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *mainWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  /* Ask the user if they're sure */
  char formatStr[] = "Are you sure you wish to delete group '%s'?";
  char messageText[strlen(formatStr) + strlen(group->groupName)];
  sprintf(messageText, formatStr, group->groupName);
  
  gint response = runConfirmationBox(mainWindow, "Delete group", messageText);
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      mainWindowDeleteSequenceGroup(mainWindow, group);
      
      /* Refresh the dialog by closing and re-opening. Only re-open if there are groups left. */
      gtk_widget_destroy(GTK_WIDGET(dialogWindow));
      if (g_list_length(mainWindowGetSequenceGroups(mainWindow)) > 0)
	{
	  showGroupSequencesDialog(mainWindow, TRUE);
	}
    }
}


/* Shows the "Group sequences" dialog. This dialog allows the user to group sequences together.
 * This tabbed dialog shows both the 'create group' and 'edit groups' dialogs in one. If the
 * 'editGroups' argument is true and groups exist, the 'Edit Groups' tab is displayed by default;
 * otherwise the 'Create Groups' tab is shown. */
void showGroupSequencesDialog(GtkWidget *mainWindow, const gboolean editGroups)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Groups", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CLOSE,
						  GTK_RESPONSE_REJECT,
						  NULL);
  
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;
  
  GtkWidget *notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(contentArea), notebook, TRUE, TRUE, 0);
  
  
  /* "CREATE GROUP" SECTION. */
  GtkWidget *section1 = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), section1, gtk_label_new("Create group"));

  GList *selectedSeqs = mainWindowGetSelectedSeqs(mainWindow);
  const gboolean seqsSelected = g_list_length(selectedSeqs) > 0;
  
  /* Radio buttons for user to specify whether to create group from current selection or from sequence names */
  GtkWidget *fromSelectionButton = gtk_radio_button_new_with_label(NULL, "From selection");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(fromSelectionButton), seqsSelected);
  gtk_box_pack_start(GTK_BOX(section1), fromSelectionButton, FALSE, FALSE, 0);

  GtkWidget *fromNameButton = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(fromSelectionButton), "From sequence name(s)");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(fromNameButton), !seqsSelected);
  gtk_box_pack_start(GTK_BOX(section1), fromNameButton, FALSE, FALSE, 0);
  
  /* Text box for user to enter sequence names to search for. Greyed out if relevant
   * radio button is not selected. Activates the dialog's default widget when enter is pressed. */
  GtkWidget *fromNameEntry = gtk_entry_new();
  gtk_widget_set_sensitive(fromNameEntry, !seqsSelected);
  gtk_entry_set_activates_default(GTK_ENTRY(fromNameEntry), TRUE);
  gtk_box_pack_start(GTK_BOX(section1), fromNameEntry, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(fromNameButton), "toggled", G_CALLBACK(onGroupSourceButtonToggled), fromNameEntry);
  
  /* Add the "add group" button here, because it is only relevant to this tab. */
  GtkWidget *addGroupButton = gtk_button_new_with_label("Create group");
  gtk_box_pack_end(GTK_BOX(section1), addGroupButton, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(addGroupButton), "clicked", G_CALLBACK(onButtonClickedAddGroup), fromNameEntry);
  
  
  /* "EDIT GROUP" SECTION. (Only relevant if some groups actually exist) */
  GList *groupList = mainWindowGetSequenceGroups(mainWindow);
  const int numRows = g_list_length(groupList) + 1; /* +1 for headers */
  const gboolean groupsExist = numRows > 1;
  
  if (groupsExist)
    {
      GtkWidget *section2 = gtk_vbox_new(FALSE, 0);
      gtk_notebook_append_page(GTK_NOTEBOOK(notebook), section2, gtk_label_new("Edit groups"));

      const int numCols = 6;
      GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
      gtk_box_pack_start(GTK_BOX(section2), GTK_WIDGET(table), FALSE, FALSE, 0);

      const int xpad = 2;
      const int ypad = 2;

      /* Add labels */
      int row = 1;
      gtk_table_attach(table, gtk_label_new("Hide"),	  2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, gtk_label_new("Highlight"), 3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, gtk_label_new("Order"),	  4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, gtk_label_new("Colour"),	  4, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, gtk_label_new("Delete"),	  6, 7, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      ++row;
      
      /* Add a set of widgets for each group */
      GList *groupItem = mainWindowGetSequenceGroups(mainWindow);
      for ( ; groupItem; groupItem = groupItem->next)
	{
	  SequenceGroup *group = (SequenceGroup*)(groupItem->data);
	  createEditGroupWidget(group, table, row, xpad, ypad);
	  ++row;
	}
      
      /* Add a button to delete all groups */
      GtkWidget *deleteGroupsButton = gtk_button_new_with_label("Delete all groups");
      gtk_box_pack_end(GTK_BOX(section2), deleteGroupsButton, FALSE, FALSE, 0);
      g_signal_connect(G_OBJECT(deleteGroupsButton), "clicked", G_CALLBACK(onButtonClickedDeleteAllGroups), NULL);
    }

  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
  
  /* If user has asked to edit groups (and some groups exist), make the second tab
   * the default and the 'close' button the default action. (Must do this after showing
   * the child widgets due to a GTK legacy whereby the notebook won't change tabs otherwise.) */
  if (editGroups && groupsExist)
    {
      gtk_notebook_next_page(GTK_NOTEBOOK(notebook));
      gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_REJECT);
    }
  else
    {
      /* Make the 'create group' button the default widget and focus it if rows
       * are selected - otherwise focus the text entry box. */
      GTK_WIDGET_SET_FLAGS(addGroupButton, GTK_CAN_DEFAULT);
      gtk_window_set_default(GTK_WINDOW(dialog), addGroupButton);

      if (seqsSelected)
	{
	  gtk_window_set_focus(GTK_WINDOW(dialog), addGroupButton);
	}
      else
	{
	  gtk_window_set_focus(GTK_WINDOW(dialog), fromNameEntry);
	}
    }
}


/***********************************************************
 *			   Settings menu                   *
 ***********************************************************/

/* Callback function called when the 'squash matches' button is toggled */
static void onSquashMatches(GtkWidget *button, gpointer data)
{
  const gboolean squash = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  GtkWidget *detailView = GTK_WIDGET(data);
  
  detailViewSquashMatches(detailView, squash);
}


/* Callback function called when the 'invert sort order' button is toggled */
static void onSortOrderToggled(GtkWidget *button, gpointer data)
{
  const gboolean invert = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  GtkWidget *detailView = GTK_WIDGET(data);
  
  detailViewSetSortInverted(detailView, invert);
}


/* Callback function called when the 'highlight differences' button is toggled */
static void onHighlightDiffsToggled(GtkWidget *button, gpointer data)
{
  const gboolean highlightDiffs = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  GtkWidget *detailView = GTK_WIDGET(data);
  
  detailViewSetHighlightDiffs(detailView, highlightDiffs);
}


/* Utility to create a check button with certain given properties, and to pack it into the parent */
static void createCheckButton(GtkWidget *parent, 
			      const char *mnemonic, 
			      const gboolean isActive, 
			      GCallback callback, 
			      GtkWidget *mainWindow)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_container_add(GTK_CONTAINER(parent), button);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), isActive);
  
  g_signal_connect(G_OBJECT(button), "toggled", callback, mainWindow);
}


/* Shows the "Settings" dialog. */
void showSettingsDialog(GtkWidget *mainWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Blixem Settings", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
  GtkWidget *vbox = createVBoxWithBorder(contentArea, borderWidth);
  
  /* Squash matches */
  createCheckButton(vbox, "_Squash matches", detailViewGetMatchesSquashed(detailView), G_CALLBACK(onSquashMatches), detailView);
  
  /* Invert sort order */
  createCheckButton(vbox, "_Invert sort order", detailViewGetSortInverted(detailView), G_CALLBACK(onSortOrderToggled), detailView);
  
  /* Highlight differences */
  createCheckButton(vbox, "_Highlight differences", detailViewGetHighlightDiffs(detailView), G_CALLBACK(onHighlightDiffsToggled), detailView);
  
  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
}


/***********************************************************
 *		    Statistics menu			   *
 ***********************************************************/

static void getStats(GtkWidget *mainWindow, GString *result, MSP *MSPlist)
{
  int numMSPs = 0;          /* count of MSPs */
  int numSeqs = 0;          /* count of sequences (excluding duplicates) */
  gint32 totalSeqSize = 0;  /* total memory used by the sequences */
  GHashTable *sequences = g_hash_table_new(NULL, NULL); /* remembers which sequences we've already seen */
  
  MSP *msp = NULL;
  for (msp = MSPlist; msp ; msp = msp->next)
    {
      ++numMSPs;
      
      if (!g_hash_table_lookup(sequences, msp->sseq))
        {
          ++numSeqs;
	  if (msp->sseq)
	    {
	      totalSeqSize += strlen(msp->sseq) * sizeof(char);
	    }
          
          g_hash_table_insert(sequences, msp->sseq, &msp->slength);
        }
    }
  
  g_hash_table_destroy(sequences);
  
  int refSeqLen = strlen(mainWindowGetRefSeq(mainWindow));
  
  /* Create the text based on the results */
  g_string_printf(result, "%s%d%s%s%d%s%s%d%s%s%d%s%s%d%s%s%d%s",
                  "Length of reference sequence\t\t= ", refSeqLen, " characters\n\n",
                  "Number of match sequences\t\t= ", numSeqs, "\n",
                  "Total memory used by sequences\t= ", totalSeqSize, " bytes\n\n",
		  "Number of MSPs\t\t\t\t\t= ", numMSPs, "\n",
                  "Size of each MSP\t\t\t\t\t= ", (int)sizeof(MSP), " bytes\n",
		  "Total memory used by MSPs\t\t= ", (int)sizeof(MSP) * numMSPs, " bytes");
}


static void showStatsDialog(GtkWidget *mainWindow, MSP *MSPlist)
{
  /* Create a dialog widget with an OK button */
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Statistics", 
                                                  GTK_WINDOW(mainWindow),
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Ensure that the dialog box (along with any children) is destroyed when the user responds. */
  g_signal_connect (dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  /* Create a text buffer containing the required text*/
  GString *displayText = g_string_sized_new(200); /* will be extended if we need more space */
  getStats(mainWindow, displayText, MSPlist);
  GtkTextBuffer *textBuffer = gtk_text_buffer_new(gtk_text_tag_table_new());
  gtk_text_buffer_set_text(GTK_TEXT_BUFFER(textBuffer), displayText->str, -1);
  g_string_free(displayText, TRUE);
  
  /* Create a text view widget and put it in the vbox area of the dialog */
  GtkWidget *textView = gtk_text_view_new_with_buffer(GTK_TEXT_BUFFER(textBuffer));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), textView, TRUE, TRUE, 0);
  
  /* Show the dialog */
  gtk_widget_show(textView);
  gtk_widget_show(dialog);
}


/***********************************************************
 *			Help menu			   *
 ***********************************************************/

/* Utility to pop up a simple confirmation dialog box with the given title and text, 
 * with just an "OK" and "Cancel" button. Blocks until the user responds, and returns
 * their response ID. */
static gint runConfirmationBox(GtkWidget *mainWindow, char *title, char *messageText)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  /* Add the message */
  GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
  GtkWidget *label = gtk_label_new(messageText);
  gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
  gtk_widget_show(label);
  
  gint response = gtk_dialog_run(GTK_DIALOG(dialog));
  
  gtk_widget_destroy(dialog);
  
  return response;
}


void displayHelp(GtkWidget *mainWindow)
{
  char *messageText = (messprintf("\
BLIXEM\n\
BLast matches In an X-windows Embedded Multiple alignment\n\
\n\
\n\
pfetch\n\
	•	Double-click a row to pfetch a sequence.\n\
\n\
\n\
MAIN MENU\n\
Right-click anywhere in the Blixem window to pop up the main menu.  The menu options are:\n\
	•	Quit: Quit Blixem.\n\
	•	Help: Display this help.\n\
	•	Print: Print the display\n\
	•	Settings: Edit settings.\n\
	•	View: Display the View dialog box. This allows you to show/hide parts of the display.\n\
	•	Create Group: Create a group of sequences.\n\
	•	Edit Groups: Edit properties for groups.\n\
	•	Deselect all: Deselect all sequences.\n\
	•	Dotter: Run Dotter and/or set Dotter parameters.\n\
	•	Feature series selection tool: ???\n\
\n\
\n\
SCROLLING\n\
	•	Middle-click in the big picture to scroll to an area.\n\
	•	Middle-click in the detail view to centre on a base.\n\
	•	Use the horizontal scrollbar at the bottom of the window to scroll the detail view.\n\
	•	Use the mouse-wheel to scroll up/down and left/right in the detail view (if your mouse has a horizontal scroll-wheel).\n\
	•	Use the mouse-wheel to scroll the display area left/right from the big picture.\n\
	•	The left/right arrow keys move the selection one nucleotide to the left/right.\n\
	•	When viewing protein matches, holding CTRL while using the left/right arrow keys moves the selection one full codon to the left/right.\n\
	•	The ',' and '.' keys scroll the screen one nucleotide to the left/right.\n\
\n\
\n\
SELECTIONS\n\
	•	You can select a sequence by left-clicking on its row in the detail view.  Selected sequences are highlighted in cyan in the big picture.\n\
	•	The name of the sequence you selected is displayed in the feedback box on the toolbar.  If multiple alignments for this sequence are currently within the display range, all of their rows will be selected.\n\
	•	You can select multiple sequences by holding down the CTRL or SHIFT keys while selecting rows.\n\
	•	You can deselect a single sequence by CTRL-clicking on its row.\n\
	•	You can deselect all sequences by right-clicking and selecting 'Deselect all', or with the SHIFT-CTRL-A keyboard shortcut.\n\
	•	You can move the selection up/down a row using the up/down arrow keys.\n\
\n\
	•	You can select a nucleotide/peptide by middle-clicking on it in the detail view.  This selects the entire column at that index, and the coordinate number on the query sequence is shown in the feedback box.  (The coordinate on the subject sequence is also shown if a subject sequence is selected.)\n\
	•	For protein matches, when a peptide is selected, the three nucleotides for that peptide (for the current reading frame) are highlighted in the header in green and red.  The current reading frame is whichever alignment list currently has the focus - click in a different list to change the reading frame.  The red highlighting indicates the specific nucleotide that is currently selected and whose coordinate is displayed in the feedback box.\n\
	•	You can move the selection to the previous/next nucleotide using the left and right arrow keys.\n\
	•	You can move the selection to the previous/next peptide by holding CTRL while using the left and right arrow keys.\n\
	•	Middle-clicking also scrolls the display to centre on the clicked coordinate.  To select a coordinate without the display scrolling, hold down CTRL as you middle-click.\n\
\n\
\n\
ZOOMING\n\
	•	Zoom in to the currently-selected region in the big picture using the +/- buttons in the top-left corner of the window, or using the CTRL '=' and CTRL '-' shortcut keys.  The 'Whole' button zooms out to show the full length of the query sequence.\n\
	•	Zoom in/out of the detail view using the +/- buttons on the toolbar, or using the '=' and '-' shortcut keys.\n\
\n\
\n\
SORTING\n\
	•	The alignments can be sorted by selecting the column you wish to sort by from the drop-down box on the toolbar.  To reverse the sort order, select the relevant option under the Settings menu.\n\
	•	The alignments can also be sorted by group by selecting the Group option from the drop-down box.  See the Groups section.\n\
\n\
\n\
GROUPS\n\
Alignments can be grouped together so that they can be sorted/highlighted/hidden etc.\n\
\n\
Creating a group from a selection:\n\
	•	Select the sequences you wish to include in the group by left-clicking their rows in the detail view.  Multiple rows can be selected by holding the CTRL or SHIFT keys while clicking.\n\
	•	Right-click and select 'Create Group', or use the SHIFT-CTRL-G shortcut key. (Note that CTRL-G will also shortcut to here if no groups currently exist.)\n\
	•	Ensure that the 'From selection' radio button is selected, and click 'Create group'.\n\
\n\
Creating a group from a sequence name:\n\
	•	Right-click and select 'Create Group', or use the SHIFT-CTRL-G shortcut key. (Or CTRL-G if no groups currently exist.)\n\
	•	Select the 'From sequence name(s)' radio button and enter the name of the sequence in the box below.  You may use the following wildcards to search for sequences: '*' for any number of characters; '?' for a single character.  For example, searching for '*X' will find all sequences whose name ends in 'X' (i.e. all exons).\n\
	•	Click 'Create group'.\n\
\n\
Editing groups:\n\
To edit a group, right-click and select 'Edit Groups', or use the CTRL-G shortcut key. You can change the following properties for a group:\n\
	•	Name: you can specify a more meaningful name to help identify the group.\n\
	•	Hide: tick this box to hide the alignments in the group from the detail view.\n\
	•	Highlight: tick this box to highlight the alignments in this group.\n\
	•	Order: when sorting by Group, alignments in a group with a lower order number will appear before those with a higher order number (or vice versa if sort order is inverted). Alignments in a group will appear before alignments that are not in a group.\n\
	•	Colour: the colour the group will be highlighted in, if 'Highlight' is enabled.  The default colour for all groups is red, so you may wish to change this if you want different groups to be highlighted in different colours.\n\
	•	To delete a single group, click on the 'Delete' button next to the group you wish to delete.\n\
	•	To delete all groups, click on the 'Delete all groups' button.\n\
\n\
\n\
DOTTER\n\
	•	To start Dotter, or to edit the Dotter parameters, right-click and select 'Dotter' or use the CTRL-D keyboard shortcut.	The Dotter settings dialog will pop up.\n\
	•	To run Dotter with the default (automatic) parameters, just hit RETURN, or click the 'Execute' button.\n\
	•	To enter manual parameters, click the 'Manual' radio button and enter the values in the 'Start' and 'End' boxes.\n\
	•	To revert to the last-saved manual parameters, click the 'Last saved' button.\n\
	•	To revert back to automatic parameters, click the 'Auto' radio button.\n\
	•	To save the parameters without running Dotter, click Save and then Cancel'.\n\
	•	To save the parameters and run Dotter, click 'Execute'.\n\
\n\
\n\
KEYBOARD SHORTCUTS\n\
	•	CTRL Q: Quit\n\
	•	CTRL H: Help\n\
	•	CTRL P: Print\n\
	•	CTRL S: Settings menu\n\
	•	CTRL V: View menu\n\
	•	SHIFT CTRL G: Create group\n\
	•	CTRL G: Edit groups (or create a group if none currently exist)\n\
	•	SHIFT CTRL A: Select all visible sequences\n\
	•	SHIFT CTRL A: Deselect all sequences\n\
	•	CTRL D: Dotter\n\
	•	Left/right arrow keys: Move the currently-selected coordinate one nucleotide to the left/right\n\
	•	CTRL and left/right arrow keys: Move the currently-selected coordinate one peptide to the left/right\n\
	•	'=' : zoom in detail view\n\
	•	'-' : zoom out detail view\n\
	•	CTRL '=' : zoom in big picture\n\
	•	CTRL '-' :  zoom out big picture\n\
	•	, (comma): scroll left one coordinate\n\
	•	. (full-stop): scroll right one coordinate\n\
	•	g: Go to coordinate\n\
\n\
\n\
SETTINGS\n\
	•	The settings menu can be accessed by right-clicking and selecting Settings, or by the shortcut CTRL-S.\n\
	•	Squash Matches: this groups multiple alignments from the same sequence together into the same row in the detail view, rather than showing them on separate rows.\n\
\n\
\n\
KEY\n\
In the detail view, the following colours and symbols have the following meanings:\n\
	•	Yellow background: query sequence\n\
	•	Cyan: identical residues\n\
	•	Violet: conserved residues\n\
	•	Grey: mismatch\n\
	•	Grey with a '.': deletion\n\
	•	Yellow vertical line: insertion\n\
	•	Thin blue vertical line: start of an exon\n\
	•	Thin dark-blue vertical line: end of an exon\n\
	•	Green background (protein matches only): the three nucleotides for the currently-selected codon. Dark green indicates the specific nucleotide that is currently displayed in the feedback box.\n\
", blixemVersion));

  /* Set a pretty big initial size */
  const int initWidth = mainWindow->allocation.width * 0.7;
  const int maxHeight = mainWindow->allocation.height * 0.7;
  
  showMessageDialog("Help", messageText, NULL, initWidth, maxHeight, TRUE, mainWindow->style->font_desc);
}


/***********************************************************
 *			  Menu actions                     *
 ***********************************************************/

/* Called when the user selects the Statistics menu option, or hits the Statistics shortcut key.
 * Pops up a dialog showing user help information */
static void onHelpMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  displayHelp(mainWindow);
}


/* Called when the user selects the View menu option, or hits the Settings shortcut key */
static void onViewMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showViewPanesDialog(mainWindow);
}


/* Called when the user selects the 'Group Sequences' menu option, or hits the relevant shortcut key */
static void onCreateGroupMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showGroupSequencesDialog(mainWindow, FALSE);
}


/* Called when the user selects the 'Groups' menu option, or hits the relevant shortcut key */
static void onEditGroupsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showGroupSequencesDialog(mainWindow, TRUE);
}


/* Called when the user selects the Settings menu option, or hits the Settings shortcut key */
static void onSettingsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showSettingsDialog(mainWindow);
}


/* Called when the user selects the Dotter menu option, or hits the Dotter shortcut key */
static void onDotterMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showDotterDialog(mainWindow);
}


/* Called when the user selects the 'Select Features' menu option, or hits the relevant shortcut key */
static void onSelectFeaturesMenu(GtkAction *action, gpointer data)
{
  selectFeatures();
}


/* Called when the user selects the 'Deselect all' menu option, or hits the relevant shortcut key */
static void onDeselectAllRows(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  mainWindowDeselectAllSeqs(mainWindow, TRUE);
}


/* Called when the user selects the Statistics menu option, or hits the Statistics shortcut key.
 * Pops up a dialog showing memory usage statistics. */
static void onStatisticsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MSP *mspList = mainWindowGetMspList(mainWindow);
  showStatsDialog(mainWindow, mspList);
}


/* Called when the user selects the Print menu option, or hits the Print shortcut key */
static void onPrintMenu(GtkAction *action, gpointer data)
{
//  GtkWidget *mainWindow = GTK_WIDGET(data);
//  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
//  
//  /* Create a print operation, using the same settings as the last print, if there was one */
//  GtkPrintOperation *print = gtk_print_operation_new();
//  
////  if (properties->printSettings != NULL)
////    {
////      gtk_print_operation_set_print_settings(print, properties->printSettings);
////    }
//  
////  g_signal_connect (print, "begin_print", G_CALLBACK (onBeginPrint), mainWindow);
////  g_signal_connect(G_OBJECT(print), "draw-page", G_CALLBACK(onDrawPage), mainWindow);
//  
//  /* Pop up the print dialog */
//  GtkPrintOperationResult printResult = gtk_print_operation_run (print, 
//								 GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
//								 GTK_WINDOW(mainWindow),
//								 NULL);
//  
//  /* If the user hit ok, remember the print settings for next time */
//  if (printResult == GTK_PRINT_OPERATION_RESULT_APPLY)
//    {
////      if (properties->printSettings != NULL)
////	{
////	  g_object_unref(properties->printSettings);
////	}
////      
////      properties->printSettings = g_object_ref(gtk_print_operation_get_print_settings(print));
//    }
//
//  g_object_unref(print);
}


/***********************************************************
 *			   Events                          *
 ***********************************************************/

/* Called after the user clicks ok in the print dialog. For now just scales the 
 * whole output to fit on a single page. */
//static void onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data)
//{
//  GtkWidget *mainWindow = GTK_WIDGET(data);
//  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
//
//  /* Set the orientation based on the stored print settings (not sure why this doesn't
//   * already get set when the settings are set in the print operation...). */
////  GtkPrintSettings *printSettings = gtk_print_operation_get_print_settings(print);
////  gtk_print_settings_set_orientation(printSettings, gtk_print_settings_get_orientation(properties->printSettings));
//
//  gtk_print_operation_set_n_pages(print, 1);
//}


static void collatePixmaps(GtkWidget *widget, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  GdkDrawable *drawable = widgetGetDrawable(widget);

  /* If this widget is visible and has a drawable set, draw it onto the main window's drawable */
  if (drawable && GTK_WIDGET_VISIBLE(widget))
    {
      /* Get the left edge of the widget in terms of the main window's coordinates */
      int xSrc = 0, ySrc = 0;
      int xDest, yDest;
      gtk_widget_translate_coordinates(widget, mainWindow, xSrc, ySrc, &xDest, &yDest);

      /* Get the coords where to draw */
      int x = xDest; /* this is a bit too simplified really but works ok for our mainly-vertical layout */
      int y = 0;
      
      if (yDest == properties->lastYCoord)
	{
	  /* Draw at the same height as the last widget */
	  y = properties->lastYStart;
	}
      else
	{
	  /* Draw at the end of the last widget, and increment the Y pos */
	  y = properties->lastYEnd;
	  
	  properties->lastYStart = properties->lastYEnd;
	  properties->lastYEnd = properties->lastYEnd + widget->allocation.height;
	  properties->lastYCoord = yDest;
	}
      
      GdkGC *gc = gdk_gc_new(widget->window);
      gdk_draw_drawable(properties->drawable, gc, drawable, xSrc, ySrc, x, y, -1, -1); /* -1 means full width/height */
    }
  
  /* If this widget is a container, recurse over its children */
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), collatePixmaps, mainWindow);
    }
}


/* Print handler - renders a specific page */
//static void onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data)
//{
//  GtkWidget *mainWindow = GTK_WIDGET(data);
//  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
//  
//  /* Create a blank white pixmap to draw on to */
//  properties->drawable = gdk_pixmap_new(mainWindow->window, mainWindow->allocation.width, mainWindow->allocation.height, -1);
//  
//  GdkGC *gc = gdk_gc_new(properties->drawable);
//  GdkColor fgColour = getGdkColor(GDK_WHITE);
//  gdk_gc_set_foreground(gc, &fgColour);
//  gdk_draw_rectangle(properties->drawable, gc, TRUE, 0, 0, mainWindow->allocation.width, mainWindow->allocation.height);
//
//  /* For each child widget that has a drawable set, draw this onto the main pixmap */
//  properties->lastYStart = 0;
//  properties->lastYEnd = 0;
//  properties->lastYCoord = -1;
//  gtk_container_foreach(GTK_CONTAINER(mainWindow), collatePixmaps, mainWindow);
//  
//  cairo_t *cr = gtk_print_context_get_cairo_context (context);
//  gdk_cairo_set_source_pixmap(cr, properties->drawable, 0, 0);
//  cairo_paint(cr);
//}


/* Mouse button handler */
static gboolean onButtonPressMainWindow(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3)
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
  }
  
  return TRUE;
}


/* Key press handler */
static gboolean onKeyPressMainWindow(GtkWidget *window, GdkEventKey *event, gpointer data)
{
  gboolean result = FALSE;
  
  const gboolean ctrlModifier = (event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK;
  const gboolean shiftModifier = (event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK;
  
  switch (event->keyval)
    {
      case GDK_Left: /* fall through */
      case GDK_Right:
	{
	  if (ctrlModifier)
	    moveSelectedDisplayIdxBy1(window, event->keyval == GDK_Left);
	  else
	    moveSelectedBaseIdxBy1(window, event->keyval == GDK_Left);
	  
	  result = TRUE;
	  break;
	}

      case GDK_comma:  /* fall through */
      case GDK_period:
	{
	  scrollDetailViewBy1(window, event->keyval == GDK_comma);
	  result = TRUE;
	  break;
	}
	
      case GDK_equal: /* fall through */
      case GDK_minus:
	{
	  const gboolean zoomIn = (event->keyval == GDK_plus || event->keyval == GDK_equal);
	  zoomMainWindow(window, zoomIn, ctrlModifier); /* if ctrl pressed, zoom big picture - else zoom detail view */
	  result = TRUE;
	  break;
	}
	
      case GDK_w: /* fall through */
      case GDK_W: /* fall through */
	{
	  if (ctrlModifier)
	    {
	      zoomWholeBigPicture(mainWindowGetBigPicture(window));
	      result = TRUE;
	    }
	  break;
	}
	
      case GDK_c: /* fall through */
      case GDK_C:
	if (ctrlModifier)
	  {
	    copySelectionToClipboard(window);
	    result = TRUE;
	  }
	break;
	
      case GDK_g: /* fall through */
      case GDK_G:
	if (!ctrlModifier)
	  {
	    goToDetailViewCoord(mainWindowGetDetailView(window), BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
	    result = TRUE;
	  }
	break;
	
      case GDK_t: /* fall through */
      case GDK_T:
	toggleStrand(mainWindowGetDetailView(window));
	result = TRUE;
	break;
	
      case GDK_1: /* fall through */
      case GDK_exclam:
	togglePaneVisibility(window, 1, ctrlModifier, shiftModifier);
	result = TRUE;
	break;

      case GDK_2:	  /* fall through */
      case GDK_quotedbl:  /* fall through */
      case GDK_at:
	togglePaneVisibility(window, 2, ctrlModifier, shiftModifier);
	result = TRUE;
	break;
	
      case GDK_3: /* fall through */
      case GDK_currency:
	togglePaneVisibility(window, 3, ctrlModifier, shiftModifier);
	result = TRUE;
	break;
    };
  
  return result;
}

/***********************************************************
 *			   Properties                      *
 ***********************************************************/

MainWindowProperties* mainWindowGetProperties(GtkWidget *widget)
{
  return widget ? (MainWindowProperties*)(g_object_get_data(G_OBJECT(widget), "MainWindowProperties")) : NULL;
}

static void onDestroyMainWindow(GtkWidget *widget)
{
  MainWindowProperties *properties = mainWindowGetProperties(widget);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "MainWindowProperties", NULL);
    }
}

static void mainWindowCreateProperties(GtkWidget *widget, 
				       GtkWidget *bigPicture, 
				       GtkWidget *detailView,
				       GtkWidget *mainmenu,
				       MSP *mspList,
				       const BlxBlastMode blastMode,
				       char *refSeq,
				       const char *refSeqName,
				       const IntRange const *refSeqRange,
				       const IntRange const *fullDisplayRange,
				       const BlxSeqType seqType,
				       char **geneticCode,
				       int numReadingFrames,
				       const gboolean gappedHsp)
{
  if (widget)
    {
      MainWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->detailView = detailView;
      properties->mainmenu = mainmenu;
      
      properties->refSeq = refSeq;
      properties->refSeqName = refSeqName ? g_strdup(refSeqName) : g_strdup("Blixem-seq");
      properties->refSeqRange.min = refSeqRange->min;
      properties->refSeqRange.max = refSeqRange->max;
      properties->fullDisplayRange.min = fullDisplayRange->min;
      properties->fullDisplayRange.max = fullDisplayRange->max;
      
      properties->mspList = mspList;
      properties->geneticCode = geneticCode;
      properties->blastMode = blastMode;
      properties->seqType = seqType;
      properties->numReadingFrames = numReadingFrames;
      
      properties->strandsToggled = FALSE;
      properties->selectedSeqs = NULL;
      properties->sequenceGroups = NULL;
      
      properties->autoDotterParams = TRUE;
      properties->dotterStart = UNSET_INT;
      properties->dotterEnd = UNSET_INT;
      properties->dotterZoom = 0;

      properties->drawable = NULL;
//      properties->printSettings = gtk_print_settings_new();
//      gtk_print_settings_set_orientation(properties->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
//      gtk_print_settings_set_quality(properties->printSettings, GTK_PRINT_QUALITY_HIGH);
      properties->lastYEnd = UNSET_INT;
      properties->lastYStart = UNSET_INT;
      properties->lastYCoord = UNSET_INT;

      g_object_set_data(G_OBJECT(widget), "MainWindowProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyMainWindow), NULL); 
    }
}

gboolean mainWindowGetStrandsToggled(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->strandsToggled : FALSE;
}

GtkWidget* mainWindowGetBigPicture(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->bigPicture : FALSE;
}

GtkWidget* mainWindowGetDetailView(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->detailView : FALSE;
}

GtkWidget* mainWindowGetMainMenu(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->mainmenu : FALSE;
}

BlxBlastMode mainWindowGetBlastMode(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->blastMode : FALSE;
}

char * mainWindowGetRefSeq(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->refSeq : NULL;
}

const char * mainWindowGetRefSeqName(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->refSeqName : NULL;
}

char** mainWindowGetGeneticCode(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->geneticCode : NULL;
}

MSP* mainWindowGetMspList(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->mspList : FALSE;
}

BlxSeqType mainWindowGetSeqType(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->seqType : FALSE;
}

IntRange* mainWindowGetFullRange(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? &properties->fullDisplayRange : NULL;
}

IntRange* mainWindowGetRefSeqRange(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? &properties->refSeqRange : NULL;
}

int mainWindowGetNumReadingFrames(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->numReadingFrames : UNSET_INT;
}

int mainWindowGetAutoDotter(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->autoDotterParams : TRUE;
}

int mainWindowGetDotterStart(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->dotterStart : UNSET_INT;
}

int mainWindowGetDotterEnd(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->dotterEnd : UNSET_INT;
}

gboolean mainWindowGetGappedHsp(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->gappedHsp : UNSET_INT;
}

/* Return the active strand - forward strand by default, reverse strand if display toggled */
static Strand mainWindowGetActiveStrand(GtkWidget *mainWindow)
{
  return mainWindowGetStrandsToggled(mainWindow) ? REVERSE_STRAND : FORWARD_STRAND;
}

/* Return the inactive strand - reverse strand by default, forward strand if display toggled */
static Strand mainWindowGetInactiveStrand(GtkWidget *mainWindow)
{
  return mainWindowGetStrandsToggled(mainWindow) ? FORWARD_STRAND : REVERSE_STRAND;
}

/* Returns a list of all the MSPs with the given sequence name */
GList *mainWindowGetSequenceMsps(GtkWidget *mainWindow, const char *seqName)
{
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  return detailViewGetSequenceMsps(detailView, seqName);
}

/* Returns the list of all sequence groups */
GList *mainWindowGetSequenceGroups(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties->sequenceGroups;
}

/* Returns the group that the given sequence belongs to, if any */
SequenceGroup *mainWindowGetSequenceGroup(GtkWidget *mainWindow, const char *seqName)
{
  SequenceGroup *result = NULL;
  
  /* Loop through all the groups until we find this sequence in one */
  GList *groupItem = mainWindowGetSequenceGroups(mainWindow);
  for ( ; groupItem; groupItem = groupItem->next)
    {
      /* Loop through all sequences in the group and compare the name */
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      GList *seqItem = group->seqList;
      
      for ( ; seqItem; seqItem = seqItem->next)
	{
	  const char *compareName = (const char*)(seqItem->data);
	  if (strcmp(compareName, seqName) == 0)
	    {
	      result = group;
	      break;
	    }
	}
    }
  
  return result;
}


/***********************************************************
 *                        Selections                       *
 ***********************************************************/

GList* mainWindowGetSelectedSeqs(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->selectedSeqs : NULL;
}


/* Get the selected sequences as a list of sequence names. The returned list
 * is formatted as comma-separated values. */
static GString* mainWindowGetSelectedSeqNames(GtkWidget *mainWindow)
{
  GList *listItem = mainWindowGetSelectedSeqs(mainWindow);
  GString *result = g_string_new_len(NULL, 50);
  gboolean first = TRUE;
  
  for ( ; listItem; listItem = listItem->next)
    {
      /* Add a separator before the name, unless it's the first one */
      if (!first)
	{
	  g_string_append(result, ", ");
	}
      else
	{
	  first = FALSE;
	}

      const char *seqName = (const char*)(listItem->data);
      g_string_append(result, seqName);
    }

  return result;
}


/* This function copys the currently-selected sequences' names to the default
 * clipboard. */
static void copySelectionToClipboard(GtkWidget *mainWindow)
{
  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
  GString *displayText = mainWindowGetSelectedSeqNames(mainWindow);
  gtk_clipboard_set_text(clipboard, displayText->str, -1);
  g_string_free(displayText, TRUE);
}


/* Update function to be called whenever the MSP selection has changed */
void mainWindowSelectionChanged(GtkWidget *mainWindow, const gboolean updateTrees)
{
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
  /* Unless requested not to, update the tree row selections to match our new selections */
  if (updateTrees)
    {
      callFuncOnAllDetailViewTrees(detailView, deselectAllRows);
      callFuncOnAllDetailViewTrees(detailView, selectRowsForSelectedSeqs);
    }
  
  /* Update the feedback box to tell the user which sequence is selected. */
  updateFeedbackBox(detailView);
  
  /* Recreate the grids (because the MSP lines will change colour with highlighting) */
  GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
  callFuncOnAllBigPictureGrids(bigPicture, widgetClearCachedDrawable);
  gtk_widget_queue_draw(bigPicture);
  
  /* Copy the selected sequence names to the "PRIMARY" clipboard, as per common
   * behaviour for X. (Not applicable for Windows.) */
#ifndef __CYGWIN__
  
  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_PRIMARY);
  GString *displayText = mainWindowGetSelectedSeqNames(mainWindow);
  gtk_clipboard_set_text(clipboard, displayText->str, -1);
  g_string_free(displayText, TRUE);
  
#endif
}


/* Call this function to select the given match sequence */
void mainWindowSelectSeq(GtkWidget *mainWindow, char *seqName, const gboolean updateTrees)
{
  if (!mainWindowIsSeqSelected(mainWindow, seqName))
    {
      MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
      properties->selectedSeqs = g_list_append(properties->selectedSeqs, seqName);
      mainWindowSelectionChanged(mainWindow, updateTrees);
    }
}


/* Call this function to deselect the given sequence */
void mainWindowDeselectSeq(GtkWidget *mainWindow, char *seqName, const gboolean updateTrees)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);

  /* See if it's in the list and, if so, get a pointer to the list element */
  GList *foundSeq = findSelectedSeqInList(properties->selectedSeqs, seqName);

  if (foundSeq)
    {
      properties->selectedSeqs = g_list_remove(properties->selectedSeqs, foundSeq->data);
      mainWindowSelectionChanged(mainWindow, updateTrees);
    }
}


/* Call this function to deselect all sequences */
void mainWindowDeselectAllSeqs(GtkWidget *mainWindow, const gboolean updateTrees)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);

  if (g_list_length(properties->selectedSeqs) > 0)
    {
      g_list_free(properties->selectedSeqs);
      properties->selectedSeqs = NULL;
      mainWindowSelectionChanged(mainWindow, updateTrees);
    }
}


/* Returns the pointer to the GList element that contains the given sequence, if the
 * sequence is in the list. */
static GList *findSelectedSeqInList(GList *list, const char *seqName)
{
  GList *result = NULL;
  
  /* We can't use g_list_find because that checks if the pointer value
   * is the same, which it might not be, because the same seq name could 
   * come from any one of several different MSPs. */
  GList *listItem = list;
  
  for ( ; listItem; listItem = listItem->next)
    {
      const char *listName = (const char*)(listItem->data);
      
      if (strcmp(listName, seqName) == 0)
	{
	  result = listItem;
	  break;
	}
    }
  
  return result;
}


/* Returns true if the given sequence is selected */
gboolean mainWindowIsSeqSelected(GtkWidget *mainWindow, const char *seqName)
{
  gboolean result = FALSE;
  
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  if (properties)
    {
      result = (findSelectedSeqInList(properties->selectedSeqs, seqName) != NULL);
    }
  
  return result;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set various properties for the main window widget */
static void setStyleProperties(GtkWidget *widget)
{
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(widget);
  const int width = gdk_screen_get_width(screen) * DEFAULT_WINDOW_WIDTH_FRACTION;
  const int height = gdk_screen_get_height(screen) * DEFAULT_WINDOW_HEIGHT_FRACTION;
  
  gtk_window_set_default_size(GTK_WINDOW(widget), width, height);
  
  gtk_container_set_border_width (GTK_CONTAINER(widget), DEFAULT_WINDOW_BORDER_WIDTH); 
  gtk_window_set_mnemonic_modifier(GTK_WINDOW(widget), GDK_MOD1_MASK); /* MOD1 is ALT on most systems */
  
  /* Set the default font size to be a bit smaller than usual */
  int origSize = pango_font_description_get_size(widget->style->font_desc) / PANGO_SCALE;
  const char *origFamily = pango_font_description_get_family(widget->style->font_desc);

  char parseString[500];
  sprintf(parseString, "gtk-font-name = \"%s %d\"", origFamily, origSize + DEFAULT_FONT_SIZE_ADJUSTMENT);
  gtk_rc_parse_string(parseString);
}


/* Create the main menu */
static GtkWidget* createMainMenu(GtkWidget *window)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  gtk_action_group_add_actions (action_group, mainMenuEntries, G_N_ELEMENTS (mainMenuEntries), window);
  
  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);
  
  GError *error = NULL;
  const char *menuDescription = userIsDeveloper() ? developerMenuDescription : standardMenuDescription;
  
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, menuDescription, -1, &error))
    {
      g_message ("building menus failed: %s", error->message);
      g_error_free (error);
      exit (EXIT_FAILURE);
    }
  
  return gtk_ui_manager_get_widget (ui_manager, "/MainMenu");
}


/* Create the main window */
GtkWidget* createMainWindow(char *refSeq, 
			    const char const *refSeqName,
			    MSP *mspList, 
			    BlxBlastMode blastMode,
			    BlxSeqType seqType, 
			    int numFrames,
			    char **geneticCode,
			    const int refSeqOffset,
			    const int startCoord1Based,
			    const SortByType sortByTypeInput,
			    const gboolean sortInverted,
			    const gboolean gappedHsp)
{
  /* If no sort type was specified, sort by ID by default */
  SortByType sortByType = (sortByTypeInput == SORTBYUNSORTED) ? SORTBYID : sortByTypeInput;

  /* Get the range of the reference sequence. If this is a DNA sequence but our
   * matches are peptide sequences, we must convert to the peptide sequence. */
  const int refSeqLen = (int)strlen(refSeq);
  IntRange refSeqRange = {refSeqOffset + 1, refSeqOffset + refSeqLen};
  IntRange fullDisplayRange = {refSeqRange.min, refSeqRange.max};
  
  if (seqType == BLXSEQ_PEPTIDE)
    {
      fullDisplayRange.min = convertDnaIdxToDisplayIdx(refSeqRange.min, seqType, 1, numFrames, FALSE, &refSeqRange, NULL);
      
      int baseNum = UNSET_INT;
      fullDisplayRange.max = convertDnaIdxToDisplayIdx(refSeqRange.max, seqType, 3, numFrames, FALSE, &refSeqRange, &baseNum);
      
      if (baseNum < numFrames)
	{
	  /* The last peptide does not have a full triplet, so cut off the range at the last full triplet */
	  fullDisplayRange.max -= 1;
	}
      
//      printf("Converted DNA sequence (%d-%d) to peptide sequence (%d-%d).\n",  refSeqRange.min, refSeqRange.max, fullDisplayRange.min, fullDisplayRange.max);
    }
  
  /* Convert the start coord (which is 1-based and on the DNA sequence) to display
   * coords (which take into account the offset and may also be peptide coords) */
  int startCoord = startCoord1Based + refSeqOffset;
  if (seqType == BLXSEQ_PEPTIDE)
    {
      startCoord = convertDnaIdxToDisplayIdx(startCoord, seqType, 1, numFrames, FALSE, &refSeqRange, NULL);
    }
  
  
  /* Create the main window */
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setStyleProperties(window);
  
  /* Create a vertical box to pack everything in */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), vbox);
  
  /* Create the main menu */
  GtkWidget *mainmenu = createMainMenu(window);
  
  /* Create the widgets. We need a single adjustment for the entire detail view, which will also be referenced
   * by the big picture view, so create it first. */
  GtkAdjustment *detailAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, /* initial value = 0 */
								      fullDisplayRange.min, /* lower value */
								      fullDisplayRange.max, /* upper value */
								      DEFAULT_SCROLL_STEP_INCREMENT, /* step increment used for mouse wheel scrolling */
								      0,   /* page increment dynamically set based on display range */
								      0)); /* page size dunamically set based on display range */
  
  GtkWidget *fwdStrandGrid = NULL, *revStrandGrid = NULL;
  
  GtkWidget *bigPicture = createBigPicture(window,
					   vbox, 
					   &fwdStrandGrid, 
					   &revStrandGrid, 
					   &fullDisplayRange);

  GtkWidget *detailView = createDetailView(window,
					   vbox, 
					   detailAdjustment, 
					   fwdStrandGrid, 
					   revStrandGrid,
					   mspList,
					   blastMode,
					   seqType,
					   numFrames,
					   refSeqName,
					   startCoord,
					   sortInverted,
					   sortByType);
  
  /* Create a custom scrollbar for scrolling the sequence column and put it at the bottom of the window */
  GtkWidget *scrollBar = createDetailViewScrollBar(detailAdjustment, detailView);
  gtk_box_pack_start(GTK_BOX(vbox), scrollBar, FALSE, TRUE, 0);
  
  /* Set required data in the main window widget */
  mainWindowCreateProperties(window, 
			     bigPicture, 
			     detailView, 
			     mainmenu,
			     mspList, 
			     blastMode, 
			     refSeq, 
			     refSeqName,
			     &refSeqRange, 
			     &fullDisplayRange, 
			     seqType, 
			     geneticCode,
			     numFrames,
			     gappedHsp);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressMainWindow), mainmenu);
  g_signal_connect(G_OBJECT(window), "key-press-event", G_CALLBACK(onKeyPressMainWindow), mainmenu);
  
  /* Add the MSP's to the trees and sort them by the initial sort mode */
  detailViewAddMspData(detailView, mspList);
  detailViewSortByType(detailView, sortByType);

  /* Initial update to set the detail view font */
  updateDetailViewFontDesc(detailView);

  /* Show the window */
  printf("Starting Blixem\n");
  gtk_widget_show_all(window);
  
  return window;
}
