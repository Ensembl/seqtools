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


typedef struct _GroupDialogData
  {
    GtkWidget *mainWindow;
    GtkWidget *searchTextBox;
    GtkWidget *searchToggleButton;
  } GroupDialogData;

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

static void			  onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
static void			  onDrawPage(GtkPrintOperation *operation, GtkPrintContext *context, gint pageNum, gpointer data);

static Strand			  mainWindowGetActiveStrand(GtkWidget *mainWindow);
static Strand			  mainWindowGetInactiveStrand(GtkWidget *mainWindow);

static GList*			  findSelectedSeqInList(GList *list, const char *seqName);
static gint			  runConfirmationBox(GtkWidget *mainWindow, char *title, char *messageText);
static void			  onButtonClickedDeleteGroup(GtkWidget *button, gpointer data);
static void			  showGroupSequencesDialog(GtkWidget *mainWindow, const gboolean createGroup);

/* Menu builders */
static const GtkActionEntry mainMenuEntries[] = {
  { "Quit",		NULL, "_Quit",		    "<control>Q",	"Quit the program",	  gtk_main_quit},
  { "Help",		NULL, "_Help",		    "<control>H",	"Display Blixem help",	  G_CALLBACK(onHelpMenu)},
  { "Print",		NULL, "_Print",		    "<control>P",	"Print",		  G_CALLBACK(onPrintMenu)},
  { "Settings",		NULL, "_Settings",	    "<control>S",	"Edit Blixem settings",	  G_CALLBACK(onSettingsMenu)},

  { "View",		NULL, "_View",		    "<control>V",	"Edit view settings",	  G_CALLBACK(onViewMenu)},
  { "CreateGroup",	NULL, "Create _Group",	    "<control>G",	"Group sequences together", G_CALLBACK(onCreateGroupMenu)},
  { "EditGroups",	NULL, "Edit Groups",	    "<shift><control>G","Groups",		  G_CALLBACK(onEditGroupsMenu)},
  { "DeselectAllRows",	NULL, "Deselect _all",	    "<shift><control>A","Deselect all",		  G_CALLBACK(onDeselectAllRows)},

  { "Dotter",		NULL, "_Dotter",	    "<control>D",	"Start Dotter",		  G_CALLBACK(onDotterMenu)},
  { "SelectFeatures",	NULL, "Feature series selection tool",	NULL,	"Feature series selection tool", G_CALLBACK(onSelectFeaturesMenu)},

  { "Statistics",	NULL, "Statistics",   NULL,			"Show memory statistics", G_CALLBACK(onStatisticsMenu)}
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
 * free'd with g_free by the caller. The given indices are 1-based. */
static gchar *copySeqSegment(const char const *inputSeq, const int idx1, const int idx2)
{
  const int minIdx = min(idx1, idx2);
  const int maxIdx = max(idx1, idx2);
  
  const int segmentLen = maxIdx - minIdx + 1;
  gchar *segment = g_malloc(sizeof(gchar) * segmentLen + 1);
  
  strncpy(segment, inputSeq + minIdx - 1, segmentLen);
  segment[segmentLen] = '\0';
  
  return segment;
}


/* Get a segment of the given sequence (which is always the DNA sequence and
 * always the forward strand) and reverse/complement it as required for the given
 * strand and reading frame (note that the sequence will only actually be reversed
 * if the 'reverse' argument is true). This function also translates it to a peptide
 * sequence if relevant (if the 'translate' flag allows it). */
gchar *getSequenceSegment(GtkWidget *mainWindow,
			  const char const *sequence,
			  const IntRange const *sequenceRange,
			  const int coord1, 
			  const int coord2,
			  const Strand strand,
			  const BlxSeqType inputSeqType,
			  const int frame,
			  const int numReadingFrames,
			  const gboolean reverse,
			  const gboolean translate)
{
  gchar *result = NULL;
  
  int qMin = min(coord1, coord2); 
  int qMax = max(coord1, coord2);

  /* If the input coords are on a peptide sequence, convert them to DNA sequence coords. */
  if (inputSeqType == BLXSEQ_PEPTIDE)
    {
      qMin = convertPeptideToDna(qMin, frame, 1, numReadingFrames);		   /* 1st base in frame */
      qMax = convertPeptideToDna(qMax, frame, numReadingFrames, numReadingFrames); /* last base in frame */
    }
  
  /* Check that the requested segment is within the sequence's range */
  if (qMin < sequenceRange->min || qMax > sequenceRange->max)
    {
      messout ( "Requested query sequence %d - %d out of available range: %d - %d\n", qMin, qMax, sequenceRange->min, sequenceRange->max);
      if (inputSeqType == BLXSEQ_PEPTIDE) messout("Input coords on peptide sequence were %d - %d", coord1, coord2);
      return NULL;
    }
  
  /* Copy the portion of interest from the reference sequence and translate as necessary */
  const BlxBlastMode mode = mainWindowGetBlastMode(mainWindow);
  
  if (mode == BLXMODE_BLASTP || mode == BLXMODE_TBLASTN)
    {
      /* Just get a straight copy of this segment from the ref seq */
      result = copySeqSegment(sequence, qMin, qMax);
      
      if (reverse)
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
	  segment = copySeqSegment(sequence, qMin, qMax);
	}
      else
	{
	  /* Get the segment of the ref seq and then complement it */
	  segment = copySeqSegment(sequence, qMin, qMax);
	  blxComplement(segment);
	  
	  if (!segment)
	    {
	      messcrash ("Error getting the reference sequence segment: Failed to complement the reference sequence for the range %d - %d.", qMin, qMax);
	    }
	}
      
      if (reverse)
	{
	  g_strreverse(segment);
	}
      
      if (mode == BLXMODE_BLASTN || !translate)
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
	      messcrash ("Error getting the reference sequence segment: Failed to translate the DNA sequence for the reference range %d - %d/\n", coord1, coord2) ;
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
      
      if (moveLeft != properties->strandsToggled)
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
      
      if (moveLeft != properties->strandsToggled)
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


/* Toggle the visibility of tree/grid panes following a number key press */
static void togglePaneVisibility(GtkWidget *mainWindow, const int number, const gboolean modifier)
{
  /* Affects grids if ctrl was pressed, trees otherwise */
  if (modifier)
    toggleGridVisibility(mainWindow, number);
  else
    toggleTreeVisibility(mainWindow, number);
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
static void showViewPanesDialog(GtkWidget *mainWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("View panes", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
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
  createVisibilityButton(bigPictureGetExonView(bigPicture), "_Exon view", bigPictureSubBox);
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
  if (gtk_widget_get_sensitive(GTK_WIDGET(entry)))
    {
      const char *searchStr = gtk_entry_get_text(entry);
      GList *resultList = NULL;

      /* Loop through all the sequences and see if the name matches the search string */
      GHashTable *seqTable = detailViewGetSeqTable(detailView);
      GList *seqNameItem = g_hash_table_get_keys(seqTable);
      
      for ( ; seqNameItem; seqNameItem = seqNameItem->next)
	{
	  char *compareName = (char *)(seqNameItem->data);
	  
	  if (wildcardSearch(compareName, searchStr))
	    {
	      resultList = g_list_append(resultList, compareName);
	    }
	}
      
      /* If we found anything, create a group */
      if (g_list_length(resultList) > 0)
	{
	  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
	  SequenceGroup *group = createSequenceGroup(properties->sequenceGroups);
	  group->seqList = resultList;
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
  showGroupSequencesDialog(mainWindow, FALSE);
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
	  showGroupSequencesDialog(mainWindow, FALSE);
	}
    }
}


/* Shows the "Group sequences" dialog. This dialog allows the user to group sequences together.
 * This tabbed dialog shows both the 'create group' and 'edit groups' dialogs in one. If the
 * 'create group' argument is true, the 'create group' tab is displayed by default; otherwise
 * the 'edit groups' tab is shown. */
static void showGroupSequencesDialog(GtkWidget *mainWindow, const gboolean createGroup)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Groups", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
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
  else if (!createGroup)
    {
      /* User has asked to edit groups but there aren't any. Warn them. */
      messout("Cannot edit groups; none exist. Create some first.");
    }

  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
  
  /* If user has asked to edit groups, make the second tab the default and the 'close'
   * button the default action. (Must do this after showing the child widgets due
   * to a GTK legacy whereby the notebook won't change tabs otherwise.) */
  if (!createGroup)
    {
      gtk_notebook_next_page(GTK_NOTEBOOK(notebook));
      gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_REJECT);
    }
  else
    {
      /* Make the add group button the default widget and focus it if rows
       * are selected - otherwise focus the text entry box. */
      gtk_widget_set_can_default(addGroupButton, TRUE);
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
  GtkWidget *mainWindow = GTK_WIDGET(data);
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
  detailViewSquashMatches(detailView, squash);
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
static void showSettingsDialog(GtkWidget *mainWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Blixem Settings", 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
  GtkWidget *vbox = createVBoxWithBorder(contentArea, borderWidth);
  
  createCheckButton(vbox, "_Squash matches", detailViewGetMatchesSquashed(detailView), G_CALLBACK(onSquashMatches), mainWindow);
    
  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
}


/***********************************************************
 *		    Statistics menu			   *
 ***********************************************************/

static void getStats(GtkWidget *mainWindow, GString *result, MSP *MSPlist)
{
  int numMSPs = 0;          //count of MSPs
  int numSeqs = 0;          //count of sequences (excluding duplicates)
  gint32 totalSeqSize = 0;  //total memory used by the sequences
  GHashTable *sequences = g_hash_table_new(NULL, NULL); //remembers which sequences we've already seen
  
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
  /* Create a modal dialog widget with an OK button */
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Statistics", 
                                                  GTK_WINDOW(mainWindow),
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Ensure that the dialog box (along with any children) is destroyed when the user responds. */
  g_signal_connect (dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  /* Create a text buffer containing the required text*/
  GString *displayText = g_string_sized_new(200); //will be extended if we need more space
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

/* Utility to pop up a simple modal dialog with the given title and text, with just an "OK" button. */
static void showModalDialog(GtkWidget *mainWindow, char *title, char *messageText)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(mainWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Add the message */
  GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
  GtkWidget *label = gtk_label_new(messageText);
  gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
  
  /* Ensure dialog is destroyed when user responds */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  gtk_widget_show_all(dialog);
}


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
				  \
				  BLIXEM - BLast matches\n\
				  In an\n\
				  X-windows\n\
				  Embedded\n\
				  Multiple alignment\n\
				  \n\
				  LEFT MOUSE BUTTON:\n\
				  Pick on boxes and sequences.\n\
				  Fetch annotation by double clicking on sequence (Requires 'efetch' to be installed.)\n\
				  \n\
				  MIDDLE MOUSE BUTTON:\n\
				  Scroll horizontally.\n\
				  \n\
				  RIGHT MOUSE BUTTON:\n\
				  Menu.  Note that the buttons Settings and Goto have their own menus.\n\
				  \n\
				  \n\
				  Keyboard shortcuts:\n\
				  \n\
				  Cntl-Q          quit application\n\
				  Cntl-P          print\n\
				  Cntl-H          help\n\
				  \n\
				  \n\
				  m/M             for mark/unmark a set of matches from the cut buffer\n\
				  \n\
				  g/G             go to the match in the cut buffer\n\
				  \n\
				  \n\
				  \n\
				  RESIDUE COLOURS:\n\
				  Yellow = Query.\n\
				  See Settings Panel for matching residues (click on Settings button).\n\
				  \n\
				  version %s\n\
				  (c) Erik Sonnhammer", blixemVersion));
  
  showModalDialog(mainWindow, "Help", messageText);
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
  showGroupSequencesDialog(mainWindow, TRUE);
}


/* Called when the user selects the 'Groups' menu option, or hits the relevant shortcut key */
static void onEditGroupsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  showGroupSequencesDialog(mainWindow, FALSE);
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
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  
  /* Create a print operation, using the same settings as the last print, if there was one */
  GtkPrintOperation *print = gtk_print_operation_new();
  
  if (properties->printSettings != NULL)
    {
      gtk_print_operation_set_print_settings(print, properties->printSettings);
    }
  
  g_signal_connect (print, "begin_print", G_CALLBACK (onBeginPrint), mainWindow);
  g_signal_connect(G_OBJECT(print), "draw-page", G_CALLBACK(onDrawPage), mainWindow);
  
  /* Pop up the print dialog */
  GtkPrintOperationResult printResult = gtk_print_operation_run (print, 
								 GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
								 GTK_WINDOW(mainWindow),
								 NULL);
  
  /* If the user hit ok, remember the print settings for next time */
  if (printResult == GTK_PRINT_OPERATION_RESULT_APPLY)
    {
      if (properties->printSettings != NULL)
	{
	  g_object_unref(properties->printSettings);
	}
      
      properties->printSettings = g_object_ref(gtk_print_operation_get_print_settings(print));
    }

  g_object_unref(print);
}


/***********************************************************
 *			   Events                          *
 ***********************************************************/

/* Called after the user clicks ok in the print dialog. For now just scales the 
 * whole output to fit on a single page. */
static void onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);

  /* Set the orientation based on the stored print settings (not sure why this doesn't
   * already get set when the settings are set in the print operation...). */
  GtkPrintSettings *printSettings = gtk_print_operation_get_print_settings(print);
  gtk_print_settings_set_orientation(printSettings, gtk_print_settings_get_orientation(properties->printSettings));

//  //gdouble scale = gtk_print_settings_get_scale(printSettings);
//  
//  gdouble windowWidth = (mainWindow->allocation.width);
//  gdouble windowHeight = (mainWindow->allocation.height);
//
//  gboolean landscape = 0;//gtk_print_settings_get_orientation(printSettings) == GTK_PAGE_ORIENTATION_LANDSCAPE;
//  gdouble pageWidth = landscape ? gtk_print_context_get_height(context) : gtk_print_context_get_width(context);
//  gdouble pageHeight = landscape ? gtk_print_context_get_width(context) : gtk_print_context_get_height(context);
//
//  gdouble hScale = (pageWidth * 100) / windowWidth;
//  gdouble vScale = (pageHeight * 100) / windowHeight;
//  gdouble scale = min(hScale, vScale);
//  scale = min(scale, 100); /* don't print bigger than 100% */
//  gtk_print_settings_set_scale(printSettings, scale);
//  gtk_print_operation_set_print_settings(print, printSettings);

//  int numWide = ceil(windowWidth / pageWidth);
//  int numHigh = ceil(windowHeight / pageHeight);

  gtk_print_operation_set_n_pages(print, 1);
}


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
static void onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  
  /* Create a blank white pixmap to draw on to */
  properties->drawable = gdk_pixmap_new(mainWindow->window, mainWindow->allocation.width, mainWindow->allocation.height, -1);
  
  GdkGC *gc = gdk_gc_new(properties->drawable);
  GdkColor fgColour = getGdkColor(GDK_WHITE);
  gdk_gc_set_foreground(gc, &fgColour);
  gdk_draw_rectangle(properties->drawable, gc, TRUE, 0, 0, mainWindow->allocation.width, mainWindow->allocation.height);

  /* For each child widget that has a drawable set, draw this onto the main pixmap */
  properties->lastYStart = 0;
  properties->lastYEnd = 0;
  properties->lastYCoord = -1;
  gtk_container_foreach(GTK_CONTAINER(mainWindow), collatePixmaps, mainWindow);
  
  cairo_t *cr = gtk_print_context_get_cairo_context (context);
  gdk_cairo_set_source_pixmap(cr, properties->drawable, 0, 0);
  cairo_paint(cr);
}


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
  
//  guint modifiers = gtk_accelerator_get_default_mod_mask();
  const gboolean ctrlModifier = (event->state & GDK_CONTROL_MASK);
  
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
	
      case GDK_g: /* fall through */
      case GDK_G:
	if (!ctrlModifier)
	  {
	    goToDetailViewCoord(mainWindowGetDetailView(window), BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
	    result = TRUE;
	  }
	break;
	
      case GDK_t:
      case GDK_T:
	toggleStrand(mainWindowGetDetailView(window));
	result = TRUE;
	break;
	
      case GDK_1:
	togglePaneVisibility(window, 1, ctrlModifier);
	result = TRUE;
	break;

      case GDK_2:
	togglePaneVisibility(window, 2, ctrlModifier);
	result = TRUE;
	break;
	
      case GDK_3:
	togglePaneVisibility(window, 3, ctrlModifier);
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
				       char *displaySeq,
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
      properties->displaySeq = displaySeq;
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
      properties->printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(properties->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
      gtk_print_settings_set_quality(properties->printSettings, GTK_PRINT_QUALITY_HIGH);
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

char* mainWindowGetDisplaySeq(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->displaySeq : NULL;
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

/* Update function to be called whenever the MSP selection has changed */
void mainWindowSelectionChanged(GtkWidget *mainWindow, const gboolean updateTrees)
{
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
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
			    int numReadingFrames,
			    char **geneticCode,
			    const int refSeqOffset,
			    const gboolean gappedHsp)
{
  /* Get the range of the reference sequence. If this is a DNA sequence but our
   * matches are peptide sequences, we must convert to the peptide sequence. */
  const int refSeqLen = (int)strlen(refSeq);
  IntRange refSeqRange = {refSeqOffset + 1, refSeqOffset + refSeqLen};
  IntRange fullDisplayRange = {refSeqRange.min, refSeqRange.max};
  char *displaySeq = refSeq;
  
  if (seqType == BLXSEQ_PEPTIDE)
    {
      displaySeq = blxTranslate(refSeq, geneticCode);
      fullDisplayRange.min = 1;
      fullDisplayRange.max = strlen(displaySeq);
      printf("Converted DNA sequence (len=%d) to peptide sequence (len=%d).\n",  refSeqLen, fullDisplayRange.max);
    }
  
  /* Create the main window */
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setStyleProperties(window);
  
  /* Create a vertical box to pack everything in */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), vbox);
  
  /* Create the main menu */
  GtkWidget *mainmenu = createMainMenu(window);
  
//  /* Create a container for our two main frames in */
//  GtkWidget *panedWidget = gtk_vpaned_new();
//  gtk_box_pack_start(GTK_BOX(vbox), panedWidget, TRUE, TRUE, 0);
  
  /* Create the widgets. We need a single adjustment for the entire detail view, which will also be referenced
   * by the big picture view, so create it first. */
  GtkAdjustment *detailAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, /* initial value = 0 */
								      0, /* lower = 0 */
								      fullDisplayRange.max - fullDisplayRange.min + 1, /* upper = lenght of ref seq */
								      DEFAULT_SCROLL_STEP_INCREMENT, /* step increment used for mouse wheel scrolling */
								      0,   /* page increment dynamically set based on display range */
								      0)); /* page size dunamically set based on display range */
  
  GtkWidget *fwdStrandGrid = NULL, *revStrandGrid = NULL;
  
  printf("Creating big picture...\n");
  GtkWidget *bigPicture = createBigPicture(window,
					   vbox, 
					   &fwdStrandGrid, 
					   &revStrandGrid, 
					   &fullDisplayRange);
  printf("Done.\n");
  
  printf("Creating detail view...\n");
  GtkWidget *detailView = createDetailView(window,
					   vbox, 
					   detailAdjustment, 
					   fwdStrandGrid, 
					   revStrandGrid,
					   mspList,
					   blastMode,
					   seqType,
					   numReadingFrames,
					   refSeqName);
  printf("Done.\n");
  
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
			     displaySeq,
			     &refSeqRange, 
			     &fullDisplayRange, 
			     seqType, 
			     geneticCode,
			     numReadingFrames,
			     gappedHsp);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressMainWindow), mainmenu);
  g_signal_connect(G_OBJECT(window), "key-press-event", G_CALLBACK(onKeyPressMainWindow), mainmenu);
  
  /* Add the MSP's to the trees */
  detailViewAddMspData(detailView, mspList);

  /* Initial update to set the detail view font */
  updateDetailViewFontDesc(detailView);

  /* Show the window */
  printf("realizing widgets...\n");
  gtk_widget_show_all(window);
  printf("Done.\n");
  
  return window;
}
