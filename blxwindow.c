/*
 *  blxWindow.c
 *  acedb
 *
 *  Created by Gemma Barson on 24/11/2009.
 *
 */

#include <SeqTools/blxwindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/blxdotter.h>
#include <SeqTools/exonview.h>
#include <SeqTools/utilities.h>
#include <gdk/gdkkeysyms.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      10   /* used to change the default border width around the blixem window */
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2   /* used to start with a smaller font than the default widget font */
#define DEFAULT_SCROLL_STEP_INCREMENT	 5    /* how many bases the scrollbar scrolls by for each increment */
#define DEFAULT_WINDOW_WIDTH_FRACTION	 0.9  /* what fraction of the screen size the blixem window width defaults to */
#define DEFAULT_WINDOW_HEIGHT_FRACTION	 0.6  /* what fraction of the screen size the blixem window height defaults to */


/* Utility struct used when comparing sequence names */
typedef struct _CompareSeqData
  {
    const char *searchStr;
    GList *matchList;
  } SeqSearchData;


/* Properties specific to the blixem window */
typedef struct _BlxWindowProperties
  {
    GtkWidget *bigPicture;	    /* The top section of the view, showing a "big picture" overview of the alignments */
    GtkWidget *detailView;	    /* The bottom section of the view, showing a detailed list of the alignments */
    GtkWidget *mainmenu;	      /* The main menu */

    BlxViewContext *blxContext;	      /* The blixem view context */

    GtkPrintSettings *printSettings;  /* Used so that we can re-use the same print settings as a previous print */
    int lastYEnd;		      /* Keeps track of where the last item ended so we can draw the next one flush to it */
    int lastYStart;		      /* Where the last item started (for drawing multiple items at same y pos) */
    int lastYCoord;		      /* Y coord of last item (so we can check if current item should be at same Y pos) */
    
  } BlxWindowProperties;



/* Local function declarations */
static BlxWindowProperties*	  blxWindowGetProperties(GtkWidget *widget);

static void			  onHelpMenu(GtkAction *action, gpointer data);
static void			  onQuit(GtkAction *action, gpointer data);
static void			  onPrintMenu(GtkAction *action, gpointer data);
static void			  onSettingsMenu(GtkAction *action, gpointer data);
static void			  onViewMenu(GtkAction *action, gpointer data);
static void			  onCreateGroupMenu(GtkAction *action, gpointer data);
static void			  onEditGroupsMenu(GtkAction *action, gpointer data);
static void			  onToggleMatchSet(GtkAction *action, gpointer data);
static void			  onDotterMenu(GtkAction *action, gpointer data);
static void			  onSelectFeaturesMenu(GtkAction *action, gpointer data);
static void			  onDeselectAllRows(GtkAction *action, gpointer data);
static void			  onStatisticsMenu(GtkAction *action, gpointer data);

static gboolean			  onKeyPressBlxWindow(GtkWidget *window, GdkEventKey *event, gpointer data);

static void			  onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
static void			  onDrawPage(GtkPrintOperation *operation, GtkPrintContext *context, gint pageNum, gpointer data);
static void			  onDestroyBlxWindow(GtkWidget *widget);

static Strand			  blxWindowGetInactiveStrand(GtkWidget *blxWindow);

static void			  onButtonClickedDeleteGroup(GtkWidget *button, gpointer data);
static void			  blxWindowGroupsChanged(GtkWidget *blxWindow);
static GtkRadioButton*		  createRadioButton(GtkBox *box, GtkRadioButton *existingButton, const char *mnemonic, const gboolean isActive, const gboolean createTextEntry, const gboolean multiline, GtkCallback callbackFunc, GtkWidget *blxWindow);
static void			  getSequencesThatMatch(gpointer listDataItem, gpointer data);
static GList*			  getSeqStructsFromText(GtkWidget *blxWindow, const char *inputText);

static void			  createCheckButton(GtkBox *box, const char *mnemonic, const gboolean isActive, GCallback callback, GtkWidget *blxWindow);
			      
/* Menu builders */
static const GtkActionEntry mainMenuEntries[] = {
  { "Quit",		NULL, "_Quit\t\t\t\tCtrl-Q",		  "<control>Q",		"Quit the program",		    G_CALLBACK(onQuit)},
  { "Help",		NULL, "_Help\t\t\t\tCtrl-H",		  "<control>H",		"Display Blixem help",		    G_CALLBACK(onHelpMenu)},
  { "Print",		NULL, "_Print\t\t\t\tCtrl-P",		  "<control>P",		"Print",			    G_CALLBACK(onPrintMenu)},
  { "Settings",		NULL, "_Settings\t\t\t\tCtrl-S",	  "<control>S",		"Edit Blixem settings",		    G_CALLBACK(onSettingsMenu)},

  { "View",		NULL, "_View\t\t\t\tv",			  "V",			"Edit view settings",		    G_CALLBACK(onViewMenu)},
  { "CreateGroup",	NULL, "Create Group\t\tShift-Ctrl-G",	  "<Shift><control>G",  "Group sequences together",	    G_CALLBACK(onCreateGroupMenu)},
  { "EditGroups",	NULL, "Edit _Groups\t\t\tCtrl-G",	  "<control>G",		"Groups",			    G_CALLBACK(onEditGroupsMenu)},
  { "ToggleMatchSet",	NULL, "Toggle _match set group\t\tg",     "G",			"Create/clear the match set group", G_CALLBACK(onToggleMatchSet)},
  { "DeselectAllRows",	NULL, "Deselect _all\t\t\tShift-Ctrl-A",  "<Shift><control>A",	"Deselect all",			    G_CALLBACK(onDeselectAllRows)},

  { "Dotter",		NULL, "_Dotter\t\t\t\tCtrl-D",		  "<control>D",		"Start Dotter",			    G_CALLBACK(onDotterMenu)},
  { "SelectFeatures",	NULL, "Feature series selection tool",	  NULL,			"Feature series selection tool",    G_CALLBACK(onSelectFeaturesMenu)},

  { "Statistics",	NULL, "Statistics",   NULL,					"Show memory statistics",	    G_CALLBACK(onStatisticsMenu)}
};


/* This defines the layout of the menu for a standard user */
static const char standardMenuDescription[] =
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
"      <menuitem action='ToggleMatchSet'/>"
"      <menuitem action='DeselectAllRows'/>"
"      <separator/>"
"      <menuitem action='Dotter'/>"
"      <menuitem action='SelectFeatures'/>"
"  </popup>"
"</ui>";


/* This defines the layout of the menu for a developer user */
static const char developerMenuDescription[] =
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
"      <menuitem action='ToggleMatchSet'/>"
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


/* Copy a segment of the given sequence (which is always the DNA sequence and
 * always the forward strand).
 *
 * The result is complemented if the reverse strand is requested, but only if
 * the allowComplement flag allows it.
 * 
 * The result is translated to a peptide sequence if the blast mode is protein 
 * matches, but only if the allowTranslate flag allows it.
 *
 * The result is reversed if the reverseResult flag is true (regardless of the
 * strand requested - this is because the caller often wants the result in the
 * opposite direction to that indicated by the strand, because the display may be
 * reversed, so we leave it up to the caller to decide).
 */
gchar *getSequenceSegment(BlxViewContext *bc,
			  const char const *dnaSequence,
			  const int coord1, 
			  const int coord2,
			  const Strand strand,
			  const BlxSeqType inputCoordType,
			  const int frame,
			  const gboolean displayRev,
			  const gboolean reverseResult,
			  const gboolean allowComplement,
			  const gboolean allowTranslate)
{
  gchar *result = NULL;

  /* Convert input coord to ref seq coords and find the min/max */
  const int qIdx1 = convertDisplayIdxToDnaIdx(coord1, inputCoordType, frame, 1, bc->numFrames, displayRev, &bc->refSeqRange);		  /* 1st base in frame */
  const int qIdx2 = convertDisplayIdxToDnaIdx(coord2, inputCoordType, frame, bc->numFrames, bc->numFrames, displayRev, &bc->refSeqRange); /* last base in frame */
  int qMin = min(qIdx1, qIdx2);
  int qMax = max(qIdx1, qIdx2);
  
  /* Check that the requested segment is within the sequence's range */
  if (qMin < bc->refSeqRange.min || qMax > bc->refSeqRange.max)
    {
      /* We might request up to 3 bases beyond the end of the range if we want to 
       * show a partial triplet at the start/end. (It's a bit tricky for the caller to
       * specify the exact bases they want here so we allow it but just limit it to 
       * the actual range so that they can't index beyond the end of the range.) Any
       * further out than one triplet is probably indicative of an error, so give a warning. */
      if (qMin < bc->refSeqRange.min - (bc->numFrames + 1) || qMax > bc->refSeqRange.max + (bc->numFrames + 1))
	{
	  if (inputCoordType == BLXSEQ_PEPTIDE)
	    printf ( "Requested query sequence %d - %d out of available range: %d - %d. Input coords on peptide sequence were %d - %d\n", qMin, qMax, bc->refSeqRange.min, bc->refSeqRange.max, coord1, coord2);
	  else
	    printf ( "Requested query sequence %d - %d out of available range: %d - %d\n", qMin, qMax, bc->refSeqRange.min, bc->refSeqRange.max);
	}
      
      if (qMax > bc->refSeqRange.max)
	qMax = bc->refSeqRange.max;
	
      if (qMin < bc->refSeqRange.min)
	qMin = bc->refSeqRange.min;
    }
  
  /* Get 0-based indices into the sequence */
  const int idx1 = qMin - bc->refSeqRange.min;
  const int idx2 = qMax - bc->refSeqRange.min;
  
  /* Copy the portion of interest from the reference sequence and translate as necessary */
  const BlxBlastMode mode = bc->blastMode;
  
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
	  
	  if (allowComplement)
	    {
	      blxComplement(segment);
	  
	      if (!segment)
		{
		  messcrash ("Error getting the reference sequence segment: Failed to complement the reference sequence for the range %d - %d.", qMin, qMax);
		}
	    }
	}
      
      if (reverseResult)
	{
	  g_strreverse(segment);
	}
      
      if (mode == BLXMODE_BLASTN || !allowTranslate)
	{
	  /* Just return the segment of DNA sequence */
	  result = segment;
	}
      else
	{
	  /* Translate the DNA sequence to a peptide sequence */
	  result = blxTranslate(segment, bc->geneticCode);
	  
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
int sequenceGetGroupOrder(GtkWidget *blxWindow, const SequenceStruct *seq)
{
  SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, seq);
  return group ? group->order : UNSET_INT;
}


/* Scroll the detail view left/right by 1 base (or by 1 page, if the modifier
 * is pressed) */
static void scrollDetailView(GtkWidget *window, const gboolean moveLeft, const gboolean modifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(window);
  
  if (moveLeft && modifier)
    scrollDetailViewLeftPage(detailView);
  else if (moveLeft)
    scrollDetailViewLeft1(detailView);
  else if (modifier)
    scrollDetailViewRightPage(detailView);
  else
    scrollDetailViewRight1(detailView);
}


/* Move the selected base index 1 base to the left/right. Moves by individual
 * DNA bases (i.e. you have to move 3 bases in order to scroll a full peptide
 * if viewing protein matches). Scrolls the detail view if necessary to keep 
 * the new base in view. */
static void moveSelectedBaseIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  BlxViewContext *blxContext = blxWindowGetContext(window);
  GtkWidget *detailView = blxWindowGetDetailView(window);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  if (detailViewProperties->selectedBaseIdx != UNSET_INT)
    {
      /* Decrement the index if moving left decrease and increment it if moving right, 
       * unless the display is toggled, in which case do the opposite */
      const int numFrames = blxContext->numFrames;
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

      IntRange *fullRange = blxWindowGetFullRange(window);
      boundsLimitValue(&newSelectedBaseIdx, fullRange);

      detailViewSetSelectedBaseIdx(detailView, newSelectedBaseIdx, detailViewProperties->selectedFrame, newBaseNum, TRUE, TRUE);
    }
}


/* Called when user pressed Home/End. If the modifier is pressed, scroll to the
*  start/end of all matches in the current selection (or all matches, if no 
* selection), or to the start/end of the entire display if the modifier is not pressed. */
static void scrollToExtremity(GtkWidget *blxWindow, const gboolean moveLeft, const gboolean modifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  if (modifier)
    {
      GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);

      if (moveLeft)
	firstMatch(detailView, selectedSeqs);
      else
	lastMatch(detailView, selectedSeqs);
    }
  else
    {
      const BlxSeqType seqType = blxWindowGetSeqType(blxWindow);
      const IntRange const *fullRange = blxWindowGetFullRange(blxWindow);

      if (moveLeft)
	setDetailViewStartIdx(detailView, fullRange->min, seqType);
      else
	setDetailViewEndIdx(detailView, fullRange->max, seqType);
    }
}


/* Jump left or right to the next/prev nearest match. Only include matches in the
 * current selection, if any rows are selected. */
static void goToMatch(GtkWidget *blxWindow, const gboolean moveLeft)
{
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (moveLeft)
    {
      prevMatch(blxWindowGetDetailView(blxWindow), selectedSeqs);
    }
  else
    {
      nextMatch(blxWindowGetDetailView(blxWindow), selectedSeqs);      
    }
}


/* Move the selected display index 1 value to the left/right. Moves by full peptides
 * if viewing protein matches. Scrolls the detail view if necessary to keep the new 
 * index in view. */
static void moveSelectedDisplayIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  GtkWidget *detailView = blxWindowGetDetailView(window);
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
      
      IntRange *fullRange = blxWindowGetFullRange(window);
      boundsLimitValue(&newSelectedBaseIdx, fullRange);
      
      detailViewSetSelectedBaseIdx(detailView, newSelectedBaseIdx, detailViewProperties->selectedFrame, detailViewProperties->selectedBaseNum, TRUE, TRUE);
    }
}


/* Zooms the display in/out. The modifiers control which section is zoomed */
static void zoomBlxWindow(GtkWidget *window, const gboolean zoomIn, const gboolean ctrl, const gboolean shift)
{
  if (ctrl)
    {
      if (shift)
	zoomWholeBigPicture(blxWindowGetBigPicture(window));
      else
	zoomBigPicture(blxWindowGetBigPicture(window), zoomIn);
    }
  else
    {
      zoomDetailView(blxWindowGetDetailView(window), zoomIn);
    }
}


/* Force a redraw of all widgets. Clears cached bitmaps etc. first */
void blxWindowRedrawAll(GtkWidget *blxWindow)
{
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  bigPictureRedrawAll(bigPicture);

  gtk_widget_queue_draw(blxWindow);
}


/* Utility to create a vbox with the given border and pack it into the given container.
 * Also put a frame around it with the given label if includeFrame is true */
static GtkWidget* createVBoxWithBorder(GtkWidget *parent, 
				       const int borderWidth,
				       const gboolean includeFrame,
				       const char *frameTitle)
{
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), borderWidth);

  if (includeFrame)
    {
      GtkWidget *frame = gtk_frame_new(frameTitle);
      gtk_container_add(GTK_CONTAINER(parent), frame);
      gtk_container_add(GTK_CONTAINER(frame), vbox);
    }
  else
    {
      gtk_container_add(GTK_CONTAINER(parent), vbox);
    }
  
  return vbox;
}

/* Utility to create an hbox with the given border and pack it into the given container */
static GtkWidget* createHBoxWithBorder(GtkWidget *parent, const int borderWidth)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(hbox), borderWidth);
  gtk_container_add(GTK_CONTAINER(parent), hbox);
  return hbox;
}


/* Utility to return true if any groups exist. Ignores the 'match set' group
 * if it doesn't have any sequences. */
static gboolean blxWindowGroupsExist(GtkWidget *blxWindow)
{
  gboolean result = FALSE;
  
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  GList *groupList = blxContext->sequenceGroups;
  
  if (g_list_length(groupList) > 1)
    {
      result = TRUE;
    }
  else if (g_list_length(groupList) == 1)
    {
      /* Only one group. If it's the match set group, check it has sequences */
      SequenceGroup *group = (SequenceGroup*)(groupList->data);
      
      if (group != blxContext->matchSetGroup || g_list_length(group->seqList) > 0)
	{
	  result = TRUE;
	}
    }
  
  return result;
}


/***********************************************************
 *			   View panes menu                 *
 ***********************************************************/

/* Toggle visibility the n'th tree. This is the active strand's frame n if displaying
 * protein matches (where we only display one strand), or the forward or reverse
 * strand tree if displaying DNA matches (where both strands are displayed). */
static void toggleTreeVisibility(GtkWidget *blxWindow, const int number)
{
  const gboolean toggled = blxWindowGetDisplayRev(blxWindow);
  const Strand activeStrand = toggled ? REVERSE_STRAND : FORWARD_STRAND;
  
  /* For protein matches, trees are always displayed in frame order (i.e. 1, 2, 3), 
   * so just use the number pressed for the frame, and the active strand for the
   * strand. */
  int frame = number;
  Strand strand = activeStrand;
  
  /* For DNA matches, the frame is always 1, but the strand depends on which number
   * was pressed: use 1 to toggle active strand, 2 for other strand */
  if (blxWindowGetSeqType(blxWindow) == BLXSEQ_DNA)
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
  
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (tree && gtk_widget_get_parent(tree))
    {
      widgetSetHidden(tree, !widgetGetHidden(tree));
    }
}


/* Toggle visibility of the active (1) or other (2) strand grid depending on the number pressed */
static void toggleGridVisibility(GtkWidget *blxWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
      const gboolean useFwdGrid = (number == 1) != blxWindowGetDisplayRev(blxWindow);
      
      GtkWidget *grid = useFwdGrid ? bigPictureGetFwdGrid(bigPicture) : bigPictureGetRevGrid(bigPicture);
      widgetSetHidden(grid, !widgetGetHidden(grid));
    }
}


/* Toggle visibility of the active (1) or other (2) strand exon view depending on the number pressed */
static void toggleExonViewVisibility(GtkWidget *blxWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
      const gboolean useFwdExonView = (number == 1) != blxWindowGetDisplayRev(blxWindow);
      
      GtkWidget *exonView = useFwdExonView ? bigPictureGetFwdExonView(bigPicture) : bigPictureGetRevExonView(bigPicture);
      widgetSetHidden(exonView, !widgetGetHidden(exonView));
    }
}


/* Toggle the visibility of tree/grid panes following a number key press */
static void togglePaneVisibility(GtkWidget *blxWindow, const int number, const gboolean modifier1, const gboolean modifier2)
{
  /* Affects big picture if modifier1 was pressed, the detail view otherwise */
  if (modifier1)
    {
      /* If modifier 2 was also pressed, affects the exon views; otherwise the grids */
      if (modifier2)
	{
	  toggleExonViewVisibility(blxWindow, number);
	}
      else
	{
	  toggleGridVisibility(blxWindow, number);
	}
    }
  else
    {
      toggleTreeVisibility(blxWindow, number);
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
  /* Some trees may have been removed from the blixem window if they are not on the active 
   * strand, so only show check boxes for those that are in the window (i.e. have a parent). 
   * Note that we get the tree container here, which consists of the tree itself plus any headers etc. */
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (gtk_widget_get_parent(tree))
    {
      const gboolean toggled = detailViewGetDisplayRev(detailView);
      gboolean isActiveStrand = ((strand == FORWARD_STRAND) != toggled);

      if (detailViewGetSeqType(detailView) == BLXSEQ_DNA)
	{
	  /* We only have 1 frame, but trees are from both strands, so distinguish between strands.
	   * Put each strand in its own frame. */
	  char text1[] = "Show _active strand";
	  char text2[] = "Show othe_r strand";

	  GtkWidget *frame = gtk_frame_new(isActiveStrand ? "Active strand" : "Other strand");
  	  gtk_container_add(GTK_CONTAINER(container), frame);
	  createVisibilityButton(tree, isActiveStrand ? text1 : text2, frame);
	}
      else
	{
	  /* All the visible trees should be in the same strand, so just distinguish by frame number.
	   * Put the strand in its own frame for consistency. */
	  char formatStr[] = "Show frame _%d";
	  char displayText[strlen(formatStr) + numDigitsInInt(frame) + 1];
	  sprintf(displayText, formatStr, frame);

	  GtkWidget *frame = gtk_frame_new(isActiveStrand ? "Active strand" : "Other strand");
	  gtk_container_add(GTK_CONTAINER(container), frame);
	  createVisibilityButton(tree, displayText, frame);
	}
    }
}


/* Callback called when the user clicks the 'bump exon view' button */
static void onBumpExonView(GtkWidget *button, gpointer data)
{
  GtkWidget *exonView = GTK_WIDGET(data);
  const gboolean expanded = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  exonViewSetExpanded(exonView, expanded);
}


/* Create the set of settings buttons to control display of an exon-view widget. */
static void createExonButtons(GtkWidget *exonView, const char *visLabel, const char *bumpLabel, GtkWidget *parent)
{
  /* Pack everything in an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(parent), hbox);
  
  /* Create a check button to control visibility of the exon view */
  createVisibilityButton(exonView, visLabel, hbox);

  /* Create a check button to control whether the exon view is expanded or compressed */
  const gboolean isBumped = exonViewGetExpanded(exonView);
  createCheckButton(GTK_BOX(hbox), bumpLabel, isBumped, G_CALLBACK(onBumpExonView), exonView);
}


/* Shows the "View panes" dialog. This dialog allows the user to show/hide certain portions of the window. */
void showViewPanesDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("View panes", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  
  /* Big picture */
  GtkWidget *bp = blxWindowGetBigPicture(blxWindow);
  GtkWidget *bpVbox = createVBoxWithBorder(contentArea, borderWidth, TRUE, "Big picture");
  
  createVisibilityButton(bp, "Show _big picture", bpVbox);
  GtkWidget *bpSubBox = createVBoxWithBorder(bpVbox, borderWidth, FALSE, NULL);
  
  GtkWidget *bpActiveStrand = createVBoxWithBorder(bpSubBox, 0, TRUE, "Active strand");
  createVisibilityButton(bigPictureGetActiveGrid(bp), "Show _grid", bpActiveStrand);
  createExonButtons(bigPictureGetActiveExonView(bp), "Show _exons    ", "_Bump exons    ", bpActiveStrand);
  
  GtkWidget *bpOtherStrand = createVBoxWithBorder(bpSubBox, 0, TRUE, "Other strand");
  createVisibilityButton(bigPictureGetInactiveGrid(bp), "Show gr_id", bpOtherStrand);
  createExonButtons(bigPictureGetInactiveExonView(bp), "Show e_xons    ", "Bum_p exons    ", bpOtherStrand);
  
  /* Detail view */
  GtkWidget *dvVbox = createVBoxWithBorder(contentArea, borderWidth, TRUE, "Alignment lists");
  GtkWidget *dv = blxWindowGetDetailView(blxWindow);
  createVisibilityButton(dv, "Show alignment _lists", dvVbox);
  
  GtkWidget *dvSubBox = createVBoxWithBorder(dvVbox, borderWidth, FALSE, NULL);
  const int numFrames = blxWindowGetNumFrames(blxWindow);
  int frame = 1;
  for ( ; frame <= numFrames; ++frame)
    {
      createTreeVisibilityButton(dv, blxWindowGetActiveStrand(blxWindow), frame, dvSubBox);
      createTreeVisibilityButton(dv, blxWindowGetInactiveStrand(blxWindow), frame, dvSubBox);
    }
  
  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
}



/***********************************************************
 *			    Find menu			   *
 ***********************************************************/

static GList* findSeqsFromName(GtkWidget *button, gpointer data)
{
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return NULL;
    }
  
  /* The text entry box was passed as the user data */
  GtkEntry *entry = GTK_ENTRY(data);
  
  if (!entry || !GTK_WIDGET_SENSITIVE(GTK_WIDGET(entry)))
    {
      messout("Could not set search string: invalid text entry box [%x]\n", entry);
      return NULL;
    }
  
  /* Extract the main blixem window from our parent window. */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
    
  /* Loop through all the sequences and see if the name matches the search string */
  const char *searchStr = gtk_entry_get_text(entry);
  GList *seqList = blxWindowGetAllMatchSeqs(blxWindow);
  SeqSearchData searchData = {searchStr, NULL};
  
  g_list_foreach(seqList, getSequencesThatMatch, &searchData);
  
  if (g_list_length(searchData.matchList) < 1)
    {
      messout("No sequences found matching text %s", searchStr);
    }
  
  return searchData.matchList;
}


/* Finds all the valid sequences blixem knows about whose names are in the given
 * text entry box (passed as the user data), and returns them in a GList of
 * SequenceStructs. Only does anything if the given radio button is active. */
static GList* findSeqsFromList(GtkWidget *button, gpointer data)
{
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return NULL;
    }
  
  /* The text entry box was passed as the user data. We should have a (multi-line) text view */
  GtkTextView *textView = GTK_TEXT_VIEW(data);
  
  if (!textView || !GTK_WIDGET_SENSITIVE(GTK_WIDGET(textView)))
    {
      messout("Could not set search string: invalid text entry box [%x]\n", textView);
      return NULL;
    }
  
  /* Get the input text from the text buffer and create the group */
  GtkTextBuffer *textBuffer = gtk_text_view_get_buffer(textView);
  
  GtkTextIter start, end;
  gtk_text_buffer_get_bounds(textBuffer, &start, &end);
  
  const char *inputText = gtk_text_buffer_get_text(textBuffer, &start, &end, TRUE);
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  GList *seqList = getSeqStructsFromText(blxWindow, inputText);
  
  if (g_list_length(seqList) < 1)
    {
      messout("No valid sequence names in buffer\n%s\n", inputText);
    }
    
  return seqList;
}


/* Callback called when requested to find sequences from a sequence name. Selects
 * the sequences and scrolls to the start of the first match in the selection */
static void onFindSeqsFromName(GtkWidget *button, gpointer data)
{
  GList *seqList = findSeqsFromName(button, data);
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      blxWindowSetSelectedSeqList(blxWindow, seqList);
      firstMatch(blxWindowGetDetailView(blxWindow), seqList);
    }
}


/* Callback called when requested to find sequences from a given list. Selects
 * the sequences ands scrolls to the start of the first match in the selection. */
static void onFindSeqsFromList(GtkWidget *button, gpointer data)
{
  GList *seqList = findSeqsFromList(button, data);
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
      
      blxWindowSetSelectedSeqList(blxWindow, seqList);
      firstMatch(blxWindowGetDetailView(blxWindow), seqList);
    }
}


void showFindDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Find sequences", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  GtkBox *contentArea = GTK_BOX(GTK_DIALOG(dialog)->vbox);
  
  GtkRadioButton *button1 = createRadioButton(contentArea, NULL, "From _name (wildcards * and ?)", TRUE, TRUE, FALSE, onFindSeqsFromName, blxWindow);
  createRadioButton(contentArea, button1, "From _list", FALSE, TRUE, TRUE, onFindSeqsFromList, blxWindow);
  
  /* Connect signals and show */
  gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(blxWindow));
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  
  gtk_widget_show_all(dialog);
}


/* Toggle the bump state. Currently only the exon view can be bumped, so this just
 * affects that. */
static void toggleBumpState(GtkWidget *blxWindow)
{
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  exonViewToggleExpanded(bigPictureGetFwdExonView(bigPicture));
  exonViewToggleExpanded(bigPictureGetRevExonView(bigPicture));
}

/***********************************************************
 *		      Group sequences menu                 *
 ***********************************************************/

/* Utility to free the given list of  strings and (if the option is true) all
 * of its data items as well. */
static void freeStringList(GList **stringList, const gboolean freeDataItems)
{
  if (stringList && *stringList)
    {
      if (freeDataItems)
	{
	  GList *item = *stringList;
	  for ( ; item; item = item->next)
	    {
	      char *strData = (char*)(item->data);
	      g_free(strData);
	      item->data = NULL;
	    }
	}
      
      g_list_free(*stringList);
      *stringList = NULL;
    }
}


/* Free the memory used by the given sequence group and its members. */
static void destroySequenceGroup(GtkWidget *blxWindow, SequenceGroup **seqGroup)
{
  if (seqGroup && *seqGroup)
    {
      /* Remove it from the list of groups */
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      blxContext->sequenceGroups = g_list_remove(blxContext->sequenceGroups, *seqGroup);
      
      /* If this is pointed to by the match-set pointer, null it */
      if (*seqGroup == blxContext->matchSetGroup)
	{
	  blxContext->matchSetGroup = NULL;
	}
      
      /* Free the memory used by the group name */
      if ((*seqGroup)->groupName)
	{
	  g_free((*seqGroup)->groupName);
	}
      
      /* Free the list of sequences */
      if ((*seqGroup)->seqList)
	{
	  freeStringList(&(*seqGroup)->seqList, (*seqGroup)->ownsSeqNames);
	}
      
      g_free(*seqGroup);
      *seqGroup = NULL;
    }
}


/* Delete a single group */
static void blxWindowDeleteSequenceGroup(GtkWidget *blxWindow, SequenceGroup *group)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext->sequenceGroups)
    {
      destroySequenceGroup(blxWindow, &group);
      blxWindowGroupsChanged(blxWindow);
    }
}


/* Delete all groups */
static void blxWindowDeleteAllSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  GList *groupItem = blxContext->sequenceGroups;
  
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      destroySequenceGroup(blxWindow, &group);
    }
  
  g_list_free(blxContext->sequenceGroups);
  blxContext->sequenceGroups = NULL;
  
  blxWindowGroupsChanged(blxWindow);
}


/* Delete all sequence structs */
static void blxWindowDeleteAllSequenceStructs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  GList *listItem = blxContext->matchSeqs;
  
  /* Delete any data the SequenceStructs own */
  for ( ; listItem; listItem = listItem->next)
    {
      SequenceStruct *seq = (SequenceStruct*)(listItem->data);
      destroySequenceStruct(seq);
    }
  
  /* Free the list itself */
  g_list_free(blxContext->matchSeqs);
  blxContext->matchSeqs = NULL;
}


/* Update function to be called whenever groups have been added or deleted,
 * or sequences have been added to or removed from a group */
static void blxWindowGroupsChanged(GtkWidget *blxWindow)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  callFuncOnAllDetailViewTrees(detailView, deselectAllRows);
  
  /* Re-sort all trees, because grouping affects sort order */
  callFuncOnAllDetailViewTrees(detailView, resortTree);
  
  /* Refilter the trees (because groups affect whether sequences are visible) */
  callFuncOnAllDetailViewTrees(detailView, refilterTree);

  /* Redraw all (because highlighting affects both big picture and detail view) */
  blxWindowRedrawAll(blxWindow);
}


/* Create a new sequence group from the given list of sequence names, with a
 * unique ID and name, and add it to the blxWindow's list of groups. The group 
 * should be destroyed with destroySequenceGroup. If ownSeqNames is true, the group
 * will take ownership of the sequence names and free them when it is destroyed. 
 * Caller can optionally provide the group name; if not provided, a default name
 * will be allocated. */
static SequenceGroup* createSequenceGroup(GtkWidget *blxWindow, GList *seqList, const gboolean ownSeqNames, const char *groupName)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Create the new group */
  SequenceGroup *group = g_malloc(sizeof(SequenceGroup));
  
  group->seqList = seqList;
  group->ownsSeqNames = ownSeqNames;
  group->hidden = FALSE;
  
  /* Find a unique ID */
  GList *lastItem = g_list_last(blxContext->sequenceGroups);
  
  if (lastItem)
    {
      SequenceGroup *lastGroup = (SequenceGroup*)(lastItem->data);
      group->groupId = lastGroup->groupId + 1;
    }
  else
    {
      group->groupId = 1;
    }

  if (groupName)
    {
      group->groupName = g_strdup(groupName);
    }
  else
    {
      /* Create a default name based on the unique ID */
      char formatStr[] = "Group%d";
      const int nameLen = strlen(formatStr) + numDigitsInInt(group->groupId);
      group->groupName = g_malloc(nameLen * sizeof(*group->groupName));
      sprintf(group->groupName, formatStr, group->groupId);
    }
  
  /* Set the order number. For simplicity, set the default order to be the same
   * as the ID number, so groups are sorted in the order they were added */
  group->order = group->groupId;

  /* Set the default highlight colour. */
  group->highlighted = TRUE;
  group->highlightColour = getGdkColor(GDK_RED);

  blxContext->sequenceGroups = g_list_append(blxContext->sequenceGroups, group);

  blxWindowGroupsChanged(blxWindow);
  
  return group;
}


/* This function sets the sequence-group-name text based on the given text entry widget */
static void onGroupNameChanged(GtkWidget *widget, gpointer data)
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
}


/* This function is called when the sequence-group-order text entry widget's
 * value has changed. It sets the new order number in the group. */
static void onGroupOrderChanged(GtkWidget *widget, gpointer data)
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
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      blxWindowGroupsChanged(blxWindow);
    }
}


/* This callback is called when the dialog settings are applied. It sets the hidden
 * status of the passed groupo based on the toggle button's state */
static void onGroupHiddenToggled(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  group->hidden = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

  /* Refilter trees and redraw all immediately show the new status */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));

  blxWindowGroupsChanged(blxWindow);
}


/* This callback is called when the toggle button for a group's "highlighted" flag is toggled.
 * It updates the group's highlighted flag according to the button's new status. */
static void onGroupHighlightedToggled(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  group->highlighted = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  /* Redraw the blixem window to immediately show the new status. The toplevel
   * parent of the button is the dialog, and the blixem window is the transient
   * parent of the dialog. */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  
  blxWindowRedrawAll(blxWindow);
}


/* Called when the user has changed the colour of a group in the 'edit groups' dialog */
static void onGroupColourChanged(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), &group->highlightColour);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), &group->highlightColour, TRUE, TRUE);
  
  /* Redraw everything in the new colours */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(button)));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  blxWindowRedrawAll(blxWindow);
}


/* This function creates a widget that allows the user to edit the group
 * pointed to by the given list item, and adds it to the table container
 * widget at the given row. */
static void createEditGroupWidget(GtkWidget *blxWindow, SequenceGroup *group, GtkTable *table, const int row, const int xpad, const int ypad)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only show the special 'match set' group if it has some sequences */
  if (group != blxContext->matchSetGroup || g_list_length(group->seqList) > 0)
    {
      /* Show the group's name in a text box that the user can edit */
      GtkWidget *nameWidget = gtk_entry_new();
      gtk_entry_set_text(GTK_ENTRY(nameWidget), group->groupName);
      gtk_entry_set_activates_default(GTK_ENTRY(nameWidget), TRUE);
      widgetSetCallbackData(nameWidget, onGroupNameChanged, group);
      
      /* Add a check box for the 'hidden' flag */
      GtkWidget *isHiddenWidget = gtk_check_button_new();
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHiddenWidget), group->hidden);
      widgetSetCallbackData(isHiddenWidget, onGroupHiddenToggled, group);

      /* Add a check box for the 'highlighted' flag */
      GtkWidget *isHighlightedWidget = gtk_check_button_new();
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHighlightedWidget), group->highlighted);
      widgetSetCallbackData(isHighlightedWidget, onGroupHighlightedToggled, group);
      
      /* Show the group's order number in an editable text box */
      GtkWidget *orderWidget = gtk_entry_new();

      char *orderStr = convertIntToString(group->order);
      gtk_entry_set_text(GTK_ENTRY(orderWidget), orderStr);
      g_free(orderStr);
      
      gtk_entry_set_activates_default(GTK_ENTRY(orderWidget), TRUE);
      gtk_widget_set_size_request(orderWidget, 30, -1);
      widgetSetCallbackData(orderWidget, onGroupOrderChanged, group);

      /* Show the group's highlight colour in a button that will also launch a colour-picker */
      GtkWidget *colourButton = gtk_color_button_new_with_color(&group->highlightColour);
      widgetSetCallbackData(colourButton, onGroupColourChanged, group);
      
      /* Create a button that will delete this group */
      GtkWidget *deleteButton = gtk_button_new_from_stock(GTK_STOCK_DELETE);
      g_signal_connect(G_OBJECT(deleteButton), "clicked", G_CALLBACK(onButtonClickedDeleteGroup), group);
      
      /* Put everything in the table */
      gtk_table_attach(table, nameWidget,		1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHiddenWidget,	2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHighlightedWidget,	3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, orderWidget,		4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, colourButton,		5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, deleteButton,		6, 7, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
    }
}


/* Called for each entry in a hash table. Compares the key of the table, which is
 * a sequence name, to the search string given in the user data. If it matches, it
 * appends the sequence name to the result list in the user data. */
static void getSequencesThatMatch(gpointer listDataItem, gpointer data)
{
  SeqSearchData *searchData = (SeqSearchData*)data;
  SequenceStruct *seq = (SequenceStruct*)listDataItem;
  
  /* Cut both down to the short versions of the variant name, in case one is the short
   * version and the other one isn't. */
  if (wildcardSearch(sequenceGetVariantName(seq), getSeqVariantName(searchData->searchStr)))
    {
      searchData->matchList = g_list_prepend(searchData->matchList, seq);
    }
}


/* If the given radio button is enabled, add a group based on the currently-
 * selected sequences. */
static void addGroupFromSelection(GtkWidget *button, gpointer data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {  
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      
      if (g_list_length(blxContext->selectedSeqs) > 0)
	{
	  GList *list = g_list_copy(blxContext->selectedSeqs); /* group takes ownership of this */
	  createSequenceGroup(blxWindow, list, FALSE, NULL);
	}
      else
	{
	  messout("Warning: cannot create group; no sequences are currently selected");
	}
    }
}


/* If the given radio button is enabled, add a group based on the search text
 * in the given text entry. */
static void addGroupFromName(GtkWidget *button, gpointer data)
{
  GList *seqList = findSeqsFromName(button, data);

  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      createSequenceGroup(blxWindow, seqList, FALSE, NULL);
    }
}


/* Utility function to take a list of names, and return a list of those that
 * that are valid sequence names from the given list of SequenceStructs. */
static GList* getSeqStructsFromNameList(GList *nameList, GList *seqList)
{
  SeqSearchData searchData = {NULL, NULL};
  
  /* Loop through all the names in the input list */
  GList *nameItem = nameList;
  for ( ; nameItem; nameItem = nameItem->next)
    {
      /* Compare this name to all names in the sequence list. If it matches,
       * add it to the result list. */
      searchData.searchStr = (const char*)(nameItem->data);
      g_list_foreach(seqList, getSequencesThatMatch, &searchData);
    }
    
  return searchData.matchList;
}


/* Utility to create a GList of SequenceStructs from a textual list of sequences.
 * Returns only valid sequences that blixem knows about. */
static GList* getSeqStructsFromText(GtkWidget *blxWindow, const char *inputText)
{
  GList *nameList = parseMatchList(inputText);
  
  /* Extract the entries from the list that are sequences that blixem knows about */
  GList *matchSeqs = blxWindowGetAllMatchSeqs(blxWindow);
  GList *seqList = getSeqStructsFromNameList(nameList, matchSeqs);
  
  /* Must free the original name list and all its data. */
  freeStringList(&nameList, TRUE);
  
  if (g_list_length(seqList) < 1)
    {
      g_list_free(seqList);
      seqList = NULL;
      messout("No valid sequence names in buffer '%s'\n", inputText);
    }
  
  return seqList;
}


/* Callback function to be used when requesting text from the clipboard to be used
 * to create the 'match set' group from the paste text */
static void createMatchSetFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data)
{
  /* Get the list of sequences to include */
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = getSeqStructsFromText(blxWindow, clipboardText);
  
  /* If a group already exists, replace its list. Otherwise create the group. */
  if (seqList)
    {
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      
      if (!blxContext->matchSetGroup)
	{
	  blxContext->matchSetGroup = createSequenceGroup(blxWindow, seqList, FALSE, "Match Set");
	}
      else
	{
	  if (blxContext->matchSetGroup->seqList)
	    {
	      g_list_free(blxContext->matchSetGroup->seqList);
	    }
	  
	  blxContext->matchSetGroup->seqList = seqList;
	}
      
      /* Reset the highlighted/hidden properties to make sure the group is initially visible */
      blxContext->matchSetGroup->highlighted = TRUE;
      blxContext->matchSetGroup->hidden = FALSE;

      blxWindowGroupsChanged(blxWindow);
    }
}


/* Callback function to be used when requesting text from the clipboard to be used
 * to find and select sequences based on the paste text. Scrolls the first selected
 * match into view. */
void findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data)
{
  /* Get the list of sequences to include */
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = getSeqStructsFromText(blxWindow, clipboardText);
  
  if (seqList)
    {
      blxWindowSetSelectedSeqList(blxWindow, seqList);
      firstMatch(blxWindowGetDetailView(blxWindow), seqList);
    }
}


/* This function toggles the match set.  That is, if the match set (a special 
 * group) exists then it deletes it; if it does not exist, then it creates it
 * from the current clipboard text (which should contain valid sequence name(s)). */
static void toggleMatchSet(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext->matchSetGroup && blxContext->matchSetGroup->seqList)
    {
      /* Clear the list of names only (don't delete the group, because we want to
       * keep any changes the user made (e.g. to the group colour etc.) for next time. */
      g_list_free(blxContext->matchSetGroup->seqList);
      blxContext->matchSetGroup->seqList = NULL;
      blxWindowGroupsChanged(blxWindow);
    }
  else
    {
      requestPrimaryClipboardText(createMatchSetFromClipboard, blxWindow);
    }
}


/* If the given radio button is enabled, add a group based on the list of sequences
 * in the given text entry. */
static void addGroupFromList(GtkWidget *button, gpointer data)
{
  GList *seqList = findSeqsFromList(button, data);
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      createSequenceGroup(blxWindow, seqList, FALSE, NULL);
    }
}


/* Called when the user has clicked the "delete all groups" button on the "group sequences" dialog. */
static void onButtonClickedDeleteAllGroups(GtkWidget *button, gpointer data)
{
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

  /* Ask the user if they're sure */
  gint response = runConfirmationBox(blxWindow, "Delete groups", "This will delete ALL groups. Are you sure?");
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      blxWindowDeleteAllSequenceGroups(blxWindow);
      
      /* Close the dialog, because there are no groups left to display. */
      gtk_widget_destroy(GTK_WIDGET(dialogWindow));
    }
}


/* Called when the user has clicked the "delete group" button on the "group sequences" dialog. */
static void onButtonClickedDeleteGroup(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  /* Ask the user if they're sure */
  char formatStr[] = "Are you sure you wish to delete group '%s'?";
  char messageText[strlen(formatStr) + strlen(group->groupName)];
  sprintf(messageText, formatStr, group->groupName);
  
  gint response = runConfirmationBox(blxWindow, "Delete group", messageText);
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      blxWindowDeleteSequenceGroup(blxWindow, group);
      
      /* Refresh the dialog by closing and re-opening. (Not ideal because it doesn't
       * remember the position). Only re-open if there are still some groups existing. */
      gtk_widget_destroy(GTK_WIDGET(dialogWindow));
      
      if (blxWindowGroupsExist(blxWindow))
	{
	  showGroupsDialog(blxWindow, TRUE);
	}
    }
}


/* Callback for when a radio button with a secondary widget (passed as the user
 * data) is toggled. It enables/disables the other widget according to whether
 * the radio button is active or not. */
static void onRadioButtonToggled(GtkWidget *button, gpointer data)
{
  GtkWidget *otherWidget = GTK_WIDGET(data);
  
  if (otherWidget)
    {
      const gboolean isActive = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
      //      gtk_widget_set_sensitive(otherWidget, isActive); 
      
      if (isActive)
	{
	  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
	  GtkWindow *mainWin = gtk_window_get_transient_for(dialogWindow);

	  gtk_window_set_focus(mainWin, otherWidget);
	}
    }
}


/* Callback when a text entry box in the groups dialog is clicked (for text that
 * is associated with a radio button). Clicking the text box activates its radio button */
static gboolean onRadioButtonTextClicked(GtkWidget *textWidget, GdkEventButton *event, gpointer data)
{
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(data), TRUE);
  
  return FALSE;
}


/* Utility to create a radio button with certain given properties, and to pack it
 * into the given container widget. Returns the radio button (so that further
 * buttons can be created in the same group by passing it as 'existingButton') */
static GtkRadioButton* createRadioButton(GtkBox *box, 
					 GtkRadioButton *existingButton,
					 const char *mnemonic, 
					 const gboolean isActive, 
					 const gboolean createTextEntry,
					 const gboolean multiline,
					 GtkCallback callbackFunc,
					 GtkWidget *blxWindow)
{
  GtkWidget *button = gtk_radio_button_new_with_mnemonic_from_widget(existingButton, mnemonic);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), isActive);
  gtk_box_pack_start(box, button, FALSE, FALSE, 0);
  
  GtkWidget *entry = NULL;
  
  if (createTextEntry && multiline)
    {
      /* Multi-line text buffer */
      GtkTextBuffer *textBuffer = gtk_text_buffer_new(gtk_text_tag_table_new());
      entry = gtk_text_view_new_with_buffer(GTK_TEXT_BUFFER(textBuffer));

      /* Specify a min height */
      const int numLines = 4;
      const int charHeight = detailViewGetCharHeight(blxWindowGetDetailView(blxWindow));
      gtk_widget_set_size_request(entry, -1, charHeight * numLines);

      GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
      gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

      GtkWidget *frame = gtk_frame_new(NULL);
      gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);

      gtk_container_add(GTK_CONTAINER(scrollWin), entry);
      gtk_container_add(GTK_CONTAINER(frame), scrollWin);
      gtk_box_pack_start(box, frame, TRUE, TRUE, 0);
    }
  else if (createTextEntry)
    {
      /* Single line text buffer */
      entry = gtk_entry_new();
      gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
      gtk_box_pack_start(box, entry, TRUE, TRUE, 0);
    }

  if (entry)
    {
      /* to do: don't want to set insensitive because want to receive clicks on text
       * box to activate it; however, it would be good to grey out the background */
//      gtk_widget_set_sensitive(entry, isActive);

      if (isActive)
	{
	  gtk_window_set_focus(GTK_WINDOW(blxWindow), entry);
	}

      g_signal_connect(G_OBJECT(button), "toggled", G_CALLBACK(onRadioButtonToggled), entry);
      g_signal_connect(G_OBJECT(entry), "button-press-event", G_CALLBACK(onRadioButtonTextClicked), button);
    }

  /* Add the callback data. This specifies what callback to use when the dialog is ok'd. */
  widgetSetCallbackData(button, callbackFunc, entry);

  return GTK_RADIO_BUTTON(button);
}


/* Called when the user chooses a different tabe on the groups dialog */
static gboolean onSwitchPageGroupsDialog(GtkNotebook *notebook, GtkNotebookPage *page, guint pageNum, gpointer data)
{
  GtkDialog *dialog = GTK_DIALOG(data);

  if (pageNum == 0)
    {
      /* For the create-groups page, set the default response to be 'accept' */
      gtk_dialog_set_default_response(dialog, GTK_RESPONSE_ACCEPT);
    }
  else
    {
      /* For other pages, set default response to be 'apply' */
      gtk_dialog_set_default_response(dialog, GTK_RESPONSE_APPLY);
    }
    
  return FALSE;
}


/* Callback called when user responds to groups dialog */
void onResponseGroupsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = FALSE;
  gboolean refresh = FALSE;
  
  /* If a notebook was passed, only call callbacks for widgets in the active tab */
  GtkNotebook *notebook = GTK_NOTEBOOK(data);
  guint pageNo = gtk_notebook_get_current_page(notebook);
  GtkWidget *page = gtk_notebook_get_nth_page(notebook, pageNo);
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      widgetCallAllCallbacks(page, NULL);
      destroy = TRUE;
      refresh = FALSE;
      break;

    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(page, NULL);
      destroy = FALSE;
      refresh = (pageNo == 0); /* if created a new group, Edit Groups section must be refreshed */
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      refresh = FALSE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
  else if (refresh)
    {
      /* Re-create the dialog, opening it on the Edit Groups page */
      GtkWindow *dialogWindow = GTK_WINDOW(dialog);
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
      gtk_widget_destroy(GTK_WIDGET(dialog));
      showGroupsDialog(blxWindow, TRUE);
    }
}


/* Shows the "Group sequences" dialog. This dialog allows the user to group sequences together.
 * This tabbed dialog shows both the 'create group' and 'edit groups' dialogs in one. If the
 * 'editGroups' argument is true and groups exist, the 'Edit Groups' tab is displayed by default;
 * otherwise the 'Create Groups' tab is shown. */
void showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Groups", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  GTK_STOCK_APPLY,
						  GTK_RESPONSE_APPLY,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  const gboolean seqsSelected = g_list_length(blxContext->selectedSeqs) > 0;

  const int numRows = g_list_length(blxContext->sequenceGroups) + 2; /* +2 for header row and delete-all button */
  const gboolean groupsExist = blxWindowGroupsExist(blxWindow);
  

  /* Create tabbed pages */
  GtkWidget *notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, TRUE, TRUE, 0);


  /* "CREATE GROUP" SECTION. */
  GtkBox *section1 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(section1), gtk_label_new("Create group"));

  GtkRadioButton *button1 = createRadioButton(section1, NULL, "From _selection", seqsSelected, FALSE, FALSE, addGroupFromSelection, blxWindow);
  createRadioButton(section1, button1, "From _name (wildcards * and ?)", !seqsSelected, TRUE, FALSE, addGroupFromName, blxWindow);
  createRadioButton(section1, button1, "From _list", FALSE, TRUE, TRUE, addGroupFromList, blxWindow);

  /* "EDIT GROUP" SECTION. */
  GtkWidget *section2 = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), section2, gtk_label_new("Edit groups"));

  const int numCols = 6;
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_box_pack_start(GTK_BOX(section2), GTK_WIDGET(table), FALSE, FALSE, 0);

  const int xpad = 2;
  const int ypad = 2;

  /* Add labels */
  int row = 1;
  gtk_table_attach(table, gtk_label_new("Group name"),	  1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Hide"),	  2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Highlight"), 3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Order"),	  4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  ++row;
  
  /* Add a set of widgets for each group */
  GList *groupItem = blxWindowGetSequenceGroups(blxWindow);
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      createEditGroupWidget(blxWindow, group, table, row, xpad, ypad);
      ++row;
    }
  
  /* Add a button to delete all groups */
  GtkWidget *deleteGroupsButton = gtk_button_new_with_label("Delete all groups");
  gtk_widget_set_size_request(deleteGroupsButton, -1, 30);
  g_signal_connect(G_OBJECT(deleteGroupsButton), "clicked", G_CALLBACK(onButtonClickedDeleteAllGroups), NULL);
//  gtk_box_pack_end(GTK_BOX(section2), deleteGroupsButton, FALSE, FALSE, 0);
  gtk_table_attach(table, deleteGroupsButton, numCols - 1, numCols + 1, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  /* Connect signals and show */
  gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(blxWindow));
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseGroupsDialog), notebook);
  g_signal_connect(notebook, "switch-page", G_CALLBACK(onSwitchPageGroupsDialog), dialog);
  
  gtk_widget_show_all(dialog);
  
  /* If user has asked to edit groups (and some groups exist), make the second tab
   * the default and the 'close' button the default action. (Must do this after showing
   * the child widgets due to a GTK legacy whereby the notebook won't change tabs otherwise.) */
  if (editGroups && groupsExist)
    {
      gtk_notebook_next_page(GTK_NOTEBOOK(notebook));
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


/* Callback function called when the 'Show SNP track' button is toggled */
static void onShowSnpTrackToggled(GtkWidget *button, gpointer data)
{
  const gboolean showSnpTrack = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  GtkWidget *detailView = GTK_WIDGET(data);
  
  detailViewSetShowSnpTrack(detailView, showSnpTrack);
}


/* Utility to create a check button with certain given properties, and to pack it into the parent */
static void createCheckButton(GtkBox *box, 
			      const char *mnemonic, 
			      const gboolean isActive, 
			      GCallback callback, 
			      GtkWidget *blxWindow)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_box_pack_start(box, button, FALSE, FALSE, 0);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), isActive);
  
  g_signal_connect(G_OBJECT(button), "toggled", callback, blxWindow);
}


/* Callback to be called when the user has entered a new column size */
static void onColumnSizeChanged(GtkWidget *widget, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)data;
  
  const gchar *newSizeText = gtk_entry_get_text(entry);
  const int newWidth = convertStringToInt(newSizeText);
  
  if (newWidth != columnInfo->width)
    {
      columnInfo->width = newWidth;

      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(widget));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
      callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns);
      resizeDetailViewHeaders(detailView);
    }
}


/* Create a set of widgets that allow columns to be resized */
static void createColumnSizeButtons(GtkWidget *parent, GtkWidget *detailView)
{
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Column sizes");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);

  /* Arrange the widgets horizontally */
  GtkWidget *hbox = createHBoxWithBorder(frame, 12);
  
  /* Loop through each column in the detail view and create a text box showing the
   * current width. */
  GList *listItem = detailViewGetColumnList(detailView);

  for ( ; listItem; listItem = listItem->next)
    {
      DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(listItem->data);
      
      /* Create a label and a text entry box, arranged vertically in a vbox */
      GtkWidget *vbox = createVBoxWithBorder(hbox, 4, FALSE, NULL);
      
      GtkWidget *label = gtk_label_new(columnInfo->title);
      gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
      
      GtkWidget *entry = gtk_entry_new();
      gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
      
      if (columnInfo->columnId == SEQUENCE_COL)
	{
	  /* The sequence column updates dynamically, so don't allow the user to edit it */
	  char displayText[] = "<dynamic>";
	  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
	  gtk_widget_set_sensitive(entry, FALSE);
	  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2); /* fudge up width a bit in case user enters longer text */
	}
      else
	{
	  char *displayText = convertIntToString(columnInfo->width);
	  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
	  g_free(displayText);

	  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
	  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

	  widgetSetCallbackData(entry, onColumnSizeChanged, (gpointer)columnInfo);
	}
    }
}


/* Callback to be called when the user has entered a new percent-ID per cell */
static void onIdPerCellChanged(GtkWidget *widget, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const int newValue = convertStringToInt(text);
  bigPictureSetIdPerCell(bigPicture, newValue);
}

 
/* Callback to be called when the user has entered a new maximum percent-ID to display */
static void onMaxPercentIdChanged(GtkWidget *widget, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const int newValue = convertStringToInt(text);
  bigPictureSetMaxPercentId(bigPicture, newValue);
}

/* Callback to be called when the user has entered a new minimum percent-ID to display */
static void onMinPercentIdChanged(GtkWidget *widget, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const int newValue = convertStringToInt(text);
  bigPictureSetMinPercentId(bigPicture, newValue);
}


/* Utility to create a text entry widget displaying the given integer value. The
 * given callback will be called when the user OK's the dialog that this widget 
 * is a child of. */
static void createTextEntryFromInt(GtkWidget *parent, 
				   const char *title, 
				   const int value, 
				   GtkCallback callbackFunc, 
				   gpointer callbackData)
{
  /* Pack label and text entry into a vbox */
  GtkWidget *vbox = createVBoxWithBorder(parent, 4, FALSE, NULL);

  GtkWidget *label = gtk_label_new(title);
  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

  GtkWidget *entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

  char *displayText = convertIntToString(value);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  g_free(displayText);

  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

  widgetSetCallbackData(entry, callbackFunc, callbackData);
}


/* Create a set of widgets to allow the user to edit grid properties */
static void createGridSettingsButtons(GtkWidget *parent, GtkWidget *bigPicture)
{
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Grid properties");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);

  /* Arrange the widgets horizontally */
  GtkWidget *hbox = createHBoxWithBorder(frame, 12);
  const IntRange const *percentIdRange = bigPictureGetPercentIdRange(bigPicture);
  
  createTextEntryFromInt(hbox, "%ID per cell", bigPictureGetIdPerCell(bigPicture), onIdPerCellChanged, bigPicture);
  createTextEntryFromInt(hbox, "Max %ID", percentIdRange->max, onMaxPercentIdChanged, bigPicture);
  createTextEntryFromInt(hbox, "Min %ID", percentIdRange->min, onMinPercentIdChanged, bigPicture);
}


/* Shows the "Settings" dialog. */
void showSettingsDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Blixem Settings", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  GTK_STOCK_APPLY,
						  GTK_RESPONSE_APPLY,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  GtkWidget *mainVBox = createVBoxWithBorder(contentArea, borderWidth, FALSE, NULL);

  /* Display options */
  GtkWidget *frame1 = gtk_frame_new("Display options");
  gtk_box_pack_start(GTK_BOX(mainVBox), frame1, FALSE, FALSE, 0);

  GtkWidget *vbox1 = createVBoxWithBorder(frame1, borderWidth, FALSE, NULL);
  
  createCheckButton(GTK_BOX(vbox1), "_Squash matches", detailViewGetMatchesSquashed(detailView), G_CALLBACK(onSquashMatches), detailView);
  createCheckButton(GTK_BOX(vbox1), "_Invert sort order", detailViewGetSortInverted(detailView), G_CALLBACK(onSortOrderToggled), detailView);
  createCheckButton(GTK_BOX(vbox1), "_Highlight differences", detailViewGetHighlightDiffs(detailView), G_CALLBACK(onHighlightDiffsToggled), detailView);
  createCheckButton(GTK_BOX(vbox1), "_Show SN_P track", detailViewGetShowSnpTrack(detailView), G_CALLBACK(onShowSnpTrackToggled), detailView);
  
  createColumnSizeButtons(mainVBox, detailView);
  createGridSettingsButtons(mainVBox, bigPicture);
  
  /* Connect signals and show */
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  gtk_widget_show_all(dialog);
}


/***********************************************************
 *		    Statistics menu			   *
 ***********************************************************/

static void getStats(GtkWidget *blxWindow, GString *result, MSP *MSPlist)
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
  
  int refSeqLen = strlen(blxWindowGetRefSeq(blxWindow));
  
  /* Create the text based on the results */
  g_string_printf(result, "%s%d%s%s%d%s%s%d%s%s%d%s%s%d%s%s%d%s",
                  "Length of reference sequence\t\t= ", refSeqLen, " characters\n\n",
                  "Number of match sequences\t\t= ", numSeqs, "\n",
                  "Total memory used by sequences\t= ", totalSeqSize, " bytes\n\n",
		  "Number of MSPs\t\t\t\t\t= ", numMSPs, "\n",
                  "Size of each MSP\t\t\t\t\t= ", (int)sizeof(MSP), " bytes\n",
		  "Total memory used by MSPs\t\t= ", (int)sizeof(MSP) * numMSPs, " bytes");
}


static void showStatsDialog(GtkWidget *blxWindow, MSP *MSPlist)
{
  /* Create a dialog widget with an OK button */
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Statistics", 
                                                  GTK_WINDOW(blxWindow),
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Ensure that the dialog box (along with any children) is destroyed when the user responds. */
  g_signal_connect (dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  /* Create a text buffer containing the required text*/
  GString *displayText = g_string_sized_new(200); /* will be extended if we need more space */
  getStats(blxWindow, displayText, MSPlist);
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
 *			About dialog			   *
 ***********************************************************/

/* Returns a string which is the name of the Blixem application. */
static char *blxGetAppName(void)
{
  return BLIXEM_TITLE ;
}

/* Returns a copyright string for the Blixem application. */
static char *blxGetCopyrightString(void)
{
  return BLIXEM_COPYRIGHT_STRING ;
}

/* Returns the Blixem website URL. */
static char *blxGetWebSiteString(void)
{
  return BLIXEM_WEBSITE_STRING ;
}

/* Returns a comments string for the Blixem application. */
static char *blxGetCommentsString(void)
{
  return BLIXEM_COMMENTS_STRING(BLIXEM_TITLE, BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE) ;
}

/* Returns a license string for the blx application. */
static char *blxGetLicenseString(void)
{
  return BLIXEM_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Blixem code. */
static char *blxGetVersionString(void)
{
  return BLIXEM_VERSION_STRING ;
}


/* Shows the 'About' dialog */
void showAboutDialog(GtkWidget *parent)
{
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (6)
  const gchar *authors[] = {BLIXEM_AUTHOR_LIST, NULL} ;

  gtk_show_about_dialog(GTK_WINDOW(parent),
			"authors", authors,
			"comments", blxGetCommentsString(), 
			"copyright", blxGetCopyrightString(),
			"license", blxGetLicenseString(),
			"name", blxGetAppName(),
			"version", blxGetVersionString(),
			"website", blxGetWebSiteString(),
			NULL) ;
#endif

  return ;
}


/***********************************************************
 *			Help menu			   *
 ***********************************************************/

void onResponseHelpDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      destroy = TRUE;
      break;
      
    case GTK_RESPONSE_HELP:
      showAboutDialog(NULL);
      destroy = FALSE;
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


void showHelpDialog(GtkWidget *blxWindow)
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
	•	Use the mouse-wheel to scroll up/down in the list of alignments.\n\
	•	Use the mouse-wheel to scroll the alignments left/right (if your mouse has a horizontal scroll-wheel).\n\
	•	Ctrl-left and Ctrl-right arrow keys scroll to the start/end of the previous/next match (limited to currently-selected sequences, if any are selected).\n\
	•	The Home/End keys scroll to the start/end of the display.\n\
	•	Ctrl-Home and Ctrl-End scroll to the start/end of the currently-selected alignments (or to the first/last alignment if none are selected).\n\
	•	The ',' (comma) and '.' (full-stop) keys scroll the display one nucleotide to the left/right.\n\
	•	You can scroll to a specific position using the Go-to button on the toolbar, or by pressing the 'p' shortcut key.\n\
\n\
\n\
SELECTIONS\n\
	•	You can select a sequence by clicking on its row in the alignment list.  Selected sequences are highlighted in cyan in the big picture.\n\
	•	You can select a sequence by clicking on it in the big picture.\n\
	•	The name of the sequence you selected is displayed in the feedback box on the toolbar.  If there are multiple alignments for the same sequence, all of them will be selected.\n\
	•	You can select multiple sequences by holding down the Ctrl or Shift keys while selecting rows.\n\
	•	You can deselect a single sequence by Ctrl-clicking on its row.\n\
	•	You can deselect all sequences by right-clicking and selecting 'Deselect all', or with the Shift-Ctrl-A keyboard shortcut.\n\
	•	You can move the selection up/down a row using the up/down arrow keys.\n\
\n\
	•	You can select a nucleotide/peptide by middle-clicking on it in the detail view.  This selects the entire column at that index, and the coordinate number on the query sequence is shown in the feedback box.  (The coordinate on the subject sequence is also shown if a subject sequence is selected.)\n\
	•	For protein matches, when a peptide is selected, the three nucleotides for that peptide (for the current reading frame) are highlighted in the header in green.  The current reading frame is whichever alignment list currently has the focus - click in a different list to change the reading frame.  The darker highlighting indicates the specific nucleotide that is currently selected (i.e. whose coordinate is displayed in the feedback box).\n\
	•	You can move the selection to the previous/next nucleotide using the left and right arrow keys.\n\
	•	You can move the selection to the previous/next peptide by holding Shift while using the left and right arrow keys.\n\
	•	You can move the selection to the start/end of the previous/next matchb by holding Ctrl while using the left and right arrow keys (limited to just the selected sequences if any are selected).\n\
	•	To select a coordinate without the display re-centering on it, hold down Ctrl as you middle-click.\n\
\n\
\n\
ZOOMING\n\
	•	Zoom in to the currently-selected region in the big picture using the +/- buttons in the top-left corner of the window, or using the Ctrl '=' and Ctrl '-' shortcut keys.  The 'Whole' button zooms out to show the full length of the query sequence.\n\
	•	Zoom in/out of the detail view using the +/- buttons on the toolbar, or using the '=' and '-' shortcut keys.\n\
\n\
\n\
COPY AND PASTE\n\
	•       When sequence(s) are selected, their names are copied to the selection buffer and can be pasted by middle-clicking.\n\
	•       To paste sequence names from the selection buffer, hit the 'f' keyboard shortcut. Blixem will select the sequences and jump to the start of the selection.\n\
        •	To copy sequence name(s) to the default clipboard, select the sequence(s) and hit Ctrl-C. Sequence names can then be pasted into other applications using Ctrl-V.\n\
	•	The clipboard text can also be pasted into Blixem using Ctrl-V. If the clipboard contains valid sequence names, those sequences will be selected and the display will jump to the start of the selection.\n\
	•	Note that text from the feedback box and some text labels (e.g. the reference sequence start/end coords) can be copied by selecting it with the mouse and then hitting Ctrl-C.\n\
	•	Text can be pasted into dialog box text entry boxes using Ctrl-V (or middle-clicking to paste from the selection buffer).\n\
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
	•	Select the sequences you wish to include in the group by left-clicking their rows in the detail view.  Multiple rows can be selected by holding the Ctrl or Shift keys while clicking.\n\
	•	Right-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Note that Ctrl-G will also shortcut to here if no groups currently exist.)\n\
	•	Ensure that the 'From selection' radio button is selected, and click 'OK'.\n\
\n\
Creating a group from a sequence name:\n\
	•	Right-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Or Ctrl-G if no groups currently exist.)\n\
	•	Select the 'From name' radio button and enter the name of the sequence in the box below.  You may use the following wildcards to search for sequences: '*' for any number of characters; '?' for a single character.  For example, searching for '*X' will find all sequences whose name ends in 'X' (i.e. all exons).\n\
	•	Click 'OK'.\n\
\n\
Creating a group from sequence name(s):\n\
	•	Right-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Or Ctrl-G if no groups currently exist.)\n\
	•	Select the 'From name(s)' radio button.\n\
        •       Enter the sequence name(s) in the text box.\n\
        •       You may use the following wildcards in a sequence name: '*' for any number of characters; '?' for a single character.  (For example, searching for '*X' will find all sequences whose name ends in 'X', i.e. all exons).\n\
        •       You may search for multiple sequence names by separating them with the following delimiters: newline, comma or semi-colon.\n\
        •       You may paste sequence names directly from ZMap: click on the feature in ZMap and then middle-click in the text box on the Groups dialog.  Grouping in Blixem works on the sequence name alone, so the feature coords will be ignored.\n\
	•	Click 'OK'.\n\
\n\
Creating a temporary 'match-set' group from the current selection:\n\
        •       You can quickly create a group from a current selection (e.g. selected features in ZMap) using the 'Toggle match set' option.\n\
        •       To create a match-set group, select the required items (e.g. in ZMap) and then select 'Toggle match set' from the right-click menu in Blixem, or hit the 'g' shortcut key.\n\
        •       To clear the match-set group, choose the 'Toggle match set' option again, or hit the 'g' shortcut key again.\n\
        •       While it exists, the match-set group can be edited like any other group, via the 'Edit Groups' dialog.\n\
        •       If you delete the match-set group from the 'Edit Groups' dialog, all settings (e.g. highlight colour) will be lost. To maintain these settings, clear the group using the 'Toggle match set' menu option (or 'g' shortcut key) instead.\n\
\n\
Editing groups:\n\
To edit a group, right-click and select 'Edit Groups', or use the Ctrl-G shortcut key. You can change the following properties for a group:\n\
	•	Name: you can specify a more meaningful name to help identify the group.\n\
	•	Hide: tick this box to hide the alignments in the alignment lists.\n\
	•	Highlight: tick this box to highlight the alignments.\n\
	•	Colour: the colour the group will be highlighted in, if 'Highlight' is enabled.  The default colour for all groups is red, so you may wish to change this if you want different groups to be highlighted in different colours.\n\
	•	Order: when sorting by Group, alignments in a group with a lower order number will appear before those with a higher order number (or vice versa if sort order is inverted). Alignments in a group will appear before alignments that are not in a group.\n\
	•	To delete a single group, click on the 'Delete' button next to the group you wish to delete.\n\
	•	To delete all groups, click on the 'Delete all groups' button.\n\
\n\
\n\
DOTTER\n\
	•	To start Dotter, or to edit the Dotter parameters, right-click and select 'Dotter' or use the Ctrl-D keyboard shortcut.	The Dotter settings dialog will pop up.\n\
	•	To run Dotter with the default (automatic) parameters, just hit RETURN, or click the 'Execute' button.\n\
	•	To enter manual parameters, click the 'Manual' radio button and enter the values in the 'Start' and 'End' boxes.\n\
	•	To revert to the last-saved manual parameters, click the 'Last saved' button.\n\
	•	To revert back to automatic parameters, click the 'Auto' radio button.\n\
	•	To save the parameters without running Dotter, click Save and then Cancel'.\n\
	•	To save the parameters and run Dotter, click 'Execute'.\n\
\n\
\n\
KEYBOARD SHORTCUTS\n\
	•	Ctrl-Q: Quit\n\
	•	Ctrl-H: Help\n\
	•	Ctrl-P: Print\n\
	•	Ctrl-S: 'Settings' menu\n\
	•	V: 'View' menu (for hiding sections of the display)\n\
	•	Shift-Ctrl-G: Create group\n\
	•	Ctrl-G: Edit groups (or create a group if none currently exist)\n\
	•	Ctrl-A: Select all sequences in the current list\n\
	•	Shift-Ctrl-A: Deselect all sequences\n\
	•	Ctrl-D: Dotter\n\
	•	Left/right arrow keys: Move the currently-selected coordinate one place to the left/right\n\
	•	Shift-Left/right: Same as left/right arrow keys, but for proteins it scrolls by a single nucleotide, rather than an entire codon.\n\
	•	Ctrl-Left/right: Scroll to the start/end of the previous/next alignment (limited to just the selected sequences, if any are selected).\n\
	•	Home/End: Scroll to the start/end of the display.\n\
	•	Ctrl-Home/End: Scroll to the first/last alignment (limited to just the selected sequences, if any are selected).\n\
	•	=: zoom in detail view\n\
	•	-: zoom out detail view\n\
	•	Ctrl-= : zoom in big picture\n\
	•	Ctrl-- :  zoom out big picture\n\
	•	Shift-Ctrl-- :  zoom out big picture to whole width\n\
	•	,: scroll left one coordinate\n\
	•	.: scroll right one coordinate\n\
	•	p: Go to position\n\
	•	t: Toggle the active strand\n\
	•	g: Toggle the 'match set' Group\n\
	•       1, 2, 3: These number keys toggle visibility of the 1st, 2nd (and 3rd, for protein matches) alignment list.\n\
	•	Ctrl-1, Ctrl-2: This toggles visibility of the 1st and 2nd big picture grid.\n\
	•       Shift-Ctrl-1, Shift-Ctrl-2: This toggles visibility of the 1st and 2nd exon views.\n\
\n\
\n\
SETTINGS\n\
	•	The settings menu can be accessed by right-clicking and selecting Settings, or by the shortcut Ctrl-S.\n\
	•	Squash Matches: this groups multiple alignments from the same sequence together into the same row in the detail view, rather than showing them on separate rows.\n\
	•	Invert Sort Order: reverse the default sort order. (Note that some columns sort ascending by default (e.g. name, start, end) and some sort descending (score and ID). This option reverses that sort order.)\n\
	•	Highlight Differences: when this option is set, matching bases are blanked out and mismatches are highlighted, making it easier to see where alignments differ from the reference sequence.\n\
	•       Column sizes: use this to change the width of the columns.\n\
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
	•	Thin blue vertical line: start boundary of an exon\n\
	•	Thin dark-blue vertical line: end boundary of an exon\n\
	•	Green background (protein matches only): the three nucleotides for the currently-selected codon. Dark green indicates the specific nucleotide that is currently displayed in the feedback box.\n\
", blixemVersion));

  /* Set a pretty big initial size */
  const int width = blxWindow->allocation.width * 0.7;
  const int maxHeight = blxWindow->allocation.height * 0.7;
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Help", 
						  NULL, 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_ABOUT,
						  GTK_RESPONSE_HELP,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  int height = maxHeight;
  GtkWidget *child = createScrollableTextView(messageText, TRUE, blxWindow->style->font_desc, &height);
  
  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseHelpDialog), NULL);
  gtk_widget_show_all(dialog);
}


/***********************************************************
 *			  Menu actions                     *
 ***********************************************************/

/* Called when the user selects the quit menu option, or hits the Quit shortcut key.
 * Pops up a dialog showing user help information */
static void onQuit(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  gtk_widget_destroy(blxWindow);
}


/* Called when the user selects the Help menu option, or hits the Help shortcut key. */
static void onHelpMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showHelpDialog(blxWindow);
}


/* Called when the user selects the View menu option, or hits the Settings shortcut key */
static void onViewMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showViewPanesDialog(blxWindow);
}


/* Called when the user selects the 'Group Sequences' menu option, or hits the relevant shortcut key */
static void onCreateGroupMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, FALSE);
}


/* Called when the user selects the 'Groups' menu option, or hits the relevant shortcut key */
static void onEditGroupsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, TRUE);
}

/* Called when the user selects the 'Toggle match set' option, or hits the relevant shortcut key */
static void onToggleMatchSet(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  toggleMatchSet(blxWindow);
}

/* Called when the user selects the Settings menu option, or hits the Settings shortcut key */
static void onSettingsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showSettingsDialog(blxWindow);
}


/* Called when the user selects the Dotter menu option, or hits the Dotter shortcut key */
static void onDotterMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showDotterDialog(blxWindow);
}


/* Called when the user selects the 'Select Features' menu option, or hits the relevant shortcut key */
static void onSelectFeaturesMenu(GtkAction *action, gpointer data)
{
  selectFeatures();
}


/* Called when the user selects the 'Deselect all' menu option, or hits the relevant shortcut key */
static void onDeselectAllRows(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  blxWindowDeselectAllSeqs(blxWindow, TRUE);
}


/* Called when the user selects the Statistics menu option, or hits the Statistics shortcut key.
 * Pops up a dialog showing memory usage statistics. */
static void onStatisticsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  MSP *mspList = blxWindowGetMspList(blxWindow);
  showStatsDialog(blxWindow, mspList);
}


/* Called when the user selects the Print menu option, or hits the Print shortcut key */
static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  
  /* Create a print operation, using the same settings as the last print, if there was one */
  GtkPrintOperation *print = gtk_print_operation_new();
  
  if (properties->printSettings != NULL)
    {
      gtk_print_operation_set_print_settings(print, properties->printSettings);
    }
  
  g_signal_connect (print, "begin_print", G_CALLBACK (onBeginPrint), blxWindow);
  g_signal_connect(G_OBJECT(print), "draw-page", G_CALLBACK(onDrawPage), blxWindow);
  
  /* Pop up the print dialog */
  GtkPrintOperationResult printResult = gtk_print_operation_run (print, 
								 GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
								 GTK_WINDOW(blxWindow),
								 NULL);
  
  /* If the user hit ok, remember the print settings for next time */
  if (printResult == GTK_PRINT_OPERATION_RESULT_APPLY)
    {
      if (properties->printSettings != NULL)
	{
	  g_object_unref(properties->printSettings);
	  properties->printSettings = NULL;
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
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);

  /* Set the orientation based on the stored print settings (not sure why this doesn't
   * already get set when the settings are set in the print operation...). */
  GtkPrintSettings *printSettings = gtk_print_operation_get_print_settings(print);
  gtk_print_settings_set_orientation(printSettings, gtk_print_settings_get_orientation(properties->printSettings));

  gtk_print_operation_set_n_pages(print, 1);
}


static void collatePixmaps(GtkWidget *widget, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  GdkDrawable *drawable = widgetGetDrawable(widget);

  /* If this widget is visible and has a drawable set, draw it onto the window's drawable */
  if (drawable && GTK_WIDGET_VISIBLE(widget))
    {
      /* Get the left edge of the widget in terms of the blixem window's coordinates */
      int xSrc = 0, ySrc = 0;
      int xDest, yDest;
      gtk_widget_translate_coordinates(widget, blxWindow, xSrc, ySrc, &xDest, &yDest);

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
      gdk_draw_drawable(widgetGetDrawable(blxWindow), gc, drawable, xSrc, ySrc, x, y, -1, -1); /* -1 means full width/height */
    }
  
  /* If this widget is a container, recurse over its children */
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), collatePixmaps, blxWindow);
    }
}


/* Print handler - renders a specific page */
static void onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  
  /* Create a blank white pixmap to draw on to */
  GdkDrawable *drawable = gdk_pixmap_new(blxWindow->window, blxWindow->allocation.width, blxWindow->allocation.height, -1);
  widgetSetDrawable(blxWindow, drawable);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor fgColour = getGdkColor(GDK_WHITE);
  gdk_gc_set_foreground(gc, &fgColour);
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, blxWindow->allocation.width, blxWindow->allocation.height);

  /* For each child widget that has a drawable set, draw this onto the main pixmap */
  properties->lastYStart = 0;
  properties->lastYEnd = 0;
  properties->lastYCoord = -1;
  gtk_container_foreach(GTK_CONTAINER(blxWindow), collatePixmaps, blxWindow);
  
  cairo_t *cr = gtk_print_context_get_cairo_context (context);
  gdk_cairo_set_source_pixmap(cr, drawable, 0, 0);
  cairo_paint(cr);
}


/* Mouse button handler */
static gboolean onButtonPressBlxWindow(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3)
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
  }
  
  return TRUE;
}


/* Key press handler */
static gboolean onKeyPressBlxWindow(GtkWidget *window, GdkEventKey *event, gpointer data)
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
	  goToMatch(window, event->keyval == GDK_Left);
	else if (shiftModifier)
	  moveSelectedBaseIdxBy1(window, event->keyval == GDK_Left);
	else
	  moveSelectedDisplayIdxBy1(window, event->keyval == GDK_Left);
	
	result = TRUE;
	break;
      }
	
      case GDK_Home:
      case GDK_End:
	scrollToExtremity(window, event->keyval == GDK_Home, ctrlModifier);
	result = TRUE;
	break;

      case GDK_comma:  /* fall through */
      case GDK_period:
	scrollDetailView(window, event->keyval == GDK_comma, ctrlModifier);
	result = TRUE;
	break;

      case GDK_underscore: /* fall through */
      case GDK_equal:	   /* fall through */
      case GDK_minus:
      {
	const gboolean zoomIn = (event->keyval == GDK_plus || event->keyval == GDK_equal);
	zoomBlxWindow(window, zoomIn, ctrlModifier, shiftModifier);
	result = TRUE;
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

      case GDK_v: /* fall through */
      case GDK_V:
	if (ctrlModifier)
	  {
	    /* Paste from the default clipboard */
	    requestDefaultClipboardText(findSeqsFromClipboard, window);
	    result = TRUE;
	  }
	break;
	
      case GDK_f: /* fall through */
      case GDK_F:
      {
	if (ctrlModifier)
	  showFindDialog(window);
	else
	  requestPrimaryClipboardText(findSeqsFromClipboard, window);
	
	result = TRUE;
	break;
      }
	
      case GDK_p: /* fall through */
      case GDK_P:
	if (!ctrlModifier)
	  {
	    goToDetailViewCoord(blxWindowGetDetailView(window), BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
	    result = TRUE;
	  }
	break;
	
      case GDK_g: /* fall through */
      case GDK_G:
	if (!ctrlModifier)
	  {	    
	    toggleMatchSet(window);
	    result = TRUE;
	  }
	break;
	
      case GDK_b: /* fall through */
      case GDK_B:
	toggleBumpState(window);
	result = TRUE;
	break;
	
      case GDK_t: /* fall through */
      case GDK_T:
	toggleStrand(blxWindowGetDetailView(window));
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

static BlxWindowProperties* blxWindowGetProperties(GtkWidget *widget)
{
  return widget ? (BlxWindowProperties*)(g_object_get_data(G_OBJECT(widget), "BlxWindowProperties")) : NULL;
}

BlxViewContext* blxWindowGetContext(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->blxContext : NULL;
}

static void onDestroyBlxWindow(GtkWidget *widget)
{
  BlxWindowProperties *properties = blxWindowGetProperties(widget);
  
  if (properties && properties->blxContext)
    {
      /* Free the list of selected sequence names (not the names themselves
       * because we don't own them). */
      if (properties->blxContext->selectedSeqs)
	{
	  g_list_free(properties->blxContext->selectedSeqs);
	  properties->blxContext->selectedSeqs = NULL;
	}
      
      blxWindowDeleteAllSequenceGroups(widget);
      blxWindowDeleteAllSequenceStructs(widget);

      if (properties->blxContext->fetchMode)
	{
	  g_free(properties->blxContext->fetchMode);
	  properties->blxContext->fetchMode = NULL;
	}
     
      /* Free the context struct itself */
      g_free(properties->blxContext);
      properties->blxContext = NULL;
    }
    
  if (properties)
    {
      if (properties->mainmenu)
	{
	  gtk_widget_destroy(properties->mainmenu);
	  properties->mainmenu = NULL;
	}
      
      /* Destroy the print settings */
      if (properties->printSettings)
	{
	  g_object_unref(properties->printSettings);
	  properties->printSettings = NULL;
	}
      
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "BlxWindowProperties", NULL);
    }
  
  /* Destroy variables created by the main program */
  blviewDestroy(widget);
  
  gtk_main_quit();
}


static BlxViewContext* blxWindowCreateContext(CommandLineOptions *options,
					      const IntRange const *refSeqRange,
					      const IntRange const *fullDisplayRange,
					      const char *paddingSeq)
{
  BlxViewContext *blxContext = g_malloc(sizeof *blxContext);
  
  blxContext->refSeq = options->refSeq;
  blxContext->refSeqName = options->refSeqName ? g_strdup(options->refSeqName) : g_strdup("Blixem-seq");
  blxContext->refSeqRange.min = refSeqRange->min;
  blxContext->refSeqRange.max = refSeqRange->max;
  blxContext->fullDisplayRange.min = fullDisplayRange->min;
  blxContext->fullDisplayRange.max = fullDisplayRange->max;
  
  blxContext->mspList = options->mspList;
  blxContext->geneticCode = options->geneticCode;
  blxContext->blastMode = options->blastMode;
  blxContext->seqType = options->seqType;
  blxContext->numFrames = options->numFrames;
  blxContext->gappedHsp = options->gappedHsp;
  blxContext->paddingSeq = paddingSeq;
  blxContext->fetchMode = g_strdup(options->fetchMode);
  blxContext->matchSeqs = NULL; /* will be initialised from MSP data */
  
  blxContext->displayRev = FALSE;
  blxContext->selectedSeqs = NULL;
  blxContext->sequenceGroups = NULL;
  blxContext->matchSetGroup = NULL;
  
  blxContext->autoDotter = TRUE;
  blxContext->dotterStart = UNSET_INT;
  blxContext->dotterEnd = UNSET_INT;
  blxContext->dotterZoom = 0;
  
  return blxContext;
}


/* Create the properties struct and initialise all values. */
static void blxWindowCreateProperties(CommandLineOptions *options,
				      BlxViewContext *blxContext,
				      GtkWidget *widget, 
				      GtkWidget *bigPicture, 
				      GtkWidget *detailView,
				      GtkWidget *mainmenu,
				      const IntRange const *refSeqRange,
				      const IntRange const *fullDisplayRange,
				      const char *paddingSeq)
{
  if (widget)
    {
      BlxWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->blxContext = blxContext;

      properties->bigPicture = bigPicture;
      properties->detailView = detailView;
      properties->mainmenu = mainmenu;

      properties->printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(properties->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
      gtk_print_settings_set_quality(properties->printSettings, GTK_PRINT_QUALITY_HIGH);
      properties->lastYEnd = UNSET_INT;
      properties->lastYStart = UNSET_INT;
      properties->lastYCoord = UNSET_INT;

      g_object_set_data(G_OBJECT(widget), "BlxWindowProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK (onDestroyBlxWindow), NULL);
    }
}

gboolean blxWindowGetDisplayRev(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->displayRev : FALSE;
}

GtkWidget* blxWindowGetBigPicture(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->bigPicture : FALSE;
}

GtkWidget* blxWindowGetDetailView(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->detailView : FALSE;
}

GtkWidget* blxWindowGetMainMenu(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->mainmenu : FALSE;
}

BlxBlastMode blxWindowGetBlastMode(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->blastMode : FALSE;
}

const char* blxWindowGetFetchMode(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->fetchMode : FALSE;
}

char * blxWindowGetRefSeq(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->refSeq : NULL;
}

const char * blxWindowGetRefSeqName(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->refSeqName : NULL;
}

char** blxWindowGetGeneticCode(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->geneticCode : NULL;
}

MSP* blxWindowGetMspList(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->mspList : NULL;
}

GList* blxWindowGetAllMatchSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->matchSeqs : NULL;
}

BlxSeqType blxWindowGetSeqType(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->seqType : FALSE;
}

IntRange* blxWindowGetFullRange(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? &blxContext->fullDisplayRange : NULL;
}

IntRange* blxWindowGetRefSeqRange(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? &blxContext->refSeqRange : NULL;
}

int blxWindowGetNumFrames(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->numFrames : UNSET_INT;
}

int blxWindowGetAutoDotter(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->autoDotter : TRUE;
}

int blxWindowGetDotterStart(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterStart : UNSET_INT;
}

int blxWindowGetDotterEnd(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterEnd : UNSET_INT;
}

int blxWindowGetDotterZoom(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterZoom : UNSET_INT;
}

gboolean blxWindowGetGappedHsp(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->gappedHsp : UNSET_INT;
}

const char* blxWindowGetPaddingSeq(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->paddingSeq : NULL;
}

/* Return the active strand - forward strand by default, reverse strand if display toggled */
Strand blxWindowGetActiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? REVERSE_STRAND : FORWARD_STRAND;
}

/* Return the inactive strand - reverse strand by default, forward strand if display toggled */
static Strand blxWindowGetInactiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? FORWARD_STRAND : REVERSE_STRAND;
}


/* Returns the list of all sequence groups */
GList *blxWindowGetSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext->sequenceGroups;
}

/* Returns the group that the given sequence belongs to, if any (assumes the sequence
 * is only in one group; otherwise it just returns the first group it finds). */
SequenceGroup *blxWindowGetSequenceGroup(GtkWidget *blxWindow, const SequenceStruct *seqToFind)
{
  SequenceGroup *result = NULL;
  
  /* Loop through all the groups until we find this sequence in one */
  GList *groupItem = blxWindowGetSequenceGroups(blxWindow);
  for ( ; groupItem; groupItem = groupItem->next)
    {
      /* See if our sequence struct is in this group's list */
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      GList *foundItem = g_list_find(group->seqList, seqToFind);
      
      if (foundItem)
	{
	  result = group;
	  break;
	}
    }
  
  return result;
}


/***********************************************************
 *                        Selections                       *
 ***********************************************************/

GList* blxWindowGetSelectedSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->selectedSeqs : NULL;
}


/* Get the selected sequences as a list of sequence names. The returned list
 * is formatted as newline-separated values. */
static GString* blxWindowGetSelectedSeqNames(GtkWidget *blxWindow)
{
  GList *listItem = blxWindowGetSelectedSeqs(blxWindow);
  GString *result = g_string_new_len(NULL, 50);
  gboolean first = TRUE;
  
  for ( ; listItem; listItem = listItem->next)
    {
      /* Add a separator before the name, unless it's the first one */
      if (!first)
	{
	  g_string_append(result, "\n");
	}
      else
	{
	  first = FALSE;
	}

      const SequenceStruct *seq = (const SequenceStruct*)(listItem->data);
      g_string_append(result, sequenceGetDisplayName(seq));
    }

  return result;
}


/* This function copys the currently-selected sequences' names to the default
 * clipboard. */
void copySelectionToClipboard(GtkWidget *blxWindow)
{
  if (g_list_length(blxWindowGetSelectedSeqs(blxWindow)) < 1)
    {
      messout("No sequences selected.");
    }
  else
    {
      GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);
      setDefaultClipboardText(displayText->str);
      g_string_free(displayText, TRUE);
    }
}


/* Update function to be called whenever the MSP selection has changed */
void blxWindowSelectionChanged(GtkWidget *blxWindow, const gboolean updateTrees)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  /* Unless requested not to, update the tree row selections to match our new selections */
  if (updateTrees)
    {
      callFuncOnAllDetailViewTrees(detailView, deselectAllRows);
      callFuncOnAllDetailViewTrees(detailView, selectRowsForSelectedSeqs);
    }
  
  /* Update the feedback box to tell the user which sequence is selected. */
  updateFeedbackBox(detailView);
  
  /* Redraw the grids and exon view (because items will change colour with selection) */
  bigPictureRedrawAll(blxWindowGetBigPicture(blxWindow));
  
  /* Copy the selected sequence names to the PRIMARY clipboard */
  GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);
  setPrimaryClipboardText(displayText->str);
  g_string_free(displayText, TRUE);
}


/* Call this function to select the given match sequence */
void blxWindowSelectSeq(GtkWidget *blxWindow, SequenceStruct *seq, const gboolean updateTrees)
{
  if (!blxWindowIsSeqSelected(blxWindow, seq))
    {
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      blxContext->selectedSeqs = g_list_append(blxContext->selectedSeqs, seq);
      blxWindowSelectionChanged(blxWindow, updateTrees);
    }
}


/* Utility function to set the list of selected sequences */
void blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList)
{
  blxWindowDeselectAllSeqs(blxWindow, FALSE);

  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->selectedSeqs = seqList;

  blxWindowSelectionChanged(blxWindow, TRUE);
}


/* Call this function to deselect the given sequence */
void blxWindowDeselectSeq(GtkWidget *blxWindow, SequenceStruct *seq, const gboolean updateTrees)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  /* See if it's in the list and, if so, get a pointer to the list element */
  GList *foundSeq = g_list_find(blxContext->selectedSeqs, seq);

  if (foundSeq)
    {
      blxContext->selectedSeqs = g_list_remove(blxContext->selectedSeqs, foundSeq->data);
      blxWindowSelectionChanged(blxWindow, updateTrees);
    }
}


/* Call this function to deselect all sequences */
void blxWindowDeselectAllSeqs(GtkWidget *blxWindow, const gboolean updateTrees)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  if (g_list_length(blxContext->selectedSeqs) > 0)
    {
      g_list_free(blxContext->selectedSeqs);
      blxContext->selectedSeqs = NULL;
      blxWindowSelectionChanged(blxWindow, updateTrees);
    }
}


/* Returns true if the given sequence is selected */
gboolean blxWindowIsSeqSelected(GtkWidget *blxWindow, const SequenceStruct *seq)
{
  GList *foundItem = NULL;
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext)
    {
      foundItem = g_list_find(blxContext->selectedSeqs, seq);
    }
  
  return (foundItem != NULL);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set various properties for the blixem window */
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
static void calcID(MSP *msp, BlxViewContext *bc)
{
  const gboolean sForward = (mspGetMatchStrand(msp) == FORWARD_STRAND);
  const gboolean qForward = (mspGetRefStrand(msp) == FORWARD_STRAND);
  
  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  msp->id = UNSET_INT;
  
  if (mspIsBlastMatch(msp) && msp->sseq && msp->sseq != bc->paddingSeq)
    {
      msp->id = 0;
      
      /* Note that this will reverse complement the ref seq if it is the reverse 
       * strand. This means that where there is no gaps array the comparison is trivial
       * as coordinates can be ignored and the two sequences just whipped through. */

      char *refSeqSegment = getSequenceSegment(bc,
					       bc->refSeq,
					       msp->qstart, 
					       msp->qend, 
					       mspGetRefStrand(msp), 
					       BLXSEQ_DNA, /* msp q coords are always on the dna sequence */
					       mspGetRefFrame(msp, bc->seqType),
					       bc->displayRev,
					       !qForward,
					       TRUE,
					       TRUE);
      
      if (!refSeqSegment)
	{
	  messout ( "calcID failed: Don't have genomic sequence %d - %d, requested for match sequence '%s' (match coords = %d - %d)\n", msp->qstart, msp->qend, msp->sname, msp->sstart, msp->send);
	  return;
	}
      
      /* We need to find the number of characters that match out of the total number */
      int numMatchingChars = 0;
      int totalNumChars = 0;
      
      if (!(msp->gaps) || arrayMax(msp->gaps) == 0)
	{
	  /* Ungapped alignments. */
	  totalNumChars = (qSeqMax - qSeqMin + 1) / bc->numFrames;
	  
	  if (bc->blastMode == BLXMODE_TBLASTN || bc->blastMode == BLXMODE_TBLASTX)
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
	  if (bc->blastMode == BLXMODE_TBLASTN)
	    {
	      printf("not implemented yet\n") ;
	    }
	  else if (bc->blastMode == BLXMODE_TBLASTX)
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
                  int q_start = qForward ? (qRangeMin - qSeqMin) / bc->numFrames : (qSeqMax - qRangeMax) / bc->numFrames;
		  
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


/* Calculate the ID, q coordinates and q frame for the given MSP and store 
 * them in the MSP struct. Returns the calculated ID (or UNSET_INT if this msp
 * is not a blast match). */
static void calcMspData(MSP *msp, BlxViewContext *bc)
{
  /* Convert the input coords (which are 1-based within the ref sequence section
   * that we're dealing with) to "real" coords (i.e. coords that the user will see). */
  int offset = bc->refSeqRange.min - 1;
  msp->qstart = msp->qstart + offset;
  msp->qend = msp->qend + offset;
  
  /* Gap coords are also 1-based, so convert those too */
  if (msp->gaps && arrayMax(msp->gaps) > 0)
    {
      int i = 0;
      for ( ; i < arrayMax(msp->gaps) ; i++)
	{
	  SMapMap *curRange = arrp(msp->gaps, i, SMapMap);
	  curRange->r1 = curRange->r1 + offset;
	  curRange->r2 = curRange->r2 + offset;
	}
    }
  
  /* Calculate the q frame that this MSP should appear in; that is, the frame in which
   * the first coord of the match will be base 1. (We can find the base number within
   * frame 1 and the required frame number is simply the same as that.) */
  /* to do: do this for exons as well; we need more info though because exons don't
   * necessarily start at base1 in their frame. */
  if (bc->seqType == BLXSEQ_PEPTIDE && mspIsBlastMatch(msp))
    {
      const gboolean reverseStrand = (mspGetRefStrand(msp) == REVERSE_STRAND);
      
      int frame = UNSET_INT;
      convertDnaIdxToDisplayIdx(msp->qstart, bc->seqType, 1, bc->numFrames, reverseStrand, &bc->refSeqRange, &frame);
      char *frameStr = convertIntToString(frame);

      if (frameStr[0] != msp->qframe[2])
	{
	  printf("Warning: calculated match frame as %c but frame in input file is %c. Sequence %s [%d - %d]\n", frameStr[0], msp->qframe[2], msp->sname, msp->sstart, msp->send);
	}  
      
      msp->qframe[2] = frameStr[0];
      g_free(frameStr);
    }
  
  /* Calculate the ID */
  if (mspIsBlastMatch(msp))
    {
      calcID(msp, bc);
    }
}


static int calculateMspData(MSP *mspList, BlxViewContext *bc)
{
  MSP *msp = mspList;
  int lowestId = UNSET_INT;
  
  for ( ; msp; msp = msp->next)
    {
      calcMspData(msp, bc);
      
      if (mspIsBlastMatch(msp) && (lowestId == UNSET_INT || msp->id < lowestId))
	{
	  lowestId = msp->id;
	}
    }
  
  return lowestId;
}


/* Create the main blixem window */
GtkWidget* createBlxWindow(CommandLineOptions *options, const char *paddingSeq)
{
  /* Get the range of the reference sequence. If this is a DNA sequence but our
   * matches are peptide sequences, we must convert to the peptide sequence. */
  const int refSeqLen = (int)strlen(options->refSeq);
  IntRange refSeqRange = {options->refSeqOffset + 1, options->refSeqOffset + refSeqLen};
  IntRange fullDisplayRange = {refSeqRange.min, refSeqRange.max};
  
  if (options->seqType == BLXSEQ_PEPTIDE)
    {
      fullDisplayRange.min = convertDnaIdxToDisplayIdx(refSeqRange.min, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);
      fullDisplayRange.max = convertDnaIdxToDisplayIdx(refSeqRange.max, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);
    }
  
  printf("Reference sequence [%d - %d], display range [%d - %d]\n", refSeqRange.min, refSeqRange.max, fullDisplayRange.min, fullDisplayRange.max);
  
  /* Convert the start coord (which is 1-based and on the DNA sequence) to display
   * coords (which take into account the offset and may also be peptide coords) */
  int startCoord = options->startCoord + options->refSeqOffset;
  if (options->seqType == BLXSEQ_PEPTIDE)
    {
      startCoord = convertDnaIdxToDisplayIdx(startCoord, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);
    }
  
  
  /* Create the main blixem window */
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
  
  BlxViewContext *blxContext = blxWindowCreateContext(options, &refSeqRange, &fullDisplayRange, paddingSeq);
  const int lowestId = calculateMspData(options->mspList, blxContext);
  
  GtkWidget *fwdStrandGrid = NULL, *revStrandGrid = NULL;
  
  GtkWidget *bigPicture = createBigPicture(window,
					   vbox, 
					   &fwdStrandGrid, 
					   &revStrandGrid,
					   options->bigPictZoom,
					   lowestId);

  GtkWidget *detailView = createDetailView(window,
					   vbox, 
					   detailAdjustment, 
					   fwdStrandGrid, 
					   revStrandGrid,
					   options->mspList,
					   options->blastMode,
					   options->seqType,
					   options->numFrames,
					   options->refSeqName,
					   startCoord,
					   options->sortInverted,
					   options->initSortMode);
  
  /* Create a custom scrollbar for scrolling the sequence column and put it at the bottom of the window */
  GtkWidget *scrollBar = createDetailViewScrollBar(detailAdjustment, detailView);
  gtk_box_pack_start(GTK_BOX(vbox), scrollBar, FALSE, TRUE, 0);
  
  /* Set required data for the blixem window */
  blxWindowCreateProperties(options,
			    blxContext,
			    window, 
			    bigPicture, 
			    detailView, 
			    mainmenu,
			    &refSeqRange, 
			    &fullDisplayRange,
			    paddingSeq);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressBlxWindow), mainmenu);
  g_signal_connect(G_OBJECT(window), "key-press-event", G_CALLBACK(onKeyPressBlxWindow), NULL);
  
  
  /* Add the MSP's to the trees and sort them by the initial sort mode. This must
   * be done after all widgets have been created, because it accesses their properties.*/
  detailViewAddMspData(detailView, options->mspList);
  detailViewSortByType(detailView, options->initSortMode);

  /* Set the detail view font (again, this accesses the widgets' properties). */
  updateDetailViewFontDesc(detailView);

  /* Calculate the number of vertical cells in the grids (again, requires properties) */
  calculateNumVCells(bigPicture);


  /* Realise the widgets */
  printf("Starting Blixem\n");
  gtk_widget_show_all(window);
  

  /* Set the initial column widths. (This must be called after the widgets are 
   * realised because it causes the scroll range to be updated, which in turn causes
   * the big picture range to be set. The widgets must be realised before this because
   * the initial big picture range depends on the detail view range, which is calculated
   * from its window's width, and this will be incorrect if it has not been realised.) */
  callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns);
  resizeDetailViewHeaders(detailView);
  
  return window;
}
