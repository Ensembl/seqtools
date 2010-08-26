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
#include <string.h>
#include <ctype.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      1   /* used to change the default border width around the blixem window */
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2   /* used to start with a smaller font than the default widget font */
#define DEFAULT_SCROLL_STEP_INCREMENT	 5    /* how many bases the scrollbar scrolls by for each increment */
#define DEFAULT_WINDOW_WIDTH_FRACTION	 0.9  /* what fraction of the screen size the blixem window width defaults to */
#define DEFAULT_WINDOW_HEIGHT_FRACTION	 0.6  /* what fraction of the screen size the blixem window height defaults to */
#define MATCH_SET_GROUP_NAME		 "Match set"

#define LOAD_DATA_TEXT                   "Load optional data"


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



#define HELP_TEXT1 "\
<b><big>BLIXEM</big></b>\n\
<b>BL</b>ast matches <b>I</b>n an <b>X</b>-windows <b>E</b>mbedded <b>M</b>ultiple alignment\n\
\n\
\n\
<span foreground=\"blue\">\
<b><big>What's new</big></b>\n\
\t•\t<b><i>Status bar</i></b>: There is now a status bar at the bottom of the window, which shows non-critical messages and warnings.\n\
\t•\t<b><i>Moused-over item feedback</i></b>: There is now an additional feedback area on the toolbar that shows information about the currently moused-over item.\n\
\t•\t<b><i>Grid scale</i></b>: The grid scale now has a finer level of granularity; you can set %ID per cell to as low as 0.1 (see the Settings dialog).\n\
</span>\
\n\
\n\
<b><big>Fetch</big></b>\n\
\t•\tDouble-click a row to fetch a sequence.\n\
\t•\tIn the Settings dialog, click the 'Load optional data' button to fetch additional information such as organism and tissue-type. This data can be viewed by enabling the corresponding column and is also displayed in the feedback area of the toolbar when a sequence is moused-over.\n\
\n\
\n\
<b><big>Main menu</big></b>\n\
Right-click anywhere in the Blixem window to pop up the main menu.  The menu options are:\n\
\n\
\t•\tQuit: Quit Blixem.\n\
\t•\tHelp: Display this help.\n\
\t•\tPrint: Print the display\n\
\t•\tSettings: Edit settings.\n\
\t•\tView: Display the View dialog box. This allows you to show/hide parts of the display.\n\
\t•\tCreate Group: Create a group of sequences.\n\
\t•\tEdit Groups: Edit properties for groups.\n\
\t•\tDeselect all: Deselect all sequences.\n\
\t•\tDotter: Run Dotter and/or set Dotter parameters.\n\
\t•\tFeature series selection tool: ???\n\
\n\
\n\
<b><big>Scrolling</big></b>\n\
\t•\tMiddle-click in the big picture to scroll to an area.\n\
\t•\tMiddle-click in the detail view to centre on a base.\n\
\t•\tUse the horizontal scrollbar at the bottom of the window to scroll the detail view.\n\
\t•\tUse the mouse-wheel to scroll up/down in the list of alignments.\n\
\t•\tUse the mouse-wheel to scroll the alignments left/right (if your mouse has a horizontal scroll-wheel).\n\
\t•\tCtrl-left and Ctrl-right arrow keys scroll to the start/end of the previous/next match (limited to currently-selected sequences, if any are selected).\n\
\t•\tThe Home/End keys scroll to the start/end of the display.\n\
\t•\tCtrl-Home and Ctrl-End scroll to the start/end of the currently-selected alignments (or to the first/last alignment if none are selected).\n\
\t•\tThe ',' (comma) and '.' (full-stop) keys scroll the display one nucleotide to the left/right.\n\
\t•\tPressing Ctrl and ',' (comma) or '.' (full-stop) scrolls the display one page to the left/right.\n\
\t•\tYou can scroll to a specific position using the Go-to button on the toolbar, or by pressing the 'p' shortcut key.\n\
\n\
\n\
<b><big>Selections</big></b>\n\
\t•\tYou can select a sequence by clicking on its row in the alignment list.  Selected sequences are highlighted in cyan in the big picture.\n\
\t•\tYou can select a sequence by clicking on it in the big picture.\n\
\t•\tThe name of the sequence you selected is displayed in the feedback box on the toolbar.  If there are multiple alignments for the same sequence, all of them will be selected.\n\
\t•\tYou can select multiple sequences by holding down the Ctrl or Shift keys while selecting rows.\n\
\t•\tYou can deselect a single sequence by Ctrl-clicking on its row.\n\
\t•\tYou can deselect all sequences by right-clicking and selecting 'Deselect all', or with the Shift-Ctrl-A keyboard shortcut.\n\
\t•\tYou can move the selection up/down a row using the up/down arrow keys.\n\
\t•\tIf the 'Show Splice Sites' option is on (see the 'Settings' dialog), splice sites will be highlighted on the reference sequence for the currently-selected alignment(s).\n\
\n\
\t•\tYou can select a nucleotide/peptide by middle-clicking on it in the detail view.  This selects the entire column at that index, and the coordinate number on the query sequence is shown in the feedback box.  (The coordinate on the subject sequence is also shown if a subject sequence is selected.)\n\
\t•\tFor protein matches, when a peptide is selected, the three nucleotides for that peptide (for the current reading frame) are highlighted in the header in green.  The current reading frame is whichever alignment list currently has the focus - click in a different list to change the reading frame.  The darker highlighting indicates the specific nucleotide that is currently selected (i.e. whose coordinate is displayed in the feedback box).\n\
\t•\tYou can move the selection to the previous/next nucleotide using the left and right arrow keys.\n\
\t•\tYou can move the selection to the previous/next peptide by holding Shift while using the left and right arrow keys.\n\
\t•\tYou can move the selection to the start/end of the previous/next matchb by holding Ctrl while using the left and right arrow keys (limited to just the selected sequences if any are selected).\n\
\t•\tTo select a coordinate without the display re-centering on it, hold down Ctrl as you middle-click.\n\
\n\
\n\
<b><big>Zooming</big></b>\n\
\t•\tZoom in to the currently-selected region in the big picture using the +/- buttons in the top-left corner of the window, or using the Ctrl '=' and Ctrl '-' shortcut keys.  The 'Whole' button zooms out to show the full length of the query sequence.\n\
\t•\tZoom in/out of the detail view using the +/- buttons on the toolbar, or using the '=' and '-' shortcut keys.\n\
\n\
\n\
<b><big>Copy and paste</big></b>\n\
\t•\tWhen sequence(s) are selected, their names are copied to the selection buffer and can be pasted by middle-clicking.\n\
\t•\tTo paste sequence names from the selection buffer, hit the 'f' keyboard shortcut. Blixem will select the sequences and jump to the start of the selection.\n\
\t•\tTo copy sequence name(s) to the default clipboard, select the sequence(s) and hit Ctrl-C. Sequence names can then be pasted into other applications using Ctrl-V.\n\
\t•\tThe clipboard text can also be pasted into Blixem using Ctrl-V. If the clipboard contains valid sequence names, those sequences will be selected and the display will jump to the start of the selection.\n\
\t•\tNote that text from the feedback box and some text labels (e.g. the reference sequence start/end coords) can be copied by selecting it with the mouse and then hitting Ctrl-C.\n\
\t•\tText can be pasted into dialog box text entry boxes using Ctrl-V (or middle-clicking to paste from the selection buffer).\n\
\n\
\n\
<b><big>Sorting</big></b>\n\
\t•\tThe alignments can be sorted by selecting the column you wish to sort by from the drop-down box on the toolbar.\n\
\t•\tYou can reverse the sort order by selecting the 'Invert sort order' option on the Settings dialog.\n\
\t•\tIf you sort by one of the optional columns (e.g. organism or tissue-type) the sort will have no effect until the data for that column is loaded. Load the data by clicking the 'Load optional data' button in the 'Columns' section of the Settings dialog.\n\
\t•\tYou can place alignments in a group and then sort by group to keep specific sequences always on top.  See the Groups section for more details.\n\
\n\
\n\
<b><big>Finding</big></b>\n\
\t•\tClick the find icon on the toolbar or press the Ctrl-F shortcut to open the find dialog.\n\
\t•\tTo search for sequences by name, use either the 'Sequence name search' or 'Sequence name list' search box.  Enter the sequence name you wish to search for and hit 'OK', 'Forward' or 'Back'.\n\
\t•\tTo search for a string of nucleotides in the reference sequence, use the 'DNA search' box.  Enter the string of nucleotides and hit 'OK', 'Forward' or 'Back'.  Note that the search is only performed on the active strand: use the 'Toggle strand' button on the toolbar or hit the 't' shortcut key to toggle which strand is active.\n\
\t•\tUsing the 'OK' button will search for the first match in the entire Blixem range.\n\
\t•\tUsing the 'Forward' or 'Back' buttons will search for the next/previous match from the current position.\n\
\t•\tAfter performing a find, you can repeat the search by hitting the F3 key (or Shift-F3 to search backwards).\n\
\t•\tYou can perform a fast find on sequence name(s) in the selection buffer by hitting the 'f' shortcut key.\n\
\n\
\n\
<b><big>Feedback about current items</big></b>\n\
There are two feedback boxes on the toolbar that display information about current items:\n\
\t•\tThe first shows info about the currently selected-alignment and/or the currently-selected index.\n\
\t•\tThe second shows info about the currently moused-over alignemnt. As a minimum this will display the sequence name. It will also show additional data such as organism and tissue type if available - if not already loaded, this data can be loaded by clicking the 'Load optional data' button on the Settings dialog.\n\
\n\
\n\
<b><big>Messages/Warnings</big></b>\n\
\t•\tThe status bar at the bottom of the window shows all messages and warnings, including non-critical messages.\n\
\t•\tCritical messages are also displayed in a popup dialog or message list.\n\
\t•\tThe user may switch between popup dialogs or the message list by selecting the 'Switch to scrolled message window' on the popup dialog or the 'Switch to popup messages' check box on the message list.\n\
\n\
\n\
<b><big>Groups</big></b>\n\
Alignments can be grouped together so that they can be sorted/highlighted/hidden etc.\n\
\n\
<b>Creating a group from a selection:</b>\n\
\t•\tSelect the sequences you wish to include in the group by left-clicking their rows in the detail view.  Multiple rows can be selected by holding the Ctrl or Shift keys while clicking.\n\
\t•\tRight-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Note that Ctrl-G will also shortcut to here if no groups currently exist.)\n\
\t•\tEnsure that the 'From selection' radio button is selected, and click 'OK'.\n\
\n\
<b>Creating a group from a sequence name:</b>\n\
\t•\tRight-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Or Ctrl-G if no groups currently exist.)\n\
\t•\tSelect the 'From name' radio button and enter the name of the sequence in the box below.  You may use the following wildcards to search for sequences: '*' for any number of characters; '?' for a single character.  For example, searching for '*X' will find all sequences whose name ends in 'X' (i.e. all exons).\n\
\t•\tClick 'OK'.\n\
\n\
<b>Creating a group from sequence name(s):</b>\n\
\t•\tRight-click and select 'Create Group', or use the Shift-Ctrl-G shortcut key. (Or Ctrl-G if no groups currently exist.)\n\
\t•\tSelect the 'From name(s)' radio button.\n\
\t•\tEnter the sequence name(s) in the text box.\n\
\t•\tYou may use the following wildcards in a sequence name: '*' for any number of characters; '?' for a single character.  (For example, searching for '*X' will find all sequences whose name ends in 'X', i.e. all exons).\n\
\t•\tYou may search for multiple sequence names by separating them with the following delimiters: newline, comma or semi-colon.\n\
\t•\tYou may paste sequence names directly from ZMap: click on the feature in ZMap and then middle-click in the text box on the Groups dialog.  Grouping in Blixem works on the sequence name alone, so the feature coords will be ignored.\n\
\t•\tClick 'OK'.\n\
\n\
<b>Creating a temporary 'match-set' group from the current selection:</b>\n\
\t•\tYou can quickly create a group from a current selection (e.g. selected features in ZMap) using the 'Toggle match set' option.\n\
\t•\tTo create a match-set group, select the required items (e.g. in ZMap) and then select 'Toggle match set' from the right-click menu in Blixem, or hit the 'g' shortcut key.\n\
\t•\tTo clear the match-set group, choose the 'Toggle match set' option again, or hit the 'g' shortcut key again.\n\
\t•\tWhile it exists, the match-set group can be edited like any other group, via the 'Edit Groups' dialog.\n\
\t•\tIf you delete the match-set group from the 'Edit Groups' dialog, all settings (e.g. highlight color) will be lost. To maintain these settings, clear the group using the 'Toggle match set' menu option (or 'g' shortcut key) instead.\n\
\n\
<b>Editing groups:</b>\n\
To edit a group, right-click and select 'Edit Groups', or use the Ctrl-G shortcut key. You can change the following properties for a group:\n\
\t•\tName: you can specify a more meaningful name to help identify the group.\n\
\t•\tHide: tick this box to hide the alignments in the alignment lists.\n\
\t•\tHighlight: tick this box to highlight the alignments.\n\
\t•\tColor: the color the group will be highlighted in, if 'Highlight' is enabled.  The default color for all groups is red, so you may wish to change this if you want different groups to be highlighted in different colors.\n\
\t•\tOrder: when sorting by Group, alignments in a group with a lower order number will appear before those with a higher order number (or vice versa if sort order is inverted). Alignments in a group will appear before alignments that are not in a group.\n\
\t•\tTo delete a single group, click on the 'Delete' button next to the group you wish to delete.\n\
\t•\tTo delete all groups, click on the 'Delete all groups' button.\n\
\n\
\n\
<b><big>Dotter</big></b>\n\
\t•\tTo start Dotter, or to edit the Dotter parameters, right-click and select 'Dotter' or use the Ctrl-D keyboard shortcut.	The Dotter settings dialog will pop up.\n\
\t•\tTo run Dotter with the default (automatic) parameters, just hit RETURN, or click the 'Execute' button.\n\
\t•\tTo enter manual parameters, click the 'Manual' radio button and enter the values in the 'Start' and 'End' boxes.\n\
\t•\tTo revert to the last-saved manual parameters, click the 'Last saved' button.\n\
\t•\tTo revert back to automatic parameters, click the 'Auto' radio button.\n\
\t•\tTo save the parameters without running Dotter, click Save and then Cancel'.\n\
\t•\tTo save the parameters and run Dotter, click 'Execute'.\n\
\n\
\n\
<b><big>Settings</big></b>\n\
The settings menu can be accessed by right-clicking and selecting Settings, or by the shortcut Ctrl-S.\n\
\t•\t<b>Squash Matches</b>:\t\tGroup multiple alignments from the same sequence together into the same row in the detail view, rather than showing them on separate rows.\n\
\t•\t<b>Show Unaligned Sequence</b>:Show any additional unaligned parts of the match sequence at the start/end of the alignment.  Specify the number of additional bases to show in 'Limit to... bases', or untick this option to show all of the unaligned sequence.  Note that this option does not work with the 'Squash Matches' option, so it will not do anything if 'Squash Matches' is on.\n\
\t•\t<b>Invert Sort Order</b>:\t\tReverse the default sort order. (Note that some columns sort ascending by default (e.g. name, start, end) and some sort descending (score and ID). This option reverses that sort order.)\n\
\t•\t<b>Highlight Differences</b>:When this option is set, matching bases are blanked out and mismatches are highlighted, making it easier to see where alignments differ from the reference sequence.\n\
\t•\t<b>Show SNP Track</b>:Shows the SNP track.\n\
\t•\t<b>Show Splice Sites</b>:Shows splice sites for the currently-selected alignment(s).  Splice sites are highlighted on the reference sequence in green (for canonical) or red (for non-canonical).  Blixem identifies GT-AG, GC-AG and AT-AC introns as canonical.\n\
\t•\t<b>Columns</b>:\t\t\tEdit the width of columns in pixels.  Set the width to 0 to hide a column. Click the 'Load optional data' button to load the data for the optional columns such as organism and tissue-type - then set the width of these columns to non-zero values to view the data. Once optional data is loaded you can also sort by it. Note that optional data is loaded on startup for DNA matches but not for protein matches, because the latter can be slow.\n\
\t•\t<b>Grid properties</b>:\t\t\tSet the maximum/minimum %ID values shown in the big picture.  Expand or contract the grid scale by adjusting '%ID per cell'.  The ID per cell can be set to as low as 0.1.  Note that setting a low ID-per-cell can result in a large number of cells and hence a very large grid, and Blixem is not currently very clever about dealing with this.\n\
\n\
\n\
<b><big>Color key</big></b>\n\
In the detail view, the following colors and symbols have the following meanings:\n\
\n\
\t<span background=\"yellow\">Yellow</span>\t\t\t\t\t\tQuery sequence\n\
\t<span background=\"cyan\">Cyan</span>\t\t\t\t\t\tIdentical match\n\
\t<span background=\"lightslateblue\">Violet</span>\t\t\t\t\t\tSimilar match\n\
\t<span background=\"grey\">Grey</span>\t\t\t\t\t\tMismatch\n\
\t<span background=\"grey\"><b> . </b></span>\t\t\t\t\t\t\tDeletion\n\
\t<span background=\"yellow\"> </span> Yellow line\t\t\t\t\tInsertion\n\
\t<span background=\"palegreen\">         </span> Empty green box\t\tCoding (CDS) exon\n\
\t<span background=\"lightcoral\">         </span> Empty red box\t\t\tNon-coding (UTR) exon\n\
\t<span background=\"royalblue\">|</span> Blue line\t\t\t\t\tStart boundary of an exon\n\
\t<span background=\"darkblue\">|</span> Dark blue line\t\t\t\tEnd boundary of an exon\n\
\t<span background=\"lightskyblue\">Blue</span> in DNA header\t\t\tNucleotides of the selected codon. Darker blue indicates the specific nucleotide displayed in the feedback box.\n\
\t<span background=\"salmon\">Red</span> in 3-frame translation\t\tSTOP codon\n\
\t<span background=\"lawngreen\">Green</span> in 3-frame translation\tMET codon\n\
\n\
\n\
<b><big>Keyboard shortcuts</big></b>\n\
\t•\t<b>Ctrl-Q</b>    Quit\n\
\t•\t<b>Ctrl-H</b>    Help\n\
\t•\t<b>Ctrl-P</b>    Print\n\
\t•\t<b>Ctrl-S</b>    'Settings' menu\n\
\n\
\t•\t<b>Ctrl-F</b>    Open the Find dialog\n\
\t•\t<b>f</b>    Perform a find on sequence name(s) in the current selection buffer\n\
\n\
\t•\t<b>Ctrl-G</b>    Edit groups (or create a group if none currently exist)\n\
\t•\t<b>Shift-Ctrl-G</b>    Create group\n\
\t•\t<b>g</b>    Toggle the 'match set' Group\n\
\n\
\t•\t<b>Shift-Ctrl-A</b>    Deselect all alignments\n\
\t•\t<b>Escape</b>    Deselect the currently-selected base index\n\
\n\
\t•\t<b>Ctrl-D</b>    Dotter\n\
\n\
\t•\t<b>Left</b>/<b>right</b> arrow keys    Move the currently-selected coordinate one place to the left/right\n\
\t•\t<b>Shift-Left</b>/<b>right</b>    Same as left/right arrow keys, but for proteins it scrolls by a single nucleotide, rather than an entire codon.\n\
\t•\t<b>Ctrl-Left</b>/<b>right</b>    Scroll to the start/end of the previous/next alignment (limited to just the selected sequences, if any are selected).\n\
\t•\t<b>Home</b>/<b>End</b>    Scroll to the start/end of the display.\n\
\t•\t<b>Ctrl-Home</b>/<b>End</b>    Scroll to the first/last alignment (limited to just the selected sequences, if any are selected).\n\
\t•\t<b>,</b> (comma)    scroll left one coordinate\n\
\t•\t<b>.</b> (period)    scroll right one coordinate\n\
\n\
\t•\t<b>=</b> (equals)    zoom in detail view\n\
\t•\t<b>-</b> (subtract)    zoom out detail view\n\
\t•\t<b>Ctrl</b>-<b>=</b>     zoom in big picture\n\
\t•\t<b>Ctrl</b>-<b>-</b>      zoom out big picture\n\
\t•\t<b>Shift</b>-<b>Ctrl</b>-<b>-</b>      zoom out big picture to whole width\n\
\n\
\t•\t<b>V</b>    'View' menu (for toggling visibility)\n\
\t•\t<b>1</b>, <b>2</b>, <b>3</b>    These number keys toggle visibility of the 1st, 2nd (and 3rd, for protein matches) alignment list.\n\
\t•\t<b>Ctrl</b>-<b>1</b>, <b>Ctrl</b>-<b>2</b>    This toggles visibility of the 1st and 2nd big picture grid.\n\
\t•\t<b>Shift</b>-<b>Ctrl</b>-<b>1</b>, <b>Shift</b>-<b>Ctrl</b>-<b>2</b>    This toggles visibility of the 1st and 2nd exon views.\n\
\n\
\t•\t<b>p</b>    Go to position\n\
\t•\t<b>t</b>    Toggle the active strand\n\
"

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
static void			  onUpdateBackgroundColor(GtkWidget *blxWindow);

static void			  onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
static void			  onDrawPage(GtkPrintOperation *operation, GtkPrintContext *context, gint pageNum, gpointer data);
static void			  onDestroyBlxWindow(GtkWidget *widget);

static BlxStrand		  blxWindowGetInactiveStrand(GtkWidget *blxWindow);

static void			  destroyBlxColor(BlxColor *blxCol);

static void			  onButtonClickedDeleteGroup(GtkWidget *button, gpointer data);
static void			  blxWindowGroupsChanged(GtkWidget *blxWindow);
static GtkRadioButton*		  createRadioButton(GtkBox *box, GtkRadioButton *existingButton, const char *mnemonic, const gboolean isActive, const gboolean createTextEntry, const gboolean multiline, BlxResponseCallback callbackFunc, GtkWidget *blxWindow);
static void			  getSequencesThatMatch(gpointer listDataItem, gpointer data);
static GList*			  getSeqStructsFromText(GtkWidget *blxWindow, const char *inputText);

static void			  createCheckButton(GtkBox *box, const char *mnemonic, const gboolean isActive, GCallback callback, gpointer data);
static void			  blxWindowSetUsePrintColors(GtkWidget *blxWindow, const gboolean usePrintColors);
static gboolean			  blxWindowGetUsePrintColors(GtkWidget *blxWindow);

static void                       blxWindowFindDnaString(GtkWidget *blxWindow, const char *inputSearchStr, const int startCoord, const gboolean searchLeft, const gboolean findAgain, GError **error);
static GList*                     findSeqsFromList(GtkWidget *blxWindow, const char *inputText, const gboolean findAgain, GError **error);
static int                        getSearchStartCoord(GtkWidget *blxWindow, const gboolean startBeginning, const gboolean searchLeft);
static GList*                     findSeqsFromName(GtkWidget *blxWindow, const char *inputText, const gboolean findAgain, GError **error);
static GtkWidget*                 dialogChildGetBlxWindow(GtkWidget *child);



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

/* Gets the dialog widget for the given dialog id. Returns null if the widget has not
 * been created yet. */
GtkWidget* getPersistentDialog(const BlxViewContext *bc, const BlxDialogId dialogId)
{
  GtkWidget *result = NULL;
  
  if (bc->dialogList[dialogId])
    {
      result = bc->dialogList[dialogId];
    }

  return result;
}

/* Add a newly-created dialog to the list of persistent dialogs. The dialog should not
 * exist in the list yet. */
void addPersistentDialog(BlxViewContext *bc, const BlxDialogId dialogId, GtkWidget *widget)
{
  if (dialogId == BLXDIALOG_NOT_PERSISTENT)
    {
      g_warning("Code error: cannot add a dialog with ID %d. Dialog will not be persistent.\n", dialogId);
    }
  else
    {
      if (bc->dialogList[dialogId])
        {
          g_warning("Creating a dialog that already exists. Old dialog will be destroyed. Dialog ID=%d.\n", dialogId);
          gtk_widget_destroy(bc->dialogList[dialogId]);
          bc->dialogList[dialogId] = NULL;
        }

      bc->dialogList[dialogId] = widget;
    }
}


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
			  const BlxStrand strand,
			  const BlxSeqType inputCoordType,
			  const int frame,
			  const gboolean displayRev,
			  const gboolean reverseResult,
			  const gboolean allowComplement,
			  const gboolean allowTranslate,
			  GError **error)
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
	    g_warning("Requested query sequence %d - %d out of available range: %d - %d. Input coords on peptide sequence were %d - %d\n", qMin, qMax, bc->refSeqRange.min, bc->refSeqRange.max, coord1, coord2);
	  else
	    g_warning("Requested query sequence %d - %d out of available range: %d - %d\n", qMin, qMax, bc->refSeqRange.min, bc->refSeqRange.max);
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
      
      if (strand == BLXSTRAND_FORWARD)
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
		  g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_SEGMENT, 
			      "Error getting sequence segment: Failed to complement the reference sequence for the range %d - %d.\n", 
			      qMin, qMax);
		  
		  g_free(segment);
		  g_free(result);
		  return NULL;
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
	      g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_SEGMENT,
			  "Error getting the sequence segment: Failed to translate the DNA sequence for the reference sequence range %d - %d\n", 
			  qMin, qMax) ;
	      
	      g_free(segment);
	      g_free(result);
	      return NULL;
	    }
	}
    }
  
  if (!result)
    {
      g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_SEGMENT, "Failed to find sequence segment for the range %d - %d\n", qMin, qMax);
    }
  
  return result;
}


/* Returns the order number of the group that this sequence belongs to, or
 * UNSET_INT if it does not belong to a group. */
int sequenceGetGroupOrder(GtkWidget *blxWindow, const BlxSequence *seq)
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


/* Move the current row selection up/down */
static gboolean moveRowSelection(GtkWidget *blxWindow, const gboolean moveUp, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  const BlxStrand activeStrand = detailViewGetSelectedStrand(detailView);
  
  GtkWidget *tree = detailViewGetTree(detailView, activeStrand, activeFrame);
  return treeMoveRowSelection(tree, moveUp, shiftModifier);
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
	{
	  zoomWholeBigPicture(blxWindowGetBigPicture(window));
	}
      else
	{
	  zoomBigPicture(blxWindowGetBigPicture(window), zoomIn);
	}
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

  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  refreshDetailViewHeaders(detailView);
  
  callFuncOnAllDetailViewTrees(detailView, widgetClearCachedDrawable, NULL);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  
  gtk_widget_queue_draw(blxWindow);
}


/* Utility to create a vbox with the given border and pack it into the given box.
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
      gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);
      gtk_container_add(GTK_CONTAINER(frame), vbox);
    }
  else
    {
      gtk_box_pack_start(GTK_BOX(parent), vbox, FALSE, FALSE, 0);
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
  const BlxStrand activeStrand = toggled ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
  
  /* For protein matches, trees are always displayed in frame order (i.e. 1, 2, 3), 
   * so just use the number pressed for the frame, and the active strand for the
   * strand. */
  int frame = number;
  BlxStrand strand = activeStrand;
  
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
	  strand = toggled ? BLXSTRAND_FORWARD : BLXSTRAND_REVERSE;
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


/* Repeat the last find operation. Searches for the next (rightwards) match unless the given
 * modifier is pressed, in which case it searches for the previous (leftwards) match */
static void findAgain(GtkWidget *blxWindow, const gboolean modifier)
{
  GError *error = NULL;

  const int startCoord = getSearchStartCoord(blxWindow, FALSE, modifier);

  /* Try the DNA search. Does nothing if last search was not a DNA search. */
  blxWindowFindDnaString(blxWindow, NULL, startCoord, modifier, TRUE, &error);

  if (error)
    {
      /* DNA search was attempted but not found. Try looping round to the beginning */
      g_error_free(error);
      error = NULL;
      const int newStart = getSearchStartCoord(blxWindow, TRUE, modifier);
      blxWindowFindDnaString(blxWindow, NULL, newStart, modifier, TRUE, &error);
    }
  
  if (!error)
    {
      /* Try the search-from-list search. Returns NULL if last search was not a list search */
      GList *seqList = findSeqsFromList(blxWindow, NULL, TRUE, &error);
      
      if (!seqList && !error)
        {
          /* Try the search-by-name search. Returns NULL if last search was not a name search. */
          seqList = findSeqsFromName(blxWindow, NULL, TRUE, &error);
        }
      
      /* If either the list or name search succeeded, select the prev/next MSP from the 
       * found sequence(s) depending on which direction we're searching. */
      if (seqList)
        {
          blxWindowSetSelectedSeqList(blxWindow, seqList);
          
          if (modifier)
            {
              prevMatch(blxWindowGetDetailView(blxWindow), seqList);
            }
          else
            {
              nextMatch(blxWindowGetDetailView(blxWindow), seqList);
            }
        }
    }
  
  if (error)
    {
      prefixError(error, "Find %s failed. ", (modifier ? "previous" : "next"));
      reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
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
static void createTreeVisibilityButton(GtkWidget *detailView, const BlxStrand strand, const int frame, GtkWidget *container)
{
  /* Some trees may have been removed from the blixem window if they are not on the active 
   * strand, so only show check boxes for those that are in the window (i.e. have a parent). 
   * Note that we get the tree container here, which consists of the tree itself plus any headers etc. */
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (gtk_widget_get_parent(tree))
    {
      const gboolean toggled = detailViewGetDisplayRev(detailView);
      gboolean isActiveStrand = ((strand == BLXSTRAND_FORWARD) != toggled);

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
	  /* All the visible trees should be in the same strand, so just distinguish by frame number. */
	  char formatStr[] = "Show frame _%d";
	  char displayText[strlen(formatStr) + numDigitsInInt(frame) + 1];
	  sprintf(displayText, formatStr, frame);

	  createVisibilityButton(tree, displayText, container);
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
void showViewPanesDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_VIEW;
  GtkWidget *dialog = getPersistentDialog(bc, dialogId);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons("View panes", 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Clear contents and re-create */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
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
  
  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}



/***********************************************************
 *			    Find menu			   *
 ***********************************************************/

static GList* findSeqsFromName(GtkWidget *blxWindow, const char *inputText, const gboolean findAgain, GError **error)
{
  static char *searchStr = NULL;
  
  if (!findAgain)
    {
      g_free(searchStr);
      searchStr = g_strdup(inputText);
    }
  
  if (!searchStr)
    {
      return NULL;
    }
  
  /* Loop through all the sequences and see if the name matches the search string */
  GList *seqList = blxWindowGetAllMatchSeqs(blxWindow);
  SeqSearchData searchData = {searchStr, NULL};
  
  g_list_foreach(seqList, getSequencesThatMatch, &searchData);
  
  if (g_list_length(searchData.matchList) < 1)
    {
      g_set_error(error, BLX_ERROR, BLX_ERROR_STRING_NOT_FOUND, "No sequences found matching text '%s'.\n", searchStr);
    }
  
  return searchData.matchList;
}


/* Utility to extract the contents of a GtkEntry and return it as a string. The result is 
 * owned by the GtkTextEntry and should not be free'd. */
static const char* getStringFromTextEntry(GtkEntry *entry)
{
  const char *result = NULL;
  
  if (!entry || !GTK_WIDGET_SENSITIVE(GTK_WIDGET(entry)))
    {
      g_warning("Could not set search string: invalid text entry box\n");
    }
  else
    {
      result = gtk_entry_get_text(entry);
    }
  
  return result;
}


/* Utility to extract the contents of a GtkTextView and return it as a string. The result is
 * owned by the GtkTextView and should not be free'd. */
static const char* getStringFromTextView(GtkTextView *textView)
{
  if (!textView || !GTK_WIDGET_SENSITIVE(GTK_WIDGET(textView)))
    {
      g_critical("Could not set search string: invalid text entry box\n");
      return NULL;
    }
  
  /* Get the input text from the text buffer and create the group */
  GtkTextBuffer *textBuffer = gtk_text_view_get_buffer(textView);
  
  GtkTextIter start, end;
  gtk_text_buffer_get_bounds(textBuffer, &start, &end);
  
  return gtk_text_buffer_get_text(textBuffer, &start, &end, TRUE);
}


/* Finds all the valid sequences blixem knows about whose names are in the given
 * text, and returns them in a GList of BlxSequences. */
static GList* findSeqsFromList(GtkWidget *blxWindow, const char *inputText, const gboolean findAgain, GError **error)
{
  static char *searchStr = NULL; /* remember last searched-for string for use with 'findAgain' option */
  
  if (!findAgain)
    {
      g_free(searchStr);
      searchStr = g_strdup(inputText);
    }
      
  if (!searchStr)
    {
      return NULL;
    }

  GList *seqList = getSeqStructsFromText(blxWindow, searchStr);

  if (g_list_length(seqList) < 1)
    {
      g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_NAME_NOT_FOUND, "No valid sequence names in search string\n%s\n", searchStr);
    }
  
  return seqList;
}


/* Callback called when requested to find sequences from a sequence name. Selects
 * the sequences and scrolls to the start of the first match in the selection */
static void onFindSeqsFromName(GtkWidget *button, const gint responseId, gpointer data)
{
  const char *inputText = NULL;
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      inputText = getStringFromTextEntry(GTK_ENTRY(data));
    }
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  GError *error = NULL;
  GList *seqList = findSeqsFromName(blxWindow, inputText, FALSE, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      blxWindowSetSelectedSeqList(blxWindow, seqList);
      
      if (responseId == BLX_RESPONSE_FORWARD)
        {
          nextMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else if (responseId == BLX_RESPONSE_BACK)
        {
          prevMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else
        {
          firstMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
    }
}


/* Callback called when requested to find sequences from a given list. Selects
 * the sequences ands scrolls to the start of the first match in the selection. */
static void onFindSeqsFromList(GtkWidget *button, const gint responseId, gpointer data)
{
  const char *inputText = NULL;
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      inputText = getStringFromTextView(GTK_TEXT_VIEW(data));
    }
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  GError *error = NULL;
  GList *seqList = findSeqsFromList(blxWindow, inputText, FALSE, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  if (seqList)
    {
      blxWindowSetSelectedSeqList(blxWindow, seqList);
      
      if (responseId == BLX_RESPONSE_FORWARD)
        {
          nextMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else if (responseId == BLX_RESPONSE_BACK)
        {
          prevMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else
        {
          firstMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
    }
}


/* Search for the given DNA string in the reference sequence. Searches for the next (rightwards)
 * value from the given start coord, unless searchLeft is true in which case it searches leftwards.
 * If findAgain is true it repeats the last DNA search. */
static void blxWindowFindDnaString(GtkWidget *blxWindow, 
                                   const char *inputSearchStr, 
                                   const int refSeqStart,
                                   const gboolean searchLeft, 
                                   const gboolean findAgain, 
                                   GError **error)
{
  /* Remember the last input string for use with findAgain */
  static char *searchStr = NULL;
  
  if (!findAgain)
    {
      /* We must copy the input string because it may not exist if/when we come to do a 'find again' */
      g_free(searchStr);
      searchStr = g_strdup(inputSearchStr);
    }
  
  const int searchStrMax = searchStr ? strlen(searchStr) - 1 : -1;
  
  if (searchStrMax < 0)
    {
      return;
    }

  const int searchStart = searchLeft ? searchStrMax : 0;
  const int searchEnd = searchLeft ? 0 : searchStrMax;
  const int searchStrIncrement = searchLeft ? -1 : 1;
  
  /* Values increase left-to-right in normal display or right-to-left in reversed display */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const gboolean searchForward = (searchLeft == bc->displayRev);
  const int refSeqIncrement = searchForward ? 1 : -1;
  
  /* We'll need to complement ref seq bases if the active strand is the reverse strand */
  const gboolean complement = (blxWindowGetActiveStrand(blxWindow) == BLXSTRAND_REVERSE);
  
  int refSeqIdx = refSeqStart;
  int searchStrIdx = searchStart;
  int matchStart = UNSET_INT;
  
  while (refSeqIdx >= bc->refSeqRange.min && refSeqIdx <= bc->refSeqRange.max && searchStrIdx >= 0 && searchStrIdx <= searchStrMax)
    {
      const char refSeqBase = getRefSeqBase(bc->refSeq, refSeqIdx, complement, &bc->refSeqRange, BLXSEQ_DNA);
      char searchStrBase = convertBaseToCorrectCase(searchStr[searchStrIdx], BLXSEQ_DNA);      
      
      if (refSeqBase == searchStrBase)
        {
          /* The base matches. If it's the first matching base, set the match-start coord (or if we're 
           * searching leftwards, then always set the match-start coord, because the start is actually 
           * the last coord that will be found). Then proceed to the next position in the search string */
          if (matchStart == UNSET_INT)
            {
              matchStart = refSeqIdx;
            }
          
          searchStrIdx += searchStrIncrement;
          refSeqIdx += refSeqIncrement;
        }
      else if (matchStart != UNSET_INT)
        {
          /* We were in a match but this base doesn't match. Reset to the start of the 
           * search string, and start looking again from one base after the place where the last
           * match started. (We need to re-check all bases from there because we're comparing
           * against a different section of the search string now.) */
          searchStrIdx = searchStart;
          refSeqIdx = matchStart + refSeqIncrement;
          matchStart = UNSET_INT;
        }
      else
        {
          refSeqIdx += refSeqIncrement;
        }
    }
  
  /* Undo the last increment, so that we have the final coords of the matching section (if found) */
  refSeqIdx -= refSeqIncrement;
  searchStrIdx -= searchStrIncrement;
  
  /* If we reached the end of the search string, then we matched the whole lot. */
  const gboolean finished = searchStrIdx == searchEnd;
  
  if (matchStart != UNSET_INT && finished)
    {
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      const int frame = 1;
      int baseNum = UNSET_INT;
      
      int result = searchLeft ? refSeqIdx : matchStart;
      result = convertDnaIdxToDisplayIdx(result, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &baseNum);
      
      detailViewSetSelectedBaseIdx(detailView, result, frame, baseNum, TRUE, FALSE);
    }
  else
    {
      g_set_error(error, BLX_ERROR, BLX_ERROR_STRING_NOT_FOUND, "The string '%s' was not found in the reference sequence searching to the %s from coord %d.\n", searchStr, (searchLeft ? "left" : "right"), refSeqStart);
    }
}


/* Get the start coord for a search. If startBeginning is false, this gets the currently-selected display
 * index (shifted by one base so that we don't start searching at the same position as a previous
 * find result) or, if no base index is selected, returns the start coord of the current display range. If
 * startBeginning is true, just start from the beginning of the reference sequence. The result is 
 * nucleotide coord on the ref sequence. */
static int getSearchStartCoord(GtkWidget *blxWindow, const gboolean startBeginning, const gboolean searchLeft)
{
  int result = UNSET_INT;
  
  const BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (startBeginning)
    {
      result = (searchLeft == bc->displayRev) ? bc->refSeqRange.min : bc->refSeqRange.max;
    }
  else  
    {
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      result = detailViewGetSelectedBaseIdx(detailView);
      
      if (result != UNSET_INT)
        {
          /* Increment by one to make sure we don't re-find a previously-found match
           * (or decrement if searching leftwards) */
           if (searchLeft)
             {
                --result;
              }
            else
              {
                ++result;
              }
        }
      else
        {
          /* The start display coord is the min coord if we're searching left and the max if searching right. */
          const IntRange const *displayRange = detailViewGetDisplayRange(detailView);
          result = searchLeft ? displayRange->max : displayRange->min;
        }

      /* Convert the display coord to a nucleotide coord */
      result = convertDisplayIdxToDnaIdx(result, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  return result;
}


/* Callback called when requested to search for a DNA string. If found, sets the currently-
 * selected base index to the coord where the matching string starts. The text entry for the
 * search string is passed as the callback data. */
static void onFindDnaString(GtkWidget *button, const gint responseId, gpointer data)
{
  /* Get the search string from the text entry. If the toggle button is not active, call
   * blxWindowFindDnaString with a NULL search string to "cancel" any previous searches
   * so that "findAgain" will not attempt to perform a DNA search). */
  const char *searchStr = NULL;

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      searchStr = getStringFromTextEntry(GTK_ENTRY(data));

      if (!searchStr || strlen(searchStr) < 1)
        {
          g_critical("DNA search failed. The search string was empty.\n");
        }
    }

  /* Search left wrt the screen if the user hit 'back' search right if 'forward' */
  const gboolean searchLeft = (responseId == BLX_RESPONSE_BACK);

  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  const gboolean startBeginning = (responseId != BLX_RESPONSE_FORWARD && responseId != BLX_RESPONSE_BACK);
  int startCoord = getSearchStartCoord(blxWindow, startBeginning, searchLeft);
  
  GError *error = NULL;
  blxWindowFindDnaString(blxWindow, searchStr, startCoord, searchLeft, FALSE, &error);
  
  if (error)
    {
      if (!startBeginning)
        {
          /* Try looping round to the beginning */
          postfixError(error, " Trying again from the %s of the range.\n", (searchLeft ? "end" : "start"));
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
          
          startCoord = getSearchStartCoord(blxWindow, TRUE, searchLeft);
          blxWindowFindDnaString(blxWindow, searchStr, startCoord, searchLeft, FALSE, &error);
        }
    }
  
  if (error)
    {
      prefixError(error, "DNA search failed. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
}


/* Show the 'Find' dialog */
void showFindDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_FIND;
  GtkWidget *dialog = getPersistentDialog(bc, dialogId);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons("Find sequences", 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_GO_BACK,
                                           BLX_RESPONSE_BACK,
                                           GTK_STOCK_GO_FORWARD,
                                           BLX_RESPONSE_FORWARD,
                                           GTK_STOCK_CLOSE,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);
  
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      GtkBox *contentArea = GTK_BOX(GTK_DIALOG(dialog)->vbox);
  
      GtkRadioButton *button1 = createRadioButton(contentArea, NULL, "Sequence _name search (wildcards * and ?)", TRUE, TRUE, FALSE, onFindSeqsFromName, blxWindow);
      createRadioButton(contentArea, button1, "_DNA search", FALSE, TRUE, FALSE, onFindDnaString, blxWindow);
      createRadioButton(contentArea, button1, "Sequence name _list search", FALSE, TRUE, TRUE, onFindSeqsFromList, blxWindow);
      
      gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(blxWindow));
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/* Show the 'Info' dialog, which displays info about the currently-selected sequence(s) */
void showInfoDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Sequence info", 
						  NULL, 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CLOSE,
						  GTK_RESPONSE_REJECT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_REJECT);

  int width = blxWindow->allocation.width * 0.7;
  int height = blxWindow->allocation.height * 0.9;
  
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Compile the message text from the selected sequence(s) */
  GString *resultStr = g_string_new("");
  const gboolean dataLoaded = blxContextGetFlag(bc, BLXFLAG_EMBL_DATA_LOADED);
  GList *seqItem = bc->selectedSeqs;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      char *seqText = blxSequenceGetInfo(blxSeq, TRUE, dataLoaded);
      g_string_append_printf(resultStr, "%s\n\n", seqText);
      g_free(seqText);
    }
  
  GtkWidget *child = createScrollableTextView(resultStr->str, TRUE, blxWindow->style->font_desc, TRUE, &height, NULL);
                             
  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  gtk_widget_show_all(dialog);
  
  g_string_free(resultStr, TRUE);
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
static void destroySequenceGroup(BlxViewContext *bc, SequenceGroup **seqGroup)
{
  if (seqGroup && *seqGroup)
    {
      /* Remove it from the list of groups */
      bc->sequenceGroups = g_list_remove(bc->sequenceGroups, *seqGroup);
      
      /* If this is pointed to by the match-set pointer, null it */
      if (*seqGroup == bc->matchSetGroup)
	{
	  bc->matchSetGroup = NULL;
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
      destroySequenceGroup(blxContext, &group);
      blxWindowGroupsChanged(blxWindow);
    }
}


/* Delete all groups */
static void blxContextDeleteAllSequenceGroups(BlxViewContext *bc)
{
  GList *groupItem = bc->sequenceGroups;
  
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      destroySequenceGroup(bc, &group);
    }
  
  g_list_free(bc->sequenceGroups);
  bc->sequenceGroups = NULL;
  
}


static void blxWindowDeleteAllSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  blxContextDeleteAllSequenceGroups(bc);
  blxWindowGroupsChanged(blxWindow);
}


/* Update function to be called whenever groups have been added or deleted,
 * or sequences have been added to or removed from a group */
static void blxWindowGroupsChanged(GtkWidget *blxWindow)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  /* Re-sort all trees, because grouping affects sort order */
  callFuncOnAllDetailViewTrees(detailView, resortTree, NULL);
  
  /* Refilter the trees (because groups affect whether sequences are visible) */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);

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
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Create the new group */
  SequenceGroup *group = g_malloc(sizeof(SequenceGroup));
  
  group->seqList = seqList;
  group->ownsSeqNames = ownSeqNames;
  group->hidden = FALSE;
  
  /* Find a unique ID */
  GList *lastItem = g_list_last(bc->sequenceGroups);
  
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

  /* Set the default highlight color. */
  group->highlighted = TRUE;

  BlxColorId colorId = (groupName && !strcmp(groupName, MATCH_SET_GROUP_NAME)) ? BLXCOLOR_MATCH_SET : BLXCOLOR_GROUP;
  GdkColor *color = getGdkColor(colorId, bc->defaultColors, FALSE, bc->usePrintColors);
  group->highlightColor = *color;

  /* Add it to the list, and update */
  bc->sequenceGroups = g_list_append(bc->sequenceGroups, group);
  blxWindowGroupsChanged(blxWindow);
  
  return group;
}


/* This function sets the sequence-group-name text based on the given text entry widget */
static void onGroupNameChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newName = gtk_entry_get_text(entry);
  
  if (!newName || strlen(newName) < 1)
    {
      g_critical("Invalid group name '%s' entered; reverting to previous group name '%s'.", newName, group->groupName);
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
static void onGroupOrderChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newOrder = gtk_entry_get_text(entry);
  
  if (!newOrder || strlen(newOrder) < 1)
    {
      g_critical("Invalid order number '%s' entered; reverting to previous order number '%d'.", newOrder, group->order);
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
static void onGroupHiddenToggled(GtkWidget *button, const gint responseId, gpointer data)
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
static void onGroupHighlightedToggled(GtkWidget *button, const gint responseId, gpointer data)
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


/* Called when the user has changed the color of a group in the 'edit groups' dialog */
static void onGroupColorChanged(GtkWidget *button, const gint responseId, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), &group->highlightColor);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), &group->highlightColor, TRUE, TRUE);
  
  /* Redraw everything in the new colors */
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

      /* Show the group's highlight color in a button that will also launch a color-picker */
      GtkWidget *colorButton = gtk_color_button_new_with_color(&group->highlightColor);
      widgetSetCallbackData(colorButton, onGroupColorChanged, group);
      
      /* Create a button that will delete this group */
      GtkWidget *deleteButton = gtk_button_new_from_stock(GTK_STOCK_DELETE);
      g_signal_connect(G_OBJECT(deleteButton), "clicked", G_CALLBACK(onButtonClickedDeleteGroup), group);
      
      /* Put everything in the table */
      gtk_table_attach(table, nameWidget,		1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHiddenWidget,	2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHighlightedWidget,	3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, orderWidget,		4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, colorButton,		5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, deleteButton,		6, 7, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
    }
}


/* Called for each entry in a hash table. Compares the key of the table, which is
 * a sequence name, to the search string given in the user data. If it matches, it
 * appends the sequence name to the result list in the user data. */
static void getSequencesThatMatch(gpointer listDataItem, gpointer data)
{
  SeqSearchData *searchData = (SeqSearchData*)data;
  BlxSequence *seq = (BlxSequence*)listDataItem;
  
  /* Try matching against the full sequence name, e.g. Em:AB123456.1 */
  gboolean found = wildcardSearch(blxSequenceGetFullName(seq), searchData->searchStr);

  if (!found)
    {
      /* Try without the prefix, e.g. AB123456.1 */
      found = wildcardSearch(blxSequenceGetVariantName(seq), searchData->searchStr);
    }
  
  if (!found)
    {
      /* Try without the variant. e.g. Em:AB123456 */
      char *seqName = g_strdup(blxSequenceGetFullName(seq));
      char *cutPoint = strchr(seqName, '.');
      
      if (cutPoint)
        {
          *cutPoint = '\0';
          found = wildcardSearch(seqName, searchData->searchStr);
        }
      
      g_free(seqName);
    }
  
  if (!found)
    {
      /* Try without the prefix or variant e.g. AB123456 */
      found = wildcardSearch(blxSequenceGetShortName(seq), searchData->searchStr);
    }

  if (!found)
    {
        const int len = strlen(searchData->searchStr);
      
      if (len > 2 &&
          (searchData->searchStr[len - 1] == 'x' || searchData->searchStr[len - 1] == 'X' ||
           searchData->searchStr[len - 1] == 'i' || searchData->searchStr[len - 1] == 'I'))
        {
    
          /* BlxSequence names don't have the 'x' or 'i' postfix for exons or introns. If the search
           * string has an 'x' or 'i' on the end, ignore it. */
          char *searchStr = g_strdup(searchData->searchStr);
          searchStr[len - 1] = '\0';
      
          found = wildcardSearch(blxSequenceGetFullName(seq), searchStr);
      
          g_free(searchStr);
        }
    }  
    
  if (found)
    {
      searchData->matchList = g_list_prepend(searchData->matchList, seq);
    }  
   
}


/* If the given radio button is enabled, add a group based on the curently-
 * selected sequences. */
static void addGroupFromSelection(GtkWidget *button, const gint responseId, gpointer data)
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
	  g_critical("Warning: cannot create group; no sequences are currently selected");
	}
    }
}


/* If the given radio button is enabled, add a group based on the search text
 * in the given text entry. */
static void addGroupFromName(GtkWidget *button, const gint responseId, gpointer data)
{
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return;
    }
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

  const char *inputText = getStringFromTextEntry(GTK_ENTRY(data));
  
  GError *error = NULL;
  GList *seqList = findSeqsFromName(blxWindow, inputText, FALSE, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      createSequenceGroup(blxWindow, seqList, FALSE, NULL);
    }
}


/* Utility function to take a list of names, and return a list of those that
 * that are valid sequence names from the given list of BlxSequences. */
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


/* Utility to create a GList of BlxSequences from a textual list of sequences.
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
      g_critical("No valid sequence names in buffer '%s'\n", inputText);
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
	  blxContext->matchSetGroup = createSequenceGroup(blxWindow, seqList, FALSE, MATCH_SET_GROUP_NAME);
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
      
      /* Refresh the groups dialog, if it happens to be open */
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
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
       * keep any changes the user made (e.g. to the group color etc.) for next time. */
      g_list_free(blxContext->matchSetGroup->seqList);
      blxContext->matchSetGroup->seqList = NULL;
      blxWindowGroupsChanged(blxWindow);

      /* Refresh the groups dialog, if it happens to be open */
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
  else
    {
      requestPrimaryClipboardText(createMatchSetFromClipboard, blxWindow);
    }
}


/* If the given radio button is enabled, add a group based on the list of sequences
 * in the given text entry. */
static void addGroupFromList(GtkWidget *button, const gint responseId, gpointer data)
{
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return;
    }
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

  /* The text entry box was passed as the user data. We should have a (multi-line) text view */
  const char *inputText = getStringFromTextView(GTK_TEXT_VIEW(data));

  GError *error = NULL;
  GList *seqList = findSeqsFromList(blxWindow, inputText, FALSE, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
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
      gtk_widget_hide_all(GTK_WIDGET(dialogWindow));
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
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
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
static gboolean onRadioButtonTextEntered(GtkWidget *textWidget, GdkEventButton *event, gpointer data)
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
					 BlxResponseCallback callbackFunc,
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
      g_signal_connect(G_OBJECT(entry), "focus-in-event", G_CALLBACK(onRadioButtonTextEntered), button);
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


/* Utility to find and return a notebook child of the given widget. Assumes there is only
 * one - if there are more it will just return the first found. Returns NULL if not found. */
static GtkNotebook* containerGetChildNotebook(GtkContainer *container)
{
  GtkNotebook *result = NULL;
  
  GList *child = gtk_container_get_children(container);
  
  for ( ; child; child = child->next)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);
      
      if (GTK_IS_NOTEBOOK(childWidget))
        {
          result = GTK_NOTEBOOK(childWidget);
          break;
        }
      else if (GTK_IS_CONTAINER(childWidget))
       {
         /* recurse */
         containerGetChildNotebook(GTK_CONTAINER(childWidget));
       }
    }
  
  return result;
}


/* Callback called when user responds to groups dialog */
void onResponseGroupsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = FALSE;
  gboolean refresh = FALSE;
  
  /* If a notebook was passed, only call callbacks for widgets in the active tab */
  GtkNotebook *notebook = containerGetChildNotebook(GTK_CONTAINER(dialog->vbox));
               
  if (!notebook)
    {
      g_warning("Expected Groups dialog to contain a notebook widget. Dialog may not refresh properly.\n");
    }

  guint pageNo = notebook ? gtk_notebook_get_current_page(notebook) : 0;
  GtkWidget *page = notebook ? gtk_notebook_get_nth_page(notebook, pageNo) : NULL;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      widgetCallAllCallbacks(page, GINT_TO_POINTER(responseId));
      destroy = TRUE;
      refresh = FALSE;
      break;

    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(page, GINT_TO_POINTER(responseId));
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
      /* Groups dialog is persistent, so hide it rather than destroying it */
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
  else if (refresh)
    {
      GtkWidget *blxWindow = GTK_WIDGET(data);
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
}


/* Shows the "Group sequences" dialog. This dialog allows the user to group sequences together.
 * This tabbed dialog shows both the 'create group' and 'edit groups' dialogs in one. If the
 * 'editGroups' argument is true and groups exist, the 'Edit Groups' tab is displayed by default;
 * otherwise the 'Create Groups' tab is shown. */
void showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_GROUPS;
  GtkWidget *dialog = getPersistentDialog(bc, dialogId);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons("Groups", 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_APPLY,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      /* Make sure we only connect the response event once */
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseGroupsDialog), blxWindow);
    }
  else
    {
      /* Refresh by deleting the dialog contents and re-creating them. */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  const gboolean seqsSelected = g_list_length(bc->selectedSeqs) > 0;

  const int numRows = g_list_length(bc->sequenceGroups) + 2; /* +2 for header row and delete-all button */
  

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
  g_signal_connect(notebook, "switch-page", G_CALLBACK(onSwitchPageGroupsDialog), dialog);
  
  gtk_widget_show_all(dialog);

  if (editGroups && notebook && blxWindowGroupsExist(blxWindow))
    {
      gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), 1); /* 'edit' page is the 2nd page */
    }
  else
    {
      gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), 0); /* 'create' page is the 1st page */
    }
  
  if (bringToFront)
    {
      /* If user has asked to edit groups (and some groups exist), make the second tab
       * the default and the 'close' button the default action. (Must do this after showing
       * the child widgets due to a GTK legacy whereby the notebook won't change tabs otherwise.) */
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/***********************************************************
 *			   Settings menu                   *
 ***********************************************************/

/* This function should be called on the child widget of a dialog box that is a transient
 * child of the main blixem window. It finds the parent dialog of the child and then finds
 * the blxWindow from the dialog. */
static GtkWidget* dialogChildGetBlxWindow(GtkWidget *child)
{
  GtkWidget *result = NULL;
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(child));
  
  if (dialogWindow)
    {
      result = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
    }
  
  return result;
}


/* Set the given flag to the given value */
void blxContextSetFlag(BlxViewContext *bc, const BlxFlag flag, const gboolean newValue)
{
  gboolean *value = &g_array_index(bc->blxFlags, gboolean, flag);
  *value = newValue;
}


/* Get the value of the given flag */
gboolean blxContextGetFlag(const BlxViewContext *bc, const BlxFlag flag)
{
  gboolean result = g_array_index(bc->blxFlags, gboolean, flag);
  
  /* Special case for the show-unaligned-sequnece option: only allow this to be true if
   * squash matches is off, because they don't work well together */
  if (result && flag == BLXFLAG_SHOW_UNALIGNED_SEQ && blxContextGetFlag(bc, BLXFLAG_SQUASH_MATCHES))
    {
      result = FALSE;
    }
  
  return result;
}


/* Updates the given flag from the given button. The passed in widget is the toggle button and
 * the data is an enum indicating which flag was toggled. Returns the new value that was set.
 * Returns the new value that was set. */
static gboolean setFlagFromButton(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  BlxFlag flag = GPOINTER_TO_INT(data);
  const gboolean newValue = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  blxContextSetFlag(bc, flag, newValue);
  
  return newValue;
}


/* This callback is called when one of the boolean flags is toggled on the settings dialog.
 * This generic function sets the flag and redraws everything; if different 
 * updates are required a custom callback function can be used instead. */
static void onToggleFlag(GtkWidget *button, gpointer data)
{
  setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  blxWindowRedrawAll(blxWindow);
}


/* Callback function called when the 'squash matches' button is toggled */
static void onSquashMatches(GtkWidget *button, gpointer data)
{
  const gboolean squash = setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  detailViewUpdateSquashMatches(detailView, squash);
}


/* Callback function called when the 'invert sort order' button is toggled */
static void onSortOrderToggled(GtkWidget *button, gpointer data)
{
  const gboolean invert = setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  detailViewUpdateSortInverted(detailView, invert);
}


/* Callback function called when the 'Show SNP track' button is toggled */
static void onShowSnpTrackToggled(GtkWidget *button, gpointer data)
{
  const gboolean showSnpTrack = setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  detailViewUpdateShowSnpTrack(detailView, showSnpTrack);
}


/* Utility to create a check button with certain given properties, and to pack it into the parent */
static void createCheckButton(GtkBox *box, 
			      const char *mnemonic, 
			      const gboolean isActive, 
			      GCallback callback, 
			      gpointer data)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_box_pack_start(box, button, FALSE, FALSE, 0);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), isActive);
  
  g_signal_connect(G_OBJECT(button), "toggled", callback, data);
}


/* Callback to be called when the user has entered a new column size */
static void onColumnSizeChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)data;
  
  const gchar *newSizeText = gtk_entry_get_text(entry);
  const int newWidth = convertStringToInt(newSizeText);
  
  if (newWidth != columnInfo->width)
    {
      columnInfo->width = newWidth;

      GtkWidget *blxWindow = dialogChildGetBlxWindow(widget);
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
      callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns, NULL);
      resizeDetailViewHeaders(detailView);
    }
}


/* Just calls gtk_widget_set_sensitive, but so that we can use it as a callback */
static void widgetSetSensitive(GtkWidget *widget, gpointer data)
{
  gboolean sensitive = GPOINTER_TO_INT(data);
  gtk_widget_set_sensitive(widget, sensitive);
}


/* Callback when the user hits the 'load embl data' button on the settings dialog */
static void onButtonClickedLoadEmblData(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Load the optional data. (Note that we don't need to re-fetch sequence data.) */
  const gboolean getSequenceData = FALSE;
  const gboolean getOptionalData = TRUE;
  GList *seqsToFetch = getSeqsToPopulate(bc->matchSeqs, getSequenceData, getOptionalData);

  GError *error = NULL;
  gboolean success = fetchSequences(seqsToFetch, bc->matchSeqs, bc->fetchMode, bc->seqType, bc->net_id, bc->port, TRUE, FALSE, bc->external, &error);
  
  g_list_free(seqsToFetch);
  
  if (error)
    {
      prefixError(error, "Error loading optional data. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  if (success)
    {
      /* Set the flag to say that the data has now been loaded */
      blxContextSetFlag(bc, BLXFLAG_EMBL_DATA_LOADED, TRUE);
      
      /* Disable the button so user can't try to load data again. */
      gtk_widget_set_sensitive(button, FALSE);
      
      /* Enable the text entry boxes for all columns. They are all in the container passed
       * as the user data. */
      GtkContainer *container = GTK_CONTAINER(data);
      gtk_container_foreach(container, widgetSetSensitive, GINT_TO_POINTER(TRUE));
      
      /* Update the flag in all columns to indicate that the data is now loaded */
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      GList *listItem = detailViewGetColumnList(detailView);
      
      for ( ; listItem; listItem = listItem->next)
        {
          DetailViewColumnInfo *columnInfo = (DetailViewColumnInfo*)(listItem->data);
          columnInfo->dataLoaded = TRUE;
        }   
      
      /* Re-sort the trees, because the new data may affect the sort order */
      callFuncOnAllDetailViewTrees(detailView, resortTree, NULL);
      detailViewRedrawAll(detailView);
    }
}


/* Create a button to allow the user to load the data for optional columns, if not already loaded */
static GtkWidget* createColumnLoadDataButton(GtkBox *box, GtkWidget *detailView)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 12);

  BlxViewContext *bc = blxWindowGetContext(detailViewGetBlxWindow(detailView));
  const gboolean dataLoaded = blxContextGetFlag(bc, BLXFLAG_EMBL_DATA_LOADED);
  
  GtkWidget *button = gtk_button_new_with_label(LOAD_DATA_TEXT);
  gtk_widget_set_sensitive(button, !dataLoaded); /* only enable if data not yet loaded */
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 12);

  GtkWidget *label = gtk_label_new("Use this button to load EMBL data for the optional columns (those greyed\nout below, if any). This may take a while if there are a lot of sequences.");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 12);

  return button;
}


/* Create a set of widgets that allow columns to be resized */
static void createColumnSizeButtons(GtkWidget *parent, GtkWidget *detailView)
{
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Columns");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);

  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);

  /* Create a button to allow the user to load the full EMBL data, if not already loaded */
  GtkWidget *button = createColumnLoadDataButton(GTK_BOX(vbox), detailView);
  
  /* Arrange the column-size boxes horizontally */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  
  /* Loop through each column in the detail view and create a text box showing the
   * current width. Compile a list of widgets that are disabled, so that we can enable
   * them if/when the user clicks the button to load their data. */
  GSList *disabledWidgets = NULL;
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
      
      if (columnInfo->columnId == BLXCOL_SEQUENCE)
	{
	  /* The sequence column updates dynamically, so don't allow the user to edit it */
	  char displayText[] = "<dynamic>";
	  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
	  gtk_widget_set_sensitive(entry, FALSE);
	  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2); /* fudge up width a bit in case user enters longer text */
	}
      else
	{
          if (!columnInfo->dataLoaded)
            {
              gtk_widget_set_sensitive(vbox, FALSE);
              disabledWidgets = g_slist_prepend(disabledWidgets, entry);
            }
          
	  char *displayText = convertIntToString(columnInfo->width);
	  gtk_entry_set_text(GTK_ENTRY(entry), displayText);

	  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
	  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

	  widgetSetCallbackData(entry, onColumnSizeChanged, (gpointer)columnInfo);
          
          g_free(displayText);
	}
    }
  
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(onButtonClickedLoadEmblData), hbox);
}


/* Callback to be called when the user has entered a new percent-ID per cell */
static void onIdPerCellChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  bigPictureSetIdPerCell(bigPicture, newValue);
}

 
/* Callback to be called when the user has entered a new maximum percent-ID to display */
static void onMaxPercentIdChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  bigPictureSetMaxPercentId(bigPicture, newValue);
}

/* Callback to be called when the user has entered a new minimum percent-ID to display */
static void onMinPercentIdChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  bigPictureSetMinPercentId(bigPicture, newValue);
}


///* Utility to create a text entry widget displaying the given integer value. The
// * given callback will be called when the user OK's the dialog that this widget 
// * is a child of. */
//static void createTextEntryFromInt(GtkWidget *parent, 
//				   const char *title, 
//				   const int value, 
//				   BlxResponseCallback callbackFunc, 
//				   gpointer callbackData)
//{
//  /* Pack label and text entry into a vbox */
//  GtkWidget *vbox = createVBoxWithBorder(parent, 4, FALSE, NULL);
//
//  GtkWidget *label = gtk_label_new(title);
//  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
//
//  GtkWidget *entry = gtk_entry_new();
//  gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
//
//  char *displayText = convertIntToString(value);
//  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
//
//  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
//  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
//
//  widgetSetCallbackData(entry, callbackFunc, callbackData);
//  
//  g_free(displayText);
//}


/* Utility to create a text entry widget displaying the given double value. The
 * given callback will be called when the user OK's the dialog that this widget 
 * is a child of. */
static void createTextEntryFromDouble(GtkWidget *parent, 
                                      const char *title, 
                                      const gdouble value, 
                                      BlxResponseCallback callbackFunc, 
                                      gpointer callbackData)
{
  /* Pack label and text entry into a vbox */
  GtkWidget *vbox = createVBoxWithBorder(parent, 4, FALSE, NULL);
  
  GtkWidget *label = gtk_label_new(title);
  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
  
  GtkWidget *entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
  
  char *displayText = convertDoubleToString(value);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  widgetSetCallbackData(entry, callbackFunc, callbackData);
  
  g_free(displayText);
}


/* Create a set of widgets to allow the user to edit grid properties */
static void createGridSettingsButtons(GtkWidget *parent, GtkWidget *bigPicture)
{
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Grid properties");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);

  /* Arrange the widgets horizontally */
  GtkWidget *hbox = createHBoxWithBorder(frame, 12);
  const DoubleRange const *percentIdRange = bigPictureGetPercentIdRange(bigPicture);
  
  createTextEntryFromDouble(hbox, "%ID per cell", bigPictureGetIdPerCell(bigPicture), onIdPerCellChanged, bigPicture);
  createTextEntryFromDouble(hbox, "Max %ID", percentIdRange->max, onMaxPercentIdChanged, bigPicture);
  createTextEntryFromDouble(hbox, "Min %ID", percentIdRange->min, onMinPercentIdChanged, bigPicture);
}


/* Callback called when user has changed a blixem color */
static void onChangeBlxColor(GtkWidget *button, const gint responseId, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  
  /* update the color */
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), color);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), color, TRUE, TRUE);
  
  /* Redraw */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  blxWindowRedrawAll(blxWindow);
}


/* Callback called when user has changed the blixem background color */
static void onChangeBackgroundColor(GtkWidget *button, const gint responseId, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  
  /* update the color */
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), color);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), color, TRUE, TRUE);
  
  /* Update */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  onUpdateBackgroundColor(blxWindow);
}


/* Create a button to allow user to change the color of the given setting */
static void createColorButton(GtkTable *table, GdkColor *color, BlxResponseCallback callbackFunc, gpointer callbackData,
			      const int row, const int column, const int xpad, const int ypad)
{
  GtkWidget *colorButton = gtk_color_button_new_with_color(color);
  widgetSetCallbackData(colorButton, callbackFunc, callbackData);

  gtk_table_attach(table, colorButton, column, column + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
}


/* Create buttons for the user to be able to change the blixem colour settings */
static void createColorButtons(GtkWidget *parent, GtkWidget *blxWindow, const int borderWidth)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* put all the colors in a table. Put the table in a scrolled window, because
   * there are likely to be many rows */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(parent), scrollWin);
  
  const int numCols = 5;
  const int numRows = BLXCOL_NUM_COLORS + 1; /* add one for header row */
  const int xpad = 2;
  const int ypad = 2;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Add a header row */
  int row = 1;
  gtk_table_attach(table, gtk_label_new("Normal   "), 2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("(selected)   "), 3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Print   "), 4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("(selected)   "), 5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* loop through all defined blixem colours */
  int colorId = BLXCOLOR_MIN + 1;
  
  for ( ; colorId < BLXCOL_NUM_COLORS; ++colorId)
    {
      ++row;
    
      BlxColor *blxCol = getBlxColor(bc->defaultColors, colorId);
      GtkWidget *label = gtk_label_new(blxCol->name);
      gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);

      /* Special callback for the background color */
      BlxResponseCallback callbackFunc = (colorId == BLXCOLOR_BACKGROUND) ? onChangeBackgroundColor : onChangeBlxColor;
      
      createColorButton(table, &blxCol->normal, callbackFunc, &blxCol->normal, row, 2, xpad, ypad);
      createColorButton(table, &blxCol->print, callbackFunc, &blxCol->print, row, 4, xpad, ypad);
      createColorButton(table, &blxCol->selected, callbackFunc, &blxCol->selected, row, 3, xpad, ypad);
      createColorButton(table, &blxCol->printSelected, callbackFunc, &blxCol->printSelected, row, 5, xpad, ypad);
    }
}


/* Called when the user toggles whether print colors should be used or not */
static void onTogglePrintColors(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  const gboolean usePrintColors = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  blxWindowSetUsePrintColors(blxWindow, usePrintColors);
}


/* Callback function called when the 'Show unaligned sequence' button is toggled */
static void onShowUnalignedSeqToggled(GtkWidget *button, gpointer data)
{
  /* Get the new value */
  const gboolean showUnalignedSeq = setFlagFromButton(button, GINT_TO_POINTER(BLXFLAG_SHOW_UNALIGNED_SEQ));
  
  /* Enable/disable the sub-options. Their widgets are all in the container passed as the data. */
  GtkWidget *subComponents = GTK_WIDGET(data);
  gtk_widget_set_sensitive(subComponents, showUnalignedSeq); 
  
  /* Get the detail view from the main window */
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  /* Update the value */
  detailViewUpdateShowUnalignedSeq(detailView, showUnalignedSeq);
}


/* Callback function called when the 'Limit unaligned bases' button is toggled */
static void onLimitUnalignedBasesToggled(GtkWidget *button, gpointer data)
{
  /* Get the new value */
  const gboolean limitUnalignedBases = setFlagFromButton(button, GINT_TO_POINTER(BLXFLAG_LIMIT_UNALIGNED_BASES));

  /* Enable/disable the sub-options. Their widgets are all in the container passed as the data. */
  GtkWidget *subComponents = GTK_WIDGET(data);
  gtk_widget_set_sensitive(subComponents, limitUnalignedBases); 
  
  /* Get the detail view from the main window */
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  /* Update the value */
  detailViewUpdateLimitUnalignedBases(detailView, limitUnalignedBases);
}


/* Callback called when the user has changed the number of additional bases to show when the
 * 'show unaligned bases' option is enabled. */
static void onSetNumUnalignedBases(GtkWidget *entry, const gint responseId, gpointer data)
{
  const char *numStr = gtk_entry_get_text(GTK_ENTRY(entry));
  int numBases = convertStringToInt(numStr);
  
  GtkWidget *detailView = GTK_WIDGET(data);
  detailViewSetNumUnalignedBases(detailView, numBases);
}


/* Create option buttons for enabling/disabling display of unaligned sequence */
static void createUnalignedSeqButtons(GtkWidget *parent, GtkWidget *detailView, BlxViewContext *bc)
{
  const gboolean showUnalignedSeq = blxContextGetFlag(bc, BLXFLAG_SHOW_UNALIGNED_SEQ);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(parent), vbox, FALSE, FALSE, 0);

  /* Create an hbox for the sub-components. Create it now so we can pass it to the main toggle
   * button callback, but don't pack it in the container till we've added the other. */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_widget_set_sensitive(hbox, showUnalignedSeq); 

  /* Main check button to enable/disable the option */
  createCheckButton(GTK_BOX(vbox), "Show _unaligned sequence (only works if Squash Matches is off)", showUnalignedSeq, G_CALLBACK(onShowUnalignedSeqToggled), hbox);

  /* Text entry box to specify the limit */
  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  widgetSetCallbackData(entry, onSetNumUnalignedBases, detailView);

  DetailViewProperties *properties = detailViewGetProperties(detailView);
  char *numStr = convertIntToString(properties->numUnalignedBases);
  
  gtk_entry_set_text(GTK_ENTRY(entry), numStr);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(numStr) + 2); /* fudge up width a bit in case user enters longer text */
  g_free(numStr);
  
  /* Check button to enable/disable setting the limit */
  GtkWidget *button = gtk_check_button_new_with_mnemonic("_Limit to ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), blxContextGetFlag(bc, BLXFLAG_LIMIT_UNALIGNED_BASES));
  g_signal_connect(G_OBJECT(button), "toggled", G_CALLBACK(onLimitUnalignedBasesToggled), entry);

  /* Pack it all in the hbox */
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("   "), FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(" additional bases"), FALSE, FALSE, 0);
}


/* Refresh the given dialog, if it is open */
void refreshDialog(const BlxDialogId dialogId, GtkWidget *blxWindow)
{
  /* This is a bit crude but does the job: if the dialog is visible, just call its
   * 'show' function to re-create its contents. Only need to do anything for persistent
   * dialogs. Note that we don't want to bring the dialog to the front, just refresh it in
   * case the user looks at it again. */
  if (blxWindow)
    {
      const BlxViewContext *bc = blxWindowGetContext(blxWindow);
      GtkWidget *dialog = getPersistentDialog(bc, dialogId);
      
      if (dialog && GTK_WIDGET_VISIBLE(dialog))
        {
           switch (dialogId)
             {
               case BLXDIALOG_SETTINGS:
                 showSettingsDialog(blxWindow, FALSE);
                 break;
                 
               case BLXDIALOG_HELP:
                 showHelpDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_FIND:
                 showFindDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_VIEW:
                 showViewPanesDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_DOTTER:
                 showDotterDialog(blxWindow, FALSE);
                 break;
                 
               case BLXDIALOG_GROUPS:
                 showGroupsDialog(blxWindow, TRUE, FALSE); /* show the 'edit' pane because we've got here by adding/deleting a group */
                 break;
                 
               default:
                 break;
             };
        }
    }
  else
    {
      g_warning("Could not refresh dialog [ID=%d]; parent window not found.\n", dialogId);
    }
}


/* Show/refresh the "Settings" dialog. */
void showSettingsDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_SETTINGS;
  GtkWidget *dialog = getPersistentDialog(bc, dialogId);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons("Blixem Settings", 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_APPLY,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Need to refresh the dialog contents, so clear and re-create content area */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);

  int borderWidth = 12;
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  /* Create separate pages for general settings and colors */
  GtkWidget *notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, TRUE, TRUE, 0);

  GtkWidget *settingsPage = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(settingsPage), gtk_label_new("General"));

  /* GENERAL PAGE */
  GtkWidget *mainVBox = createVBoxWithBorder(settingsPage, borderWidth, FALSE, NULL);

  /* Display options */
  GtkWidget *vbox1 = createVBoxWithBorder(mainVBox, borderWidth, TRUE, "Display options");
  
  createCheckButton(GTK_BOX(vbox1), "_Squash matches", blxContextGetFlag(bc, BLXFLAG_SQUASH_MATCHES), G_CALLBACK(onSquashMatches), GINT_TO_POINTER(BLXFLAG_SQUASH_MATCHES));
  createUnalignedSeqButtons(vbox1, detailView, bc);
  createCheckButton(GTK_BOX(vbox1), "_Invert sort order", blxContextGetFlag(bc, BLXFLAG_INVERT_SORT), G_CALLBACK(onSortOrderToggled), GINT_TO_POINTER(BLXFLAG_INVERT_SORT));
  createCheckButton(GTK_BOX(vbox1), "_Highlight differences", blxContextGetFlag(bc, BLXFLAG_HIGHLIGHT_DIFFS), G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_HIGHLIGHT_DIFFS));
  createCheckButton(GTK_BOX(vbox1), "Show SN_P track", blxContextGetFlag(bc, BLXFLAG_SHOW_SNP_TRACK), G_CALLBACK(onShowSnpTrackToggled), GINT_TO_POINTER(BLXFLAG_SHOW_SNP_TRACK));
  createCheckButton(GTK_BOX(vbox1), "Show Sp_lice Sites for selected seqs", blxContextGetFlag(bc, BLXFLAG_SHOW_SPLICE_SITES), G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_SPLICE_SITES));
  
  GtkWidget *pfetchBox = createVBoxWithBorder(mainVBox, borderWidth, TRUE, "Fetch mode");
  createPfetchDropDownBox(GTK_BOX(pfetchBox), blxWindow);
  
  createColumnSizeButtons(mainVBox, detailView);
  createGridSettingsButtons(mainVBox, bigPicture);

  /* APPEARANCE PAGE */
  GtkWidget *appearancePage = gtk_vbox_new(FALSE, borderWidth);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(appearancePage), gtk_label_new("Appearance"));

  const gboolean usePrintColours = blxWindowGetUsePrintColors(blxWindow);
  createCheckButton(GTK_BOX(appearancePage), "Use _print colours", usePrintColours, G_CALLBACK(onTogglePrintColors), blxWindow);
  createColorButtons(appearancePage, blxWindow, borderWidth);

  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/***********************************************************
 *		    Statistics menu			   *
 ***********************************************************/

static void getStats(GtkWidget *blxWindow, GString *result, MSP *MSPlist)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Compile data for sequences */
  int totalNumSeqs = 0;             /* total number of sequences */
  int numValidSeqs = 0;             /* how many sequences have sequence data filled in */

  gint32 seqDataSize = 0;           /* total memory used by the sequence data */
  gint32 seqStructSize = 0;         /* total memory used by the sequence structs */

  GList *seqItem = bc->matchSeqs;
  for ( ; seqItem; seqItem = seqItem->next)
    {
      ++totalNumSeqs;
      seqStructSize += sizeof(BlxSequence);
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      
      if (blxSeq->sequence)
        {
          ++numValidSeqs;
          seqDataSize += blxSequenceGetLength(blxSeq) * sizeof(char);
        }
    }
  

  /* Compile data for MSPs */
  int numMSPs = 0;                /* total number of MSPs */
  gint32 mspStructSize = 0;       /* total memory used by the MSP structs */

  MSP *msp = NULL;
  for (msp = MSPlist; msp ; msp = msp->next)
    {
      ++numMSPs;
      mspStructSize += sizeof(MSP);
    }

  
  /* Other data */
  int refSeqLen = strlen(blxWindowGetRefSeq(blxWindow));

  
  /* Create the text based on the results */
  g_string_printf(result, "%s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s",
                  "Length of reference sequence\t\t\t\t\t\t\t= ", refSeqLen, " characters\n\n",
                  "Total number of match sequences\t\t\t\t\t\t= ", totalNumSeqs, "\n",
                  "Number of match sequences containing sequence data\t= ", numValidSeqs, "\n",
                  "Total memory used by sequence data\t\t\t\t\t= ", seqDataSize, " bytes\n\n",
                  "Size of each sequence struct\t\t\t\t\t\t\t= ", (int)sizeof(BlxSequence), " bytes\n",
                  "Total memory used by sequence structs\t\t\t\t\t= ", seqStructSize, " bytes\n\n",
		  "Number of MSPs\t\t\t\t\t\t\t\t\t\t= ", numMSPs, "\n",
                  "Size of each MSP\t\t\t\t\t\t\t\t\t\t= ", (int)sizeof(MSP), " bytes\n",
		  "Total memory used by MSP structs\t\t\t\t\t\t= ", (int)sizeof(MSP) * numMSPs, " bytes");
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
      /* If it's a persistent dialog, just hide it, otherwise destroy it */
      const gboolean isPersistent = GPOINTER_TO_INT(data);
      
      if (isPersistent)
        {
          gtk_widget_hide_all(GTK_WIDGET(dialog));
        }
      else
        {
          gtk_widget_destroy(GTK_WIDGET(dialog));
        }
    }
}


void showHelpDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_HELP;
  GtkWidget *dialog = getPersistentDialog(bc, dialogId);
  
  if (!dialog)
    {
      /* Create the dialog */
      dialog = gtk_dialog_new_with_buttons("Help", 
                                           NULL, 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_ABOUT,
                                           GTK_RESPONSE_HELP,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

      /* Set a pretty big initial size. */
      const int width = blxWindow->allocation.width * 0.7;
      int height = blxWindow->allocation.height * 0.9;
      char *messageText = blxprintf(HELP_TEXT1, blixemVersion);

      GtkWidget *child = createScrollableTextView(messageText, TRUE, blxWindow->style->font_desc, TRUE, &height, NULL);
      
      g_free(messageText);

      gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
      gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);
  
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseHelpDialog), GINT_TO_POINTER(TRUE));
    }
      
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
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
  showHelpDialog(blxWindow, TRUE);
}


/* Called when the user selects the View menu option, or hits the Settings shortcut key */
static void onViewMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showViewPanesDialog(blxWindow, TRUE);
}


/* Called when the user selects the 'Group Sequences' menu option, or hits the relevant shortcut key */
static void onCreateGroupMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, FALSE, TRUE);
}


/* Called when the user selects the 'Groups' menu option, or hits the relevant shortcut key */
static void onEditGroupsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, TRUE, TRUE);
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
  showSettingsDialog(blxWindow, TRUE);
}


/* Called when the user selects the Dotter menu option, or hits the Dotter shortcut key */
static void onDotterMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showDotterDialog(blxWindow, TRUE);
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
  blxWindowDeselectAllSeqs(blxWindow);
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
  
  GdkColor fgColor;
  GError *error = NULL;
  getColorFromString(BLX_WHITE, &fgColor, &error);

  gdk_gc_set_foreground(gc, &fgColor);
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


/* Handlers for specific key presses */
static gboolean onKeyPressEscape(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  /* Reset the selected base index. Leave the selected frame as it is, though. */
  GtkWidget *detailView = blxWindowGetDetailView(window);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  detailViewSetSelectedBaseIdx(detailView, UNSET_INT, properties->selectedFrame, UNSET_INT, TRUE, TRUE);
  return TRUE;
}

static gboolean onKeyPressLeftRight(GtkWidget *window, const gboolean left, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  if (ctrlModifier)
    {
      goToMatch(window, left);
    }
  else if (shiftModifier)
    {
      moveSelectedBaseIdxBy1(window, left);
    }
  else
    {
      moveSelectedDisplayIdxBy1(window, left);
    }
  
  return TRUE;
}

static gboolean onKeyPressUpDown(GtkWidget *window, const gboolean up, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = moveRowSelection(window, up, ctrlModifier, shiftModifier);
  return result;
}

static gboolean onKeyPressHomeEnd(GtkWidget *window, const gboolean home, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollToExtremity(window, home, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressComma(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollDetailView(window, TRUE, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressPeriod(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollDetailView(window, FALSE, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressPlusMinus(GtkWidget *window, const gboolean zoomIn, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  zoomBlxWindow(window, zoomIn, ctrlModifier, shiftModifier);
  return TRUE;
}

static gboolean onKeyPressC(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (ctrlModifier)
    {
      copySelectionToClipboard(window);
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressV(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (ctrlModifier)
    {
      /* Paste from the default clipboard */
      requestDefaultClipboardText(findSeqsFromClipboard, window);
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressF(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  if (ctrlModifier)
    {
      showFindDialog(window, TRUE);
    }
  else
    {
      requestPrimaryClipboardText(findSeqsFromClipboard, window);
    }

  return TRUE;
}

static gboolean onKeyPressP(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (!ctrlModifier)
    {
      goToDetailViewCoord(blxWindowGetDetailView(window), BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressG(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (!ctrlModifier)
    {	    
      toggleMatchSet(window);
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressB(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  toggleBumpState(window);
  
  /* Refresh the view dialog, if it is open */
  refreshDialog(BLXDIALOG_VIEW, window);
  
  return TRUE;
}

static gboolean onKeyPressT(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  toggleStrand(blxWindowGetDetailView(window));
  return TRUE;
}

static gboolean onKeyPressI(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  showInfoDialog(window);
  return TRUE;
}

static gboolean onKeyPressNumber(GtkWidget *window, const int number, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  togglePaneVisibility(window, number, ctrlModifier, shiftModifier);
  
  /* Refresh the view dialog, if it is open */
  refreshDialog(BLXDIALOG_VIEW, window);
  
  return TRUE;
}

static gboolean onKeyPressF3(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  findAgain(window, shiftModifier);
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
      case GDK_Escape:	    result = onKeyPressEscape(window, ctrlModifier, shiftModifier);	      break;
      
      case GDK_Left:	    result = onKeyPressLeftRight(window, TRUE, ctrlModifier, shiftModifier);  break;
      case GDK_Right:	    result = onKeyPressLeftRight(window, FALSE, ctrlModifier, shiftModifier); break;
      
      case GDK_Up:	    result = onKeyPressUpDown(window, TRUE, ctrlModifier, shiftModifier);     break;
      case GDK_Down:	    result = onKeyPressUpDown(window, FALSE, ctrlModifier, shiftModifier);    break;
	
      case GDK_Home:	    result = onKeyPressHomeEnd(window, TRUE, ctrlModifier, shiftModifier);    break;
      case GDK_End:	    result = onKeyPressHomeEnd(window, FALSE, ctrlModifier, shiftModifier);   break;

      case GDK_comma:	    result = onKeyPressComma(window, ctrlModifier, shiftModifier);	      break;
      case GDK_period:	    result = onKeyPressPeriod(window, ctrlModifier, shiftModifier);	      break;

      case GDK_equal:	    /* fall through */
      case GDK_plus:	    result = onKeyPressPlusMinus(window, TRUE, ctrlModifier, shiftModifier);  break;
      
      case GDK_minus:	    /* fall through */
      case GDK_underscore:  result = onKeyPressPlusMinus(window, FALSE, ctrlModifier, shiftModifier); break;
      
      case GDK_F3:	    result = onKeyPressF3(window, ctrlModifier, shiftModifier);		      break;

      case GDK_c:	    /* fall through */
      case GDK_C:	    result = onKeyPressC(window, ctrlModifier, shiftModifier);		      break;

      case GDK_v:	    /* fall through */
      case GDK_V:	    result = onKeyPressV(window, ctrlModifier, shiftModifier);		      break;
	
      case GDK_f:	    /* fall through */
      case GDK_F:	    result = onKeyPressF(window, ctrlModifier, shiftModifier);		      break;

      case GDK_g:	    /* fall through */
      case GDK_G:	    result = onKeyPressG(window, ctrlModifier, shiftModifier);		      break;

      case GDK_p:	    /* fall through */
      case GDK_P:	    result = onKeyPressP(window, ctrlModifier, shiftModifier);		      break;
		
      case GDK_b:	    /* fall through */
      case GDK_B:	    result = onKeyPressB(window, ctrlModifier, shiftModifier);		      break;
	
      case GDK_t:	    /* fall through */
      case GDK_T:	    result = onKeyPressT(window, ctrlModifier, shiftModifier);		      break;

      case GDK_i:	    /* fall through */
      case GDK_I:	    result = onKeyPressI(window, ctrlModifier, shiftModifier);		      break;

      case GDK_1:	    /* fall through */
      case GDK_exclam:	    result = onKeyPressNumber(window, 1, ctrlModifier, shiftModifier);	      break;

      case GDK_2:	    /* fall through */
      case GDK_quotedbl:    /* fall through */
      case GDK_at:	    result = onKeyPressNumber(window, 2, ctrlModifier, shiftModifier);	      break;

      case GDK_3:	    /* fall through */
      case GDK_currency:    result = onKeyPressNumber(window, 3, ctrlModifier, shiftModifier);	      break;
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


static void destroyBlxContext(BlxViewContext **bc)
{
  if (bc && *bc)
    {
      /* Free the list of selected sequence names (not the names themselves
       * because we don't own them). */
      if ((*bc)->selectedSeqs)
	{
	  g_list_free((*bc)->selectedSeqs);
	  (*bc)->selectedSeqs = NULL;
	}
      
      blxContextDeleteAllSequenceGroups(*bc);
      
      if ((*bc)->fetchMode)
	{
	  g_free((*bc)->fetchMode);
	  (*bc)->fetchMode = NULL;
	}
      
      if ((*bc)->defaultColors)
	{
	  BlxColorId i = BLXCOLOR_MIN + 1;
	  for (; i < BLXCOL_NUM_COLORS; ++i)
	    {
	      BlxColor *blxColor = &g_array_index((*bc)->defaultColors, BlxColor, i);
	      destroyBlxColor(blxColor);
	    }

	  g_array_free((*bc)->defaultColors, TRUE);
	  (*bc)->defaultColors = NULL;
	}
    
      if ((*bc)->blxFlags)
	{
	  g_array_free((*bc)->blxFlags, FALSE);
	}
      
      destroyMspList(&((*bc)->mspList));
      destroyBlxSequenceList(&((*bc)->matchSeqs));
      blxDestroyGffTypeList(&((*bc)->supportedTypes));
      
      /* Free the context struct itself */
      g_free((*bc));
      *bc = NULL;
    }
}


static void onDestroyBlxWindow(GtkWidget *widget)
{
  BlxWindowProperties *properties = blxWindowGetProperties(widget);
  
  if (properties)
    {
      destroyBlxContext(&properties->blxContext);

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

  /* Reset any globals */
  blviewResetGlobals();
  
  gtk_main_quit();
}


/* Free all memory used by a blixem color */
static void destroyBlxColor(BlxColor *blxColor)
{
  if (blxColor)
    {
      g_free(blxColor->name);
      g_free(blxColor->desc);
    }
}


/* Create a Blixem color. The result should be destroyed with destroyBlxColor. Returns
 * NULL if there was any problem. If "selected" state colors are not passed, this function
 * calculates a slightly darker shade of the normal color to use for when it is selected.
 * Inserts the new color into the given list */
static void createBlxColor(BlxViewContext *bc,
			   BlxColorId colorId,
			   const char *name, 
			   const char *desc, 
			   const char *normalCol, 
			   const char *printCol,
			   const char *normalColSelected,
			   const char *printColSelected)
{
  BlxColor *result = &g_array_index(bc->defaultColors, BlxColor, colorId);
				    
  result->transparent = FALSE;
  GError *error = NULL;
  
  if (!normalCol)
    {
      result->transparent = TRUE;
    }
  else if (getColorFromString(normalCol, &result->normal, &error)) 
    {
      /* find a "selected" version of it, if not passed one */
      if (normalColSelected)
	{
	  if (!getColorFromString(normalColSelected, &result->selected, &error))
	    {
              prefixError(error, "Error getting 'selected' color: using normal color instead. ");
              reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
	      result->selected = result->normal;
	    }
	}
      else
	{
	  getSelectionColor(&result->normal, &result->selected); 
	}
      
      /* Parse the print color */
      if (getColorFromString(printCol, &result->print, &error))
	{
	  /* find a "selected" version of it, if not passed one */
	  if (printColSelected)
	    {
	      if (!getColorFromString(printColSelected, &result->printSelected, &error))
		{
		  prefixError(error, "Error getting print color for selected items: using normal print color instead. ");
                  reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
		  result->printSelected = result->print;
		}
	    }
	  else
	    {
	      getSelectionColor(&result->print, &result->printSelected); 
	    }
	}
      else
	{
	  /* Error parsing the print color: use the normal color but give a warning */
          prefixError(error, "Error getting print colors: using normal colors instead. ");
          reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
	  result->print = result->normal;
	  result->printSelected = result->selected;
	}
    }
  else
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      g_free(result);
      result = NULL;
    }

  if (result)
    {
      /* Set the other properties */
      result->name = g_strdup(name);
      result->desc = g_strdup(desc);
    }
}


/* This basically does what gdk_color_to_string does, but that function isn't available in
 * older versions of GDK... */
static char* convertColorToString(GdkColor *color)
{
//  char *result = gdk_color_to_string(&widget->style->bg[GTK_STATE_NORMAL]);
//  return result;

  /* Need to make sure the color is allocated (otherwise the 'pixel' field may be zero) */
  gboolean failures[1];
  gdk_colormap_alloc_colors(gdk_colormap_get_system(), color, 1, TRUE, TRUE, failures);

  const int hexLen = 8; /* to fit a string of the format '#ffffff', plus the terminating character */
  
  char *result = g_malloc(sizeof(char) * hexLen);
  sprintf(result, "#%x", color->pixel);
  
  return result;
}


/* Create the colors that blixem will use for various specific purposes */
static void createBlxColors(BlxViewContext *bc, GtkWidget *widget)
{
  /* Initialise the array with empty BlxColor structs */
  bc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), BLXCOL_NUM_COLORS);
  int i = BLXCOLOR_MIN + 1;
  
  for ( ; i < BLXCOL_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(bc->defaultColors, *blxColor);
    }
  
  /* Get the default background color of our widgets (i.e. that inherited from the theme).
   * Convert it to a string so we can use the same creation function as the other colors */
  char *defaultBgColorStr = convertColorToString(&widget->style->bg[GTK_STATE_NORMAL]);
  createBlxColor(bc, BLXCOLOR_BACKGROUND, "Background", "Background color", defaultBgColorStr, BLX_WHITE, "#bdbdbd", NULL);
  
  /* reference sequence */
  createBlxColor(bc, BLXCOLOR_REF_SEQ, "Reference sequence", "Default background color for the reference sequence", BLX_YELLOW, BLX_LIGHT_GREY, "#d0d000", NULL);
  
  /* matches */
  createBlxColor(bc, BLXCOLOR_MATCH, "Exact match", "Exact match", "#00ffe5", BLX_LIGHT_GREY, "#00c3b0", NULL);
  createBlxColor(bc, BLXCOLOR_CONS, "Conserved match", "Conserved match", "#78b4f0", BLX_LIGHT_GREY, "#5c98d5", NULL);
  createBlxColor(bc, BLXCOLOR_MISMATCH, "Mismatch", "Mismatch", "#cacaca", BLX_WHITE, "#989898", NULL);
  createBlxColor(bc, BLXCOLOR_INSERTION, "Insertion", "Insertion", BLX_YELLOW, BLX_DARK_GREY, NULL, NULL);
  
  /* exons */
  createBlxColor(bc, BLXCOLOR_EXON_START, "Exon start", "Exon start boundary", BLX_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_EXON_END, "Exon end", "Exon end boundary", BLX_DARK_BLUE, BLX_GREY, NULL, NULL);

  createBlxColor(bc, BLXCOLOR_EXON_FILL, "Exon fill color", "Exon fill color in big picture", BLX_YELLOW, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_EXON_LINE, "Exon line color", "Exon line color in big picture", BLX_BLUE, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_CDS_FILL, "CDS fill color", "Coding section fill color in big picture", BLX_PALE_GREEN, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_CDS_LINE, "CDS line color", "Coding section line color in big picture", BLX_DARK_GREEN, BLX_GREY, BLX_VERY_DARK_GREEN, NULL);
  createBlxColor(bc, BLXCOLOR_UTR_FILL, "Exon fill color (UTR)", "Untranslated region fill color in big picture", BLX_LIGHT_RED, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_UTR_LINE, "Exon line color (UTR)", "Untranslated region line color in big picture", BLX_DARK_RED, BLX_GREY, BLX_VERY_DARK_RED, NULL);
  
  /* codons */
  createBlxColor(bc, BLXCOLOR_CODON, "Codon nucleotides", "Codon nucleotides", BLX_LIGHT_SKY_BLUE, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_MET, "MET codons", "MET codons", BLX_LAWN_GREEN, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_STOP, "STOP codons", "MET codons", BLX_LIGHT_SALMON, BLX_LIGHT_GREY, NULL, NULL);
  
  /* SNPs */
  createBlxColor(bc, BLXCOLOR_SNP, "SNPs", "SNPs", BLX_ORANGE, BLX_GREY, NULL, NULL);

  /* Big Picture */
  createBlxColor(bc, BLXCOLOR_GRID_LINE, "Grid lines", "Big Picture grid lines", BLX_YELLOW, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_GRID_TEXT, "Grid text", "Big Picture grid text", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_HIGHLIGHT_BOX, "Highlight box", "Highlight box in the big picture", "#c7c7c7", BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_PREVIEW_BOX, "Preview box", "Preview box in the big picture", "#b1b1b1", BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_MSP_LINE, "Big picture match line", "Color of the lines representing matches in the Big Picture", BLX_BLACK, BLX_BLACK, BLX_CYAN, BLX_GREY);

  /* groups */
  createBlxColor(bc, BLXCOLOR_GROUP, "Default group color", "Default highlight color for a new group", BLX_ORANGE_RED, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_MATCH_SET, "Default match set color", "Default color for the match set group (applies only when it is created for the first time or after being deleted)", BLX_RED, BLX_LIGHT_GREY, NULL, NULL);
  
  /* misc */
  createBlxColor(bc, BLXCOLOR_UNALIGNED_SEQ, "Unaligned sequence", "Addition sequence in the match that is not part of the alignment", defaultBgColorStr, BLX_WHITE, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_CANONICAL, "Canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are canonical", BLX_GREEN, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_NON_CANONICAL, "Non-canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are not canonical", BLX_RED, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_POLYA_TAIL, "polyA tail", "polyA tail", BLX_RED, BLX_GREY, NULL, NULL);
  createBlxColor(bc, BLXCOLOR_TREE_GRID_LINES, "Tree grid lines", "Tree grid lines", BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY);
  createBlxColor(bc, BLXCOLOR_CLIP_MARKER, "Clipped-match indicator", "Marker to indicate a match has been clipped to the display range", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  
  g_free(defaultBgColorStr);
}


static BlxViewContext* blxWindowCreateContext(CommandLineOptions *options,
					      const IntRange const *refSeqRange,
					      const IntRange const *fullDisplayRange,
					      const char *paddingSeq,
                                              GList *seqList,
                                              GSList *supportedTypes,
					      GtkWidget *widget,
					      GtkWidget *statusBar,
                                              const char *net_id,
                                              int port,
                                              const gboolean External)
{
  BlxViewContext *blxContext = g_malloc(sizeof *blxContext);
  
  blxContext->statusBar = statusBar;
  
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
  blxContext->matchSeqs = seqList;
  blxContext->supportedTypes = supportedTypes;
  
  blxContext->displayRev = FALSE;
  blxContext->net_id = net_id;
  blxContext->port = port;
  blxContext->external = External;
  
  blxContext->selectedSeqs = NULL;
  blxContext->sequenceGroups = NULL;
  blxContext->matchSetGroup = NULL;
  
  blxContext->autoDotter = TRUE;
  blxContext->dotterStart = UNSET_INT;
  blxContext->dotterEnd = UNSET_INT;
  blxContext->dotterZoom = 0;
  
  blxContext->defaultColors = NULL;
  blxContext->usePrintColors = FALSE;
  
  createBlxColors(blxContext, widget);
  
  blxContext->blxFlags = g_array_sized_new(FALSE, TRUE, sizeof(gboolean), BLXFLAG_NUM_FLAGS);
  
  /* Initialise all the flags to false */
  int flag = BLXFLAG_MIN + 1;
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      blxContextSetFlag(blxContext, flag, FALSE);
    }
  
  /* Set any specific flags that we want initialised to TRUE */
  blxContextSetFlag(blxContext, BLXFLAG_LIMIT_UNALIGNED_BASES, TRUE);
  blxContextSetFlag(blxContext, BLXFLAG_EMBL_DATA_LOADED, options->parseFullEmblInfo);
  
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < BLXDIALOG_NUM_DIALOGS; ++dialogId)
    {
      blxContext->dialogList[dialogId] = NULL;
    }
  
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
BlxStrand blxWindowGetActiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
}

/* Return the inactive strand - reverse strand by default, forward strand if display toggled */
static BlxStrand blxWindowGetInactiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? BLXSTRAND_FORWARD : BLXSTRAND_REVERSE;
}


/* Returns the list of all sequence groups */
GList *blxWindowGetSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext->sequenceGroups;
}

/* Returns the group that the given sequence belongs to, if any (assumes the sequence
 * is only in one group; otherwise it just returns the first group it finds). */
SequenceGroup *blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind)
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


static gboolean blxWindowGetUsePrintColors(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  return bc->usePrintColors;
}


/* Recursively set the background color for the given widget and all its children */
static void setWidgetBackgroundColor(GtkWidget *widget, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  GtkWidget *blxWindow = gtk_widget_get_toplevel(widget);
  
  blxWindow->style->bg[GTK_STATE_NORMAL] = *color;
  gtk_widget_modify_bg(widget, GTK_STATE_NORMAL, color);
  gtk_widget_modify_base(widget, GTK_STATE_NORMAL, color);
  
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), setWidgetBackgroundColor, data);
    }
}


/* This should be called whenever the background color has changed */
static void onUpdateBackgroundColor(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  GdkColor *defaultBgColor = getGdkColor(BLXCOLOR_BACKGROUND, bc->defaultColors, FALSE, bc->usePrintColors);
  setWidgetBackgroundColor(blxWindow, defaultBgColor);
  
  blxWindowRedrawAll(blxWindow);
}


/* This sets the 'use print colors' flag and then updates the display */
static void blxWindowSetUsePrintColors(GtkWidget *blxWindow, const gboolean usePrintColors)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  bc->usePrintColors = usePrintColors;
  onUpdateBackgroundColor(blxWindow);
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

      const BlxSequence *seq = (const BlxSequence*)(listItem->data);
      g_string_append(result, blxSequenceGetDisplayName(seq));
    }

  return result;
}


/* This function copys the currently-selected sequences' names to the default
 * clipboard. */
void copySelectionToClipboard(GtkWidget *blxWindow)
{
  if (g_list_length(blxWindowGetSelectedSeqs(blxWindow)) < 1)
    {
      g_critical("No sequences selected.");
    }
  else
    {
      GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);
      setDefaultClipboardText(displayText->str);
      g_string_free(displayText, TRUE);
    }
}


/* Update function to be called whenever the MSP selection has changed */
void blxWindowSelectionChanged(GtkWidget *blxWindow)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);

  /* Re-filter the detail-view trees, because selecting an MSP can cause other MSPs to 
   * become visible (i.e. polyA tails). */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  
  /* Redraw */
  updateFeedbackBox(detailView);
  blxWindowRedrawAll(blxWindow);
  
  /* Copy the selected sequence names to the PRIMARY clipboard */
  GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);
  setPrimaryClipboardText(displayText->str);
  g_string_free(displayText, TRUE);
  
  /* Refresh the dotter dialog, if it happens to be open */
  refreshDialog(BLXDIALOG_DOTTER, blxWindow);
}


/* Call this function to select the given match sequence */
void blxWindowSelectSeq(GtkWidget *blxWindow, BlxSequence *seq)
{
  if (!blxWindowIsSeqSelected(blxWindow, seq))
    {
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      blxContext->selectedSeqs = g_list_append(blxContext->selectedSeqs, seq);
      blxWindowSelectionChanged(blxWindow);
    }
}


/* Utility function to set the list of selected sequences */
void blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList)
{
  blxWindowDeselectAllSeqs(blxWindow);

  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->selectedSeqs = seqList;

  blxWindowSelectionChanged(blxWindow);
}


/* Call this function to deselect the given sequence */
void blxWindowDeselectSeq(GtkWidget *blxWindow, BlxSequence *seq)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  /* See if it's in the list and, if so, get a pointer to the list element */
  GList *foundSeq = g_list_find(blxContext->selectedSeqs, seq);

  if (foundSeq)
    {
      blxContext->selectedSeqs = g_list_remove(blxContext->selectedSeqs, foundSeq->data);
      blxWindowSelectionChanged(blxWindow);
    }
}


/* Call this function to deselect all sequences */
void blxWindowDeselectAllSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  if (g_list_length(blxContext->selectedSeqs) > 0)
    {
      g_list_free(blxContext->selectedSeqs);
      blxContext->selectedSeqs = NULL;
      blxWindowSelectionChanged(blxWindow);
    }
}


/* Returns true if the given sequence is selected */
gboolean blxWindowIsSeqSelected(GtkWidget *blxWindow, const BlxSequence *seq)
{
  GList *foundItem = NULL;
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext)
    {
      foundItem = g_list_find(blxContext->selectedSeqs, seq);
    }
  
  return (foundItem != NULL);
}


/* Set the given sequence as selected or unselected, depending on the given argument */
void blxWindowSetSeqSelected(GtkWidget *blxWindow, BlxSequence *seq, const gboolean selected)
{
  if (selected)
    {
      blxWindowSelectSeq(blxWindow, seq);
    }
  else
    {
      blxWindowDeselectSeq(blxWindow, seq);
    }
}


BlxSequence* blxWindowGetLastSelectedSeq(GtkWidget *blxWindow)
{
  BlxSequence *result = NULL;
  
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (g_list_length(selectedSeqs) > 0)
    {
      /* Get the last-selected sequence */
      GList *lastItem = g_list_last(selectedSeqs);
      result = (BlxSequence*)(lastItem->data);
    }
  
  return result;
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
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
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
  const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  
  msp->id = 0.0;
  
  if (mspIsBlastMatch(msp))
    {
      /* If there is no sequence data, leave the ID as zero */
      const char *matchSeq = mspGetMatchSeq(msp);

      if (matchSeq)
        {
          /* Note that this will reverse complement the ref seq if it is the reverse 
           * strand. This means that where there is no gaps array the comparison is trivial
           * as coordinates can be ignored and the two sequences just whipped through. */
          GError *error = NULL;
          
          char *refSeqSegment = getSequenceSegment(bc,
                                                   bc->refSeq,
                                                   msp->qRange.min, 
                                                   msp->qRange.max, 
                                                   mspGetRefStrand(msp), 
                                                   BLXSEQ_DNA, /* msp q coords are always on the dna sequence */
                                                   mspGetRefFrame(msp, bc->seqType),
                                                   bc->displayRev,
                                                   !qForward,
                                                   TRUE,
                                                   TRUE,
                                                   &error);
          
          if (!refSeqSegment)
            {
	      prefixError(error, "Failed to calculate ID for sequence '%s' (match coords = %d - %d). ", mspGetSName(msp), msp->sRange.min, msp->sRange.max);
              reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
              return;
            }
          
          /* We need to find the number of characters that match out of the total number */
          int numMatchingChars = 0;
          int totalNumChars = 0;
          const int numGaps = msp->gaps ? g_slist_length(msp->gaps) : 0;
          
          if (numGaps == 0)
            {
              /* Ungapped alignments. */
              totalNumChars = (msp->qRange.max - msp->qRange.min + 1) / bc->numFrames;
              
              if (bc->blastMode == BLXMODE_TBLASTN || bc->blastMode == BLXMODE_TBLASTX)
                {
                  int i = 0;
                  for ( ; i < totalNumChars; i++)
                    {
                      if (toupper(matchSeq[i]) == toupper(refSeqSegment[i]))
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
                      int sIndex = sForward ? msp->sRange.min + i - 1 : msp->sRange.max - i - 1;
                      if (toupper(matchSeq[sIndex]) == toupper(refSeqSegment[i]))
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
                  GSList *rangeItem = msp->gaps;
                  
                  for ( ; rangeItem; rangeItem = rangeItem->next)
                    {
                      CoordRange *range = (CoordRange*)(rangeItem->data);
                      
                      int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
                      getCoordRangeExtents(range, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
                      
                      totalNumChars += sRangeMax - sRangeMin + 1;
                      
                      /* Note that refSeqSegment is just the section of the ref seq relating to this msp.
                       * We need to translate the first coord in the range (which is in terms of the full
                       * reference sequence) into coords in the cut-down ref sequence. */
                      int q_start = qForward ? (qRangeMin - msp->qRange.min) / bc->numFrames : (msp->qRange.max - qRangeMax) / bc->numFrames;
                      
                      /* We can index sseq directly (but we need to adjust by 1 for zero-indexing). We'll loop forwards
                       * through sseq if we have the forward strand or backwards if we have the reverse strand,
                       * so start from the lower or upper end accordingly. */
                      int s_start = sForward ? sRangeMin - 1 : sRangeMax - 1 ;
                      
                      int sIdx = s_start, qIdx = q_start ;
                      while ((sForward && sIdx < sRangeMax) || (!sForward && sIdx >= sRangeMin - 1))
                        {
                          if (toupper(matchSeq[sIdx]) == toupper(refSeqSegment[qIdx]))
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
          
          msp->id = (100.0 * numMatchingChars / totalNumChars);
          
          g_free(refSeqSegment);
        }
    }
  
  return ;
}


/* Calculate the reference sequence reading frame that the given MSP belongs to, if not
 * already set. Requires either the phase to be set, or the frame to already be set; otherwise
 * assumes a phase of 0 and gives a warning. */
static void calcReadingFrame(MSP *msp, const BlxViewContext *bc)
{
  /* Always calculate frame if the phase is known, because the old code that used to pass the
   * reading frame seemed to occasionally pass an incorrect reading frame. */
  if (msp->phase != UNSET_INT || msp->qFrame == UNSET_INT)
    {
      if (msp->phase == UNSET_INT)
        {
          g_warning("Phase is not specified for MSP '%s' (q=%d-%d; s=%d-%d) - assuming phase 0.\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
          msp->phase = 0;
        }
      
      /* Get the first coord of the first complete codon. This is the start coord (5' end) of the match
       * plus the phase (if non-zero), which initially gets stored in the qFrame field in the MSP... */
      const int coord = mspGetQStart(msp) + msp->phase;

      /* Find the reading frame that this coord belongs in. This is the same as the base number within
       * reading frame 1. */
      int frame = UNSET_INT;
      const gboolean invertCoords = (mspGetRefStrand(msp) == BLXSTRAND_REVERSE);

      convertDnaIdxToDisplayIdx(coord, bc->seqType, 1, bc->numFrames, invertCoords, &bc->refSeqRange, &frame);
      
      if (frame != UNSET_INT)
        {
          if (msp->qFrame != UNSET_INT && msp->qFrame != frame && bc->seqType == BLXSEQ_PEPTIDE)
            {
              g_warning("Calculated reading frame '%d' differs from parsed reading frame '%d'; using calculated frame (%s).\n", frame, msp->qFrame, mspGetSName(msp));
            }
            
          msp->qFrame = frame;
        }
    }

  if (msp->qFrame == UNSET_INT)
    {
      g_warning("Reading frame is not set for MSP '%s' (q=%d-%d; s=%d-%d) - setting to 1.\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
      msp->qFrame = 1;
    }

  /* Messy, for backwards compatibility... set the frame number in the qframe string */
  char *frameStr = convertIntToString(msp->qFrame);
  msp->qframe[2] = frameStr[0];
  g_free(frameStr);
}


/* Calculate the ID and q frame for the given MSP and store 
 * them in the MSP struct. Returns the calculated ID (or UNSET_INT if this msp
 * is not a blast match). */
static void calcMspData(MSP *msp, BlxViewContext *bc)
{  
  /* Calculate the ref seq reading frame this match should be shown against */
  calcReadingFrame(msp, bc);

  /* Calculate the ID */
  if (mspIsBlastMatch(msp))
    {
      calcID(msp, bc);
    }
}


static gdouble calculateMspData(MSP *mspList, BlxViewContext *bc)
{
  MSP *msp = mspList;
  gdouble lowestId = -1.0;
  
  for ( ; msp; msp = msp->next)
    {
      calcMspData(msp, bc);
      
      if (mspIsBlastMatch(msp) && (lowestId == -1.0 || msp->id < lowestId))
	{
	  lowestId = msp->id;
	}
    }
  
  return lowestId;
}


/* Create the main blixem window */
GtkWidget* createBlxWindow(CommandLineOptions *options, 
                           const char *paddingSeq, 
                           GList *seqList, 
                           GSList *supportedTypes,
                           const char *net_id, 
                           int port,
                           const gboolean External)
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
  
  printf("Reference sequence [%d - %d], display range [%d - %d]\n", 
	 refSeqRange.min, refSeqRange.max, fullDisplayRange.min, fullDisplayRange.max);
  
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

  /* Set the message handlers again and pass the window, now we know it */
  g_log_set_default_handler(defaultMessageHandler, window);
  g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL, popupMessageHandler, window);

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
  
  /* Create a status bar */
  GtkWidget *statusBar = gtk_statusbar_new();
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), TRUE);
  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");

  BlxViewContext *blxContext = blxWindowCreateContext(options, 
                                                      &refSeqRange, 
                                                      &fullDisplayRange, 
                                                      paddingSeq, 
                                                      seqList, 
                                                      supportedTypes,
                                                      window, 
						      statusBar,
                                                      net_id, 
                                                      port, 
                                                      External);
  
  const gdouble lowestId = calculateMspData(options->mspList, blxContext);
  
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
					   options->initSortColumn,
                                           options->parseFullEmblInfo);
  
  /* Create a custom scrollbar for scrolling the sequence column and put it at the bottom of the window */
  GtkWidget *scrollBar = createDetailViewScrollBar(detailAdjustment, detailView);
  gtk_box_pack_start(GTK_BOX(vbox), scrollBar, FALSE, TRUE, 0);

  /* Put the statusbar at the bottom */
  gtk_box_pack_start(GTK_BOX(vbox), statusBar, FALSE, TRUE, 0);

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
  detailViewSetSortColumn(detailView, options->initSortColumn);

  /* Set the detail view font (again, this accesses the widgets' properties). */
  updateDetailViewFontDesc(detailView);

  /* Calculate the number of vertical cells in the grids (again, requires properties) */
  calculateNumVCells(bigPicture);


  /* Realise the widgets */
  printf("Starting Blixem\n");
  gtk_widget_show_all(window);

  /* If the options don't say to show the reverse strand grid, hide it now. (This must be done
   * after showing the widgets, or it will get shown again in show_all.) */
  if (!options->bigPictRev)
    {
      widgetSetHidden(revStrandGrid, TRUE);
    }

  /* Set the initial column widths. (This must be called after the widgets are 
   * realised because it causes the scroll range to be updated, which in turn causes
   * the big picture range to be set. The widgets must be realised before this because
   * the initial big picture range depends on the detail view range, which is calculated
   * from its window's width, and this will be incorrect if it has not been realised.) */
  callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns, NULL);
  resizeDetailViewHeaders(detailView);

  return window;
}
