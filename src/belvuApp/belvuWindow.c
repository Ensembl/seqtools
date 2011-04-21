/*  File: belvuWindow.c
 *  Author: Gemma Barson, 2011-04-11
 *  Copyright (c) 2009 - 2010 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: see belvuWindow.h
 *----------------------------------------------------------------------------
 */

#include <belvuApp/belvuWindow.h>
#include <belvuApp/belvuAlignment.h>
#include <belvuApp/belvu_.h>
#include <gtk/gtk.h>
#include <string.h>


#define DEFAULT_WINDOW_BORDER_WIDTH      1    /* used to change the default border width around the blixem window */
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2   /* used to start with a smaller font than the default widget font */
#define MAIN_BELVU_WINDOW_NAME           "BelvuWindow"
#define WRAPPED_BELVU_WINDOW_NAME        "WrappedBelvuWindow"
#define DEFAULT_BELVU_WINDOW_WIDTH_FRACTION     0.95   /* default width of belvu window (as fraction of screen width) */
#define DEFAULT_BELVU_WINDOW_HEIGHT_FRACTION    0.45   /* default height of belvu window (as fraction of screen height) */
#define DEFAULT_WRAP_WINDOW_WIDTH_FRACTION      0.6    /* default width of wrap window (as fraction of screen width) */
#define DEFAULT_WRAP_WINDOW_HEIGHT_FRACTION     0.85   /* default height of wrap window (as fraction of screen height) */

/* Properties specific to the belvu window */
typedef struct _BelvuWindowProperties
  {
    BelvuContext *bc;                   /* The belvu context */
    GtkActionGroup *actionGroup;        /* Holds the menu and toolbar actions */
    GtkWidget *statusBar;		/* Message bar at the bottom of the main window */
    GtkWidget *belvuAlignment;		/* The widget that draws the alignments */
  
    GdkCursor *defaultCursor;		/* default cursor */
    GdkCursor *removeSeqsCursor;	/* cursor to use when removing sequences */
  } BelvuWindowProperties;




/* Local function declarations */
static void                      onCloseMenu(GtkAction *action, gpointer data);
static void                      onQuitMenu(GtkAction *action, gpointer data);
static void                      onHelpMenu(GtkAction *action, gpointer data);
static void                      onAboutMenu(GtkAction *action, gpointer data);
static void                      onPrintMenu(GtkAction *action, gpointer data);
static void                      onWrapMenu(GtkAction *action, gpointer data);
static void                      onMakeTreeMenu(GtkAction *action, gpointer data);
static void                      onTreeOptsMenu(GtkAction *action, gpointer data);
static void                      onConsPlotMenu(GtkAction *action, gpointer data);
static void                      onSaveMenu(GtkAction *action, gpointer data);
static void                      onSaveAsMenu(GtkAction *action, gpointer data);
static void                      onOutputMenu(GtkAction *action, gpointer data);
static void                      onCompareMenu(GtkAction *action, gpointer data);
static void                      onCleanUpMenu(GtkAction *action, gpointer data);

static void                      onrmPickedMenu(GtkAction *action, gpointer data);
static void                      onRemoveSeqsMenu(GtkAction *action, gpointer data);
static void                      onCancelRemoveSeqs(GtkAction *action, gpointer data);
static void                      onrmGappySeqsPromptMenu(GtkAction *action, gpointer data);
static void                      onrmPartialMenu(GtkAction *action, gpointer data);
static void                      onmkNonRedundantPromptMenu(GtkAction *action, gpointer data);
static void                      onrmOutliersMenu(GtkAction *action, gpointer data);
static void                      onrmScoreMenu(GtkAction *action, gpointer data);
static void                      onrmColumnPromptMenu(GtkAction *action, gpointer data);
static void                      onrmColumnLeftMenu(GtkAction *action, gpointer data);
static void                      onrmColumnRightMenu(GtkAction *action, gpointer data);
static void                      onrmColumnCutoffMenu(GtkAction *action, gpointer data);
static void                      onrmEmptyColumnsInteractMenu(GtkAction *action, gpointer data);
static void                      onrmEmptyColumnsToggleMenu(GtkAction *action, gpointer data);
static void                      onreadLabelsMenu(GtkAction *action, gpointer data);
static void                      onselectGapsMenu(GtkAction *action, gpointer data);
static void                      onhideMenu(GtkAction *action, gpointer data);
static void                      onunhideMenu(GtkAction *action, gpointer data);

static void                      onToggleSchemeType(GtkRadioAction *action, GtkRadioAction *current, gpointer data);
static void                      onToggleResidueScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data);
static void                      onToggleConsScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data);
static void                      onToggleSortOrder(GtkRadioAction *action, GtkRadioAction *current, gpointer data);

static void                      ontogglePaletteMenu(GtkAction *action, gpointer data);
static void                      oncolorByResIdMenu(GtkAction *action, gpointer data);
static void                      onsaveColorCodesMenu(GtkAction *action, gpointer data);
static void                      onloadColorCodesMenu(GtkAction *action, gpointer data);
static void                      onignoreGapsMenu(GtkAction *action, gpointer data);
static void                      onprintColorsMenu(GtkAction *action, gpointer data);
static void                      onmarkupMenu(GtkAction *action, gpointer data);
static void                      oncolorRectMenu(GtkAction *action, gpointer data);
static void                      onlowercaseMenu(GtkAction *action, gpointer data);
static void                      oneditColorCodesMenu(GtkAction *action, gpointer data);

static void                      showHelpDialog();
static void			 showAboutDialog(GtkWidget *parent);
static void                      showWrapDialog(GtkWidget *belvuWindow);
static void                      showWrapWindow(GtkWidget *belvuWindow, const int linelen, const gchar *title);
static void                      getWrappedWindowDrawingArea(GtkWidget *window, gpointer data);

static void			 startRemovingSequences(GtkWidget *belvuWindow);
static void			 endRemovingSequences(GtkWidget *belvuWindow);

static gboolean                  onButtonPressBelvu(GtkWidget *window, GdkEventButton *event, gpointer data);

static BelvuWindowProperties*    belvuWindowGetProperties(GtkWidget *widget);

/***********************************************************
 *                      Menus and Toolbar                  *
 ***********************************************************/

#define rmEmptyColumnsStr "Automatically remove all-gap columns after sequence removals"

#define colorSimStr            "Average similarity by Blosum62"
#define colorIdStr             "Percent identity only"
#define colorIdSimStr          "Percent identity + Blosum62"
#define colorRectStr           "Display colors (faster without)"
#define printColorsStr         "Use gray shades (for printing)"
#define ignoreGapsStr          "Ignore gaps in conservation calculation"
#define thresholdStr           "Only colour residues above %id threshold"


/* Define the menu actions for standard menu entries */
static const GtkActionEntry menuEntries[] = {
  { "FileMenuAction",  NULL, "_File"},
  { "EditMenuAction",  NULL, "_Edit"},
  { "ColorMenuAction", NULL, "_Color"},
  { "SortMenuAction",  NULL, "_Sort"},
  { "HelpMenuAction",  NULL, "_Help"},

  { "ByResidueMenuAction", NULL, "Choose color scheme"},
  { "ByConsMenuAction",    NULL, "Choose color scheme"},

  { "CancelRemove", NULL,	      "Cancel remove sequences", "Escape",  "Cancel 'remove sequences' mode",       G_CALLBACK(onCancelRemoveSeqs)},

  { "Close",	GTK_STOCK_CLOSE,      "_Close",               "<control>W", "Close",                            G_CALLBACK(onCloseMenu)},
  { "Quit",	GTK_STOCK_QUIT,       "_Quit",                "<control>Q", "Quit  Ctrl+Q",                     G_CALLBACK(onQuitMenu)},
  { "Help",	GTK_STOCK_HELP,       "_Help",                "<control>H", "Display help  Ctrl+H",             G_CALLBACK(onHelpMenu)},
  { "About",	GTK_STOCK_ABOUT,      "A_bout",               NULL,         "About",                            G_CALLBACK(onAboutMenu)},
  { "Print",	GTK_STOCK_PRINT,      "_Print",               "<control>P", "Print  Ctrl+P",                    G_CALLBACK(onPrintMenu)},
  { "Wrap",	NULL,                 "_Wrap for printing",   NULL,         "Wrap alignments for printing",     G_CALLBACK(onWrapMenu)},
  { "MakeTree",	NULL,                 "_Make tree",           NULL,         "Make tree from current alignment", G_CALLBACK(onMakeTreeMenu)},
  { "TreeOpts",	GTK_STOCK_PROPERTIES, "_Tree options",        NULL,         "Tree options",                     G_CALLBACK(onTreeOptsMenu)},
  { "ConsPlot",	NULL,                 "Conservation p_lot",   NULL,         "Plot conservation profile",        G_CALLBACK(onConsPlotMenu)},
  { "Save",	GTK_STOCK_SAVE,       "_Save",                "<control>S", "Save alignment",                   G_CALLBACK(onSaveMenu)},
  { "SaveAs",	GTK_STOCK_SAVE_AS,    "Save _as...",          NULL,         "Save alignment as",                G_CALLBACK(onSaveAsMenu)},
  { "Output",	NULL,                 "_Output score/coords", NULL,         "Output current alignment's score and coords",  G_CALLBACK(onOutputMenu)},
  { "Compare",	NULL,                 "Co_mpare all",         NULL,         "Compare all vs all, list identity",            G_CALLBACK(onCompareMenu)},
  { "CleanUp",	GTK_STOCK_CLOSE,      "Clean _up windows",    NULL,         "Clean up windows",                 G_CALLBACK(onCleanUpMenu)},

  {"rmPicked",               NULL,    "Remove highlighted line",  NULL, "Remove highlighted line",  G_CALLBACK(onrmPickedMenu)},
  {"rmMany",		     NULL,    "Remove many sequences",    NULL, "Remove many sequences",    G_CALLBACK(onRemoveSeqsMenu)},
  {"rmGappySeqsPrompt",      NULL,    "Remove gappy sequences",   NULL, "Remove gappy sequences",   G_CALLBACK(onrmGappySeqsPromptMenu)},
  {"rmPartial",              NULL,    "Remove partial sequences", NULL, "Remove partial sequences", G_CALLBACK(onrmPartialMenu)},
  {"mkNonRedundantPrompt",   NULL,    "Make non-redundant",       NULL, "Make non-redundant",       G_CALLBACK(onmkNonRedundantPromptMenu)},
  {"rmOutliers",             NULL,    "Remove outliers",          NULL, "Remove outliers",          G_CALLBACK(onrmOutliersMenu)},
  {"rmScore",                NULL,    "Remove sequences below given score",                    NULL, "Remove sequences below given score",                    G_CALLBACK(onrmScoreMenu)},
  {"rmColumnPrompt",         NULL,    "Remove columns",           NULL, "Remove columns",           G_CALLBACK(onrmColumnPromptMenu)},
  {"rmColumnLeft",           NULL,    "<- Remove columns to the left of cursor (inclusive)",   NULL, "<- Remove columns to the left of cursor (inclusive)",   G_CALLBACK(onrmColumnLeftMenu)},
  {"rmColumnRight",          NULL,    "Remove columns to the right of cursor (inclusive) -> ", NULL, "Remove columns to the right of cursor (inclusive) -> ", G_CALLBACK(onrmColumnRightMenu)},
  {"rmColumnCutoff",         NULL,    "Remove columns according to conservation",              NULL, "Remove columns according to conservation",              G_CALLBACK(onrmColumnCutoffMenu)},
  {"rmEmptyColumnsInteract", NULL,    "Remove gappy columns",     NULL, "Remove gappy columns",   G_CALLBACK(onrmEmptyColumnsInteractMenu)},
  {"rmEmptyColumnsToggle",   NULL,    rmEmptyColumnsStr,          NULL, rmEmptyColumnsStr,        G_CALLBACK(onrmEmptyColumnsToggleMenu)},
  {"readLabels",             NULL,    "Read labels of highlighted sequence and spread them",   NULL, "Read labels of highlighted sequence and spread them",   G_CALLBACK(onreadLabelsMenu)},
  {"selectGaps",             NULL,    "Select gap character",     NULL, "Select gap character",   G_CALLBACK(onselectGapsMenu)},
  {"hide",                   NULL,    "Hide highlighted line",    NULL, "Hide highlighted line",  G_CALLBACK(onhideMenu)},
  {"unhide",                 NULL,    "Unhide all hidden lines",  NULL,"Unhide all hidden lines", G_CALLBACK(onunhideMenu)},

  {"togglePalette",          NULL, "Toggle last palette",               "T",  "Last palette",                    G_CALLBACK(ontogglePaletteMenu)},
  {"saveColorCodes",       NULL, "Save colour scheme",                NULL, "Save current colour scheme",        G_CALLBACK(onsaveColorCodesMenu)},
  {"loadColorCodes",       NULL, "Load colour scheme",                NULL, "Read colour scheme from file",      G_CALLBACK(onloadColorCodesMenu)},
  {"editColorCodes",       NULL, "Edit colour scheme",                NULL, "Open window to edit colour scheme", G_CALLBACK(oneditColorCodesMenu)},
};

/* Define the menu actions for toggle menu entries */
static const GtkToggleActionEntry toggleMenuEntries[] = {
{"colorByResId",         NULL, thresholdStr,                        NULL, thresholdStr,                        G_CALLBACK(oncolorByResIdMenu), FALSE},
{"ignoreGaps",           NULL, ignoreGapsStr,                       NULL, ignoreGapsStr,                       G_CALLBACK(onignoreGapsMenu), FALSE},
{"printColors",          NULL, printColorsStr,                      NULL, printColorsStr,                      G_CALLBACK(onprintColorsMenu), FALSE},
{"markup",               NULL, "Exclude highlighted from calculations", NULL, "Exclude highlighted from calculations", G_CALLBACK(onmarkupMenu), FALSE},
{"colorRect",            NULL, colorRectStr,                        NULL, colorRectStr,                        G_CALLBACK(oncolorRectMenu), TRUE},
{"lowercase",            NULL, "Highlight lowercase characters",    NULL, "Highlight lowercase characters",    G_CALLBACK(onlowercaseMenu), FALSE},
};

/* Define the menu actions for radio-button menu entries */
static const GtkRadioActionEntry schemeMenuEntries[] = {
{"ColorByResidue",       NULL, "By _residue",                       NULL, "Color by residue",                  BELVU_SCHEME_TYPE_RESIDUE},
{"ColorByCons",          NULL, "By _conservation",                  NULL, "Color by conservation",             BELVU_SCHEME_TYPE_CONS}
};

static const GtkRadioActionEntry residueSchemeMenuEntries[] = {
{"colorSchemeStandard",  NULL, "Erik's",                            NULL, "Erik's",                            BELVU_SCHEME_ERIK},
{"colorSchemeGibson",    NULL, "Toby's",                            NULL, "Toby's",                            BELVU_SCHEME_GIBSON},
{"colorSchemeCys",       NULL, "Cys/Gly/Pro",                       NULL, "Cys/Gly/Pro",                       BELVU_SCHEME_CGP},
{"colorSchemeEmpty",     NULL, "Clean slate",                       NULL, "Clean slate",                       BELVU_SCHEME_NONE}
};

static const GtkRadioActionEntry consSchemeMenuEntries[] = {
{"colorSim",             NULL, colorSimStr,                         NULL, colorSimStr,                         BELVU_SCHEME_BLOSUM},
{"colorId",              NULL, colorIdStr,                          NULL, colorIdStr,                          BELVU_SCHEME_ID},
{"colorIdSim",           NULL, colorIdSimStr,                       NULL, colorIdSimStr,                       BELVU_SCHEME_ID_BLOSUM}
};

static const GtkRadioActionEntry sortMenuEntries[] = {
{"scoreSort",            NULL, "by score",                     NULL, "Sort by score",                     BELVU_SORT_SCORE},
{"alphaSort",            NULL, "alphabetically",               NULL, "Sort alphabetically",               BELVU_SORT_ALPHA},
{"organismSort",         NULL, "by swissprot organism",        NULL, "Sort by swissprot organism",        BELVU_SORT_ORGANISM},
{"treeSort",             NULL, "by tree order",                NULL, "Sort by tree order",                BELVU_SORT_TREE},
{"simSort",              NULL, "by similarity to highlighted sequence", NULL, "Sort by similarity to highlighted sequence", BELVU_SORT_SIM},
{"idSort",               NULL, "by identity to highlighted sequence",   NULL, "Sort by identity to highlighted sequence",   BELVU_SORT_ID}
};


/* Define the menu layout */
static const char standardMenuDescription[] =
"<ui>"
/* MAIN MENU BAR */
"  <accelerator action='CancelRemove'/>"
"  <menubar name='MenuBar'>"
     /* File menu */
"    <menu action='FileMenuAction'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Wrap'/>"
"      <menuitem action='Print'/>"
"      <separator/>"
"      <menuitem action='MakeTree'/>"
"      <menuitem action='TreeOpts'/>"
"      <separator/>"
"      <menuitem action='ConsPlot'/>"
"      <separator/>"
"      <menuitem action='Save'/>"
"      <menuitem action='SaveAs'/>"
"      <menuitem action='Output'/>"
"      <separator/>"
"      <menuitem action='Compare'/>"
"      <menuitem action='CleanUp'/>"
"    </menu>"
    /* Edit menu */
"    <menu action='EditMenuAction'>"
"      <menuitem action='rmPicked'/>"
"      <menuitem action='rmMany'/>"
"      <menuitem action='rmGappySeqsPrompt'/>"
"      <menuitem action='rmPartial'/>"
"      <menuitem action='mkNonRedundantPrompt'/>"
"      <menuitem action='rmOutliers'/>"
"      <menuitem action='rmScore'/>"
"      <separator/>"
"      <menuitem action='rmColumnPrompt'/>"
"      <menuitem action='rmColumnLeft'/>"
"      <menuitem action='rmColumnRight'/>"
"      <menuitem action='rmColumnCutoff'/>"
"      <menuitem action='rmEmptyColumnsInteract'/>"
"      <menuitem action='rmEmptyColumnsToggle'/>"
"      <separator/>"
"      <menuitem action='readLabels'/>"
"      <menuitem action='selectGaps'/>"
"      <menuitem action='hide'/>"
"      <menuitem action='unhide'/>"
"    </menu>"
    /* Color menu */
"    <menu action='ColorMenuAction'>"
"      <menuitem action='ColorByResidue'/>"
"      <menu action='ByResidueMenuAction'>"
"        <menuitem action='colorSchemeStandard'/>"
"        <menuitem action='colorSchemeGibson'/>"
"        <menuitem action='colorSchemeCys'/>"
"        <menuitem action='colorSchemeEmpty'/>"
"      </menu>"
"      <menuitem action='colorByResId'/>"
"      <separator/>"
"      <menuitem action='ColorByCons'/>"
"      <menu action='ByConsMenuAction'>"
"        <menuitem action='colorSim'/>"
"        <menuitem action='colorId'/>"
"        <menuitem action='colorIdSim'/>"
"      </menu>"
"      <menuitem action='ignoreGaps'/>"
"      <menuitem action='printColors'/>"
"      <separator/>"
"      <menuitem action='togglePalette'/>"
"      <menuitem action='markup'/>"
"      <menuitem action='colorRect'/>"
"      <menuitem action='lowercase'/>"
"      <separator/>"
"      <menuitem action='editColorCodes'/>"
"      <menuitem action='saveColorCodes'/>"
"      <menuitem action='loadColorCodes'/>"
"    </menu>"
    /* Sort menu */
"    <menu action='SortMenuAction'>"
"      <menuitem action='scoreSort'/>"
"      <menuitem action='alphaSort'/>"
"      <menuitem action='organismSort'/>"
"      <menuitem action='treeSort'/>"
"      <menuitem action='simSort'/>"
"      <menuitem action='idSort'/>"
"    </menu>"
/* Help menu */
"    <menu action='HelpMenuAction'>"
"      <menuitem action='Help'/>"
"      <menuitem action='About'/>"
"    </menu>"
"  </menubar>"
/* Main context menu */
"  <popup name='ContextMenu' accelerators='true'>"
"    <menuitem action='Quit'/>"
"    <menuitem action='Wrap'/>"
"    <menuitem action='Print'/>"
"    <separator/>"
"    <menuitem action='MakeTree'/>"
"    <menuitem action='TreeOpts'/>"
"    <separator/>"
"    <menuitem action='ConsPlot'/>"
"    <separator/>"
"    <menuitem action='Save'/>"
"    <menuitem action='SaveAs'/>"
"    <menuitem action='Output'/>"
"    <separator/>"
"    <menuitem action='Compare'/>"
"    <menuitem action='CleanUp'/>"
"  </popup>"
/* Wrapped-alignments window context menu */
"  <popup name='WrapContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <menuitem action='Wrap'/>"
"  </popup>"
/* Tree context menu */
"  <popup name='TreeContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <menuitem action='TreeOpts'/>"
//"    <menuitem action='SaveTree'/>"
//"    <menuitem action='FindOrthogs'/>"
//"    <menuitem action='ShowOrgs'/>"
"  </popup>"
/* Conservation-plot context menu */
"  <popup name='PlotContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
//"    <menuitem action='PlotOpts'/>"
"  </popup>"
/* Toolbar */
"  <toolbar name='Toolbar'>"
"    <toolitem action='Help'/>"
"    <toolitem action='About'/>"
"  </toolbar>"
"</ui>";


static void enableMenuAction(GtkActionGroup *action_group, const char *actionName, const gboolean enable)
{
  GtkAction *action = gtk_action_group_get_action(action_group, actionName);
  gtk_action_set_sensitive(action, enable);
}


/* Certain actions need to be enabled/disabled depending on certain properties */
static void greyOutInvalidActions(BelvuContext *bc, GtkActionGroup *action_group)
{
  enableMenuAction(action_group, "Output", bc->displayScores);
  enableMenuAction(action_group, "rmScore", bc->displayScores);
  enableMenuAction(action_group, "scoreSort", bc->displayScores);
  
  enableMenuAction(action_group, "colorByResId", !bc->color_by_conserv);

  enableMenuAction(action_group, "ignoreGaps", bc->color_by_conserv);
  enableMenuAction(action_group, "printColors", bc->color_by_conserv);
}


/* Utility function to create the UI manager for the menus */
static GtkUIManager* createUiManager(GtkWidget *window, 
                                     BelvuContext *bc, 
                                     GtkActionGroup **actionGroupOut)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  gtk_action_group_add_actions(action_group, menuEntries, G_N_ELEMENTS(menuEntries), window);
  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS(toggleMenuEntries), window);

  BelvuSchemeType schemeType = bc->color_by_conserv ? BELVU_SCHEME_TYPE_CONS : BELVU_SCHEME_TYPE_RESIDUE;
  
  gtk_action_group_add_radio_actions(action_group, schemeMenuEntries, G_N_ELEMENTS(schemeMenuEntries), schemeType, G_CALLBACK(onToggleSchemeType), window);
  gtk_action_group_add_radio_actions(action_group, residueSchemeMenuEntries, G_N_ELEMENTS(residueSchemeMenuEntries), BELVU_SCHEME_ERIK, G_CALLBACK(onToggleResidueScheme), window);
  gtk_action_group_add_radio_actions(action_group, consSchemeMenuEntries, G_N_ELEMENTS(consSchemeMenuEntries), BELVU_SCHEME_BLOSUM, G_CALLBACK(onToggleConsScheme), window);
  gtk_action_group_add_radio_actions(action_group, sortMenuEntries, G_N_ELEMENTS(sortMenuEntries), BELVU_SORT_TREE, G_CALLBACK(onToggleSortOrder), window);

  greyOutInvalidActions(bc, action_group);

  GtkUIManager *ui_manager = gtk_ui_manager_new();
  gtk_ui_manager_insert_action_group(ui_manager, action_group, 0);
  gtk_ui_manager_set_add_tearoffs(ui_manager, TRUE);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group(ui_manager);
  gtk_window_add_accel_group(GTK_WINDOW(window), accel_group);
  
  if (actionGroupOut)
    *actionGroupOut = action_group;
  
  return ui_manager;
}


/* Create a menu */
static GtkWidget* createBelvuMenu(GtkWidget *window, 
                                  const char *menuDescription, 
                                  const char *path, 
                                  GtkUIManager *ui_manager)
{
  GError *error = NULL;
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, menuDescription, -1, &error))
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  GtkWidget *menu = gtk_ui_manager_get_widget (ui_manager, path);
  
  return menu;
}


/* The following functions implement the menu actions */

/* FILE MENU ACTIONS */
static void onCloseMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  gtk_widget_destroy(window);
}

static void onQuitMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  gtk_widget_destroy(belvuWindow);
}

static void onHelpMenu(GtkAction *action, gpointer data)
{
  showHelpDialog();
}

static void onAboutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showAboutDialog(belvuWindow);
}

static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);

  static GtkPageSetup *pageSetup = NULL;
  static GtkPrintSettings *printSettings = NULL;
  
  if (!pageSetup)
    {
      pageSetup = gtk_page_setup_new();
      gtk_page_setup_set_orientation(pageSetup, GTK_PAGE_ORIENTATION_PORTRAIT);
    }
  
  if (!printSettings)
    {
      printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(printSettings, GTK_PAGE_ORIENTATION_PORTRAIT);
      gtk_print_settings_set_quality(printSettings, GTK_PRINT_QUALITY_HIGH);
      gtk_print_settings_set_resolution(printSettings, DEFAULT_PRINT_RESOLUTION);
    }
  
  /* If this is the wrapped-alignment window, just print the main drawing area */
  const char *name = gtk_widget_get_name(window);
  GtkWidget *widgetToPrint = NULL;
  
  if (stringsEqual(name, WRAPPED_BELVU_WINDOW_NAME, TRUE))
    {
      getWrappedWindowDrawingArea(window, &widgetToPrint);
    }
  
  if (widgetToPrint)
    blxPrintWidget(widgetToPrint, window, &printSettings, &pageSetup, TRUE, PRINT_FIT_WIDTH);
  else
    blxPrintWidget(window, window, &printSettings, &pageSetup, FALSE, PRINT_FIT_BOTH);
}

static void onWrapMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showWrapDialog(belvuWindow);
}

static void onMakeTreeMenu(GtkAction *action, gpointer data)
{
}

static void onTreeOptsMenu(GtkAction *action, gpointer data)
{
}

static void onConsPlotMenu(GtkAction *action, gpointer data)
{
}

static void onSaveMenu(GtkAction *action, gpointer data)
{
}

static void onSaveAsMenu(GtkAction *action, gpointer data)
{
}

static void onOutputMenu(GtkAction *action, gpointer data)
{
}

static void onCompareMenu(GtkAction *action, gpointer data)
{
}

static void onCleanUpMenu(GtkAction *action, gpointer data)
{
}

/* EDIT MENU ACTIONS */
static void onrmPickedMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  if (!bc->highlightedAln) 
    {
      g_critical("Please select a sequence to remove.\n");
      return;
    }

  const int idx = bc->highlightedAln->nr - 1;
  g_array_remove_index(bc->alignArr, idx);
  arrayOrder(bc->alignArr);
  
  bc->saved = FALSE;
  
  g_message("Removed %s/%d-%d.  %d sequences left.\n\n", bc->highlightedAln->name, bc->highlightedAln->start, bc->highlightedAln->end, bc->alignArr->len);

  bc->highlightedAln = NULL;
  
  rmFinaliseGapRemoval(bc);
  
  updateOnVScrollSizeChaged(properties->belvuAlignment);
}

static void onRemoveSeqsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  startRemovingSequences(belvuWindow);
}

static void onCancelRemoveSeqs(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  endRemovingSequences(belvuWindow);
}

static void onrmGappySeqsPromptMenu(GtkAction *action, gpointer data)
{
}

static void onrmPartialMenu(GtkAction *action, gpointer data)
{
}

static void onmkNonRedundantPromptMenu(GtkAction *action, gpointer data)
{
}

static void onrmOutliersMenu(GtkAction *action, gpointer data)
{
}

static void onrmScoreMenu(GtkAction *action, gpointer data)
{
}

static void onrmColumnPromptMenu(GtkAction *action, gpointer data)
{
}

static void onrmColumnLeftMenu(GtkAction *action, gpointer data)
{
}

static void onrmColumnRightMenu(GtkAction *action, gpointer data)
{
}

static void onrmColumnCutoffMenu(GtkAction *action, gpointer data)
{
}

static void onrmEmptyColumnsInteractMenu(GtkAction *action, gpointer data)
{
}

static void onrmEmptyColumnsToggleMenu(GtkAction *action, gpointer data)
{
}

static void onreadLabelsMenu(GtkAction *action, gpointer data)
{
}

static void onselectGapsMenu(GtkAction *action, gpointer data)
{
}

static void onhideMenu(GtkAction *action, gpointer data)
{
}

static void onunhideMenu(GtkAction *action, gpointer data)
{
}

/* COLOR MENU ACTIONS */
static void onToggleSchemeType(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  properties->bc->color_by_conserv = !properties->bc->color_by_conserv;

  greyOutInvalidActions(properties->bc, properties->actionGroup);
}

static void onToggleResidueScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
}

static void onToggleConsScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
}

static void onToggleSortOrder(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
}

static void ontogglePaletteMenu(GtkAction *action, gpointer data)
{
}

static void oncolorByResIdMenu(GtkAction *action, gpointer data)
{
}

static void onsaveColorCodesMenu(GtkAction *action, gpointer data)
{
}

static void onloadColorCodesMenu(GtkAction *action, gpointer data)
{
}

static void onignoreGapsMenu(GtkAction *action, gpointer data)
{
}

static void onprintColorsMenu(GtkAction *action, gpointer data)
{
}

static void onmarkupMenu(GtkAction *action, gpointer data)
{
}

static void oncolorRectMenu(GtkAction *action, gpointer data)
{
}

static void onlowercaseMenu(GtkAction *action, gpointer data)
{
}

static void oneditColorCodesMenu(GtkAction *action, gpointer data)
{
}

/* SORT MENU ACTIONS */


/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuWindowProperties* belvuWindowGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuWindowProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuWindowProperties")) : NULL;
}

static void onDestroyBelvuWindow(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (properties)
    {
      destroyBelvuContext(&properties->bc);
  
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuWindow), "BelvuWindowProperties", NULL);
    }

  gtk_main_quit();  
}


/* Create the properties struct and initialise all values. */
static void belvuWindowCreateProperties(GtkWidget *belvuWindow, 
					BelvuContext *bc, 
					GtkActionGroup *actionGroup,
					GtkWidget *statusBar,
					GtkWidget *belvuAlignment)
{
  if (belvuWindow)
    {
      BelvuWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->actionGroup = actionGroup;
      properties->statusBar = statusBar;
      properties->belvuAlignment = belvuAlignment;
    
      properties->defaultCursor = NULL; /* get from gdkwindow once it is shown */
      properties->removeSeqsCursor = gdk_cursor_new(GDK_PIRATE);
      
      g_object_set_data(G_OBJECT(belvuWindow), "BelvuWindowProperties", properties);
      g_signal_connect(G_OBJECT(belvuWindow), "destroy", G_CALLBACK (onDestroyBelvuWindow), NULL);
    }
}


/***********************************************************
 *		      Remove sequences			   *
 ***********************************************************/

/* This should be called when the 'removing sequences' option has changed */
static void updateSequenceRemovalMode(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (properties->bc->removingSeqs)
    {
      gdk_window_set_cursor(belvuWindow->window, properties->removeSeqsCursor);
      g_message("Double-click on sequences to remove; press Esc or right-click to cancel.\n");
    }
  else
    {
      gdk_window_set_cursor(belvuWindow->window, properties->defaultCursor);
      g_message("Finished removing sequences.\n");
    }
}


static void startRemovingSequences(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  properties->bc->removingSeqs = TRUE;
  updateSequenceRemovalMode(belvuWindow);
}

static void endRemovingSequences(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  properties->bc->removingSeqs = FALSE;
  updateSequenceRemovalMode(belvuWindow);
}


/***********************************************************
 *			About dialog			   *
 ***********************************************************/

/* Returns a string which is the name of the Blixem application. */
static char *belvuGetAppName(void)
{
  return BELVU_TITLE ;
}

/* Returns a copyright string for the Blixem application. */
static char *belvuGetCopyrightString(void)
{
  return BELVU_COPYRIGHT_STRING ;
}

/* Returns the Blixem website URL. */
static char *belvuGetWebSiteString(void)
{
  return BELVU_WEBSITE_STRING ;
}

/* Returns a comments string for the Blixem application. */
static char *belvuGetCommentsString(void)
{
  return BELVU_COMMENTS_STRING() ;
}

/* Returns a license string for the belvu application. */
static char *belvuGetLicenseString(void)
{
  return BELVU_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Blixem code. */
static char *belvuGetVersionString(void)
{
  return BELVU_VERSION_STRING ;
}


/* Shows the 'About' dialog */
static void showAboutDialog(GtkWidget *parent)
{
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (6)
  const gchar *authors[] = {AUTHOR_LIST, NULL} ;
  
  gtk_show_about_dialog(GTK_WINDOW(parent),
			"authors", authors,
			"comments", belvuGetCommentsString(), 
			"copyright", belvuGetCopyrightString(),
			"license", belvuGetLicenseString(),
			"name", belvuGetAppName(),
			"version", belvuGetVersionString(),
			"website", belvuGetWebSiteString(),
			NULL) ;
#endif
  
  return ;
}

/***********************************************************
 *                      Help dialog                        *
 ***********************************************************/

static void showHelpDialog()
{
  GError *error = NULL;

  /* The docs should live in /share/doc/seqtools/, in the same parent
   * directory that our executable's 'bin' directory is in. Open the 'quick
   * start' page. */
  char rel_path[100] = "../share/doc/seqtools/belvu_quick_start.html";

  /* Find the executable's path */
  char *exe = NULL;
  gboolean ok = findCommand(g_get_prgname(), &exe);

  if (ok)
    {
      /* Get the executable's directory */
      char *dir = g_path_get_dirname(exe);
      
      ok = dir != NULL;
      
      if (ok)
        {
          /* Get the path to the html page */
          char *path = blxprintf("%s/%s", dir, rel_path);
          
          ok = path != NULL;
          
          if (ok)
            {
              g_message("Opening help page '%s'\n", path);
              seqtoolsLaunchWebBrowser(path, &error);
              g_free(path);
            }

          g_free(dir);
        }
    }
  
  if (!ok)
    {
      if (error)
        reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      else
        g_critical("Could not find help documentation: %s\n", rel_path);
    }
}


/***********************************************************
 *                         Wrap window                     *
 ***********************************************************/

/* Create a text entry with a label.  'labelText' gives the label text and 
 * defaultInput gives the default input to show in the text entry (may be null).
 * Adds the result to table, if given, and returns the text entry widget */
static GtkWidget* createTextEntryWithLabel(const char *labelText, 
                                           const char *defaultInput, 
                                           GtkTable *table,
                                           const int col,
                                           const int row)
{
  const int xpad = 2;
  const int ypad = 2;
  
  /* Create the label in the given column */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1, 0);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create the entry in the next column */
  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

  gtk_table_attach(table, entry, col + 1, col + 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  if (defaultInput)
    {
      gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);
      const int defaultLen = min(strlen(defaultInput) * 8, 500);
      gtk_widget_set_size_request(entry, defaultLen, -1);
    }
  
  return entry;
}


/* This shows a dialog that asks the user for settings for the wrap-alignment
 * view and then opens the wrap-alignment window on ok. */
static void showWrapDialog(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Belvu - wrap alignment", 
                                                  GTK_WINDOW(belvuWindow), 
                                                  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  /* Create a text entry for the line width and title */
  GtkWidget *table = gtk_table_new(2, 2, FALSE);
  gtk_box_pack_start(GTK_BOX(contentArea), table, TRUE, TRUE, 0);
  
  GtkWidget *widthEntry = createTextEntryWithLabel("Line width", "80", GTK_TABLE(table), 0, 0);
  GtkWidget *titleEntry = createTextEntryWithLabel("Title", properties->bc->Title, GTK_TABLE(table), 0, 1);
  
  gtk_widget_show_all(dialog);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const gchar *inputText = gtk_entry_get_text(GTK_ENTRY(widthEntry));
      const int linelen = convertStringToInt(inputText);
      
      const gchar *title = gtk_entry_get_text(GTK_ENTRY(titleEntry));
      
      showWrapWindow(belvuWindow, linelen, title);
    }
  
  gtk_widget_destroy(dialog);
}


/* Utility function to return the main drawing area of the wrapped-alignment window */
static void getWrappedWindowDrawingArea(GtkWidget *widget, gpointer data)
{
  GtkWidget **result = (GtkWidget**)data;

  if (*result != NULL)
    {
      return;
    }
  else if (GTK_IS_LAYOUT(widget))
    {
      *result = widget;
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), getWrappedWindowDrawingArea, result);
    }
}


static void setWrapWindowStyleProperties(GtkWidget *window)
{
  gtk_widget_set_name(window, WRAPPED_BELVU_WINDOW_NAME);
  
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);

  const int width = screenWidth * DEFAULT_WRAP_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_WRAP_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


static void showWrapWindow(GtkWidget *belvuWindow, const int linelen, const gchar *title)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  /* Create the window */
  GtkWidget *wrapWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setWrapWindowStyleProperties(wrapWindow);
  
  /* Create the context menu and set a callback to show it */
  GtkUIManager *uiManager = createUiManager(wrapWindow, properties->bc, NULL);
  GtkWidget *contextmenu = createBelvuMenu(wrapWindow, standardMenuDescription, "/WrapContextMenu", uiManager);
  
  gtk_widget_add_events(wrapWindow, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(wrapWindow), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(wrapWindow), vbox);
  
  /* Add the alignment section */
  GtkWidget *wrappedAlignment = createBelvuAlignment(properties->bc, title, linelen);
  gtk_box_pack_start(GTK_BOX(vbox), wrappedAlignment, TRUE, TRUE, 0);
  
  gtk_widget_show_all(wrapWindow);
  gtk_window_present(GTK_WINDOW(wrapWindow));
}

/***********************************************************
 *                         Events                      *
 ***********************************************************/

/* Mouse button handler */
static gboolean onButtonPressBelvu(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      GtkMenu *contextMenu = GTK_MENU(data);
      gtk_menu_popup (contextMenu, NULL, NULL, NULL, NULL, event->button, event->time);
      handled = TRUE;
    }
  
  return handled;
}



/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Create the colors that belvu will use for various specific purposes */
static void createBelvuColors(BelvuContext *bc, GtkWidget *widget)
{
  /* Initialise the array with empty BlxColor structs */
  bc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), BELCOLOR_NUM_COLORS);
  int i = BELCOLOR_MIN + 1;
  
  for ( ; i < BELCOLOR_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(bc->defaultColors, *blxColor);
    }
  
  createBlxColor(bc->defaultColors, BELCOLOR_BACKGROUND, "Background", "Background color", BLX_WHITE, BLX_WHITE, "#bdbdbd", NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_ALIGN_TEXT, "Text color for alignments", "Text color for alignments", BLX_BLACK, BLX_BLACK, NULL, NULL);
}


/* Set various properties for the main belvu window components */
static void setStyleProperties(GtkWidget *window, GtkToolbar *toolbar)
{
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_BELVU_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_BELVU_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  gtk_container_set_border_width (GTK_CONTAINER(window), DEFAULT_WINDOW_BORDER_WIDTH); 
  gtk_window_set_mnemonic_modifier(GTK_WINDOW(window), GDK_MOD1_MASK); /* MOD1 is ALT on most systems */
  
  /* Set the default font size to be a bit smaller than usual */
  int origSize = pango_font_description_get_size(window->style->font_desc) / PANGO_SCALE;
  const char *origFamily = pango_font_description_get_family(window->style->font_desc);

  char parseString[500];
  sprintf(parseString, "gtk-font-name = \"%s %d\"", origFamily, origSize + DEFAULT_FONT_SIZE_ADJUSTMENT);
  gtk_rc_parse_string(parseString);


  /* Set toolbar style properties */
  gtk_toolbar_set_style(toolbar, GTK_TOOLBAR_ICONS);
  gtk_toolbar_set_icon_size(toolbar, GTK_ICON_SIZE_SMALL_TOOLBAR);
}


gboolean createBelvuWindow(BelvuContext *bc, BlxMessageData *msgData)
{
  gboolean ok = TRUE;
  
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_widget_set_name(window, MAIN_BELVU_WINDOW_NAME);
  
  createBelvuColors(bc, window);
  
  /* Create the status bar */
  GtkWidget *statusBar = gtk_statusbar_new();
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), TRUE);
  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), FALSE);

  /* Set the window and statusbar in the message handler data, now that we know them */
  msgData->parent = GTK_WINDOW(window);
  msgData->statusBar = GTK_STATUSBAR(statusBar);

  /* Create the menu and toolbar */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(window, bc, &actionGroup);
  GtkWidget *menubar = createBelvuMenu(window, standardMenuDescription, "/MenuBar", uiManager);
  GtkWidget *contextmenu = createBelvuMenu(window, standardMenuDescription, "/ContextMenu", uiManager);
  GtkWidget *toolbar = createBelvuMenu(window, standardMenuDescription, "/Toolbar", uiManager);

  /* Set the style properties */
  setStyleProperties(window, GTK_TOOLBAR(toolbar));

  /* Create the alignment section */
  GtkWidget *belvuAlignment = createBelvuAlignment(bc, NULL, UNSET_INT);
  
  /* We'll put everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), GTK_WIDGET(vbox));

  gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), toolbar, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), belvuAlignment, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), statusBar, FALSE, FALSE, 0);

  /* Connect signals */
  gtk_widget_add_events(window, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
//  graphRegister(PICK, boxPick);
//  graphRegister(MIDDLE_DOWN, middleDown);
//  graphRegister(RESIZE, belvuRedraw);
//  graphRegister(KEYBOARD, keyboard);
//  graphRegister(DESTROY, belvuDestroy) ;
//  

  belvuWindowCreateProperties(window, bc,  actionGroup, statusBar, belvuAlignment);
  
  gtk_widget_show_all(window);
  
  /* Set the default cursor (can only get the window's cursor after window is shown) */
  BelvuWindowProperties *properties = belvuWindowGetProperties(window);
  properties->defaultCursor = gdk_window_get_cursor(window->window);

  return ok;
}
