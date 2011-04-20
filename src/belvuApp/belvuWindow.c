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
    BelvuContext *bc;	              /* The belvu context */
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

static void                      showHelpDialog();
static void                      showWrapDialog(GtkWidget *belvuWindow);
static void                      showWrapWindow(GtkWidget *belvuWindow, const int linelen, const gchar *title);
static void                      getWrappedWindowDrawingArea(GtkWidget *window, gpointer data);

static gboolean                  onButtonPressBelvu(GtkWidget *window, GdkEventButton *event, gpointer data);

static BelvuWindowProperties*    belvuWindowGetProperties(GtkWidget *widget);

/***********************************************************
 *                      Menus and Toolbar                  *
 ***********************************************************/

/* Define the menu actions */
static const GtkActionEntry menuEntries[] = {
  { "FileMenuAction",  NULL, "_File"},
  { "EditMenuAction",  NULL, "_Edit"},
  { "ColorMenuAction", NULL, "_Color"},
  { "SortMenuAction",  NULL, "_Sort"},

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
  { "CleanUp",	GTK_STOCK_CLOSE,      "Clean _up windows",    NULL,         "Clean up windows",                 G_CALLBACK(onCleanUpMenu)}

};


/* Define the menu layout */
static const char standardMenuDescription[] =
"<ui>"
"  <menubar name='MenuBar'>"
"    <menu action='FileMenuAction'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
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
"    <menu action='EditMenuAction'>"
"    </menu>"
"    <menu action='ColorMenuAction'>"
"    </menu>"
"    <menu action='SortMenuAction'>"
"    </menu>"
"  </menubar>"
"  <popup name='ContextMenu' accelerators='true'>" /* main window context menu */
"    <menuitem action='Quit'/>"
"    <menuitem action='Help'/>"
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
"  <popup name='WrapContextMenu' accelerators='true'>" /* wrap-window context menu */
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <menuitem action='Wrap'/>"
"  </popup>"
"  <toolbar name='Toolbar'>"
"    <toolitem action='Help'/>"
"    <toolitem action='About'/>"
"  </toolbar>"
"</ui>";



/* Utility function to create the UI manager for the menus */
static GtkUIManager* createUiManager(GtkWidget *window)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  gtk_action_group_add_actions(action_group, menuEntries, G_N_ELEMENTS(menuEntries), window);
//  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS (toggleMenuEntries), window);
//  gtk_action_group_add_radio_actions(action_group, radioMenuEntries, G_N_ELEMENTS (radioMenuEntries), hspMode, G_CALLBACK(onToggleHspMode), window);
  
  GtkUIManager *ui_manager = gtk_ui_manager_new();
  gtk_ui_manager_insert_action_group(ui_manager, action_group, 0);
  gtk_ui_manager_set_add_tearoffs(ui_manager, TRUE);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group(ui_manager);
  gtk_window_add_accel_group(GTK_WINDOW(window), accel_group);
  
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
static void belvuWindowCreateProperties(GtkWidget *belvuWindow, BelvuContext *bc)
{
  if (belvuWindow)
    {
      BelvuWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;

      g_object_set_data(G_OBJECT(belvuWindow), "BelvuWindowProperties", properties);
      g_signal_connect(G_OBJECT(belvuWindow), "destroy", G_CALLBACK (onDestroyBelvuWindow), NULL);
    }
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
  /* Create the window */
  GtkWidget *wrapWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setWrapWindowStyleProperties(wrapWindow);
  
  /* Create the context menu and set a callback to show it */
  GtkUIManager *uiManager = createUiManager(wrapWindow);
  GtkWidget *contextmenu = createBelvuMenu(wrapWindow, standardMenuDescription, "/WrapContextMenu", uiManager);
  
  gtk_widget_add_events(wrapWindow, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(wrapWindow), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(wrapWindow), vbox);
  
  /* Add the alignment section */
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
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
  GtkUIManager *uiManager = createUiManager(window);
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

  belvuWindowCreateProperties(window, bc);
  
  gtk_widget_show_all(window);
  

  return ok;
}
