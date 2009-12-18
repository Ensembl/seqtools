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
#include <SeqTools/bigpicturemspline.h>
#include <SeqTools/utilities.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      4
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2


/* Local function declarations */
static void			  blxShowStats(GtkAction *action, gpointer data);
static const char const*	  mainWindowGetRefSeq(GtkWidget *mainWindow);
static MSP*			  mainWindowGetMspList(GtkWidget *mainWindow);


/* Menu builders */
static const GtkActionEntry mainMenuEntries[] = {
  { "Quit",	  NULL, "_Quit",	"<control>Q",	"Quit the program",	  gtk_main_quit},
  { "Help",	  NULL, "_Help",	"<control>H",	"Display help",		  blxHelp},
  { "Print",	  NULL, "_Print",	"<control>P",	"Print",		  NULL},
  { "Settings",	  NULL, "_Settings",	NULL,		"Change settings",	  NULL},
  { "Dotter",	  NULL, "_Dotter",	NULL,		"Start Dotter",		  NULL},
  { "Statistics", NULL, "_Statistics",	NULL,		"Show memory statistics", G_CALLBACK(blxShowStats)}
};


static const char *mainMenuDescription =
"<ui>"
"  <popup name='MainMenu'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
"      <menuitem action='Print'/>"
"      <menuitem action='Settings'/>"
"      <separator/>"
"      <menuitem action='Dotter'/>"
"      <separator/>"
"      <menuitem action='Statistics'/>"
"  </popup>"
"</ui>";


/***********************************************************
 *			   Menu Utilities                  *
 ***********************************************************/

static void showModalDialog(char *title, char *messageText)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  NULL, 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Add the message */
  //  GtkWidget *contentArea = gtk_dialog_get_content_area(GTK_DIALOG(dialog)); //not in pre 2.14 versions
  GtkWidget *vbox = GTK_DIALOG(dialog)->vbox;
  GtkWidget *label = gtk_label_new(messageText);
  gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);
  
  /* Ensure dialog is destroyed when user responds */
  g_signal_connect_swapped(dialog, "response", G_CALLBACK(gtk_widget_destroy), dialog);
  
  gtk_widget_show_all(dialog);
}


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
                                                  GTK_WINDOW(blixemWindow),
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  /* Ensure that the dialog box (along with any children) is destroyed when the user responds. */
  g_signal_connect_swapped (dialog, "response", G_CALLBACK(gtk_widget_destroy), dialog);
  
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
 *			  Menu actions                     *
 ***********************************************************/

/* Pop up a dialog reporting on various statistics about the process */
static void blxShowStats(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MSP *mspList = mainWindowGetMspList(mainWindow);
  showStatsDialog(mainWindow, mspList);
}

/* Pop up a dialog showing help information */
void blxHelp(void)
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
  
  showModalDialog("Help", messageText);
}


/***********************************************************
 *			   Events                          *
 ***********************************************************/

static gboolean onButtonPressMainWindow(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3)
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
  }
  
  return FALSE;
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
				       MSP *mspList,
				       const BlxBlastMode blastMode,
				       const char const *refSeq)
{
  if (widget)
    {
      MainWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->detailView = detailView;
      properties->mspList = mspList;
      properties->blastMode = blastMode;
      properties->refSeq = refSeq;
      properties->strandsToggled = FALSE;
      
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

BlxBlastMode mainWindowGetBlastMode(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->blastMode : FALSE;
}

static const char const* mainWindowGetRefSeq(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->refSeq : NULL;
}

static MSP* mainWindowGetMspList(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->mspList : FALSE;
}

/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

static void setStyleProperties(GtkWidget *widget)
{
  gtk_container_set_border_width (GTK_CONTAINER(widget), DEFAULT_WINDOW_BORDER_WIDTH); 

  /* Set the default font size to be a bit smaller than usual */
  int origSize = pango_font_description_get_size(widget->style->font_desc) / PANGO_SCALE;
  const char *origFamily = pango_font_description_get_family(widget->style->font_desc);

  char parseString[500];
  sprintf(parseString, "gtk-font-name = \"%s %d\"", origFamily, origSize + DEFAULT_FONT_SIZE_ADJUSTMENT);
  gtk_rc_parse_string(parseString);
}


static GtkWidget* createMainMenu(GtkWidget *window)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  gtk_action_group_add_actions (action_group, mainMenuEntries, G_N_ELEMENTS (mainMenuEntries), window);
  
  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);
  
  GError *error = NULL;
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, mainMenuDescription, -1, &error))
    {
      g_message ("building menus failed: %s", error->message);
      g_error_free (error);
      exit (EXIT_FAILURE);
    }
  
  return gtk_ui_manager_get_widget (ui_manager, "/MainMenu");
}


GtkWidget* createMainWindow(char *refSeq, 
			    MSP *mspList, 
			    BlxBlastMode blastMode,
			    BlxSeqType seqType, 
			    int numReadingFrames)
{
  int refSeqEnd = strlen(refSeq) ;
  IntRange refSeqRange = {1, refSeqEnd};
  IntRange bigPictureDisplayRange  = {1, refSeqEnd};
  
  /* Create the main window */
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setStyleProperties(window);
  
  /* Create a vertical box to pack everything in */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), vbox);
  
  /* Create the main menu */
  GtkWidget *mainmenu = createMainMenu(window);
  
  /* Create a vertical-paned widget to put our two main frames in */
  GtkWidget *panedWidget = gtk_vpaned_new();
  gtk_box_pack_start(GTK_BOX(vbox), panedWidget, TRUE, TRUE, 0);
  
  /* Create the widgets. We need a single adjustment for the entire detail view, which will also be referenced
   * by the big picture view, so create it first. */
  GtkAdjustment *detailAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, refSeqRange.max - refSeqRange.min + 1, 1, 0, 0));
  
  GtkWidget *fwdStrandGrid = NULL, *revStrandGrid = NULL, *header = NULL;
  
  printf("Creating big picture...\n");
  GtkWidget *bigPicture = createBigPicture(window,
					   panedWidget, 
					   &header, 
					   &fwdStrandGrid, 
					   &revStrandGrid, 
					   &bigPictureDisplayRange, 
					   &refSeqRange);
  printf("Done.\n");
  
  printf("Creating detail view...\n");
  GtkWidget *detailView = createDetailView(window,
					   panedWidget, 
					   detailAdjustment, 
					   fwdStrandGrid, 
					   revStrandGrid,
					   mspList,
					   refSeq,
					   blastMode,
					   seqType,
					   numReadingFrames,
					   &refSeqRange);
  printf("Done.\n");
  
  /* Create a custom scrollbar for scrolling the sequence column and put it at the bottom of the window */
  GtkWidget *scrollBar = createDetailViewScrollBar(detailAdjustment, detailView);
  gtk_box_pack_start(GTK_BOX(vbox), scrollBar, FALSE, TRUE, 0);
  
  /* Set required data in the main window widget */
  mainWindowCreateProperties(window, bigPicture, detailView, mspList, blastMode, refSeq);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressMainWindow), mainmenu);
  
  /* Show the window */
  printf("realizing widgets...\n");
  gtk_widget_show_all(window);
  printf("Done.\n");
  
  return window;
}
