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
#include <SeqTools/utilities.h>
#include <gdk/gdkkeysyms.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      4
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2
#define DEFAULT_SCROLL_STEP_INCREMENT	 5

/* Local function declarations */
static void			  onShowStatistics(GtkAction *action, gpointer data);
static void			  onPrint(GtkAction *action, gpointer data);
static void			  onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
static void			  onDrawPage(GtkPrintOperation *operation, GtkPrintContext *context, gint pageNum, gpointer data);

static MSP*			  mainWindowGetMspList(GtkWidget *mainWindow);


/* Menu builders */
static const GtkActionEntry mainMenuEntries[] = {
  { "Quit",	      NULL, "_Quit",	    "<control>Q",	"Quit the program",	  gtk_main_quit},
  { "Help",	      NULL, "_Help",	    "<control>H",	"Display help",		  onDisplayHelp},
  { "Print",	      NULL, "_Print",	    "<control>P",	"Print",		  G_CALLBACK(onPrint)},
  { "Settings",	      NULL, "_Settings",    "<control>S",	"Settings",		  NULL},
  { "Dotter",	      NULL, "_Dotter",	    NULL,		"Start Dotter",		  NULL},
  { "Statistics",     NULL, "S_tatistics",  "<control>T",	"Show memory statistics", G_CALLBACK(onShowStatistics)}
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
 * sequence if relevant. */
gchar *getSequenceSegment(GtkWidget *mainWindow, 
			const char const *sequence,
			const IntRange const *sequenceRange,
		        const int coord1, 
		        const int coord2,
		        const Strand strand,
			const BlxSeqType inputSeqType,
		        const int frame,
			const int numReadingFrames,
		        const gboolean reverse)
{
  gchar *result = NULL;
  
  int qMin = min(coord1, coord2); 
  int qMax = max(coord1, coord2);

  /* If the input coords are on a peptide sequence, convert them to DNA sequence coords. */
  if (inputSeqType == BLXSEQ_PEPTIDE)
    {
      qMin = convertPeptideToDna(qMin, frame, numReadingFrames);			/* 1st base in frame */
      qMax = convertPeptideToDna(qMax, frame, numReadingFrames) + numReadingFrames - 1; /* last base in frame */
    }
  
  /* Check that the requested segment is within the sequence's range */
  if (qMin < sequenceRange->min || qMax > sequenceRange->max)
    {
      messout ( "Requested query sequence %d - %d out of available range: %d - %d\n", qMin, qMax, sequenceRange->min, sequenceRange->max);
      if (inputSeqType == BLXSEQ_PEPTIDE) messout("Input coords on peptide sequence were %d - %d", coord1, coord2);
      return NULL;
    }
  
  /* Copy the portion of interest from the reference sequence and translate as necessary */
  BlxBlastMode mode = mainWindowGetBlastMode(mainWindow);
			    
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
      
      if (mode == BLXMODE_BLASTN) 
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


/* Called when the user selects the Print menu option, or hits the Print shortcut key */
static void onPrint(GtkAction *action, gpointer data)
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


/* Pop up a dialog reporting on various statistics about the process */
static void onShowStatistics(GtkAction *action, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  MSP *mspList = mainWindowGetMspList(mainWindow);
  showStatsDialog(mainWindow, mspList);
}


/* Pop up a dialog showing help information */
void onDisplayHelp()
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


/* Print handler - renders a specific page */
static void onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data)
{
  GtkWidget *mainWindow = GTK_WIDGET(data);
  
  gboolean toggled = mainWindowGetStrandsToggled(mainWindow);
  GtkWidget *bigPicture = mainWindowGetBigPicture(mainWindow);
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);

  GtkWidget *grid1 = toggled ? properties->revStrandGrid : properties->fwdStrandGrid;
  GtkWidget *grid2 = toggled ? properties->fwdStrandGrid : properties->revStrandGrid;
  GdkDrawable *grid1pixmap = gridGetProperties(grid1)->drawable;
  GdkDrawable *grid2pixmap = gridGetProperties(grid2)->drawable;

  /* Create a blank white pixmap */
  GdkDrawable *pixmap = gdk_pixmap_new(mainWindow->window, mainWindow->allocation.width, mainWindow->allocation.height, -1);
  GdkGC *gc = gdk_gc_new(pixmap);
  GdkColor fgColour = getGdkColor(GDK_WHITE);
  gdk_gc_set_foreground(gc, &fgColour);
  gdk_draw_rectangle(pixmap, gc, TRUE, 0, 0, mainWindow->allocation.width, mainWindow->allocation.height);

  /* Draw the grids onto the main pixmap */
  gdk_draw_drawable(pixmap, gc, grid1pixmap, 0, 0, 0, 0, -1, -1);
  gdk_draw_drawable(pixmap, gc, grid2pixmap, 0, 0, 0, grid1->allocation.height, -1, -1);
  
  cairo_t *cr = gtk_print_context_get_cairo_context (context);
  gdk_cairo_set_source_pixmap(cr, pixmap, 0, 0);
  cairo_paint(cr);
}


static gboolean onButtonPressMainWindow(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3)
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
  }
  
  return FALSE;
}


static gboolean keyPressLeftRight(GtkWidget *window, GdkEventKey *event)
{
  MainWindowProperties *properties = mainWindowGetProperties(window);
  GtkWidget *detailView = properties->detailView;
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  IntRange *fullRange = mainWindowGetFullRange(window);
  IntRange *displayRange = &detailViewProperties->displayRange;
  
  if (detailViewProperties->selectedBaseIdx != UNSET_INT)
    {
      /* Move the selection left/right by 1 base */
      if (event->keyval == GDK_Left)
	{
	  detailViewProperties->selectedBaseIdx -= 1;
	}
      else if (event->keyval == GDK_Right)
	{
	  detailViewProperties->selectedBaseIdx += 1;
	}
      
      /* Limit it to within the reference sequence range */
      if (detailViewProperties->selectedBaseIdx < fullRange->min)
	{
	  detailViewSetSelectedBaseIdx(detailView, fullRange->min);
	}
      else if (detailViewProperties->selectedBaseIdx > fullRange->max)
	{
	  detailViewSetSelectedBaseIdx(detailView, fullRange->max);
	}
      
      /* If we've moved outside the current display range, scroll by 1 base.
       * (This should probably also jump to the the selected base if it was previously
       * out of view - at the moment it just scrolls by 1, which is usually not enough 
       * to bring it into view.) */
      if (detailViewProperties->selectedBaseIdx > displayRange->max)
	{
	  scrollDetailViewRight1(detailView);
	}
      else if (detailViewProperties->selectedBaseIdx < displayRange->min)
	{
	  scrollDetailViewLeft1(detailView);
	}
    }
  
  return TRUE;
}


static gboolean onKeyPressMainWindow(GtkWidget *window, GdkEventKey *event, gpointer data)
{
  gboolean result = FALSE;
  
  if (event->keyval == GDK_Left || event->keyval == GDK_Right)
    {
      keyPressLeftRight(window, event);
      result = TRUE;
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
				       char *refSeq,
				       char *displaySeq,
				       const IntRange const *refSeqRange,
				       const IntRange const *fullDisplayRange,
				       const BlxSeqType seqType,
				       char **geneticCode,
				       int numReadingFrames)
{
  if (widget)
    {
      MainWindowProperties *properties = g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->detailView = detailView;
      properties->mspList = mspList;
      properties->blastMode = blastMode;
      properties->refSeq = refSeq;
      properties->displaySeq = displaySeq;
      properties->refSeqRange.min = refSeqRange->min;
      properties->refSeqRange.max = refSeqRange->max;
      properties->fullDisplayRange.min = fullDisplayRange->min;
      properties->fullDisplayRange.max = fullDisplayRange->max;
      properties->strandsToggled = FALSE;
      properties->seqType = seqType;
      properties->geneticCode = geneticCode;
      properties->numReadingFrames = numReadingFrames;
      properties->selectedMsps = NULL;
      
      properties->printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(properties->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
      
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

char * mainWindowGetRefSeq(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->refSeq : NULL;
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

static MSP* mainWindowGetMspList(GtkWidget *mainWindow)
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

GList* mainWindowGetSelectedMsps(GtkWidget *mainWindow)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return properties ? properties->selectedMsps : NULL;
}


/* Update function to be called whenever the MSP selection has changed */
static void selectionChanged(GtkWidget *mainWindow, const gboolean updateTrees)
{
  GtkWidget *detailView = mainWindowGetDetailView(mainWindow);
  
  if (updateTrees)
    {
      callFuncOnAllDetailViewTrees(detailView, deselectAllRows);
      callFuncOnAllDetailViewTrees(detailView, selectRowsForSelectedMsps);
    }
  
  /* Update the feedback box to tell the user which sequence is selected. */
  updateFeedbackBox(detailView);
  
  /* Redraw the grids */
  gtk_widget_queue_draw(mainWindowGetBigPicture(mainWindow));
}


void mainWindowSelectMsp(GtkWidget *mainWindow, MSP *msp, const gboolean updateTrees)
{
  if (!mainWindowIsMspSelected(mainWindow, msp))
    {
      MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
      properties->selectedMsps = g_list_prepend(properties->selectedMsps, msp);
      selectionChanged(mainWindow, updateTrees);
    }
}

void mainWindowDeselectMsp(GtkWidget *mainWindow, MSP *msp, const gboolean updateTrees)
{
  if (mainWindowIsMspSelected(mainWindow, msp))
    {
      MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
      properties->selectedMsps = g_list_remove(properties->selectedMsps, msp);
      selectionChanged(mainWindow, updateTrees);
    }
}

void mainWindowDeselectAllMsps(GtkWidget *mainWindow, const gboolean updateTrees)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);

  if (g_list_length(properties->selectedMsps) > 0)
    {
      g_list_free(properties->selectedMsps);
      properties->selectedMsps = NULL;
      selectionChanged(mainWindow, updateTrees);
    }
}

gboolean mainWindowIsMspSelected(GtkWidget *mainWindow, MSP *msp)
{
  MainWindowProperties *properties = mainWindowGetProperties(mainWindow);
  return (g_list_find(properties->selectedMsps, msp) != NULL);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set various properties for the main window widget */
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
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, mainMenuDescription, -1, &error))
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
			    const int refSeqOffset)
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
  
  /* Create a vertical-paned widget to put our two main frames in */
  GtkWidget *panedWidget = gtk_vpaned_new();
  gtk_box_pack_start(GTK_BOX(vbox), panedWidget, TRUE, TRUE, 0);
  
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
					   panedWidget, 
					   &fwdStrandGrid, 
					   &revStrandGrid, 
					   &fullDisplayRange);
  printf("Done.\n");
  
  printf("Creating detail view...\n");
  GtkWidget *detailView = createDetailView(window,
					   panedWidget, 
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
			     mspList, 
			     blastMode, 
			     refSeq, 
			     displaySeq,
			     &refSeqRange, 
			     &fullDisplayRange, 
			     seqType, 
			     geneticCode,
			     numReadingFrames);
  
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
