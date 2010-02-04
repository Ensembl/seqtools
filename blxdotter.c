/*
 *  blxdotter.c
 *  acedb
 *
 *  Created by Gemma Barson on 03/02/2010.
 *
 */

#include <SeqTools/blxdotter.h>
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/utilities.h>
#include <SeqTools/dotter.h>


typedef struct _DotterDialogData
  {
    GtkWidget *dialog;
    GtkWidget *blxWindow;
    
    GtkWidget *startCoordWidget;
    GtkWidget *endCoordWidget;
    
    GtkWidget *autoButton;
    GtkWidget *manualButton;
    GtkWidget *lastSavedButton;
    GtkWidget *fullRangeButton;
    
    char *autoStart;
    char *autoEnd;
    char *lastSavedStart;
    char *lastSavedEnd;
    char *fullRangeStart;
    char *fullRangeEnd;
  } DotterDialogData;

/* Local function declarations */
static gboolean	      smartDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotter_start_out, int *dotter_end_out);
static char*	      fetchSeqRaw(char *seqname);
static char*	      fetchSequence(char *seqname, char *fetch_prog);
static void	      blxCallDotter(GtkWidget *blxWindow, const gboolean hspsOnly);
static char*	      getDotterSSeq(GtkWidget *blxWindow);


/*******************************************************************
 *                      Dotter settings dialog                     *
 *******************************************************************/

static void dotterDialogSaveSettings(DotterDialogData *dialogData)
{
  MainWindowProperties *properties = mainWindowGetProperties(dialogData->blxWindow);
  GtkEntry *startEntry = GTK_ENTRY(dialogData->startCoordWidget);
  GtkEntry *endEntry = GTK_ENTRY(dialogData->endCoordWidget);
  
  properties->autoDotterParams = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton));
  
  if (!properties->autoDotterParams)
    {
      /* Save the manual parameters entered */
      properties->dotterStart = atoi(gtk_entry_get_text(startEntry));
      properties->dotterEnd = atoi(gtk_entry_get_text(endEntry));
    }  
}


/* Called when the user has hit a response button on the dotter settings dialog */
static void onResponseDotterDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)(data);
  gboolean destroy = TRUE;
  
  switch (responseId)
    {
      case GTK_RESPONSE_ACCEPT:
	dotterDialogSaveSettings(dialogData);
	blxCallDotter(dialogData->blxWindow, FALSE);
	break;
	
      case GTK_RESPONSE_APPLY:
	dotterDialogSaveSettings(dialogData);
	destroy = FALSE;
	break;
	
      default:
	break;
    };

  if (destroy)
    {
      g_free(dialogData->autoStart);
      g_free(dialogData->autoEnd);
      g_free(dialogData->lastSavedStart);
      g_free(dialogData->lastSavedEnd);
      g_free(dialogData);

      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/* Called when the 'full range' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end of the full ref seq range */
static void onFullRangeButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;

  gtk_entry_set_text(GTK_ENTRY(dialogData->startCoordWidget), dialogData->fullRangeStart);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endCoordWidget), dialogData->fullRangeEnd);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the 'last saved' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end coords that were last saved */
static void onLastSavedButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startCoordWidget), dialogData->lastSavedStart);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endCoordWidget), dialogData->lastSavedEnd);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the auto/manual radio button is toggled. */
static void onRadioButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  
  GtkEntry *startEntry = GTK_ENTRY(dialogData->startCoordWidget);
  GtkEntry *endEntry = GTK_ENTRY(dialogData->endCoordWidget);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton)))
    {
      gtk_entry_set_text(startEntry, dialogData->autoStart);
      gtk_entry_set_text(endEntry, dialogData->autoEnd);
      
      gtk_widget_set_sensitive(GTK_WIDGET(startEntry), FALSE);
      gtk_widget_set_sensitive(GTK_WIDGET(endEntry), FALSE);
    }
  else
    {
      /* Manual coords. Leave values as they are but unlock the boxes so they can be edited. */
      gtk_widget_set_sensitive(GTK_WIDGET(startEntry), TRUE);
      gtk_widget_set_sensitive(GTK_WIDGET(endEntry), TRUE);
    }
}


/* Converts the given integer to a string. The result must be free'd with g_free */
static char* convertIntToString(const int value)
{
  char result[numDigitsInInt(value) + 1];
  sprintf(result, "%d", value);
  return g_strdup(result);
}


/* Pop up a dialog to allow the user to edit dotter parameters and launch dotter */
void showDotterDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Dotter settings", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_EXECUTE,
						  GTK_RESPONSE_ACCEPT,
						  GTK_STOCK_SAVE,
						  GTK_RESPONSE_APPLY,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  NULL);
  
  GtkContainer *contentArea = GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(dialog)));
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  int spacing = 4;

  const IntRange const *refSeqRange = mainWindowGetRefSeqRange(blxWindow);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(blxWindow);

  /* Find the various possibilities for start/end coords */
  char *lastSavedStart = convertIntToString(mainWindowGetDotterStart(blxWindow));
  char *lastSavedEnd = convertIntToString(mainWindowGetDotterEnd(blxWindow));
  char *fullRangeStart = convertIntToString(rightToLeft ? refSeqRange->max : refSeqRange->min);
  char *fullRangeEnd = convertIntToString(rightToLeft ? refSeqRange->min : refSeqRange->max);

  int autoStartCoord = UNSET_INT, autoEndCoord = UNSET_INT;
  if (getDotterSSeq(blxWindow))
    {
      smartDotterRange(blxWindow, getDotterSSeq(blxWindow), &autoStartCoord, &autoEndCoord);
    }

  char *autoStart = autoStartCoord != UNSET_INT ? convertIntToString(autoStartCoord) : fullRangeStart;
  char *autoEnd = autoEndCoord != UNSET_INT ? convertIntToString(autoEndCoord) : fullRangeEnd;
  
  /* Create a container for the child widgets */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_container_add(contentArea, GTK_WIDGET(hbox));

  /* Radio buttons for auto/manual params */
  GtkBox *vbox3 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox3), FALSE, FALSE, spacing);
  
  GtkWidget *autoButton = gtk_radio_button_new_with_mnemonic(NULL, "_Auto");
  GtkWidget *manualButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(autoButton), "_Manual");
  gtk_box_pack_start(vbox3, autoButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox3, manualButton, FALSE, FALSE, spacing);

  /* Buttons to populate the entry boxes with the last-saved or full-range values */
  GtkWidget *lastSavedButton = gtk_button_new_with_mnemonic("_Last saved ->");
  GtkWidget *fullRangeButton = gtk_button_new_with_mnemonic("_Full range ->");
  gtk_box_pack_start(vbox3, lastSavedButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox3, fullRangeButton, FALSE, FALSE, spacing);

  /* Disable last-saved button if no saved values exist */
  if (mainWindowGetDotterStart(blxWindow) == UNSET_INT)
    {
      gtk_widget_set_sensitive(lastSavedButton, FALSE);
    }

  
  /* coord labels */
  GtkBox *vbox1 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox1), FALSE, FALSE, spacing);

  GtkWidget *label1 = gtk_label_new("<i>Start:</i>");
  gtk_label_set_use_markup(GTK_LABEL(label1), TRUE);
  gtk_box_pack_start(vbox1, label1, FALSE, FALSE, spacing);

  GtkWidget *label2 = gtk_label_new("<i>End:</i>");
  gtk_label_set_use_markup(GTK_LABEL(label2), TRUE);
  gtk_box_pack_start(vbox1, label2, FALSE, FALSE, spacing);
  
  /* coord entry boxes */
  GtkBox *vbox2 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox2), FALSE, FALSE, spacing);
  
  GtkWidget *startCoordWidget = gtk_entry_new();
  gtk_box_pack_start(vbox2, startCoordWidget, FALSE, FALSE, spacing);
  gtk_entry_set_activates_default(GTK_ENTRY(startCoordWidget), TRUE);
  
  GtkWidget *endCoordWidget = gtk_entry_new();
  gtk_box_pack_start(vbox2, endCoordWidget, FALSE, FALSE, spacing);
  gtk_entry_set_activates_default(GTK_ENTRY(endCoordWidget), TRUE);
  
  /* Set the initial state of the toggle buttons and entry widgets */
  if(mainWindowGetAutoDotter(blxWindow))
    {
      gtk_entry_set_text(GTK_ENTRY(startCoordWidget), autoStart);
      gtk_entry_set_text(GTK_ENTRY(endCoordWidget), autoEnd);
      gtk_widget_set_sensitive(startCoordWidget, FALSE);
      gtk_widget_set_sensitive(endCoordWidget, FALSE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(autoButton), TRUE);
    }
  else
    {
      gtk_entry_set_text(GTK_ENTRY(startCoordWidget), lastSavedStart);
      gtk_entry_set_text(GTK_ENTRY(endCoordWidget), lastSavedEnd);
      gtk_widget_set_sensitive(startCoordWidget, TRUE);
      gtk_widget_set_sensitive(endCoordWidget, TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manualButton), TRUE);
    }
  
  /* Connect signals and show dialog */
  DotterDialogData *data = g_malloc(sizeof(DotterDialogData));
  data->dialog = dialog;
  data->blxWindow = blxWindow;
  data->startCoordWidget = startCoordWidget;
  data->endCoordWidget = endCoordWidget;
  data->autoButton = autoButton;
  data->manualButton = manualButton;
  data->lastSavedButton = lastSavedButton;
  data->fullRangeButton = fullRangeButton;
  data->autoStart = autoStart;
  data->autoEnd = autoEnd;
  data->lastSavedStart = lastSavedStart;
  data->lastSavedEnd = lastSavedEnd;
  data->fullRangeStart = fullRangeStart;
  data->fullRangeEnd = fullRangeEnd;
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDotterDialog), data);
  g_signal_connect(autoButton, "toggled", G_CALLBACK(onRadioButtonToggled), data);
  g_signal_connect(manualButton, "toggled", G_CALLBACK(onRadioButtonToggled), data);
  g_signal_connect(lastSavedButton, "clicked", G_CALLBACK(onLastSavedButtonClicked), data);
  g_signal_connect(fullRangeButton, "clicked", G_CALLBACK(onFullRangeButtonClicked), data);

  gtk_widget_show_all(dialog);
 
}


/*******************************************************************
 *		      Dotter utility functions                     *
 *******************************************************************/

/* Convert a blixem blast mode to a char for passing to dotter */
char getDotterMode(const BlxBlastMode blastMode)
{
  char type = ' ';
  
  if (blastMode == BLXMODE_BLASTP || blastMode == BLXMODE_TBLASTN)
    {
      type = 'P';
    }
  else if (blastMode == BLXMODE_BLASTX)
    {
      type = 'X';
    }
  else if (blastMode == BLXMODE_BLASTN || blastMode == BLXMODE_TBLASTX)
    {
      type = 'N';
    }
  
  return type;
}


/* Get the start/end coords. If the autoDotterParams flag is set, calculate coords
 * automatically - otherwise use the stored manual coords */
void getDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotterStart, int *dotterEnd)
{
  const gboolean autoDotterParams = mainWindowGetAutoDotter(blxWindow);
  
  if (!autoDotterParams)
    {
      /* Use manual coords */
      *dotterStart = mainWindowGetDotterStart(blxWindow);
      *dotterEnd = mainWindowGetDotterEnd(blxWindow);
      
      if (*dotterStart == UNSET_INT || *dotterEnd == UNSET_INT)
	{
	  messerror("Manual dotter parameters were requested but one or more coord is not set. Values are start=%d, end=%d. Attempting to calculate parameters.");
	}
    }
  
  if (*dotterStart == UNSET_INT || *dotterEnd == UNSET_INT)
    {
      if (!smartDotterRange(blxWindow, dotterSSeq, dotterStart, dotterEnd))
	{
	  messout("Could not start dotter - error calculating coordinate range.");
	  return;
	}
    }
}


/* Utility to fetch the selected match sequence or get it from the selected MSP.
 * This function assumes that if multiple MSPs are selected, that they are all for 
 * the same match sequence. returns null if no MSPs are selected */
static char* getDotterSSeq(GtkWidget *blxWindow)
{
  char *dotterSSeq = NULL;
  
  /* Get the sequence name from the first selected MSP. */
  GList *selectedMspList = mainWindowGetSelectedMsps(blxWindow);
   
  if (g_list_length(selectedMspList) > 0)
    {
      /* If we're in seqbl mode, only part of the sequence is in in the MSP. */
      const BlxBlastMode blastMode = mainWindowGetBlastMode(blxWindow);
      if (blastMode != BLXMODE_TBLASTN)
	{
	  MSP *msp = (MSP*)(selectedMspList->data);
	  dotterSSeq = fetchSeqRaw(msp->sname);
	}

      if (!dotterSSeq)
	{
	  /* Check if sequence is passed from acedb */
	  if (blastMode != BLXMODE_TBLASTX)
	    {
	      printf("Looking for sequence stored internally ... ");
	      
	      GList *listItem = selectedMspList;
	      for ( ; listItem ; listItem = listItem->next)
		{
		  const MSP *msp = (MSP*)listItem->data;
		  
		  //	      if (msp->sseq != padseq) //to do: implement this, if still required.
		  {
		    dotterSSeq = g_strdup(msp->sseq);
		    break;
		  }
		}
	      
	      if (!dotterSSeq) printf("not ");
	      printf("found\n");
	    }
	}

      if (!dotterSSeq)
	{
	  printf("Can't fetch subject sequence for dotter - aborting\n");
	  messout("Can't fetch subject sequence for dotter - aborting\n");
	}
      else if (strchr(dotterSSeq, '-') || blastMode == BLXMODE_TBLASTN)
	{
	  messout("Note: the sequence passed to dotter is incomplete");
	}
    }
  
  return dotterSSeq;
}


//static int smartDotterRangeSelf(void)
//{
//    int len, mid;
//
//    len = 2000;
//    mid = dispstart + plusmin*displen/2;
//
//    dotterStart = mid - plusmin*len/2;
//    dotterEnd = mid + plusmin*len/2;
//
//    /* Keep it within bounds */
//    if (dotterStart < 1) dotterStart = 1;
//    if (dotterStart > qlen) dotterStart = qlen;
//    if (dotterEnd > qlen) dotterEnd = qlen;
//    if (dotterEnd < 1) dotterEnd = 1;

//    return 1;
//}


/* Attempts to set the range of dotter in some sort of sensible way. The problem is that
 * hits can occur over a much wider range than the user is looking at, so the function
 * attempts to find the range of hits that corresponds to what the user can see.
 * Returns TRUE if it managed to find sequences and set a sensible range, FALSE otherwise.
 * NOTE: This function assumes that only a single sequence can be selected at any one time. */
static gboolean smartDotterRange(GtkWidget *blxWindow,
				 const char *dotterSSeq, 
				 int *dotter_start_out, 
				 int *dotter_end_out)
{
  gboolean result = FALSE;

  GtkWidget *bigPicture = mainWindowGetBigPicture(blxWindow);
  const IntRange const *bigPicRange = bigPictureGetDisplayRange(bigPicture);
  const BlxBlastMode blastMode = mainWindowGetBlastMode(blxWindow);
  const BlxSeqType seqType = mainWindowGetSeqType(blxWindow);
  const int numFrames = mainWindowGetNumReadingFrames(blxWindow);
  const gboolean rightToLeft = mainWindowGetStrandsToggled(blxWindow);
  char activeStrand = (rightToLeft ? '-' : '+') ;

  
  /* Loop through all selected MSPs. We'll estimate the wanted query region from 
   * the extent of the HSP's that are completely within view. */
  GList *listItem = mainWindowGetSelectedMsps(blxWindow);
  int qMin = UNSET_INT, qMax = UNSET_INT;
  const char *selectedSeqName = NULL;
  
  for ( ; listItem; listItem = listItem->next)
    {
      const MSP *msp = (MSP*)listItem->data;
      const int minMspCoord = min(msp->displayStart, msp->displayEnd);
      const int maxMspCoord = max(msp->displayStart, msp->displayEnd);

      /* Check if the MSP is in a visible tree row and is entirely within the big picture range */
      if ((msp->qframe[1] == activeStrand || (blastMode == BLXMODE_BLASTN)) &&
	  (minMspCoord >= bigPicRange->min && maxMspCoord <= bigPicRange->max))
	{
	  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
	  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
	  
	  if (!selectedSeqName)
	    {
	      selectedSeqName = msp->sname;
	    }

	  /* Extrapolate qMin backwards to the start of the match sequence (i.e. where
	   * s==0) and qMax forwards to the end of the match sequence (i.e. where s==sLength). */
	  int distToSMin = sSeqMin - 1;
	  int distToSMax = 200; /* default amount if sequence not found or if mode is tblastn */
	  
	  if (blastMode != BLXMODE_TBLASTN && dotterSSeq)
	    {
	      distToSMax = strlen(dotterSSeq) - sSeqMax;
	    }

	  /* If the match sequence is a peptide sequence, convert the number of peptide
	   * coords we want to traverse to the equivalent number of DNA coords */
	  if (seqType == BLXSEQ_PEPTIDE)
	    {
	      distToSMin *= numFrames;
	      distToSMax *= numFrames;
	    }

	  /* If the strands are in opposite directions, the low end of the ref 
	   * sequence corresponds to the high of the match sequence, and vice versa. */
	  const gboolean sameDirection = (msp->qframe[1] == msp->sframe[1]);
	  qSeqMin -= sameDirection ? distToSMin : distToSMax;
	  qSeqMax += sameDirection ? distToSMax : distToSMin;

	  /* Set the results */
	  if (qMin == UNSET_INT || qSeqMin < qMin)
	    {
	      qMin = qSeqMin;
	    }
	  
	  if (qMax == UNSET_INT || qSeqMax > qMax)
	    {
	      qMax = qSeqMax;
	    }
	}
    }


  if (qMin == UNSET_INT)
    {
      messout("Could not find any matches on the '%c' strand to %s.", activeStrand, selectedSeqName);
      result = FALSE;
    }
  else
    {
      /* Due to gaps, we might miss the ends - add some more */
      int extend = 0.1 * (qMax - qMin) ;
      qMin -= extend ;
      qMax += extend ;

      if (blastMode == BLXMODE_BLASTX || blastMode == BLXMODE_TBLASTX)
	{
	  /* If sstart and send weren't in the end exons, we'll miss those - add some more */
	  extend = 0.2 * (qMax - qMin) ;
	  qMin -= extend ;
	  qMax += extend ;
	}

      /* Keep it within bounds */
      const IntRange const *refSeqRange = mainWindowGetRefSeqRange(blxWindow);
      boundsLimitValue(&qMin, refSeqRange);
      boundsLimitValue(&qMax, refSeqRange);

      /* Apply min and max limits:  min 500 residues, max 10 Mb dots */
      int numDnaCoords = qMax - qMin;
      int midCoord = qMin + numDnaCoords/2;

      if (numDnaCoords < 500)
	{
	  numDnaCoords = 500;
	}

      const int numPeptideCoords = (seqType == BLXSEQ_PEPTIDE) ? numDnaCoords / 3 : numDnaCoords;
      if (numDnaCoords * numPeptideCoords > 1e7)
	{
	  numDnaCoords = 1e7 / numPeptideCoords;
	}

      qMin = midCoord - (numDnaCoords / 2) ;
      qMax = midCoord + (numDnaCoords / 2) ;

      /* Bounds check again */
      boundsLimitValue(&qMin, refSeqRange);
      boundsLimitValue(&qMax, refSeqRange);

      /* Return the start/end. The values start low and end high in normal 
       * left-to-right display, or vice-versa if the display is reversed. */
      *dotter_start_out = rightToLeft ? qMax : qMin;
      *dotter_end_out = rightToLeft ? qMin : qMax;

      result = TRUE;
    }

  return result;
}


/* Get a sequence entry using either efetch or pfetch. */
static char *fetchSeqRaw(char *seqname)
{
  char *result = NULL ;

  if (!*seqname)
    {
      messout ( "Nameless sequence - skipping Efetch\n");
    }
  else
    {
      char *fetch_prog = blxGetFetchProg();
      char *seq_buf = fetchSequence(seqname, fetch_prog);
      
      if (seq_buf)
	{
	  result = g_strdup(seq_buf);
	}
    }
  
  return result ;
}


/* Common routine to call efetch or pfetch to retrieve a sequence entry. */
static char *fetchSequence(char *seqname, char *fetch_prog)
{
  char *result = NULL;
  
  static char *fetchstr = NULL ;

  /* I have no idea why blixem trys to open this file, it doesn't work unless you are in
   * a directory that you have write access to anyway....I await an answer from Eric.
   * Meanwhile at least now we record the error rather than crashing later in this
   * routine when we can't open the file. */
  FILE *outputFile = fopen("myoutput", "w");
  
  if (!outputFile)
    {
      messerror("%s", "Cannot open blixem sequence output file: \"myoutput\"") ;
    }

  if (!strcmp(fetch_prog, "pfetch"))
    {
      /* --client gives logging information to pfetch server,
       * -q  Sequence only output (one line) */
      fetchstr = hprintf(0, "%s --client=acedb_%s_%s -q '%s' &",
			 fetch_prog, getSystemName(), getLogin(TRUE), seqname) ;
    }
  else
    {
      fetchstr = hprintf(0, "%s -q '%s'", fetch_prog, seqname) ;
    }

  printf("%sing %s...\n", fetch_prog, seqname);

  /* Try and get the sequence, if we overrun the buffer then we need to try again. */
  FILE *pipe = (FILE*)popen(fetchstr, "r");
  
  if (!pipe)
    {
      messcrash("Failed to open pipe %s\n", fetchstr);
    }
  
  /* Create a string to read the file into. We don't know the length we need, so
   * use a GString, which will grow automatically. */
  GString *resultString = g_string_sized_new(200);
  
  while (!feof(pipe))
    {
      unsigned char inputChar = fgetc(pipe);
      g_string_append_c(resultString, inputChar);
    }

  if (ferror(pipe))
    {
      g_string_free(resultString, TRUE);
    }
  else
    {
      result = resultString->str;
      g_string_free(resultString, FALSE);

      if (outputFile)
	{
	  fprintf(outputFile, "%s", result);
	}
    }

  /* Clean up */
  pclose(pipe);

  if (fetchstr)
    {
      messfree(fetchstr);
    }

  if (outputFile)
    {
      fclose(outputFile);
    }
  
  return result;
}


/*******************************************************************
 *		      Functions to call dotter                     *
 *******************************************************************/

static void blxCallDotter(GtkWidget *blxWindow, const gboolean hspsOnly)
{
  GList *selectedMsps = mainWindowGetSelectedMsps(blxWindow);
  if (g_list_length(selectedMsps) < 1)
    {
      messout("Select a sequence first");
      return;
    }
  
  /* Get the match sequence. Blixem uses g_malloc consistently now to allocate 
   * strings but unfortunately dotter will free this string with messfree, so 
   * we need to copy the result into a string allocated with messalloc. */
  char *dotterSSeqTemp = getDotterSSeq(blxWindow);
  char *dotterSSeq = messalloc(strlen(dotterSSeqTemp) + 1);
  strcpy(dotterSSeq, dotterSSeqTemp);
  g_free(dotterSSeqTemp);
  
  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT;
  getDotterRange(blxWindow, dotterSSeq, &dotterStart, &dotterEnd);
  
  /* Get the reference sequence name */
  const char *dotterQName = mainWindowGetRefSeqName(blxWindow);
  
  /* Get the section of reference sequence that we're interested in */
  const MSP *msp = (MSP*)selectedMsps->data;
  const Strand strand = msp->qframe[1] == '+' ? FORWARD_STRAND : REVERSE_STRAND;
  
  char frameStr[2];
  frameStr[0] = msp->qframe[2];
  frameStr[1] = 0;
  int frame = atoi(frameStr);
  
  char *querySeqSegmentTemp = getSequenceSegment(blxWindow, 
						 mainWindowGetRefSeq(blxWindow),
						 mainWindowGetRefSeqRange(blxWindow),
						 dotterStart,
						 dotterEnd, 
						 strand,
						 BLXSEQ_DNA, /* calculated dotter coords are always in terms of DNA seq */
						 frame,
						 mainWindowGetNumReadingFrames(blxWindow),
						 TRUE,  /* allow sequence to be reversed if reverse strand */
						 FALSE); /* don't allow translation to a peptide seq */
  
  if (!querySeqSegmentTemp)
    {
      messerror("Cannot start dotter - failed to get query sequence.");
      return;
    }
  
  
  /* Again, dotter will free this with messfree, so we need to pass a string allocated with messalloc */
  char *querySeqSegment = messalloc(strlen(querySeqSegmentTemp) + 1);
  strcpy(querySeqSegment, querySeqSegmentTemp);
  g_free(querySeqSegmentTemp);
  
  
  /* Get the match sequence name (chopping off the letters before the colon, if there is one). */
  char *dotterSName = strchr(msp->sname, ':');
  if (dotterSName)
    {
      dotterSName++;
    }
  else
    {
      dotterSName = msp->sname;
    }
  
  /* Get the offset */
  const IntRange const *refSeqRange = mainWindowGetRefSeqRange(blxWindow);
  int offset = min(dotterStart, dotterEnd) - 1 + refSeqRange->min - 1;
  
  /* Get the options */
  static char opts[] = "     ";
  opts[0] = mainWindowGetStrandsToggled(blxWindow) ? 'R' : ' ';
  opts[1] = hspsOnly ? 'H' : ' ';
  opts[2] = mainWindowGetGappedHsp(blxWindow) ? 'G' : ' ';
  
  /* Get the mode */
  char type = getDotterMode(mainWindowGetBlastMode(blxWindow));
  
  /* Get the list of all MSPs */
  MSP *mspList = mainWindowGetMspList(blxWindow);
  
  int dotterZoom = 0; /* to do: implement this */
  
  
  printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);
  
  printf("  query sequence: name -  %s, offset - %d\n"
	 "subject sequence: name -  %s, offset - %d\n", dotterQName, refSeqRange->min - 1, dotterSName, 0);
  
  dotter(type, opts, dotterQName, querySeqSegment, offset, dotterSName, dotterSSeq, 0,
	 0, 0, NULL, NULL, NULL, 0.0, dotterZoom, mspList, refSeqRange->min - 1, 0, 0);
}


//static void callDotterSelf(const MSP const *msplist)
//{//
//    static char opts[] = "     ";
//    char type, *queryseq;
//    int  offset;
//
//    type = ' ';
//    if (blastp || tblastn || tblastx)
//      type = 'P';
//    else if (blastx || blastn)
//      type = 'N';
//
//    if (smartDotter)
//	if (!smartDotterRangeSelf())
//	  return;
//
//    if (!*dotterqname)
//      {
//	if (!*qname_G)
//	  strcpy(dotterqname, "Blixem-seq");
//	else
//	  strncpy(dotterqname, qname_G, LONG_NAMESIZE);
//	dotterqname[LONG_NAMESIZE] = 0;
//      }
//
//    /* Get query sequence */
//    /* Can't do reversed strand since Dotter can't reverse vertical scale */
//    /* Avoid translating queryseq by pretending to be blastn - very sneaky and dangerous */
//    if (blastx || tblastx)
//      blastn = 1;
//
//    if (!(queryseq = getqseq(min(dotterStart, dotterEnd), max(dotterStart, dotterEnd), q)))
//      return;
//
//    if (blastx || tblastx)
//      blastn = 0;
//
//    dottersseq = g_malloc(strlen(queryseq)+1);
//    strcpy(dottersseq, queryseq);
//
//    offset = min(dotterStart, dotterEnd)-1 + qoffset;
//
//    printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);
//
//    dotter(type, opts, dotterqname, queryseq, offset, dotterqname, dottersseq, offset,
//	   0, 0, NULL, NULL, NULL, 0.0, dotterZoom, msplist, qoffset, 0, 0);
//}
