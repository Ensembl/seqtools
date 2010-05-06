/*
 *  blxdotter.c
 *  acedb
 *
 *  Created by Gemma Barson on 03/02/2010.
 *
 */

#include <SeqTools/blxdotter.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/bigpicture.h>
#include <SeqTools/utilities.h>
#include <SeqTools/dotter.h>

typedef struct _DotterDialogData
  {
    GtkWidget *blxWindow;
    
    GtkWidget *autoButton;
    GtkWidget *manualButton;

    GtkWidget *startEntry;
    GtkWidget *endEntry;
    GtkWidget *zoomEntry;
  } DotterDialogData;


/* Local function declarations */
static gboolean	      smartDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotter_start_out, int *dotter_end_out);
static char*	      fetchSeqRaw(const char *seqname, const char *fetchMode);
static char*	      fetchSequence(const char *seqname, char *fetch_prog);
static char*	      getDotterSSeq(GtkWidget *blxWindow);


/*******************************************************************
 *                      Dotter settings dialog                     *
 *******************************************************************/

/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the dotter mode according to the toggle state of the 
 * "auto" button. */
static void onSaveDotterMode(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->autoDotter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the start parameter that was entered (if manual dotter
 * params are being used). */
static void onSaveDotterStart(GtkWidget *entry, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  /* Only save the parameter if we are using manual parameters */
  if (!blxContext->autoDotter)
    {
      blxContext->dotterStart = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
    }  
}

static void onSaveDotterEnd(GtkWidget *entry, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!blxContext->autoDotter)
    {
      blxContext->dotterEnd = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
    }  
}

static void onSaveDotterZoom(GtkWidget *entry, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!blxContext->autoDotter)
    {
      blxContext->dotterZoom = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
    }  
}


/* Called when the 'last saved' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end coords that were last saved */
static void onLastSavedButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  char *startString = convertIntToString(bc->dotterStart);
  char *endString = convertIntToString(bc->dotterEnd);
  char *zoomString = convertIntToString(bc->dotterZoom);
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the 'full range' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end of the full ref seq range */
static void onFullRangeButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  char *startString = convertIntToString(bc->displayRev ? bc->refSeqRange.max : bc->refSeqRange.min);
  char *endString = convertIntToString(bc->displayRev ? bc->refSeqRange.min : bc->refSeqRange.max);
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the 'big picture range' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end of the big picture's current display range */
static void onBpRangeButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(dialogData->blxWindow);
  const IntRange const *displayRange = bigPictureGetDisplayRange(bigPicture);
  
  char *startString = convertIntToString(bc->displayRev ? displayRange->max : displayRange->min);
  char *endString = convertIntToString(bc->displayRev ? displayRange->min : displayRange->max);
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the user has hit a response button on the dotter settings dialog */
static void onResponseDotterDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  DotterDialogData *dialogData = (DotterDialogData*)(data);
  
  switch (responseId)
    {
      case GTK_RESPONSE_ACCEPT:
	widgetCallAllCallbacks(GTK_WIDGET(dialog), NULL);
	destroy = callDotter(dialogData->blxWindow, FALSE); /* keep dialog open if not successful */
	break;
	
      case GTK_RESPONSE_APPLY:
	widgetCallAllCallbacks(GTK_WIDGET(dialog), NULL);
	destroy = FALSE;
	break;
	
      default:
	break;
    };

  if (destroy)
    {
      g_free(dialogData);
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/* Called when the auto/manual radio button is toggled. */
static void onRadioButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton)))
    {
      /* Recalculate auto start/end in case user has selected a different sequence */
      int autoStart = UNSET_INT, autoEnd = UNSET_INT;
      smartDotterRange(dialogData->blxWindow, getDotterSSeq(dialogData->blxWindow), &autoStart, &autoEnd);

      if (autoStart == UNSET_INT)
	autoStart = bc->displayRev ? bc->refSeqRange.max : bc->refSeqRange.min;
      
      if (autoEnd == UNSET_INT)
	autoEnd = bc->displayRev ? bc->refSeqRange.min : bc->refSeqRange.max;
      
      char *startString = convertIntToString(autoStart);
      char *endString = convertIntToString(autoEnd);
      char *zoomString = convertIntToString(bc->dotterZoom);

      /* Display the new values */
      gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
      
      g_free(startString);
      g_free(endString);
      g_free(zoomString);
      
      /* Lock out the entry boxes so they cannot be edited */
      gtk_widget_set_sensitive(dialogData->startEntry, FALSE);
      gtk_widget_set_sensitive(dialogData->endEntry, FALSE);
      gtk_widget_set_sensitive(dialogData->zoomEntry, FALSE);
    }
  else
    {
      /* Manual coords. Leave values as they are but unlock the boxes so they can be edited. */
      gtk_widget_set_sensitive(dialogData->startEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->endEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->zoomEntry, TRUE);
    }
}


static GtkWidget* createTextEntry(GtkTable *table, 
				  int col, 
				  int row, 
				  const int xpad, 
				  const int ypad, 
				  char *title,
				  GtkCallback callbackFunc,
				  GtkWidget *blxWindow,
				  const int initValue)
{
  /* Create a label in the given column */
  GtkWidget *label = gtk_label_new(title);
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create the text entry widget in the next column in the same row */
  ++col;
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  /* Set the initial text */
  char *initText = convertIntToString(initValue);
  gtk_entry_set_text(GTK_ENTRY(entry), initText);
  g_free(initText);
  
  /* Add the callback data. This specifies what callback to use when the user hits OK or Apply on the dialog. */
  widgetSetCallbackData(entry, callbackFunc, blxWindow);
  
  return entry;
}


/* Pop up a dialog to allow the user to edit dotter parameters and launch dotter */
void showDotterDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Dotter settings", 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  GTK_STOCK_SAVE,
						  GTK_RESPONSE_APPLY,
						  GTK_STOCK_EXECUTE,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  GtkContainer *contentArea = GTK_CONTAINER(GTK_DIALOG(dialog)->vbox);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const int spacing = 4;

  /* Create a container for the child widgets */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_container_add(contentArea, GTK_WIDGET(hbox));

  /* Radio buttons for auto/manual params */
  GtkBox *vbox3 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox3), FALSE, FALSE, spacing);
  
  GtkWidget *autoButton = gtk_radio_button_new_with_mnemonic(NULL, "_Auto");
  gtk_box_pack_start(vbox3, autoButton, FALSE, FALSE, spacing);
  widgetSetCallbackData(autoButton, onSaveDotterMode, blxWindow);
		   
  GtkWidget *manualButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(autoButton), "_Manual");
  gtk_box_pack_start(vbox3, manualButton, FALSE, FALSE, spacing);

  /* Buttons that the user can click to populate the parameter boxes with certain values */
  GtkWidget *lastSavedButton = gtk_button_new_with_mnemonic("_Last saved ->");
  GtkWidget *fullRangeButton = gtk_button_new_with_mnemonic("_Full range ->");
  GtkWidget *bpRangeButton = gtk_button_new_with_mnemonic("_Big picture range ->");

  gtk_box_pack_start(vbox3, lastSavedButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox3, fullRangeButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox3, bpRangeButton, FALSE, FALSE, spacing);

  /* Disable last-saved button if no saved values exist */
  if (bc->dotterStart == UNSET_INT)
    {
      gtk_widget_set_sensitive(lastSavedButton, FALSE);
    }

  /* Create the text entry boxes for the dotter parameters */
  GtkTable *table = GTK_TABLE(gtk_table_new(3, 2, FALSE));
  gtk_box_pack_start(hbox, GTK_WIDGET(table), FALSE, FALSE, spacing);
  int xpad = 4, ypad = 4;
  
  GtkWidget *startEntry = createTextEntry(table, 1, 1, xpad, ypad, "<i>Start:</i>", onSaveDotterStart, blxWindow, bc->dotterStart);
  GtkWidget *endEntry = createTextEntry(table, 1, 2, xpad, ypad, "<i>End:</i>", onSaveDotterEnd, blxWindow, bc->dotterEnd);
  GtkWidget *zoomEntry = createTextEntry(table, 1, 3, xpad, ypad, "<i>Zoom:</i>", onSaveDotterZoom, blxWindow, bc->dotterZoom);

  /* Create the data struct we need to pass to the toggled callback, and connect signals */
  DotterDialogData *dialogData = g_malloc(sizeof(DotterDialogData));
  dialogData->blxWindow = blxWindow;
  dialogData->autoButton = autoButton;
  dialogData->manualButton = manualButton;
  dialogData->startEntry = startEntry;
  dialogData->endEntry = endEntry;
  dialogData->zoomEntry = zoomEntry;
  
  /* There is an issue if the user selects a different sequence while the dotter dialog
   * is still open: the auto range does not update automatically for the new sequence. To 
   * mitigate this, connect the 'clicked' signal as well as the toggle signal, so that they can
   * click on the 'auto' toggle button and have it refresh, even if that button is already selected.*/
  g_signal_connect(G_OBJECT(autoButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(manualButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(autoButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(manualButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  
  g_signal_connect(G_OBJECT(lastSavedButton), "clicked", G_CALLBACK(onLastSavedButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(fullRangeButton), "clicked", G_CALLBACK(onFullRangeButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(bpRangeButton), "clicked", G_CALLBACK(onBpRangeButtonClicked), dialogData);

  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDotterDialog), dialogData);
  
  /* Set the initial state of the toggle buttons and entry widgets */
  if (bc->autoDotter)
    {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manualButton), TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(autoButton), TRUE);
    }
  else
    {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(autoButton), TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manualButton), TRUE);
    }
  
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


/* Get the start/end coords. If the autoDotter flag is set, calculate coords
 * automatically - otherwise use the stored manual coords */
void getDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotterStart, int *dotterEnd, int *dotterZoom)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (!bc->autoDotter)
    {
      /* Use manual coords */
      *dotterStart = bc->dotterStart;
      *dotterEnd = bc->dotterEnd;
      *dotterZoom = bc->dotterZoom;
      
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
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the selected sequence name */
  if (g_list_length(bc->selectedSeqs) > 0)
    {
      const char *seqName = (const char*)(bc->selectedSeqs->data);
      GList *mspList = blxWindowGetSequenceMsps(blxWindow, seqName);

      /* If we're in seqbl mode, only part of the sequence is in the MSP. */
      const BlxBlastMode blastMode = bc->blastMode;
      if (blastMode != BLXMODE_TBLASTN)
	{
	  const char *fetchMode = bc->fetchMode;
	  dotterSSeq = fetchSeqRaw(seqName, fetchMode);
	  
	  /* If the match is on the reverse s strand, we need to modify it, because
	   * dotter does not currently handle it. */
	  if (dotterSSeq && g_list_length(mspList) > 0)
	    {
	      const MSP *msp = (const MSP*)(mspList->data);
	      const gboolean displayRev = bc->displayRev;
	      const gboolean qForward = (mspGetRefStrand(msp) == FORWARD_STRAND);
	      
	      if (qForward && displayRev)
		{
		  /* This should be consolidated with the code below */
		  blxComplement(dotterSSeq);
		  //g_strreverse(dotterSSeq);
		}
	    }
	}

      if (!dotterSSeq)
	{
	  /* Check if sequence is passed from acedb */
	  if (blastMode != BLXMODE_TBLASTX)
	    {
	      printf("Looking for sequence stored internally ... ");
	      
	      /* Loop through all MSPs in the selected sequence */
	      GList *mspListItem = mspList;
	      
	      for ( ; mspListItem ; mspListItem = mspListItem->next)
		{
		  MSP *msp = (MSP*)(mspListItem->data);
		  
		  if (msp->sseq != bc->paddingSeq)
		    {
		      dotterSSeq = g_strdup(msp->sseq);
		      break;
		    }
		}
	      
	      if (!dotterSSeq) printf("not ");
	      printf("found\n");
	      
	      /* If the match is on the reverse s strand, we need to modify it, because
	       * dotter does not currently handle it. */
	      if (dotterSSeq && g_list_length(mspList) > 0)
		{
		  const MSP *msp = (const MSP*)(mspList->data);
		  const gboolean displayRev = bc->displayRev;
		  const gboolean sForward = (mspGetMatchStrand(msp) == FORWARD_STRAND);
		  const gboolean qForward = (mspGetRefStrand(msp) == FORWARD_STRAND);
		  const gboolean sameDirection = (qForward == sForward);
		  
		  if (bc->seqType == BLXSEQ_DNA && ((sameDirection && displayRev) || (!sForward && !displayRev)))
		    {
		      /* Complementing the match sequence here maintains existing dotter behaviour.
		       * However, I think this is wrong - the match shows agains the wrong strand
		       * in the alignment view in dotter. I think we should reverse it here and
		       * not complement it; however, then the dot view shows the match in the 
		       * wrong place, because it doesn't currently handle reverse s coords. */
		      blxComplement(dotterSSeq);
		      //g_strreverse(dotterSSeq);
		    }
		}
	    }
	}

      if (dotterSSeq && (strchr(dotterSSeq, SEQUENCE_CHAR_PAD) || blastMode == BLXMODE_TBLASTN))
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
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Check that a sequence is selected */
  GList *selectedSeqs = bc->selectedSeqs;
  if (g_list_length(selectedSeqs) < 1)
    {
      return result;
    }

  const IntRange const *bigPicRange = bigPictureGetDisplayRange(bc->bigPicture);
  char activeStrand = (bc->displayRev ? '-' : '+') ;

  /* Loop through all MSPs in the selected sequence. We'll estimate the wanted
   * query region from the extent of the HSP's that are completely within view. */
  const char *selectedSeqName = (const char*)(selectedSeqs->data);
  int qMin = UNSET_INT, qMax = UNSET_INT;
  
  GList *mspListItem = blxWindowGetSequenceMsps(blxWindow, selectedSeqName);
  for ( ; mspListItem ; mspListItem = mspListItem->next)
    {
      const MSP *msp = (MSP*)(mspListItem->data);
      const int qFrame = mspGetRefFrame(msp, bc->seqType);
      
      /* Get the msp start/end in terms of display coords, and find the min/max */
      int base1, base2;
      const int coord1 = convertDnaIdxToDisplayIdx(msp->qstart, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &base1);
      const int coord2 = convertDnaIdxToDisplayIdx(msp->qend, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &base2);
      const int minMspCoord = min(coord1, coord2);
      const int maxMspCoord = max(coord1, coord2);

      /* Check if the MSP is in a visible tree row and is entirely within the big picture range */
      if ((msp->qframe[1] == activeStrand || (bc->blastMode == BLXMODE_BLASTN)) &&
	  (minMspCoord >= bigPicRange->min && maxMspCoord <= bigPicRange->max))
	{
	  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
	  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
	  
	  /* Extrapolate qMin backwards to the start of the match sequence (i.e. where
	   * s==0) and qMax forwards to the end of the match sequence (i.e. where s==sLength). */
	  int distToSMin = sSeqMin - 1;
	  int distToSMax = 200; /* default amount if sequence not found or if mode is tblastn */
	  
	  if (bc->blastMode != BLXMODE_TBLASTN && dotterSSeq)
	    {
	      distToSMax = strlen(dotterSSeq) - sSeqMax;
	    }

	  /* If the match sequence is a peptide sequence, convert the number of peptide
	   * coords we want to traverse to the equivalent number of DNA coords */
	  if (bc->seqType == BLXSEQ_PEPTIDE)
	    {
	      distToSMin *= bc->numFrames;
	      distToSMax *= bc->numFrames;
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

      if (bc->blastMode == BLXMODE_BLASTX || bc->blastMode == BLXMODE_TBLASTX)
	{
	  /* If sstart and send weren't in the end exons, we'll miss those - add some more */
	  extend = 0.2 * (qMax - qMin) ;
	  qMin -= extend ;
	  qMax += extend ;
	}

      /* Keep it within bounds */
      boundsLimitValue(&qMin, &bc->refSeqRange);
      boundsLimitValue(&qMax, &bc->refSeqRange);

      /* Apply min and max limits:  min 500 residues, max 10 Mb dots */
      int numDnaCoords = qMax - qMin;
      int midCoord = qMin + numDnaCoords/2;

      if (numDnaCoords < 500)
	{
	  numDnaCoords = 500;
	}

      const int numPeptideCoords = (bc->seqType == BLXSEQ_PEPTIDE) ? numDnaCoords / bc->numFrames : numDnaCoords;
      if (numDnaCoords * numPeptideCoords > 1e7)
	{
	  numDnaCoords = 1e7 / numPeptideCoords;
	}

      qMin = midCoord - (numDnaCoords / 2) ;
      qMax = midCoord + (numDnaCoords / 2) ;

      /* Bounds check again */
      boundsLimitValue(&qMin, &bc->refSeqRange);
      boundsLimitValue(&qMax, &bc->refSeqRange);

      /* Return the start/end. The values start low and end high in normal 
       * left-to-right display, or vice-versa if the display is reversed. */
      *dotter_start_out = bc->displayRev ? qMax : qMin;
      *dotter_end_out = bc->displayRev ? qMin : qMax;

      result = TRUE;
    }

  return result;
}


/* Get a sequence entry using either efetch or pfetch. */
static char *fetchSeqRaw(const char *seqname, const char *fetchMode)
{
  char *result = NULL ;

  if (!*seqname)
    {
      messout ( "Nameless sequence - skipping Efetch\n");
    }
  else
    {
      char *fetch_prog = blxGetFetchProg(fetchMode);
      char *seq_buf = fetchSequence(seqname, fetch_prog);
      
      if (seq_buf)
	{
	  result = g_strdup(seq_buf);
	}
    }
  
  return result ;
}


/* Common routine to call efetch or pfetch to retrieve a sequence entry. */
static char *fetchSequence(const char *seqname, char *fetch_prog)
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
      
      if (!feof(pipe))
	{
	  g_string_append_c(resultString, inputChar);
	}
    }

  if (ferror(pipe) || resultString->len < 1)
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

/* Call dotter. Returns true if dotter was called; false if we quit trying. */
gboolean callDotter(GtkWidget *blxWindow, const gboolean hspsOnly)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  const int numSeqsSelected = g_list_length(bc->selectedSeqs);
  
  if (numSeqsSelected < 1)
    {
      messout("Select a sequence first");
      return FALSE;
    }
  else if (numSeqsSelected > 1)
    {
      messout("Dotter cannot be called on multiple sequences. Select a single sequence and try again.");
      return FALSE;
    }
  
  /* Check this sequence is a valid blast match (just check the first MSP;
   * they must all the same type if they have the same seq name) */
  const char *selectedSeqName = (const char*)(bc->selectedSeqs->data);
  GList *selectedMsps = blxWindowGetSequenceMsps(blxWindow, selectedSeqName);
  const MSP *firstMsp = (const MSP*)(selectedMsps->data);

  if (!mspIsBlastMatch(firstMsp))
    {
      messout("Select a valid match sequence first");
      return FALSE;
    }
  
  /* Get the match sequence. Blixem uses g_malloc consistently now to allocate 
   * strings but unfortunately dotter will free this string with messfree, so 
   * we need to copy the result into a string allocated with messalloc. */
  char *dotterSSeqTemp = getDotterSSeq(blxWindow);
  if (!dotterSSeqTemp)
    {
      printf("Error starting Dotter: failed to fetch subject sequence\n");
      messout("Error starting Dotter: failed to fetch subject sequence\n");
      return FALSE;
    }
  
  char *dotterSSeq = messalloc(strlen(dotterSSeqTemp) + 1);
  strcpy(dotterSSeq, dotterSSeqTemp);
  g_free(dotterSSeqTemp);
  
  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT, dotterZoom = 0;
  getDotterRange(blxWindow, dotterSSeq, &dotterStart, &dotterEnd, &dotterZoom);
  
  /* Get the reference sequence name */
  const char *dotterQName = bc->refSeqName;
  
  /* Get the section of reference sequence that we're interested in */
  const Strand strand = bc->seqType == BLXSEQ_DNA ? mspGetRefStrand(firstMsp) : blxWindowGetActiveStrand(blxWindow);
  const int frame = mspGetRefFrame(firstMsp, bc->seqType);
  
  char *querySeqSegmentTemp = getSequenceSegment(bc,
						 bc->refSeq,
						 dotterStart,
						 dotterEnd, 
						 strand,
						 BLXSEQ_DNA,	  /* calculated dotter coords are always in terms of DNA seq */
						 frame,
						 FALSE,		  /* input coords are always left-to-right, even if display reversed */
						 bc->displayRev,  /* whether to reverse */
						 bc->displayRev,  /* whether to allow rev strands to be complemented */
						 FALSE);	  /* don't allow translation to a peptide seq */
  
  if (!querySeqSegmentTemp)
    {
      messerror("Cannot start dotter - failed to get query sequence.");
      return FALSE;
    }
  
  
  /* Again, dotter will free this with messfree, so we need to pass a string allocated with messalloc */
  char *querySeqSegment = messalloc(strlen(querySeqSegmentTemp) + 1);
  strcpy(querySeqSegment, querySeqSegmentTemp);
  g_free(querySeqSegmentTemp);
  
  
  /* Get the match sequence name (chopping off the letters before the colon, if there is one). */
  char *dotterSName = strchr(firstMsp->sname, ':');
  if (dotterSName)
    {
      dotterSName++;
    }
  else
    {
      dotterSName = firstMsp->sname;
    }
  
  const int offset = min(dotterStart, dotterEnd) - 1;
  
  /* Get the options */
  static char opts[] = "     ";
  opts[0] = bc->displayRev ? 'R' : ' ';
  opts[1] = hspsOnly ? 'H' : ' ';
  opts[2] = bc->gappedHsp ? 'G' : ' ';
  
  /* Get the mode */
  char type = getDotterMode(bc->blastMode);
  
  /* Get the list of all MSPs */
  printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);
  
  printf("  query sequence: name -  %s, offset - %d\n"
	 "subject sequence: name -  %s, offset - %d\n", dotterQName, offset, dotterSName, 0);
  
  dotter(type, opts, dotterQName, querySeqSegment, offset, dotterSName, dotterSSeq, 0,
	 0, 0, NULL, NULL, NULL, 0.0, dotterZoom, bc->mspList, bc->refSeqRange.min - 1, 0, 0);
  
  return TRUE;
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
