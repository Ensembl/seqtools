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
#include <SeqTools/detailview.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blxmsp.h>
#include <SeqTools/dotter.h>
#include <string.h>
#include <stdlib.h>
#include <regular.h> /* for getSystemName etc. */

#define DEFAULT_DOTTER_RANGE_SELF	2000

typedef struct _DotterDialogData
  {
    GtkWidget *blxWindow;           /* pointer to the main blixem window */
    
    GtkWidget *autoButton;          /* the radio button on the dialog for automatic dotter parameters */
    GtkWidget *manualButton;        /* the radio button on the dialog for manual dotter parameters */

    GtkWidget *startEntry;          /* the text entry box on the dialog for the start coord */
    GtkWidget *endEntry;            /* the text entry box on the dialog for the end coord */
    GtkWidget *zoomEntry;           /* the text entry box on the dialog for the zoom value */
    
    gboolean callOnSelf;            /* whether to call dotter on the query seq versus itself */
    gboolean hspsOnly;              /* whether to call dotter on HSPs only */
    
    char *dotterSSeq;               /* the match sequence to call dotter on */
  } DotterDialogData;


/* Error codes and domain */
#define BLX_DOTTER_ERROR g_quark_from_string("Dotter")

typedef enum {
  BLX_DOTTER_ERROR_NO_SEQS,	      /* no sequences are selected */
  BLX_DOTTER_ERROR_TOO_MANY_SEQS,     /* too many sequences are selected */
  BLX_DOTTER_ERROR_INVALID_SEQ,       /* selected sequence(s) are invalid */
  BLX_DOTTER_ERROR_NOT_FOUND,	      /* failed to find the sequence */
  BLX_DOTTER_ERROR_NO_REF_SEQ,	      /* failed to find the query sequence segment */
  BLX_DOTTER_ERROR_INTERNAL_SEQ,      /* using internally-stored sequence (because fetch failed) */
  BLX_DOTTER_ERROR_NO_MATCHES,        /* there are no matches on the requested sequence */
  BLX_DOTTER_ERROR_NO_SEQ_DATA,       /* the match sequence has no sequence data (e.g. if could not pfetch it) */
  BLX_DOTTER_ERROR_INVALID_STRAND,    /* there are no match sequences on the correct strand */
  BLX_DOTTER_ERROR_NO_EXE             /* no dotter executable found */
} BlxDotterError;


/* Local function declarations */
static gboolean	      getDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, const gboolean callOnSelf, const gboolean autoDotter, int *dotterStart, int *dotterEnd, int *dotterZoom, GError **error);
static gboolean	      smartDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotter_start_out, int *dotter_end_out, GError **error);
static gboolean	      smartDotterRangeSelf(GtkWidget *blxWindow, int *dotter_start_out, int *dotter_end_out, GError **error);
static char*	      fetchSeqRaw(const char *seqname, const char *fetchMode);
static char*	      fetchSequence(const char *seqname, char *fetch_prog);
static gboolean	      callDotterSelf(GtkWidget *blxWindow, GError **error);


/*******************************************************************
 *                      Dotter settings dialog                     *
 *******************************************************************/

/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the dotter mode according to the toggle state of the 
 * "auto" button. */
static gboolean onSaveDotterMode(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->autoDotter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the start parameter that was entered (if manual dotter
 * params are being used). */
static gboolean onSaveDotterStart(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Only save the parameter if we are using manual parameters */
  if (!bc->autoDotter)
    {
      const int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      
      if (valueWithinRange(newVal, &bc->refSeqRange))
        {
          bc->dotterStart = newVal;
        }
      else
        {
          result = FALSE;
          g_critical("Start value %d is outside reference sequence range %d -> %d. Value not saved.\n", newVal, bc->refSeqRange.min, bc->refSeqRange.max);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterEnd(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!bc->autoDotter)
    {
      const int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      
      if (valueWithinRange(newVal, &bc->refSeqRange))
        {
          bc->dotterEnd = newVal;
        }
      else
        {
          result = FALSE;
          g_critical("End value %d is outside reference sequence range %d -> %d. Value not saved.\n", newVal, bc->refSeqRange.min, bc->refSeqRange.max);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterZoom(GtkWidget *entry, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!blxContext->autoDotter)
    {
      blxContext->dotterZoom = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
    }  
  
  return TRUE;
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

  const int qStart = convertDisplayIdxToDnaIdx(displayRange->min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int qEnd = convertDisplayIdxToDnaIdx(displayRange->max, bc->seqType, 1, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  
  char *startString = convertIntToString(qStart);
  char *endString = convertIntToString(qEnd);
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the "call on self" button is toggled on the dotter dialog */
static void onSelfButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  dialogData->callOnSelf = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  /* Recalculate auto start/end */
  int autoStart = UNSET_INT, autoEnd = UNSET_INT;
  const gboolean autoDotter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton));

  getDotterRange(dialogData->blxWindow, dialogData->dotterSSeq, dialogData->callOnSelf, autoDotter, &autoStart, &autoEnd, NULL, NULL);

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
}


/* Called when the "call on self" button is toggled on the dotter dialog */
static void onHspsButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  dialogData->hspsOnly = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
}


/* Called when the user has hit a response button on the dotter settings dialog */
static void onResponseDotterDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  DotterDialogData *dialogData = (DotterDialogData*)(data);
  
  GError *error = NULL;
  
  switch (responseId)
    {
      case GTK_RESPONSE_ACCEPT:
        /* Only continue to call dotter if saving the values is successful */
	if (widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId)))
          {
            if (dialogData->callOnSelf)
              {
                destroy = callDotterSelf(dialogData->blxWindow, &error);
              }
            else
              {
                destroy = callDotter(dialogData->blxWindow, dialogData->hspsOnly, dialogData->dotterSSeq, &error);
              }
          }
        else
          {
            destroy = FALSE; /* there was an error, so leave the dialog open */
          }
        
	break;
	
      case GTK_RESPONSE_APPLY:
	widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
	destroy = FALSE;
	break;
	
      default:
	break;
    };

  /* If any errors were found, report them */
  if (error)
    {
      prefixError(error, "Could not start Dotter. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  if (destroy)
    {
      /* The dialog is persistent so hide it rather than destroying it. */
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
}


static void onDestroyDotterDialog(GtkWidget *dialog, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  
  if (dialogData)
    {
      if (dialogData->dotterSSeq)
        {
          g_free(dialogData->dotterSSeq);
          dialogData->dotterSSeq = NULL;
        }
      
      g_free(dialogData);
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
      getDotterRange(dialogData->blxWindow, dialogData->dotterSSeq, dialogData->callOnSelf, TRUE, &autoStart, &autoEnd, NULL, NULL);

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
				  BlxResponseCallback callbackFunc,
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


/* Utility to get the title for the dotter dialog. Uses the selected sequence name if a single
 * sequence is selected, or shows <no sequences> or <multiple sequences>. The result should be
 * free'd with g_free. */
static const char* getDotterTitle(const BlxViewContext *bc)
{
  const char *result = NULL;
  
  GString *resultStr = g_string_new("Dotter sequence: ");
  
  const int numSeqs = g_list_length(bc->selectedSeqs);
  
  if (numSeqs == 1)
    {
      const BlxSequence *blxSeq = (const BlxSequence*)(bc->selectedSeqs->data);
      g_string_append(resultStr, blxSequenceGetDisplayName(blxSeq));
    }
  else if (numSeqs < 1)
    {
      g_string_append(resultStr, "<no sequences selected>");
    }
  else
    {
      g_string_append(resultStr, "<multiple sequences>");
    }
  
  result = resultStr->str;
  g_string_free(resultStr, FALSE);
  
  return result;
}


/* Pop up a dialog to allow the user to edit dotter parameters and launch dotter */
void showDotterDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_DOTTER;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  static DotterDialogData *dialogData = NULL;
  
  const char *title = getDotterTitle(bc);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons(title, 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_SAVE,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_EXECUTE,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      /* Create the dialog data struct first time round, but re-populate it each time. Create 
       * a destructor function that will free the struct. */
      dialogData = g_malloc(sizeof(DotterDialogData));
      dialogData->blxWindow = NULL;
      dialogData->autoButton = NULL;
      dialogData->manualButton = NULL;
      dialogData->startEntry = NULL;
      dialogData->endEntry = NULL;
      dialogData->zoomEntry = NULL;
      dialogData->callOnSelf = FALSE;
      dialogData->hspsOnly = FALSE;
      dialogData->dotterSSeq = NULL;

      g_signal_connect(G_OBJECT(dialog), "destroy", G_CALLBACK(onDestroyDotterDialog), dialogData);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDotterDialog), dialogData);
    }
  else
    {
      /* Refresh the dialog by clearing its contents an re-creating it */
      dialogClearContentArea(GTK_DIALOG(dialog));
      
      gtk_window_set_title(GTK_WINDOW(dialog), title);
    }
  
  GtkContainer *contentArea = GTK_CONTAINER(GTK_DIALOG(dialog)->vbox);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  const int spacing = 4;

  /* Create a container for the child widgets */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_container_add(contentArea, GTK_WIDGET(hbox));

  /* Radio buttons for auto/manual params */
  GtkBox *vbox1 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox1), FALSE, FALSE, spacing);
  
  GtkWidget *autoButton = gtk_radio_button_new_with_mnemonic(NULL, "_Auto");
  gtk_box_pack_start(vbox1, autoButton, FALSE, FALSE, spacing);
  widgetSetCallbackData(autoButton, onSaveDotterMode, blxWindow);
		   
  GtkWidget *manualButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(autoButton), "_Manual");
  gtk_box_pack_start(vbox1, manualButton, FALSE, FALSE, spacing);

  /* Buttons that the user can click to populate the parameter boxes with certain values */
  GtkWidget *lastSavedButton = gtk_button_new_with_mnemonic("_Last saved ->");
  GtkWidget *fullRangeButton = gtk_button_new_with_mnemonic("_Full range ->");
  GtkWidget *bpRangeButton = gtk_button_new_with_mnemonic("_Big picture range ->");

  gtk_box_pack_start(vbox1, lastSavedButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox1, fullRangeButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox1, bpRangeButton, FALSE, FALSE, spacing);

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

  /* Create check buttons for the various options */
  GtkBox *vbox3 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox3), FALSE, FALSE, spacing);
  
  GtkWidget *selfButton = gtk_check_button_new_with_mnemonic("Call on _self");
  gtk_box_pack_start(vbox3, selfButton, FALSE, FALSE, 0);
  GtkWidget *hspsButton = gtk_check_button_new_with_mnemonic("_HSPs only");
  gtk_box_pack_start(vbox3, hspsButton, FALSE, FALSE, 0);
  
  /* Create the data struct we need to pass to the toggled callback, and connect signals */
  dialogData->blxWindow = blxWindow;
  dialogData->autoButton = autoButton;
  dialogData->manualButton = manualButton;
  dialogData->startEntry = startEntry;
  dialogData->endEntry = endEntry;
  dialogData->zoomEntry = zoomEntry;
  dialogData->callOnSelf = FALSE;
  dialogData->hspsOnly = FALSE;
  
  if (dialogData->dotterSSeq)
    {
      g_free(dialogData->dotterSSeq);
      dialogData->dotterSSeq = NULL;
    }

  dialogData->dotterSSeq = getDotterSSeq(blxWindow, NULL);
  
  /* There is an issue if the user selects a different sequence while the dotter dialog
   * is still open: the auto range does not update automatically for the new sequence. To 
   * mitigate this, connect the 'clicked' signal as well as the toggle signal, so that they can
   * click on the 'auto' toggle button and have it refresh, even if that button is already selected.*/
//  g_signal_connect(G_OBJECT(autoButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
//  g_signal_connect(G_OBJECT(manualButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(autoButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(manualButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  
  g_signal_connect(G_OBJECT(lastSavedButton), "clicked", G_CALLBACK(onLastSavedButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(fullRangeButton), "clicked", G_CALLBACK(onFullRangeButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(bpRangeButton), "clicked", G_CALLBACK(onBpRangeButtonClicked), dialogData);

  g_signal_connect(G_OBJECT(selfButton), "toggled", G_CALLBACK(onSelfButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(hspsButton), "toggled", G_CALLBACK(onHspsButtonToggled), dialogData);
  
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
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
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


/* Get the start/end coords. If the passed autoDotter flag is true, calculate coords
 * automatically - otherwise use the stored manual coords */
static gboolean getDotterRange(GtkWidget *blxWindow, 
			       const char *dotterSSeq, 
			       const gboolean callOnSelf,
                               const gboolean autoDotter,
			       int *dotterStart, 
			       int *dotterEnd, 
			       int *dotterZoom, 
			       GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */

  gboolean success = TRUE;
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (!autoDotter)
    {
      /* Use manual coords */
      if (dotterStart) *dotterStart = bc->dotterStart;
      if (dotterEnd) *dotterEnd = bc->dotterEnd;
      if (dotterZoom) *dotterZoom = bc->dotterZoom;
      
      if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
	{
	  g_debug("Manual dotter parameters were requested but one or more coord is not set (start=%d, end=%d). Calculating automatic parameters instead.\n",
		  *dotterStart, *dotterEnd);
	}
    }
  
  if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
    {
      GError *tmpError = NULL;
      
      if (callOnSelf)
	success = smartDotterRangeSelf(blxWindow, dotterStart, dotterEnd, &tmpError);
      else
	success = smartDotterRange(blxWindow, dotterSSeq, dotterStart, dotterEnd, &tmpError);

      if (error && tmpError)
	{
	  g_propagate_error(error, tmpError);
	  prefixError(*error, "Failed to calculate dotter coordinates. ");
	}
    }
  
  return success;
}


/* Utility to fetch the selected match sequence or get it from the selected MSP.
 * This function assumes that if multiple MSPs are selected, that they are all for 
 * the same match sequence. Returns null if no MSPs are selected, with details of the error
 * in 'error'.  If the sequence was found but there were warnings, it returns non-null with
 * the warnings in 'error'. The return value should be free'd with g_free */
char* getDotterSSeq(GtkWidget *blxWindow, GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */
  
  char *dotterSSeq = NULL;
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the selected sequence name */
  if (g_list_length(bc->selectedSeqs) < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return dotterSSeq;
    }
  
  const BlxSequence *blxSeq = (const BlxSequence*)(bc->selectedSeqs->data);

  /* If we're in seqbl mode, only part of the sequence is in the MSP. */
  const BlxBlastMode blastMode = bc->blastMode;
  if (blastMode != BLXMODE_TBLASTN)
    {
      const char *fetchMode = bc->fetchMode;
      dotterSSeq = fetchSeqRaw(blxSequenceGetFullName(blxSeq), fetchMode);
    }

  if (!dotterSSeq && blastMode != BLXMODE_TBLASTX)
    {
      /* Check if sequence is stored internally (i.e. it was passed from acedb) */
      g_debug("Looking for sequence stored internally... ");
    
      dotterSSeq = g_strdup(blxSequenceGetSeq(blxSeq));
      
      if (!dotterSSeq)
	{
	  g_debug("not found.\n");
	  g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NOT_FOUND, "Failed to find sequence for '%s'.\n", blxSequenceGetFullName(blxSeq));
	  return FALSE;
	}

      g_debug("found.\n");
      
      /* Dotter expects the passed sequence to be forwards and un-complemented but if this is
       * the reverse strand the sequence will be complemented, so un-complement it. */
      if (blxSeq->strand == BLXSTRAND_REVERSE)
        {
          blxComplement(dotterSSeq);
        }
    }

  if (dotterSSeq && (strchr(dotterSSeq, SEQUENCE_CHAR_PAD) || blastMode == BLXMODE_TBLASTN))
    {
      g_warning("The sequence for '%s' is incomplete.\n", blxSequenceGetFullName(blxSeq));
    }
  
  return dotterSSeq;
}


/* Attempts to set the range of dotter in a sensible way, when calling dotter on the reference
 * sequence versus itself. */
static gboolean smartDotterRangeSelf(GtkWidget *blxWindow,
				     int *dotter_start_out,
				     int *dotter_end_out,
				     GError **error)
{
  /* We'll just use a large-ish range centred on the current display range */
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  const IntRange const *displayRange = detailViewGetDisplayRange(detailView);

  const int mid = getRangeCentre(displayRange);
  const int offset = DEFAULT_DOTTER_RANGE_SELF / 2;

  *dotter_start_out = mid - offset;
  *dotter_end_out = *dotter_start_out + DEFAULT_DOTTER_RANGE_SELF;

  const BlxViewContext *bc = blxWindowGetContext(blxWindow);
  boundsLimitValue(dotter_start_out, &bc->refSeqRange);
  boundsLimitValue(dotter_end_out, &bc->refSeqRange);
  
  return TRUE;
}


/* Attempts to set the range of dotter in some sort of sensible way. The problem is that
 * hits can occur over a much wider range than the user is looking at, so the function
 * attempts to find the range of hits that corresponds to what the user can see.
 * Returns TRUE if it managed to find sequences and set a sensible range, FALSE otherwise.
 * NOTE: This function assumes that only a single sequence can be selected at any one time. */
static gboolean smartDotterRange(GtkWidget *blxWindow,
				 const char *dotterSSeq, 
				 int *dotter_start_out, 
				 int *dotter_end_out,
				 GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */

  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Check that a sequence is selected */
  GList *selectedSeqs = bc->selectedSeqs;
  if (g_list_length(selectedSeqs) < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }

  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  const IntRange const *bigPicRange = bigPictureGetDisplayRange(bigPicture);
  const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD) ;

  /* Loop through all MSPs in the selected sequence. We'll estimate the wanted
   * query region from the extent of the HSP's that are completely within view. */
  const BlxSequence *selectedSeq = (const BlxSequence*)(selectedSeqs->data);
  int qMin = UNSET_INT, qMax = UNSET_INT;
  GList *mspListItem = selectedSeq->mspList;  
  
  for ( ; mspListItem ; mspListItem = mspListItem->next)
    {
      const MSP *msp = (MSP*)(mspListItem->data);
      const int qFrame = mspGetRefFrame(msp, bc->seqType);
      
      /* Get the msp start/end in terms of display coords, and find the min/max */
      int base1, base2;
      const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &base1);
      const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, qFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &base2);
      const int minMspCoord = min(coord1, coord2);
      const int maxMspCoord = max(coord1, coord2);

      /* Check if the MSP is in a visible tree row and is entirely within the big picture range */
      if ((msp->qStrand == activeStrand || (bc->blastMode == BLXMODE_BLASTN)) &&
	  (minMspCoord >= bigPicRange->min && maxMspCoord <= bigPicRange->max))
	{
	  int qSeqMin = msp->qRange.min;
	  int qSeqMax = msp->qRange.max;
	  int sSeqMin = msp->sRange.min;
	  int sSeqMax = msp->sRange.max;
	  
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
	  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
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
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_MATCHES, 
		  "Could not find any matches on the '%c' strand of the selected sequence '%s'.", 
		  activeStrand, blxSequenceGetFullName(selectedSeq));

      return FALSE;
    }
  
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
  
  return TRUE;
}


/* Get a sequence entry using either efetch or pfetch. */
static char *fetchSeqRaw(const char *seqname, const char *fetchMode)
{
  char *result = NULL ;

  if (!*seqname)
    {
      g_warning( "Nameless sequence - skipping %s\n", fetchMode);
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
      g_warning("%s", "Cannot open blixem sequence output file: 'myoutput'.\n") ;
    }

  if (!strcmp(fetch_prog, "pfetch"))
    {
      /* --client gives logging information to pfetch server,
       * -q  Sequence only output (one line) */
      fetchstr = hprintf(0, "%s --client=%s_%s_%s -q '%s' &",
			 fetch_prog, g_get_prgname(), getSystemName(), getLogin(TRUE), seqname) ;
    }
  else
    {
      fetchstr = hprintf(0, "%s -q '%s'", fetch_prog, seqname) ;
    }

  g_debug("%sing %s...\n", fetch_prog, seqname);

  /* Try and get the sequence, if we overrun the buffer then we need to try again. */
  FILE *pipe = (FILE*)popen(fetchstr, "r");
  
  if (!pipe)
    {
      g_critical("Failed to open pipe %s\n", fetchstr);
      return result;
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

/* Find an executable and return its complete pathname.
 */
static int findCommand (char *command, char **retp)
{
#if !defined(NO_POPEN)
  static char retstr[1025] ;
  char *path, file[1025], retval;
  int found=0;
  
  /* Don't use csh - fails if the path is not set in .cshrc * /
   if (access(csh, X_OK)) {
   g_critical("Could not find %s", csh);
   return 0;
   }
   if (!(pipe = (FILE *)popen(messprintf("%s -cf \"which %s\"", csh, command), "r"))) {
   return 0;
   }
   
   while (!feof(pipe))
   fgets(retval, 1024, pipe);
   retval[1024] = 0;
   pclose(pipe);
   
   if (cp = strchr(retval, '\n')) *cp = 0;
   if (retp) *retp = retval;
   
   / * Check if whatever "which" returned is an existing and executable file * /
   if (!access(retval, F_OK) && !access(retval, X_OK))
   return 1;
   else
   return 0;
   */
  
  path = g_malloc(strlen(getenv("PATH"))+1);
  /* Don't free 'path' since it changes later on - never mind, 
   we're only calling it once */
  
  strcpy(path, getenv("PATH"));
  path = strtok(path, ":");
  while (path) {
    strcpy(file, path);
    strcat(file,"/");
    strcat(file, command);
    if (!access(file, F_OK) && !access(file, X_OK)) {
      found = 1;
      break;
    }
    
    path = strtok(0, ":");
  }
  
  if (found) {
    strcpy(retstr, file);
    retval = 1;
  }
  else {
    strcpy(retstr, "Can't find executable 'dotter' in path.\n");
    retval = 0;
  }
  
  if (retp) *retp = retstr;
  return retval;
  
#endif
}


/* Call dotter as an external process */
gboolean callDotterExternal(BlxViewContext *bc,
                            int dotterZoom, 
                            IntRange *refSeqRange,
                            char *refSeqSegment,
                            const char *dotterSName,
                            char *dotterSSeq,
                            const gboolean hspsOnly,
                            char *Xoptions,
                            GError **error)
{
#if !defined(NO_POPEN)

  const int qlen = strlen(refSeqSegment);
  const int slen = strlen(dotterSSeq);
  
  boundsLimitRange(refSeqRange, &bc->refSeqRange, FALSE);

  int xstart = refSeqRange->min;
  int xend = refSeqRange->max;
  int ystart = 1;
  int yend = slen;
  
  static char *dotterBinary = NULL;
  
  /* Open pipe to new dotterBinary */
  if (!dotterBinary) 
    { 
      printf("Looking for Dotter ...\n");
      
      if (!findCommand("dotter", &(dotterBinary))) 
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_EXE, "No dotter executable found in path.\n$PATH=%s\n", getenv("PATH"));
          dotterBinary = NULL;
          return FALSE;
        }
    }
  
  g_debug("Calling %s with region: %d,%d - %d,%d\n", dotterBinary, xstart, ystart, xend, yend);
  fflush(stdout);
  
  char *pipeText = blxprintf("/bin/csh -cf \"%s -z %d -q %d -s %d%s%s -S '%s' %d '%s' %d %s %s\"", 
                             dotterBinary, 
                             dotterZoom, 
                             xstart - 1, 
                             ystart - 1, 
                             bc->displayRev ? " -r" : "",
                             hspsOnly ? " -H" : "",
                             bc->refSeqName, 
                             qlen, 
                             dotterSName, 
                             slen,
                             dotterBinary,
                             (Xoptions ? Xoptions : ""));
  
  g_debug("Sending to pipe: %s\n", pipeText);
  
  FILE *pipe = (FILE *)popen(pipeText, "w");
  
  fwrite(refSeqSegment, 1, qlen, pipe);
  fwrite(dotterSSeq, 1, slen, pipe);
  
  /* Pass on features */
  GList *seqItem = bc->matchSeqs;
  for ( ; seqItem; seqItem = seqItem->next) 
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      writeBlxSequenceToOutput(pipe, blxSeq, refSeqRange);
    }
  
  fprintf(pipe, "%c\n", EOF);
  fflush(pipe);
#endif
  
  return TRUE;
}


/* Call dotter. Returns true if dotter was called; false if we quit trying. */
gboolean callDotter(GtkWidget *blxWindow, const gboolean hspsOnly, char *dotterSSeqIn, GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const int numSeqsSelected = g_list_length(bc->selectedSeqs);
  
  if (numSeqsSelected < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }
  else if (numSeqsSelected > 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_TOO_MANY_SEQS, "Cannot run dotter on multiple sequences.\n");
      return FALSE;
    }
  
  /* Check this sequence is a valid blast match (just check the first MSP;
   * they must all the same type if they have the same seq name) */
  const BlxSequence *selectedSeq = (const BlxSequence*)(bc->selectedSeqs->data);
  const MSP *firstMsp = (const MSP*)(selectedSeq->mspList->data);

  if (!mspIsBlastMatch(firstMsp))
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_INVALID_SEQ, "You must select a valid match sequence first.\n");
      return FALSE;
    }
  
  /* We will display the active strand as the main strand in dotter */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);

  if (bc->seqType == BLXSEQ_PEPTIDE)
    {
      GList *mspItem = selectedSeq->mspList;
      gboolean found = FALSE;
      
      for ( ; mspItem && !found; mspItem = mspItem->next)
        {
          const MSP const *msp = (const MSP const *)(mspItem->data);
          found = (msp->qStrand == qStrand);
        }
      
      if (!found)
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_INVALID_STRAND, "You must select a match sequence on the active strand, or toggle strands first.\n");
          return FALSE;
        }
    }
  
  /* Make a copy of the match sequence, because dotter takes ownership of this. */
  char *dotterSSeq = NULL;
  if (dotterSSeqIn)
    {
      dotterSSeq = g_strdup(dotterSSeqIn);
    }
  else
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQ_DATA, "No sequence data for this sequence.\n");
      return FALSE;
    }
  
  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT, dotterZoom = 0;
  GError *tmpError = NULL;
  
  if (!getDotterRange(blxWindow, dotterSSeq, FALSE, bc->autoDotter, &dotterStart, &dotterEnd, &dotterZoom, &tmpError))
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const int frame = mspGetRefFrame(firstMsp, bc->seqType);
  IntRange dotterRange;
  intrangeSetValues(&dotterRange, dotterStart, dotterEnd);

  char *querySeqSegment = getSequenceSegment(bc->refSeq,
                                             &dotterRange,
                                             BLXSTRAND_FORWARD,   /* always pass forward strand to dotter */
                                             BLXSEQ_DNA,	  /* calculated dotter coords are always nucleotide coords */
                                             BLXSEQ_DNA,          /* required sequence is in nucleotide coords */
                                             frame,
					     bc->numFrames,
					     &bc->refSeqRange,
					     bc->blastMode,
					     bc->geneticCode,
                                             FALSE,		  /* input coords are always left-to-right, even if display reversed */
                                             FALSE,               /* always pass forward strand to dotter */
                                             FALSE,               /* always pass forward strand to dotter */
                                             &tmpError);
  
  if (!querySeqSegment)
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
    }
  
  /* Get the match sequence name (chopping off the letters before the colon, if there is one). */
  const char *dotterSName = strchr(mspGetSName(firstMsp), ':');
  if (dotterSName)
    {
      dotterSName++;
    }
  else
    {
      dotterSName = mspGetSName(firstMsp);
    }
  
  const int offset = dotterRange.min - 1;
  
  /* Get the options */
//  static char opts[] = "     ";
//  opts[0] = bc->displayRev ? 'R' : ' ';
//  opts[1] = hspsOnly ? 'H' : ' ';
//  opts[2] = bc->gappedHsp ? 'G' : ' ';
  
  /* Get the mode */
//  char type = getDotterMode(bc->blastMode);
  
  /* Get the list of all MSPs */
  printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);
  
  printf("  query sequence: name -  %s, offset - %d\n"
	 "subject sequence: name -  %s, offset - %d\n", bc->refSeqName, offset, dotterSName, 0);

  return callDotterExternal(bc, dotterZoom, &dotterRange, querySeqSegment, dotterSName, dotterSSeq, hspsOnly, NULL, error);
  
//  dotter(type, opts, bc->refSeqName, querySeqSegment, offset, qStrand, dotterSName, dotterSSeq, 0,
//	 selectedSeq->strand, 0, 0, NULL, NULL, NULL, 0.0, dotterZoom, bc->mspList, bc->matchSeqs, 0, 0, 0);
  
//  return TRUE;
}


/* Call dotter on the reference sequence versus itself.
 *
 * The follow notes on this are from http://sonnhammer.sbc.su.se/Dotter.html:
 *
 * When looking for overlaps between many sequences, for instance when assembling contigs, it can
 * be uselful to make a dotplot of all sequences vs. each other. This way any overlap will show up 
 * as a diagonal in the corner of a subsequence dotplot. Dotter has a built-in mechanism for this. 
 * To run Dotter on many sequences at once, concatenate the sequence files (in fasta format). Then
 * run dotter on the concatenated sequence file against itself, and green partitioning lines will
 * appear between the sequences. At each partitioning line, the name of the following sequence is 
 * printed. These lines can be turned on and off with the button "Draw lines a segment ends" in
 * the "Feature series selection tool", which is launched from the main menu. 
 */
static gboolean callDotterSelf(GtkWidget *blxWindow, GError **error)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the auto range, if requested */
  int dotterStart = UNSET_INT;
  int dotterEnd = UNSET_INT;
  int dotterZoom = 0;
  
  GError *tmpError = NULL;
  if (!getDotterRange(blxWindow, NULL, TRUE, bc->autoDotter, &dotterStart, &dotterEnd, &dotterZoom, &tmpError))
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  
  /* Set the options string (no options are required) */
//  static char opts[] = "     ";

  /* Set the type */
  char type = ' ';
  
  if (bc->blastMode == BLXMODE_BLASTP || bc->blastMode == BLXMODE_TBLASTN || bc->blastMode == BLXMODE_TBLASTX)
    {
      type = 'P';
    }
  else if (bc->blastMode == BLXMODE_BLASTX || bc->blastMode == BLXMODE_BLASTN)
    {
      type = 'N';
    }

  /* Get the section of reference sequence that we're interested in */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);
  const int frame = 1;
  IntRange dotterRange;
  intrangeSetValues(&dotterRange, dotterStart, dotterEnd);
    
  char *querySeqSegment = getSequenceSegment(bc->refSeq,
                                             &dotterRange,
                                             qStrand,
                                             BLXSEQ_DNA,	  /* calculated dotter coords are always in nucleotide coords */
                                             BLXSEQ_DNA,      /* required sequence is in nucleotide coords */
                                             frame,
                                             bc->numFrames,
                                             &bc->refSeqRange,
                                             bc->blastMode,
                                             bc->geneticCode,
                                             FALSE,		  /* input coords are always left-to-right, even if display reversed */
                                             bc->displayRev,  /* whether to reverse */
                                             bc->displayRev,  /* whether to allow rev strands to be complemented */
                                             &tmpError);
  
  if (!querySeqSegment)
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
    }
  
  /* Make a copy of the reference sequence segment to pass as the match sequence */
  char *dotterSSeq = messalloc(strlen(querySeqSegment) + 1);
  strcpy(dotterSSeq, querySeqSegment);

//  int offset = dotterRange.min - 1;

  printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);

  callDotterExternal(bc, dotterZoom, &dotterRange, querySeqSegment, bc->refSeqName, dotterSSeq, FALSE, NULL, error);
  
//  dotter(type, opts, bc->refSeqName, querySeqSegment, offset, BLXSTRAND_FORWARD, bc->refSeqName, dotterSSeq, offset,
//	 BLXSTRAND_FORWARD, 0, 0, NULL, NULL, NULL, 0.0, dotterZoom, bc->mspList, bc->matchSeqs, 0, 0, 0);
  
  return TRUE;
}



