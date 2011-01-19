/*  File: blxdotter.c
 *  Author: Gemma Barson, 2010-02-03
 *  Copyright (c) 2010 Genome Research Ltd
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
 * Description: See blxdotter.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/blxdotter.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/detailview.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

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
  BLX_DOTTER_ERROR_NO_EXE,            /* no dotter executable found */
  BLX_DOTTER_ERROR_PIPE,	      /* error opening pipe */
  BLX_DOTTER_ERROR_FORK		      /* error forking process */
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
  
  GString *resultStr = g_string_new("Blixem - Dotter sequence: ");
  
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

  GError *tmpError = NULL;
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
      /* Calculate automatic coords */
      if (callOnSelf)
	success = smartDotterRangeSelf(blxWindow, dotterStart, dotterEnd, &tmpError);
      else
	success = smartDotterRange(blxWindow, dotterSSeq, dotterStart, dotterEnd, &tmpError);
    }

  if (success && !tmpError)
    {
      /* Check that there are valid MSPs within the dotter range. Set a warning if not. */
      gboolean found = FALSE;
      const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
      const int qMin = min(*dotterStart, *dotterEnd);
      const int qMax = max(*dotterStart, *dotterEnd);

      GList *seqItem = bc->selectedSeqs;
      for ( ; seqItem; seqItem = seqItem->next)
        {
          const BlxSequence *seq = (const BlxSequence*)(seqItem->data);
          GList *mspItem = seq->mspList;  
          
          for ( ; mspItem ; mspItem = mspItem->next)
            {
              const MSP *msp = (MSP*)(mspItem->data);
              
              if ((msp->qStrand == activeStrand || bc->blastMode == BLXMODE_BLASTN) &&
                  msp->qRange.min <= qMax && msp->qRange.max >= qMin)
                {
                  found = TRUE;
                  break;
                }
            }
        }
      
      if (!found)
        {
          g_set_error(&tmpError, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_MATCHES,
            "There were no matches for the selected sequence(s) within the selected dotter range.\nThis may be because the range was trimmed - try zooming in on the region of interest.\n");
        }
    }

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
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
 * Returns true but sets the error if there is a warning.
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
      
      IntRange mspDisplayRange;
      intrangeSetValues(&mspDisplayRange, coord1, coord2);

      /* Check if the MSP is in a visible strand is entirely within the big picture range */
      if ((msp->qStrand == activeStrand || (bc->blastMode == BLXMODE_BLASTN)) &&
	  (mspDisplayRange.min >= bigPicRange->min && mspDisplayRange.max <= bigPicRange->max))
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

  if (qMin == UNSET_INT && qMax == UNSET_INT)
    {
      /* No alignments found. Give a warning, and use the big picture range. */
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_MATCHES, 
                  "There were no matches for the selected sequence(s) within the big picture range.\nZoom out to ensure alignments lie entirely within the big picture range.");
      
      qMin = bigPicRange->min;
      qMax = bigPicRange->max;
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

  /* Keep it within bounds. */
  boundsLimitValue(&qMin, &bc->refSeqRange);
  boundsLimitValue(&qMax, &bc->refSeqRange);

  /* Apply min and max limits:  min 500 residues, max 10 Mb dots. */
  int qLen = qMax - qMin;
  int midCoord = qMin + qLen/2;

  if (qLen < 500)
    {
      qLen = 500;
    }

  const int sLen = blxSequenceGetLength(selectedSeq);
  if (qLen * sLen > 1e7)
    {
      qLen = 1e7 / sLen;
    }

  qMin = midCoord - (qLen / 2) ;
  qMax = midCoord + (qLen / 2) ;

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
      fetchstr = blxprintf("%s --client=%s_%s_%s -q '%s' &",
			   fetch_prog, g_get_prgname(), g_get_host_name(), g_get_user_name(), seqname) ;
    }
  else
    {
      fetchstr = blxprintf("%s -q '%s'", fetch_prog, seqname) ;
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
      g_free(fetchstr);
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
static gboolean findCommand (char *command, char **resultOut)
{
  gboolean found = FALSE;

#if !defined(NO_POPEN)
  static char result[1025];
  char fileName[1025];
  FILE *file = NULL;
    
  char *pathEnv = g_malloc(strlen(getenv("PATH"))+1);
  strcpy(pathEnv, getenv("PATH"));

  /* Try each path in the enviroment var */
  char *path = strtok(pathEnv, ":");
  
  while (path) 
    {
      strcpy(fileName, path);
      strcat(fileName,"/");
      strcat(fileName, command);

      /* Check if the file exists in this path */
      file = fopen(fileName, "r");
      if (file)  //!access(fileName, F_OK) && !access(fileName, X_OK)) 
	{
	  found = TRUE;
	  fclose(file);
	  break;
	}
      
      path = strtok(0, ":");
    }
  
  g_free(pathEnv);
  
  if (found) 
    {
      strcpy(result, fileName);
      found = TRUE;
    }
  else 
    {
      strcpy(result, "Can't find executable 'dotter' in path.\n");
    }
    
  if (resultOut)
    {
      *resultOut = result;
    }
#else
  strcpy(result, "Can't open executable 'dotter' - popen command is not defined.\n");
#endif

  return found;
}


/* This actually executes the dotter child process */
static void callDotterChildProcess(char *dotterBinary, 
				   const int dotterZoom,
				   const int qstart,
				   const int sstart,
				   const int qlen,
				   const int slen,
				   const gboolean hspsOnly,
				   const BlxStrand sStrand,
				   int *pipes, 
				   BlxViewContext *bc,
				   const char *dotterSName,
				   char *Xoptions)
{
  DEBUG_OUT("callDotterChildProcess\n");

  /* Create the argument list */
  GSList *argList = NULL;
  argList = g_slist_append(argList, g_strdup(dotterBinary));
  argList = g_slist_append(argList, g_strdup("-z"));
  argList = g_slist_append(argList, convertIntToString(dotterZoom));
  argList = g_slist_append(argList, g_strdup("-q"));
  argList = g_slist_append(argList, convertIntToString(qstart));
  argList = g_slist_append(argList, g_strdup("-s"));
  argList = g_slist_append(argList, convertIntToString(sstart));
  
  if (bc->displayRev)		    argList = g_slist_append(argList, g_strdup("-r"));
  if (sStrand == BLXSTRAND_REVERSE) argList = g_slist_append(argList, g_strdup("-v"));
  if (hspsOnly)			    argList = g_slist_append(argList, g_strdup("-H"));
  
  argList = g_slist_append(argList, g_strdup("-S"));
  argList = g_slist_append(argList, g_strdup(bc->refSeqName));
  argList = g_slist_append(argList, convertIntToString(qlen));
  argList = g_slist_append(argList, g_strdup(dotterSName));
  argList = g_slist_append(argList, convertIntToString(slen));
  argList = g_slist_append(argList, g_strdup(dotterBinary));
  argList = g_slist_append(argList, g_strdup(Xoptions));

  if (Xoptions)
    argList = g_slist_append(argList, NULL); /* add null on end, if Xoptions wasn't already null */

  /* Convert the list to an array */
  DEBUG_OUT("Args = ");
  char *args[g_slist_length(argList)];
  GSList *item = argList;
  int i = 0;
  
  for ( ; item; item = item->next)
    {
      char *arg = (char*)(item->data);
      args[i] = arg;
      ++i;
      DEBUG_OUT(", %s", arg  );
    }
  DEBUG_OUT("\n");
    
  DEBUG_OUT("Executing dotter\n");
  dup2(pipes[0], 0);
  close(pipes[0]);
  execv(args[0], args);
  
  exit(1);
}


/* Call dotter as an external process */
gboolean callDotterExternal(BlxViewContext *bc,
                            int dotterZoom, 
                            IntRange *refSeqRange,
                            char *refSeqSegment,
                            const char *dotterSName,
                            char *dotterSSeq,
			    const BlxStrand sStrand,
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
      g_debug("Looking for Dotter ...\n");
      
      if (!findCommand("dotter", &(dotterBinary))) 
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_EXE, "No dotter executable found in path.\n$PATH=%s\n", getenv("PATH"));
          dotterBinary = NULL;
          return FALSE;
        }
    }
  
  g_debug("Calling %s with region: %d,%d - %d,%d\n", dotterBinary, xstart, ystart, xend, yend);
  fflush(stdout);

  /* Pass on the features via pipes. */
  int pipes[2];
  if (pipe (pipes))
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_PIPE, "Error opening pipe to Dotter.\n");
      return FALSE;
    }

  /* Create the child process */
  pid_t pid = fork();
  bc->spawnedProcesses = g_slist_append(bc->spawnedProcesses, GINT_TO_POINTER(pid));
  
  if (pid < 0)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_FORK, "Error forking process for Dotter.\n");
      return FALSE;
    }
  else if (pid == 0)
    {
      /* Child process. Execute dotter. First, close the write end of the pipe */
      close(pipes[1]);

      DEBUG_OUT("Calling dotter child process\n");
      callDotterChildProcess(dotterBinary, dotterZoom, xstart - 1, ystart - 1, qlen, slen, 
			     hspsOnly, sStrand, pipes, bc, dotterSName, Xoptions);
    }
  else
    {
      /* Parent process. Pipe sequences and MSPs to dotter. First, close the read end of the pipe. */
      close(pipes[0]);

      g_debug("Spawned process %d\n", pid);
      DEBUG_OUT("Piping sequences to dotter...\n");
      
      fflush(stdout);
      FILE *pipe = fdopen(pipes[1], "w");

      /* Pass the sequences */
      fwrite(refSeqSegment, 1, qlen, pipe);
      fwrite(dotterSSeq, 1, slen, pipe);
      
      DEBUG_OUT("...done\n");
      
      /* Pass the features */
      DEBUG_OUT("Piping features to dotter...\n");
      GList *seqItem = bc->matchSeqs;
      for ( ; seqItem; seqItem = seqItem->next) 
	{
	  BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
	  writeBlxSequenceToOutput(pipe, blxSeq, refSeqRange);
	}
      
      fprintf(pipe, "%c\n", EOF);
      fflush(pipe);
      
      DEBUG_OUT("...done\n");
    }
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
  GError *rangeError = NULL;
  
  gboolean ok = getDotterRange(blxWindow, dotterSSeq, FALSE, bc->autoDotter, &dotterStart, &dotterEnd, &dotterZoom, &rangeError);

  if (!ok)
    {
      prefixError(rangeError, "Error calculating dotter range. ");
      g_propagate_error(error, rangeError);
      return FALSE;
    }
  else if (ok && rangeError)
    {
      /* There was a warning when calculating the range. Ask the user if they want to continue. */
      prefixError(rangeError, "Warning: ");
      postfixError(rangeError, "\nContinue?");

      ok = (runConfirmationBox(blxWindow, "Blixem - Warning", rangeError->message) == GTK_RESPONSE_ACCEPT);
      g_error_free(rangeError);
      rangeError = NULL;
      
      if (!ok)
	return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const int frame = mspGetRefFrame(firstMsp, bc->seqType);
  IntRange dotterRange;
  intrangeSetValues(&dotterRange, dotterStart, dotterEnd);
  GError *seqError = NULL;

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
                                             &seqError);
  
  if (!querySeqSegment)
    {
      g_propagate_error(error, seqError);
      return FALSE;
    }
  else
    {
      /* If there was an error set but the sequence was still returned then it's a non-critical warning */
      reportAndClearIfError(&seqError, G_LOG_LEVEL_WARNING);
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
  
  /* Get the list of all MSPs */
  g_message("Calling dotter on match '%s' with reference sequence region: %d -> %d\n", dotterSName, dotterStart, dotterEnd);
  
  g_debug("reference sequence: name =  %s, offset = %d\n"
          "    match sequence: name =  %s, offset = %d\n", 
          bc->refSeqName, offset, dotterSName, 0);

  return callDotterExternal(bc, dotterZoom, &dotterRange, querySeqSegment, dotterSName, dotterSSeq, selectedSeq->strand, hspsOnly, NULL, error);
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
  char *dotterSSeq = g_malloc(strlen(querySeqSegment) + 1);
  strcpy(dotterSSeq, querySeqSegment);

  g_message("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);

  callDotterExternal(bc, dotterZoom, &dotterRange, querySeqSegment, bc->refSeqName, dotterSSeq, BLXSTRAND_FORWARD, FALSE, NULL, error);

  return TRUE;
}

