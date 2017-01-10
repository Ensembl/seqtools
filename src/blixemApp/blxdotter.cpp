/*  File: blxdotter.c
 *  Author: Gemma Barson, 2010-02-03
 *  Copyright (c) 2006-2017 Genome Research Ltd
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

#include <blixemApp/blxdotter.hpp>
#include <blixemApp/blxcontext.hpp>
#include <blixemApp/blxwindow.hpp>
#include <blixemApp/bigpicture.hpp>
#include <blixemApp/detailview.hpp>
#include <seqtoolsUtils/utilities.hpp>
#include <seqtoolsUtils/blxmsp.hpp>
#include <gbtools/gbtools.hpp>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>

using namespace std;


#define DEFAULT_DOTTER_RANGE_SELF	2000
#define DEFAULT_WINDOW_WIDTH            500
#define MAX_WINDOW_WIDTH_FRACTION       0.3

typedef struct _DotterDialogData
  {
    GtkWidget *blxWindow;           /* pointer to the main blixem window */
    GtkWidget *dialog;
    GtkWidget *notebook;
    guint notebookPage;             /* current notebook page */

    GtkWidget *autoButton;          /* radio button for automatic dotter parameters */
    GtkWidget *manualButton;        /* radio button for manual dotter parameters */
    GtkWidget *transcriptButton;    /* radio button for dottering against transcript rather than
                                     * reference sequence range */

    GtkWidget *startEntry;          /* the text entry box on the dialog for the start coord */
    GtkWidget *endEntry;            /* the text entry box on the dialog for the end coord */
    GtkWidget *zoomEntry;           /* the text entry box on the dialog for the zoom value */
    
    DotterMatchType matchType;      /* whether to call dotter on the selected match, an adhoc seq,
                                     * or the query seq versus itself */
    DotterRefType refType;          /* whether to use the ref seq or a transcript seq */
    gboolean hspsOnly;              /* whether to call dotter on HSPs only */
    gboolean sleep;                 /* whether to sleep dotter on startup */
    
    GtkWidget *selectedButton;      /* the radio button for the use-selected-sequence option */
    GtkWidget *adhocButton;         /* the radio button for the use-adhoc-sequence option */
    GtkWidget *selfButton;          /* the radio button for the use-self option */

    GtkWidget *adhocSeqText;       /* text entry buffer for pasting sequence to dotter against into */
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
  BLX_DOTTER_ERROR_FORK,	      /* error forking process */
  BLX_DOTTER_ERROR_NEGATED_COORD,     /* warning that an original coord was not in range so we negated it */
  BLX_DOTTER_ERROR_OUT_OF_RANGE,      /* dotter coord was out of range */
  BLX_DOTTER_ERROR_NO_TRANSCRIPT      /* tried to dotter vs a transcript but none (or multiple) selected */
} BlxDotterError;


/* Local function declarations */
static gboolean	      getDotterRange(GtkWidget *blxWindow, DotterMatchType matchType, const DotterRefType refType, int *dotterStart, int *dotterEnd, int *dotterZoom, GError **error);
static gboolean	      smartDotterRange(GtkWidget *blxWindow, int *dotter_start_out, int *dotter_end_out, GError **error);
static gboolean	      callDotterOnSelf(DotterDialogData *dialogData, GError **error);
static gboolean	      callDotterOnAdhocSeq(DotterDialogData *dialogData, GError **error);
static char*          getSelectedSequenceDNA(const BlxSequence *selectedSeq, const BlxContext *bc, GError **error);
static void           textGetSeqDetails(const char *text, char **sequence, char **sequenceName);
static char*          getDotterTitle(const BlxContext *bc, const DotterMatchType matchType, const DotterRefType refType, const char *adhocSeq); 
static char*          getDotterTitleAdhocSeq(const BlxContext *bc, const char *adhocSeq, const DotterRefType refType);
static const char*    getDotterRefSeqName(const BlxContext *bc, const gboolean transcript);
static gboolean       boundsCheckDotterCoord(int *coordIn, BlxContext *bc, GError **error);
static int            getDisplayCoord(const int coordIn, BlxContext *bc);
static char*          textViewGetText(GtkWidget *textView);


/*******************************************************************
 *                       Response callbacks                        *
 *******************************************************************/

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
            switch (dialogData->matchType)
              {
                default: /* fall through */
                case BLXDOTTER_MATCH_SELECTED:
                  {
                    /* The 'transcript' flag indicates that we're dottering vs the selected
                     * transcript rather than the reference sequence*/
                    destroy = callDotterOnSelectedSeqs(dialogData->blxWindow, 
                                                       dialogData->hspsOnly, 
                                                       dialogData->sleep,
                                                       dialogData->refType,
                                                       &error);
                    break;
                  }
                case BLXDOTTER_MATCH_ADHOC:
                  destroy = callDotterOnAdhocSeq(dialogData, &error);
                  break;
                case BLXDOTTER_MATCH_SELF:
                  destroy = callDotterOnSelf(dialogData, &error);
                  break;
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


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the dotter mode according to the toggle state of the 
 * "auto" button. */
static gboolean onSaveDotterMode(GtkWidget *button, const gint responseId, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  GtkWidget *blxWindow = dialogData->blxWindow;
  BlxContext *bc = blxWindowGetContext(blxWindow);
  
  bc->dotterRefType = dialogData->refType;

  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the sequence we should dotter against to be the selected
 * sequence if the button is active. */
static gboolean onSaveDotterSelectedSeq(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxContext *blxContext = blxWindowGetContext(blxWindow);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    blxContext->dotterMatchType = BLXDOTTER_MATCH_SELECTED;

  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the sequence we should dotter against to be the adhoc
 * sequence if the button is active. */
static gboolean onSaveDotterAdhoc(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxContext *blxContext = blxWindowGetContext(blxWindow);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    blxContext->dotterMatchType = BLXDOTTER_MATCH_ADHOC;

  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the adhoc sequence text. */
static gboolean onSaveDotterAdhocSeq(GtkWidget *textView, const gint responseId, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *blxContext = blxWindowGetContext(dialogData->blxWindow);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton)))
    {
      if (blxContext->dotterAdhocSeq)
        {
          g_free(blxContext->dotterAdhocSeq);
          blxContext->dotterAdhocSeq = NULL;
        }

      GtkTextBuffer *textBuffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textView));

      GtkTextIter start, end;
      gtk_text_buffer_get_bounds(textBuffer, &start, &end);

      blxContext->dotterAdhocSeq = gtk_text_buffer_get_text(textBuffer, &start, &end, TRUE);
    }

  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the sequence we should dotter against to be the 
 * reference sequence if the button is active. */
static gboolean onSaveDotterSelf(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxContext *blxContext = blxWindowGetContext(blxWindow);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    blxContext->dotterMatchType = BLXDOTTER_MATCH_SELF;

  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the state of the 'HSPs only' button. */
static gboolean onSaveDotterHsps(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->dotterHsps = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the state of the 'sleep' button. */
static gboolean onSaveDotterSleep(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->dotterSleep = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the start parameter that was entered (if manual dotter
 * params are being used). */
static gboolean onSaveDotterStart(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;

  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);

  /* Only save the parameter if we are using manual parameters */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->manualButton)))
    {
      int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      GError *error = NULL;
      
      result = boundsCheckDotterCoord(&newVal, bc, &error);
      
      if (result)
        {
          bc->dotterStart = newVal;
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }
      else
        {
          postfixError(error, "Value not saved.");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterEnd(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->manualButton)))
    {
      int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));

      GError *error = NULL;
      
      result = boundsCheckDotterCoord(&newVal, bc, &error);
      
      if (result)
        {
          bc->dotterEnd = newVal;
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }
      else
        {
          postfixError(error, "Value not saved.");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterZoom(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *blxContext = blxWindowGetContext(dialogData->blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->manualButton)))
    {
      const int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      
      if (newVal < 0)
        {
          result = FALSE;
          g_critical("Zoom must be greater than zero (or equal to zero to use auto zoom)");
        }
      else
        {
          blxContext->dotterZoom = newVal;
        }
    }  
  
  return result;
}


/*******************************************************************
 *                            Events                               *
 *******************************************************************/

/* Called when the user changes the page on the notebook. Caches the page so we can reinstate it
 * when we refresh the dialog */
static gboolean onChangeCurrentPage(GtkNotebook *notebook, GtkWidget *page, gint pageNum, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;

  dialogData->notebookPage = pageNum;

  return TRUE;
}


/* Called when the 'last saved' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end coords that were last saved */
static void onLastSavedButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  char *startString = convertIntToString(getDisplayCoord(bc->dotterStart, bc));
  char *endString = convertIntToString(getDisplayCoord(bc->dotterEnd, bc));
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
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  const int startCoord = (bc->displayRev ? bc->refSeqRange.max() : bc->refSeqRange.min());
  const int endCoord = (bc->displayRev ? bc->refSeqRange.min() : bc->refSeqRange.max());
  
  char *startString = convertIntToString(getDisplayCoord(startCoord, bc));
  char *endString = convertIntToString(getDisplayCoord(endCoord, bc));
  
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
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(dialogData->blxWindow);
  const IntRange* const displayRange = bigPictureGetDisplayRange(bigPicture);

  int qStart = convertDisplayIdxToDnaIdx(displayRange->min(), bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  int qEnd = convertDisplayIdxToDnaIdx(displayRange->max(), bc->seqType, bc->numFrames, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange);

  boundsLimitValue(&qStart, &bc->refSeqRange);
  boundsLimitValue(&qEnd, &bc->refSeqRange);
  
  char *startString = convertIntToString(getDisplayCoord(qStart, bc));
  char *endString = convertIntToString(getDisplayCoord(qEnd, bc));
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the user enters in the adhoc-sequence text entry box. Sets the adhoc-sequence
 * toggle button to be the active one and updates the dialog title. */
static gboolean onAdhocSeqEntered(GtkWidget *widget, GdkEvent *event, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);

  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton), TRUE);

  char *adhocSeq = textViewGetText(dialogData->adhocSeqText);
  char *title = getDotterTitleAdhocSeq(bc, adhocSeq, dialogData->refType);

  gtk_window_set_title(GTK_WINDOW(dialogData->dialog), title);

  g_free(title);
  g_free(adhocSeq);

  /* Need to force a redraw of the dialog to get the new title */
  /* gb10: queue_draw doesn't work, not sure how we would do this */
  //gtk_widget_queue_draw(dialogData->dialog);
  //
  //while (gtk_events_pending())
  //  gtk_main_iteration();

  return FALSE; /* allow other handlers to continue */
}


/* Called when the auto/manual/transcript radio button is toggled to select where to get the
 * ref seq coords from. */
static void onRefTypeToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  char *title = getDotterTitle(bc, dialogData->matchType, dialogData->refType, bc->dotterAdhocSeq);
  gtk_window_set_title(GTK_WINDOW(dialogData->dialog), title);
  g_free(title);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->manualButton)))
    {
      dialogData->refType = BLXDOTTER_REF_MANUAL;

      /* Manual coords. Leave values as they are but unlock the boxes so they can be edited. */
      gtk_widget_set_sensitive(dialogData->startEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->endEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->zoomEntry, TRUE);
    }
  else
    {
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->transcriptButton)))
        dialogData->refType = BLXDOTTER_REF_TRANSCRIPT;
      else
        dialogData->refType = BLXDOTTER_REF_AUTO;

      /* Recalculate auto start/end in case user has selected a different sequence */
      int autoStart = UNSET_INT, autoEnd = UNSET_INT;
      getDotterRange(dialogData->blxWindow, dialogData->matchType, dialogData->refType,
                     &autoStart, &autoEnd, NULL, NULL);

      if (autoStart == UNSET_INT && autoEnd == UNSET_INT)
        {
          autoStart = bc->displayRev ? bc->refSeqRange.max() : bc->refSeqRange.min();
          autoEnd = bc->displayRev ? bc->refSeqRange.min() : bc->refSeqRange.max();
        }
      
      char *startString = convertIntToString(getDisplayCoord(autoStart, bc));
      char *endString = convertIntToString(getDisplayCoord(autoEnd, bc));
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
}


/* Called when the user toggles what type of match sequence to dotter against */
static void onMatchTypeToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->selfButton)))
    {
      dialogData->matchType = BLXDOTTER_MATCH_SELF;
    }
  else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton)))
    {
      dialogData->matchType = BLXDOTTER_MATCH_ADHOC;
      gtk_widget_grab_focus(dialogData->adhocSeqText);
    }
  else
    {
      dialogData->matchType = BLXDOTTER_MATCH_SELECTED;
    }

  /* Update the title of the dialog box to reflect the sequence we're dottering */
  char *adhocSeq = textViewGetText(dialogData->adhocSeqText);
  char *title = getDotterTitle(bc, dialogData->matchType, dialogData->refType, adhocSeq);

  gtk_window_set_title(GTK_WINDOW(dialogData->dialog), title);

  g_free(title);
  g_free(adhocSeq);

  /* If using auto coords, recalculate them */ 
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton)))
    {
      int autoStart = UNSET_INT, autoEnd = UNSET_INT;
      getDotterRange(dialogData->blxWindow, dialogData->matchType, dialogData->refType, &autoStart, &autoEnd, NULL, NULL);

      if (autoStart == UNSET_INT)
        autoStart = bc->displayRev ? bc->refSeqRange.max() : bc->refSeqRange.min();
  
      if (autoEnd == UNSET_INT)
        autoEnd = bc->displayRev ? bc->refSeqRange.min() : bc->refSeqRange.max();
  
      char *startString = convertIntToString(getDisplayCoord(autoStart, bc));
      char *endString = convertIntToString(getDisplayCoord(autoEnd, bc));
      char *zoomString = convertIntToString(bc->dotterZoom);

      /* Display the new values */
      gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
  
      g_free(startString);
      g_free(endString);
      g_free(zoomString);
    }
}


/* Called when the "call on self" button is toggled on the dotter dialog */
static void onHspsButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  dialogData->hspsOnly = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
}


/* Called when the "sleep" button is toggled on the dotter dialog */
static void onSleepButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  dialogData->sleep = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
}


static void onDestroyDotterDialog(GtkWidget *dialog, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  
  if (dialogData)
    {
      delete dialogData;
    }
}


/*******************************************************************
 *                      Creating the dialog                        *
 *******************************************************************/

static GtkWidget* createTextEntry(GtkTable *table, 
				  int col, 
				  int row, 
				  const int xpad, 
				  const int ypad, 
				  const char *title,
				  BlxResponseCallback callbackFunc,
                                  gpointer callbackData,
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
  widgetSetCallbackData(entry, callbackFunc, callbackData);
  
  return entry;
}


/* Set the initial size of the dotter window. We add a bit of width if possible because the title
 * can be quite long and it updates dynamically (and may get longer) with user-selected options. */
static void dotterDialogSetDefaultSize(GtkWidget *dialog, GtkWidget *blxWindow)
{
  /* We'll set the default window width to a preferred default, but make it smaller
   * if it exceeds a maximum percentage of the screen size. Just use the height 
   * allocated for the widgets. */
  int maxWidth = 300;
  gbtools::GUIGetTrueMonitorSizeFraction(dialog, MAX_WINDOW_WIDTH_FRACTION, 1, &maxWidth, NULL);
  
  const int width = min(DEFAULT_WINDOW_WIDTH, maxWidth);

  gtk_window_set_default_size(GTK_WINDOW(dialog), width, -1);
}


static void createCoordsTab(DotterDialogData *dialogData, const int spacing)
{
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);

  /* Put everything in a table */
  int numRows = 5;
  int numCols = 3;
  int row = 0;
  int col = 0;
  int xpad = spacing;
  int ypad = spacing;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));

  /* Append the table as a new tab to the notebook */
  gtk_notebook_append_page(GTK_NOTEBOOK(dialogData->notebook), GTK_WIDGET(table), gtk_label_new("Reference sequence"));

  /* Radio buttons for auto/manual params. */
  dialogData->autoButton = gtk_radio_button_new_with_mnemonic(NULL, "_Auto");
  widgetSetCallbackData(dialogData->autoButton, onSaveDotterMode, dialogData);
  gtk_widget_set_tooltip_text(dialogData->autoButton, "Dotter against the visible big picture range");
  gtk_table_attach(table, dialogData->autoButton, 0, 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  dialogData->manualButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(dialogData->autoButton), "_Manual");
  widgetSetCallbackData(dialogData->manualButton, onSaveDotterMode, dialogData);
  gtk_widget_set_tooltip_text(dialogData->manualButton, "Dotter against manually-entered coords. Once saved, the same coords will be used until you change them manually again.");
  gtk_table_attach(table, dialogData->manualButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  dialogData->transcriptButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(dialogData->autoButton), "_Transcript");
  widgetSetCallbackData(dialogData->transcriptButton, onSaveDotterMode, dialogData);
  gtk_widget_set_tooltip_text(dialogData->transcriptButton, "Dotter against the currently selected transcript Note that if you are dottering against a selected sequence then you need to select both the sequence and the transcript. Hold Ctrl while you left-click to select multiple features.");
  gtk_table_attach(table, dialogData->transcriptButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  /* Buttons that the user can click to populate the parameter boxes with certain values */
  GtkWidget *lastSavedButton = gtk_button_new_with_mnemonic("_Last saved ->");
  gtk_widget_set_tooltip_text(lastSavedButton, "Put the last-saved coords in the manual coords boxes");
  gtk_table_attach(table, lastSavedButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  GtkWidget *fullRangeButton = gtk_button_new_with_mnemonic("_Full range ->");
  gtk_widget_set_tooltip_text(fullRangeButton, "Put the full range in the manual coords boxes");
  gtk_table_attach(table, fullRangeButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  GtkWidget *bpRangeButton = gtk_button_new_with_mnemonic("_Big picture range ->");
  gtk_widget_set_tooltip_text(bpRangeButton, "Put the big picture range in the manual coords boxes");
  gtk_table_attach(table, bpRangeButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  /* Disable last-saved button if no saved values exist */
  if (bc->dotterStart == UNSET_INT)
    {
      gtk_widget_set_sensitive(lastSavedButton, FALSE);
    }

  /* Create the text entry boxes for the dotter parameters */
  ++col;
  row = 0;

  dialogData->startEntry = createTextEntry(table, col, row, xpad, ypad, "<i>Start:</i>", onSaveDotterStart, 
                                           dialogData, getDisplayCoord(bc->dotterStart, bc));
  ++row;
  dialogData->endEntry = createTextEntry(table, col, row, xpad, ypad, "<i>End:</i>", onSaveDotterEnd, 
                                         dialogData, getDisplayCoord(bc->dotterEnd, bc));
  ++row;
  dialogData->zoomEntry = createTextEntry(table, col, row, xpad, ypad, "<i>Zoom:</i>", onSaveDotterZoom, 
                                          dialogData, bc->dotterZoom);
  ++row;

  gtk_widget_set_tooltip_text(dialogData->zoomEntry, "The level of zoom to open Dotter with (higher values zoom in)");

  /* There is an issue if the user selects a different sequence while the dotter dialog
   * is still open: the auto range does not update automatically for the new sequence. To 
   * mitigate this, connect the 'clicked' signal so that they can
   * click on the 'auto' toggle button and have it refresh, even if that button is already selected.*/
  g_signal_connect(G_OBJECT(dialogData->autoButton),      "clicked", G_CALLBACK(onRefTypeToggled), dialogData);
  g_signal_connect(G_OBJECT(dialogData->manualButton),    "clicked", G_CALLBACK(onRefTypeToggled), dialogData);
  g_signal_connect(G_OBJECT(dialogData->transcriptButton),"clicked", G_CALLBACK(onRefTypeToggled), dialogData);

  g_signal_connect(G_OBJECT(lastSavedButton), "clicked", G_CALLBACK(onLastSavedButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(fullRangeButton), "clicked", G_CALLBACK(onFullRangeButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(bpRangeButton),   "clicked", G_CALLBACK(onBpRangeButtonClicked), dialogData);

  /* Set the initial state of the toggle buttons and entry widgets */
  switch (dialogData->refType)
    {
      case BLXDOTTER_REF_AUTO:
        /* unset and then set to force callback to be called to set the auto coords */
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->transcriptButton), TRUE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->autoButton), TRUE);
        break;
      case BLXDOTTER_REF_MANUAL:
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
        break;
      case BLXDOTTER_REF_TRANSCRIPT:
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->transcriptButton), TRUE);
        break;
    }
}


static void createSequenceTab(DotterDialogData *dialogData, const int spacing)
{
  GtkWidget *blxWindow = dialogData->blxWindow;
  BlxContext *bc = blxWindowGetContext(blxWindow);

  /* Put everything in a table */
  int numRows = 5;
  int numCols = 3;
  int row = 0;
  int col = 0;
  int xpad = spacing;
  int ypad = spacing;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));

  /* Append the box as a new tab to the notebook */
  gtk_notebook_append_page(GTK_NOTEBOOK(dialogData->notebook), GTK_WIDGET(table), gtk_label_new("Match sequence"));

  /* Create radio buttons to choose whether to dotter against the selected match sequence, a
   * manually entered adhoc sequence, or the reference sequence against itself. */
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, spacing));
  gtk_table_attach(table, GTK_WIDGET(vbox), col, col + 1, row, row + 2, GTK_FILL, GTK_FILL, xpad, ypad);

  dialogData->selectedButton = gtk_radio_button_new_with_mnemonic(NULL, "_Selected match(es)");
  gtk_widget_set_tooltip_text(dialogData->selectedButton, "Dotter against the selected match sequence");
  widgetSetCallbackData(dialogData->selectedButton, onSaveDotterSelectedSeq, blxWindow);
  gtk_box_pack_start(vbox, dialogData->selectedButton, FALSE, FALSE, spacing);

  dialogData->selfButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(dialogData->selectedButton), "Sel_f");
  gtk_widget_set_tooltip_text(dialogData->selfButton, "Dotter the reference sequence against itself");
  widgetSetCallbackData(dialogData->selfButton, onSaveDotterSelf, blxWindow);
  gtk_box_pack_start(vbox, dialogData->selfButton, FALSE, FALSE, spacing);

  ++col;
  row = 0;
  dialogData->adhocButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(dialogData->selectedButton), "_Enter adhoc sequence");
  gtk_widget_set_tooltip_text(dialogData->adhocButton, "Dotter against a manually entered sequence. This can be in raw or FASTA format.");
  widgetSetCallbackData(dialogData->adhocButton, onSaveDotterAdhoc, blxWindow);
  gtk_table_attach(table, dialogData->adhocButton, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;

  GtkTextBuffer *textBuffer = gtk_text_buffer_new(gtk_text_tag_table_new());
  if (bc->dotterAdhocSeq)
    gtk_text_buffer_set_text(textBuffer, bc->dotterAdhocSeq, -1);

  dialogData->adhocSeqText = gtk_text_view_new_with_buffer(GTK_TEXT_BUFFER(textBuffer));
  widgetSetCallbackData(dialogData->adhocSeqText, onSaveDotterAdhocSeq, dialogData);

  const int numLines = 4;
  const int charHeight = getTextHeight(dialogData->adhocSeqText, "A");
  gtk_widget_set_size_request(dialogData->adhocSeqText, -1, roundNearest(charHeight * numLines));

  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  GtkWidget *frame = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
  gtk_container_add(GTK_CONTAINER(scrollWin), dialogData->adhocSeqText);
  gtk_container_add(GTK_CONTAINER(frame), scrollWin);

  gtk_table_attach(table, frame, col, col + 1, row, row + 1, 
                   (GtkAttachOptions)(GTK_FILL | GTK_EXPAND), (GtkAttachOptions)(GTK_FILL | GTK_EXPAND), xpad, ypad);
  ++row;


  /* Connect signals */
  g_signal_connect(G_OBJECT(dialogData->selectedButton),  "clicked", G_CALLBACK(onMatchTypeToggled), dialogData);
  g_signal_connect(G_OBJECT(dialogData->adhocButton),    "clicked", G_CALLBACK(onMatchTypeToggled), dialogData);
  g_signal_connect(G_OBJECT(dialogData->selfButton),      "clicked", G_CALLBACK(onMatchTypeToggled), dialogData);

  /* Add a callback to activate the adhocButton toggle button when the user enters some text... */
  g_signal_connect(G_OBJECT(dialogData->adhocSeqText), "focus-in-event", G_CALLBACK(onAdhocSeqEntered), dialogData);

  /* Set the initial state of the toggle buttons (unset it and then set it to
   * force the callback to be called) */
  switch (dialogData->matchType)
    {
      default: /* fall through */
      case BLXDOTTER_MATCH_SELECTED:
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selectedButton), TRUE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selfButton), FALSE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton), FALSE);
        break;
     
      case BLXDOTTER_MATCH_SELF:
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selectedButton), FALSE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selfButton), TRUE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton), FALSE);
        break;

      case BLXDOTTER_MATCH_ADHOC:
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selectedButton), FALSE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->selfButton), FALSE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->adhocButton), TRUE);
        break;
    }
}


static void createOptionsTab(DotterDialogData *dialogData, const int spacing)
{
  BlxContext *bc = blxWindowGetContext(dialogData->blxWindow);

  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 0));

  /* Append the box as a new tab to the notebook */
  gtk_notebook_append_page(GTK_NOTEBOOK(dialogData->notebook), GTK_WIDGET(vbox), gtk_label_new("Options"));

  /* Add a tick box to set hsps-only */
  GtkWidget *hspsButton = gtk_check_button_new_with_mnemonic("_HSPs only");
  gtk_widget_set_tooltip_text(hspsButton, "Start dotter showing high-scoring pairs only");
  gtk_box_pack_start(vbox, hspsButton, FALSE, FALSE, spacing);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hspsButton), bc->dotterHsps);
  widgetSetCallbackData(hspsButton, onSaveDotterHsps, dialogData->blxWindow);

  /* Add a tick box to set dotter sleep on startup (for debugging) */
  GtkWidget *sleepButton = gtk_check_button_new_with_mnemonic("_Sleep on startup");
  gtk_widget_set_tooltip_text(sleepButton, "Sleep dotter for a short period on startup (for debugging)");
  gtk_box_pack_start(vbox, sleepButton, FALSE, FALSE, spacing);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sleepButton), bc->dotterSleep);
  widgetSetCallbackData(sleepButton, onSaveDotterSleep, dialogData->blxWindow);
  
  g_signal_connect(G_OBJECT(hspsButton), "toggled", G_CALLBACK(onHspsButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(sleepButton), "toggled", G_CALLBACK(onSleepButtonToggled), dialogData);
}


/* First time round, create the dotter dialog. Subsequent calls re-use the same dialog but clear
 * the contents ready to re-populate it. If resetData is true, reset entries to last-saved values */
static GtkWidget* getOrCreateDotterDialog(GtkWidget *blxWindow, 
                                          const gboolean resetValues,
                                          DotterDialogData **dialogData_inout)
{
  BlxContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_DOTTER;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  DotterDialogData *dialogData = *dialogData_inout;
  gboolean createContent = FALSE;
  
  if (!dialog)
    {
      createContent = TRUE;

      dialog = gtk_dialog_new_with_buttons(NULL, 
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
      dotterDialogSetDefaultSize(dialog, blxWindow);

      /* Create the dialog data struct first time round, but re-populate it each time. Create 
       * a destructor function that will free the struct. */
      *dialogData_inout = new DotterDialogData;
      dialogData = *dialogData_inout;
      dialogData->blxWindow = blxWindow;
      dialogData->notebookPage = 0;
      dialogData->dialog = dialog;
      dialogData->matchType = bc->dotterMatchType;
      dialogData->refType = bc->dotterRefType;
      dialogData->hspsOnly = bc->dotterHsps;
      dialogData->sleep = bc->dotterSleep;
      dialogData->autoButton = NULL;
      dialogData->manualButton = NULL;
      dialogData->transcriptButton = NULL;
      dialogData->startEntry = NULL;
      dialogData->endEntry = NULL;
      dialogData->zoomEntry = NULL;
      dialogData->adhocSeqText = NULL;

      g_signal_connect(G_OBJECT(dialog), "destroy", G_CALLBACK(onDestroyDotterDialog), dialogData);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDotterDialog), dialogData);
    }
  else if (resetValues)
    {
      /* Refresh the dialog by clearing its contents an re-creating it */
      createContent = TRUE;
      dialogClearContentArea(GTK_DIALOG(dialog));      

      /* Reset to last-saved values */
      dialogData->matchType = bc->dotterMatchType;
      dialogData->refType = bc->dotterRefType;
      dialogData->hspsOnly = bc->dotterHsps;
      dialogData->sleep = bc->dotterSleep;
    }
  
  if (createContent)
    {
      /* Create the dialog content */
      GtkContainer *contentArea = GTK_CONTAINER(GTK_DIALOG(dialog)->vbox);
      gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
      const int spacing = 4;

      /* Create the dialog content */
      dialogData->notebook = gtk_notebook_new();
      gtk_container_add(contentArea, dialogData->notebook);
      g_signal_connect(G_OBJECT(dialogData->notebook), "switch-page", G_CALLBACK(onChangeCurrentPage), dialogData);
  
      createCoordsTab(dialogData, spacing);
      createSequenceTab(dialogData, spacing);
      createOptionsTab(dialogData, spacing);
    }

  char *adhocSeq = textViewGetText(dialogData->adhocSeqText);
  char *title = getDotterTitle(bc, dialogData->matchType, dialogData->refType, adhocSeq);

  gtk_window_set_title(GTK_WINDOW(dialog), title);

  g_free(title);
  g_free(adhocSeq);

  return dialog;
}


/* Pop up a dialog to allow the user to edit dotter parameters and launch dotter */
void showDotterDialog(GtkWidget *blxWindow, const gboolean resetValues)
{
  static DotterDialogData *dialogData = NULL;

  GtkWidget *dialog = getOrCreateDotterDialog(blxWindow, resetValues, &dialogData);
  
  gtk_widget_show_all(dialog);

  onRefTypeToggled(NULL, dialogData);

  if (resetValues)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/*******************************************************************
 *		      Title utility functions                      *
 *******************************************************************/

/* Utility to get the title for the dotter dialog. Uses the selected sequence name if a single
 * sequence is selected, or shows <no sequences> or <multiple sequences>. The result should be
 * free'd with g_free. */
static char* getDotterTitleSelectedSeq(const BlxContext *bc, const DotterRefType refType)
{
  char *result = NULL;
  GString *resultStr = g_string_new(blxGetTitlePrefix(bc));

  /* Get ref seq name */
  const gboolean transcript = (refType == BLXDOTTER_REF_TRANSCRIPT);
  const char *refSeqName = getDotterRefSeqName(bc, transcript);
  g_string_append_printf(resultStr, "Dotter %s vs ", refSeqName);

  /* Get selected match seq name */
  GList *selectedSeqs = bc->getSelectedSeqsByType(BLXSEQUENCE_MATCH);
  const int numSeqs = g_list_length(selectedSeqs);
  
  if (numSeqs == 1)
    {
      const BlxSequence *blxSeq = (const BlxSequence*)(selectedSeqs->data);
      g_string_append(resultStr, blxSequenceGetName(blxSeq));
    }
  else if (numSeqs < 1)
    {
      g_string_append(resultStr, "<no sequences selected>");
    }
  else
    {
      g_string_append_printf(resultStr, "<%d sequences>", numSeqs);
    }
  
  result = g_string_free(resultStr, FALSE);
  
  return result;
}


/* Get the dotter title for the adhoc sequence text */
static char *getDotterTitleAdhocSeq(const BlxContext *bc, const char *adhocSeq, const DotterRefType refType)
{
  char *result = NULL;
  GString *resultStr = g_string_new(blxGetTitlePrefix(bc));

  /* Get ref seq name */
  const gboolean transcript = (refType == BLXDOTTER_REF_TRANSCRIPT);
  const char *refSeqName = getDotterRefSeqName(bc, transcript);
  g_string_append_printf(resultStr, "Dotter %s vs ", refSeqName);

  if (adhocSeq)
    {
      char *seqName = NULL;
      textGetSeqDetails(adhocSeq, NULL, &seqName);

      if (seqName)
        {
          g_string_append(resultStr, seqName);
          g_free(seqName);
        }
      else
        {
          g_string_append(resultStr, "<adhoc sequence>");
        }
    }
  else
    {
      g_string_append(resultStr, "<adhoc sequence>");
    }

  result = g_string_free(resultStr, FALSE);
  return result;
}


/* Get the dotter title for the ref sequence vs itself */
static char *getDotterTitleSelf(const BlxContext *bc, const DotterRefType refType)
{
  char *result = NULL;
  GString *resultStr = g_string_new(blxGetTitlePrefix(bc));

  /* Get ref seq name */
  const gboolean transcript = (refType == BLXDOTTER_REF_TRANSCRIPT);
  const char *refSeqName = getDotterRefSeqName(bc, transcript);
  g_string_append_printf(resultStr, "%s vs %s", refSeqName, refSeqName);

  result = g_string_free(resultStr, FALSE);
  return result;
}


/* Get the title for the dotter dialog. */
static char *getDotterTitle(const BlxContext *bc, 
                            const DotterMatchType matchType, 
                            const DotterRefType refType,
                            const char *adhocSeq)
{
  char *result = NULL;

  switch (matchType)
    {
      default:
      case BLXDOTTER_MATCH_SELECTED: 
        result = getDotterTitleSelectedSeq(bc, refType);
        break;

      case BLXDOTTER_MATCH_ADHOC: 
        result = getDotterTitleAdhocSeq(bc, adhocSeq, refType);
        break;

      case BLXDOTTER_MATCH_SELF: 
        result = getDotterTitleSelf(bc, refType);
        break;
    }

  return result;
}


/*******************************************************************
 *		      Dotter utility functions                     *
 *******************************************************************/

/* Get the display coord version of the given coord (i.e. negated if the
 * 'negate coords' option is enabled and the display is reversed). Only 
 * applicable to reference sequence coords. */
static int getDisplayCoord(const int coordIn, BlxContext *bc)
{
  int result = coordIn;
  
  if (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS])
    {
      result *= -1;
    }
  
  return result;
}


/* Check whether the given dotter coord is within the blixem ref seq range.
 * The coord should be as the user sees it, i.e. negated if the 'negate coords'
 * option is enabled and the display is reversed - in this case in input arg
 * will be updated to reflect the real coord. 
 * Returns false and sets the error if not in range. Also sets the error but 
 * returns true if we succeeded but with a warning. */
static gboolean boundsCheckDotterCoord(int *coordIn, BlxContext *bc, GError **error)
{
  gboolean ok = TRUE;
  
  const gboolean negate = bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS];
  
  int coord = *coordIn;
  
  if (negate)
    {
      /* Display coords are negated - un-negate to get the real coord */
      coord *= -1;
    }
  
  if (!valueWithinRange(coord, &bc->refSeqRange))
    {
      /* Try negating it in case the user missed the minus sign off. */
      coord *= -1;
      
      if (valueWithinRange(coord, &bc->refSeqRange))
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NEGATED_COORD,
                      "Negated coord '%d' because original coord was not in range.\n", *coordIn);
        }
      else
        {
          ok = FALSE;
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_OUT_OF_RANGE,
                      "Coord '%d' is outside reference sequence range [%d -> %d].\n", 
                      *coordIn, 
                      (negate ? bc->refSeqRange.max() * -1 : bc->refSeqRange.min()),
                      (negate ? bc->refSeqRange.min() * -1 : bc->refSeqRange.max()));
        }
    }
  
  *coordIn = coord;
  
  return ok;
}


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



/* Get the start/end coords. If the passed refType flag is auto, calculate coords
 * automatically - otherwise use the stored manual coords. If refType is transcript
 * then don't set the coords (they get set later) */
static gboolean getDotterRange(GtkWidget *blxWindow, 
                               const DotterMatchType matchType,
                               const DotterRefType refType,
			       int *dotterStart, 
			       int *dotterEnd, 
			       int *dotterZoom, 
			       GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */

  GError *tmpError = NULL;
  gboolean success = TRUE;
  
  BlxContext *bc = blxWindowGetContext(blxWindow);
  

  if (refType == BLXDOTTER_REF_TRANSCRIPT)
    {
      BlxSequence *transcriptSeq = bc->getSelectedTranscript(NULL);

      if (transcriptSeq)
        {
          if (dotterStart)
            *dotterStart = transcriptSeq->qRangeFwd.min();

          if (dotterEnd)
            *dotterEnd = transcriptSeq->qRangeFwd.max();
        }
    }
  else if (refType == BLXDOTTER_REF_MANUAL)
    {
      /* Use manual coords */
      if (dotterStart) 
        *dotterStart = bc->dotterStart;

      if (dotterEnd) 
        *dotterEnd = bc->dotterEnd;

      if (dotterZoom) 
        *dotterZoom = bc->dotterZoom;
      
      if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
	{
	  g_debug("Manual dotter parameters were requested but one or more coord is not set (start=%d, end=%d). Calculating automatic parameters instead.\n",
		  *dotterStart, *dotterEnd);
	}
    }
  
  if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
    {
      /* Calculate automatic coords */
      success = smartDotterRange(blxWindow, dotterStart, dotterEnd, &tmpError);
    }

  if (success && !tmpError && matchType == BLXDOTTER_MATCH_SELECTED)
    {
      /* Check that there are valid MSPs within the dotter range. Set a warning if not. */
      gboolean found = FALSE;
      const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
      const int qMin = min(*dotterStart, *dotterEnd);
      const int qMax = max(*dotterStart, *dotterEnd);

      GList *selectedSeqs = bc->getSelectedSeqsByType(BLXSEQUENCE_MATCH);
      GList *seqItem = selectedSeqs;
      
      for ( ; seqItem; seqItem = seqItem->next)
        {
          const BlxSequence *seq = (const BlxSequence*)(seqItem->data);
          GList *mspItem = seq->mspList;  
          
          for ( ; mspItem ; mspItem = mspItem->next)
            {
              const MSP *msp = (MSP*)(mspItem->data);
              
              if ((msp->qStrand == activeStrand || bc->blastMode == BLXMODE_BLASTN) &&
                  msp->qRange.min() <= qMax && msp->qRange.max() >= qMin)
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


/* Attempts to set the range of dotter in a sensible way, when calling dotter on the reference
 * sequence versus itself. For now it just uses the big picture range. It's different to setting
 * the big picture range using manual coords because it gets updated whenever the user scrolls,
 * whereas manual coords don't change until the user hits the button again.
 * gb10: this is a much simpler replacement for the ifdef'd out smartDotterRange. We could do
 * something much smarter but at the moment it doesn't seem to be necessary. */
static gboolean smartDotterRange(GtkWidget *blxWindow,
                                 int *dotter_start_out,
                                 int *dotter_end_out,
                                 GError **error)
{
  /* We'll just use a large-ish range centred on the current display range */
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  BlxContext *bc = blxWindowGetContext(blxWindow);
  const IntRange* const displayRange = bigPictureGetDisplayRange(bigPicture);
  
  /* Convert to DNA coords */
  int start = convertDisplayIdxToDnaIdx(displayRange->min(), bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  int end = convertDisplayIdxToDnaIdx(displayRange->max(), bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  
  if (dotter_start_out)
    *dotter_start_out = start;
  
  if (dotter_end_out)
    *dotter_end_out = end;

  boundsLimitValue(dotter_start_out, &bc->refSeqRange);
  boundsLimitValue(dotter_end_out, &bc->refSeqRange);
  
  return TRUE;
}


#ifdef USE_OLD_AUTO_DOTTER_COORDS
/* gb10: Users never liked this way of calculating auto dotter coords because it can return a range
 * that is nowhere near the visible range in the big picture. I've changed auto coords to just
 * use the visible big picture range instead for now. I've left the old code here because it's _almost_
 * good - it just needs to do something more sensible when clipping the range - so this might be
 * useful in future if we can improve it.  */

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

  BlxContext *bc = blxWindowGetContext(blxWindow);

  /* Check that a sequence is selected */
  GList *selectedSeqs = bc->getSelectedSeqsByType(BLXSEQUENCE_MATCH);
  if (g_list_length(selectedSeqs) < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }

  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  const IntRange* const bigPicRange = bigPictureGetDisplayRange(bigPicture);
  const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD) ;

  /* Loop through all MSPs in the selected sequence. We'll estimate the wanted
   * query region from the extent of the HSP's that are completely within view. */
  const BlxSequence *selectedSeq = (const BlxSequence*)(selectedSeqs->data);
  int qMin = UNSET_INT, qMax = UNSET_INT;
  GList *mspListItem = selectedSeq->mspList;  
  
  for ( ; mspListItem ; mspListItem = mspListItem->next)
    {
      const MSP *msp = (MSP*)(mspListItem->data);
      
      /* Get the msp start/end in terms of display coords, and find the min/max */
      const IntRange* const mspDisplayRange = mspGetDisplayRange(msp);

      /* Check if the MSP is in a visible strand is entirely within the big picture range */
      if ((msp->qStrand == activeStrand || (bc->blastMode == BLXMODE_BLASTN)) &&
	  rangesOverlap(mspDisplayRange, bigPicRange))
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
      
      qMin = convertDisplayIdxToDnaIdx(bigPicRange->min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
      qMax = convertDisplayIdxToDnaIdx(bigPicRange->max, bc->seqType, 1, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange);
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
#endif


/* Utility to fetch the selected match sequence's DNA, or get it from the selected MSP.
 * This function assumes that if multiple MSPs are selected, that they are all for 
 * the same match sequence. Returns null if no MSPs are selected, with details of the error
 * in 'error'.  If the sequence was found but there were warnings, it returns non-null with
 * the warnings in 'error'. The return value should be free'd with g_free */
static char* getSelectedSequenceDNA(const BlxSequence *blxSeq, const BlxContext *bc, GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */
  
  char *dotterSSeq = NULL;
  
  /* If we're in seqbl mode, only part of the sequence is stored
   * internally, so try to fetch the full sequence.
   * gb10: I don't think this is applicable any more (or even if it
   * will work if partial sequences are stored). If we do need to do
   * a fetch here then we will need to look for a fetch method that
   * returns the fasta sequence (rather than the embl entry). */
  const BlxBlastMode blastMode = bc->blastMode;

  if (!dotterSSeq && blastMode != BLXMODE_TBLASTX)
    {
      /* Check if sequence is stored internally (i.e. it was passed from acedb) */
      g_debug("Looking for sequence stored internally... ");
    
      dotterSSeq = g_strdup(blxSequenceGetSequence(blxSeq));
      
      if (!dotterSSeq)
	{
	  g_debug("not found.\n");
	  g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NOT_FOUND, "Failed to find sequence for '%s'.\n", blxSequenceGetName(blxSeq));
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
      g_warning("The sequence for '%s' is incomplete.\n", blxSequenceGetName(blxSeq));
    }
  
  return dotterSSeq;
}


static const char* getDotterRefSeqName(const BlxContext *bc, const gboolean transcript)
{
  const char *result = "<no transcript selected>";
  g_return_val_if_fail(bc, result);

  if (transcript)
    {
      /* Get the sequence for a transcript */
      int num_found = 0;
      BlxSequence *transcriptSeq = bc->getSelectedTranscript(&num_found);
      
      if (transcriptSeq)
        result = blxSequenceGetName(transcriptSeq);
      else if (num_found > 0)
        result = "<multiple transcripts>";
    }
  else
    {
      result = bc->refSeqName;
    }

  return result;
}

/* Get the section of reference sequence to dotter against. Also sets the dotterRange if getting
 * the sequence for a transcript. */
static char* getDotterRefSeq(GtkWidget *blxWindow,
                             IntRange *dotterRange, 
                             const int frame, 
                             const gboolean transcript,
                             const char** refSeqName_out,
                             GError **error)
{
  char *result = NULL;
  g_return_val_if_fail(blxWindow && dotterRange, result);

  BlxContext *bc = blxWindowGetContext(blxWindow);
  GError *tmpError = NULL;

  if (transcript)
    {
      /* Get the sequence for a transcript */
      BlxSequence *transcriptSeq = blxWindowGetSelectedTranscript(blxWindow);
      
      if (transcriptSeq)
        {
          result = blxSequenceGetSplicedSequence(transcriptSeq, bc->refSeq, &bc->refSeqRange, &tmpError);
          
          if (refSeqName_out)
            *refSeqName_out = blxSequenceGetName(transcriptSeq);
        }
      else
        {
          if (refSeqName_out)
            *refSeqName_out = "<no transcript selected>";

          g_set_error(&tmpError, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_TRANSCRIPT,
                      "Please select one (and only one) transcript to dotter against.\n\nUse Ctrl to select multiple features i.e. the transcript and the sequence(s) you want to dotter against it.");
        }
      
      /* Set the range of coords of the result. */
      if (result)
        dotterRange->set(1, strlen(result));
    }
  else
    {
      /* Get the sequence for a range */
      result = getSequenceSegment(bc->refSeq,
                                  dotterRange,
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

      if (refSeqName_out)
        *refSeqName_out = bc->refSeqName;
    }

  if (tmpError)
    g_propagate_error(error, tmpError);

  return result;
}


/* Utility to return the texxt contents of a GtkTextView. The result must be free'd by
 * the caller with g_free */
static char *textViewGetText(GtkWidget *textView)
{
  char *result = NULL;
  g_return_val_if_fail(textView && GTK_IS_TEXT_VIEW(textView), result);

  GtkTextBuffer *textBuffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textView));
  GtkTextIter start, end;
  gtk_text_buffer_get_bounds(textBuffer, &start, &end);  

  result = gtk_text_buffer_get_text(textBuffer, &start, &end, TRUE);
  return result;
}


/* Extract the sequence DNA from the text. Also sets the sequence name if 
 * the text is in fasta format; otherwise just expects the raw sequence. This also removes any
 * newlines or whitespace from the sequence */
static void textGetSeqDetails(const char *text, char **sequence, char **sequenceName)
{
  if (!text || !(*text) || !(sequence || sequenceName))
    return ;

  GString *sequenceStr = NULL;
  GString *nameStr = NULL;
  gboolean parsingHeader = FALSE;
  gboolean parsingName = FALSE;
  const char *cp = text;
      
  for ( ; cp && *cp; ++cp)
    {
      if (*cp == '>')
        {
          /* Start of FASTA header */
          parsingName = TRUE;
          parsingHeader = TRUE;
        }
      else if (isWhitespaceChar(*cp))
        {
          /* If we were parsing the name then stop because it ends at the first space. Ignore
           * other whitespace chars. */
          parsingName = FALSE;
        }
      else if (isNewlineChar(*cp))
        {
          /* If we were parsing the header then stop because it ends at the first newline.
           * Ignore other newline chars. */
          parsingHeader = FALSE;
        }
      else if (parsingHeader)
        {
          /* If parsing the name section of the header append it to the name (if name was 
           * requested). Ignore the rest of the header. */
          if (parsingName && sequenceName)
            {
              if (!nameStr)
                nameStr = g_string_new(NULL);

              g_string_append_c(nameStr, *cp);
            }
        }
      else
        {
          /* Must now be parsing the sequence. If the sequence wasn't requested then quit. */
          if (!sequence)
            break;

          if (!sequenceStr)
            sequenceStr = g_string_new(NULL);

          g_string_append_c(sequenceStr, *cp);
        }
    }
      
  if (sequenceStr && sequence)
    *sequence = g_string_free(sequenceStr, FALSE);
  else if (sequenceStr)
    g_string_free(sequenceStr, TRUE);

  if (nameStr && sequenceName)
    *sequenceName = g_string_free(nameStr, FALSE);
  else if (nameStr)
    g_string_free(nameStr, TRUE);
}


/* Extract the sequence DNA from the text in the given widget. Also sets the sequence name if 
 * the text is in fasta format; otherwise just expects the raw sequence. This also removes any
 * newlines or whitespace from the sequence */
static void textViewGetSeqDetails(GtkWidget *textView, char **sequence, char **sequenceName)
{
  char *text = textViewGetText(textView);
  textGetSeqDetails(text, sequence, sequenceName);
  g_free(text);
}


/*******************************************************************
 *		      Functions to call dotter                     *
 *******************************************************************/

/* This actually executes the dotter child process */
static GIOChannel* callDotterChildProcess(GtkWidget *blxWindow,
                                          const char *dotterBinary, 
                                          const int dotterZoom,
                                          const gboolean hspsOnly,
                                          const gboolean sleep,
                                          const char *seq1Name,
                                          const IntRange* const seq1Range,
                                          const BlxStrand seq1Strand,
                                          const gboolean seq1DisplayRev,
                                          const char *seq2Name,
                                          const IntRange* const seq2Range,
                                          const BlxStrand seq2Strand,
                                          const gboolean seq2DisplayRev,
                                          BlxContext *bc,
                                          GPid *childPid,
                                          GError **error)
{
  DEBUG_OUT("callDotterChildProcess\n");

  GIOChannel *ioChannel = NULL;
  char *seq1OffsetStr = convertIntToString(seq1Range->min() - 1);
  char *seq2OffsetStr = convertIntToString(seq2Range->min() - 1);
  char *seq1LenStr = convertIntToString(seq1Range->length());
  char *seq2LenStr = convertIntToString(seq2Range->length());
  
  /* Create the argument list - start with any options we want to pass */
  GSList *argList = NULL;
  argList = g_slist_append(argList, g_strdup(dotterBinary));
  argList = g_slist_append(argList, g_strdup("-z"));
  argList = g_slist_append(argList, convertIntToString(dotterZoom));
  argList = g_slist_append(argList, g_strdup("-q"));
  argList = g_slist_append(argList, seq1OffsetStr);
  argList = g_slist_append(argList, g_strdup("-s"));
  argList = g_slist_append(argList, seq2OffsetStr);
  argList = g_slist_append(argList, g_strdup("--horizontal-type"));
  argList = g_slist_append(argList, g_strdup("d"));
  argList = g_slist_append(argList, g_strdup("--vertical-type"));

  if (bc->seqType == BLXSEQ_PEPTIDE)
    argList = g_slist_append(argList, g_strdup("p"));
  else
    argList = g_slist_append(argList, g_strdup("d"));

  if (bc->flags[BLXFLAG_ABBREV_TITLE])
    argList = g_slist_append(argList, g_strdup("--abbrev-title-on"));
  else
    argList = g_slist_append(argList, g_strdup("--abbrev-title-off"));

  if (bc->windowColor)
    argList = g_slist_append(argList, g_strdup_printf("--session_colour=%s", bc->windowColor));
  
  if (seq1Strand == BLXSTRAND_REVERSE)      argList = g_slist_append(argList, g_strdup("-r"));
  if (seq2Strand == BLXSTRAND_REVERSE)	    argList = g_slist_append(argList, g_strdup("-v"));
  if (seq1DisplayRev)			    argList = g_slist_append(argList, g_strdup("--reverse-h-display"));
  if (seq2DisplayRev)			    argList = g_slist_append(argList, g_strdup("--reverse-v-display"));
  if (hspsOnly)				    argList = g_slist_append(argList, g_strdup("-H"));
  if (sleep)				    argList = g_slist_append(argList, g_strdup("--sleep=15"));
  if (bc->flags[BLXFLAG_NEGATE_COORDS])	    argList = g_slist_append(argList, g_strdup("-N"));

  /* now tell Dotter that we're calling it internally from another SeqTools
   * program, so that it knows to expect piped data */
  argList = g_slist_append(argList, g_strdup("-S"));
  
  /* Now pass the required arguments. These must be in the correct order. */
  argList = g_slist_append(argList, g_strdup(seq1Name));
  argList = g_slist_append(argList, seq1LenStr);
  argList = g_slist_append(argList, g_strdup(seq2Name));
  argList = g_slist_append(argList, seq2LenStr);
  argList = g_slist_append(argList, g_strdup(dotterBinary));

  /* Now pass the screen number - we want to start dotter on the same screen as blixem's main
   * window */
  if (blxWindow && GTK_IS_WINDOW(blxWindow))
    {
      GdkScreen *screen = gtk_window_get_screen(GTK_WINDOW(blxWindow));
      if (screen)
        {
          char *screenStr = convertIntToString(gdk_screen_get_number(screen));
          argList = g_slist_append(argList, g_strdup("--screen"));
          argList = g_slist_append(argList, screenStr);
        }
    }

  /* Terminate list with null */
  argList = g_slist_append(argList, NULL);

  /* Convert the list to an array */
  DEBUG_OUT("Args = ");
  char *argv[g_slist_length(argList)];
  GSList *item = argList;
  int i = 0;
  
  for ( ; item; item = item->next)
    {
      char *arg = (char*)(item->data);
      argv[i] = arg;
      ++i;
      DEBUG_OUT(", %s", arg  );
    }
  DEBUG_OUT("\n");
    
  DEBUG_OUT("Executing dotter\n");

  int standard_input = 0;

  gboolean ok = g_spawn_async_with_pipes(NULL, //inherit parent' working directory
                                         argv,
                                         NULL, //inherit parent's environment
                                         (GSpawnFlags)0,
                                         NULL, 
                                         NULL,
                                         childPid,
                                         &standard_input,
                                         NULL,
                                         NULL,
                                         error);

  if (ok)
    ioChannel = g_io_channel_unix_new(standard_input);

  return ioChannel;
}


/* Call dotter as an external process */
gboolean callDotterExternal(GtkWidget *blxWindow,
                            BlxContext *bc,
                            int dotterZoom, 
                            const gboolean hspsOnly,
                            const gboolean sleep,
                            const char *seq1Name,
                            IntRange *seq1Range,
                            char *seq1,
                            const BlxStrand seq1Strand,
			    const gboolean seq1DisplayRev,
                            const char *seq2Name,
                            IntRange *seq2Range,
                            char *seq2,
			    const BlxStrand seq2Strand,
			    const gboolean seq2DisplayRev,
                            const gboolean clipRange,
                            const IntRange* const refSeqRange,
                            const BlxSequence *transcriptSeq,
                            GError **error)
{
  if (clipRange)
    seq1Range->boundsLimit(&bc->refSeqRange, FALSE);

  static char *dotterBinary = NULL;

  /* Open pipe to new dotterBinary */
  if (!dotterBinary) 
    { 
      g_debug("Looking for Dotter ...\n");
      dotterBinary = g_find_program_in_path("dotter");

      if (!dotterBinary)
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_EXE, "No dotter executable found in path.\n$PATH=%s\n", getenv("PATH"));
          dotterBinary = NULL;
          return FALSE;
        }
    }
  
  g_debug("Calling %s with region: %d,%d - %d,%d\n", dotterBinary, seq1Range->min(), seq2Range->min(), seq1Range->max(), seq2Range->max());

  /* Create the child process */
  GPid childPid = 0;
  gsize bytes_written = 0;
  GError *tmpError = NULL;
  
  GIOChannel *ioChannel = callDotterChildProcess(blxWindow, dotterBinary, dotterZoom, hspsOnly, sleep,
                                                 seq1Name, seq1Range, seq1Strand, seq1DisplayRev,
                                                 seq2Name, seq2Range, seq2Strand, seq2DisplayRev,
                                                 bc, &childPid, &tmpError);

  if (ioChannel)
    {
      g_debug("Spawned process %d\n", childPid);
      bc->spawnedProcesses = g_slist_append(bc->spawnedProcesses, GINT_TO_POINTER((int)childPid));
    }
  else
    {
      if (tmpError)
        prefixError(tmpError, "Error creating child process for Dotter: ");
      else
        g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_EXE, "Error creating child process for Dotter.\n");
    }

  if (!tmpError)
    {
      /* Pass the sequences */
      DEBUG_OUT("Piping sequences to dotter...\n");

      g_io_channel_write_chars(ioChannel, seq1, seq1Range->length(), &bytes_written, &tmpError);

      if (!tmpError)
        g_io_channel_write_chars(ioChannel, seq2, seq2Range->length(), &bytes_written, &tmpError);

      DEBUG_OUT("...done\n");
    }

  if (!tmpError)
    {
      /* Pass the features */
      DEBUG_OUT("Piping features to dotter...\n");

      if (transcriptSeq)
        {
          writeTranscriptToOutput(ioChannel, transcriptSeq, seq1Range, refSeqRange, &tmpError);
        }
      else
        {
          GList *seqItem = bc->matchSeqs;
          for ( ; seqItem; seqItem = seqItem->next) 
            {
              BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
              writeBlxSequenceToOutput(ioChannel, blxSeq, seq1Range, seq2Range, &tmpError);
            }
        }
      
      DEBUG_OUT("...done\n");
    }
      
  /* Close the channel */
  g_io_channel_shutdown(ioChannel, TRUE, NULL);

  if (tmpError)
    g_propagate_error(error, tmpError);

  return TRUE;
}


/* Call dotter on the given sequence. Returns true if dotter was called; 
 * false if we quit trying. */
gboolean callDotterOnSelectedSeq(GtkWidget *blxWindow, 
                                 const gboolean hspsOnly,
                                 const gboolean sleep,
                                 const DotterRefType refType,
                                 const BlxSequence *selectedSeq,
                                 GError **error)
{
  BlxContext *bc = blxWindowGetContext(blxWindow);

  /* Check this sequence is a valid blast match (just check the first MSP;
   * they must all the same type if they have the same seq name) */
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
          const MSP* const msp = (const MSP*)(mspItem->data);
          found = (msp->qStrand == qStrand);
        }
      
      if (!found)
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_INVALID_STRAND, "You must select a match sequence on the active strand, or toggle strands first.\n");
          return FALSE;
        }
    }
  
  /* Make a copy of the match sequence, because dotter takes ownership of this. */
  char *dotterSSeq = getSelectedSequenceDNA(selectedSeq, bc, NULL);

  if (!dotterSSeq)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQ_DATA, "No sequence data for this sequence.\n");
      return FALSE;
    }
  
  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT, dotterZoom = 0;
  GError *rangeError = NULL;
  
  gboolean ok = getDotterRange(blxWindow, BLXDOTTER_MATCH_SELECTED, refType, &dotterStart, &dotterEnd, &dotterZoom, &rangeError);

  if (!ok)
    {
      prefixError(rangeError, "Error calculating dotter range. ");
      g_propagate_error(error, rangeError);
      return FALSE;
    }
  else if (ok && rangeError && error) /* if error is null, don't issue warnings */
    {
      /* There was a warning when calculating the range. Ask the user if they want to continue. */
      prefixError(rangeError, "Warning: ");
      postfixError(rangeError, "\nContinue?");

      char *title = g_strdup_printf("%sWarning", blxGetTitlePrefix(bc));
      ok = (runConfirmationBox(blxWindow, title, rangeError->message) == GTK_RESPONSE_ACCEPT);
      
      g_free(title);
      g_error_free(rangeError);
      rangeError = NULL;
      
      if (!ok)
	return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const int frame = mspGetRefFrame(firstMsp, bc->seqType);
  IntRange dotterRange;
  dotterRange.set(dotterStart, dotterEnd);
  GError *seqError = NULL;
  const char *refSeqName = NULL;
  const gboolean transcript = (refType == BLXDOTTER_REF_TRANSCRIPT);

  char *refSeqSegment = getDotterRefSeq(blxWindow, &dotterRange, frame, transcript, &refSeqName, &seqError);
  
  if (!refSeqSegment)
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

  IntRange sRange = {1, (int)strlen(dotterSSeq)};
  
  const int offset = dotterRange.min() - 1;
  const BlxStrand refSeqStrand = blxWindowGetActiveStrand(blxWindow);
  
  if (transcript)
    g_message("Calling dotter on %s vs %s\n", refSeqName, dotterSName);
  else
    g_message("Calling dotter on %s [%d,%d] vs %s\n", refSeqName, dotterStart, dotterEnd, dotterSName);
  
  g_debug("reference sequence: name =  %s, offset = %d\n"
          "    match sequence: name =  %s, offset = %d\n", 
          refSeqName, offset, dotterSName, 0);
  
  const gboolean revHozScale = (refSeqStrand == BLXSTRAND_REVERSE);
  const gboolean revVertScale = FALSE; /* don't rev match seq scale, because it would show in dotter with -ve coords, but blixem always shows +ve coords */
  const gboolean clipRange = !transcript; /* don't clip the range if using a transcript range */

  /* If dottering vs a specific transcript, get it now */
  const BlxSequence *transcriptSeq = NULL;
  if (refType == BLXDOTTER_REF_TRANSCRIPT)
    transcriptSeq = blxWindowGetSelectedTranscript(blxWindow);

  return callDotterExternal(blxWindow, bc, dotterZoom, hspsOnly, sleep,
                            refSeqName, &dotterRange, refSeqSegment, refSeqStrand, revHozScale,
                            dotterSName, &sRange, dotterSSeq, selectedSeq->strand, revVertScale,
                            clipRange, &bc->refSeqRange, transcriptSeq, error);
}


/* Call dotter on the currently-selected sequence(s). Returns true if all dotters were; 
 * successful; false if there was a problem. */
gboolean callDotterOnSelectedSeqs(GtkWidget *blxWindow, 
                                  const gboolean hspsOnly,
                                  const gboolean startup_sleep,
                                  const DotterRefType refType,
                                  GError **error)
{
  gboolean result = FALSE;
  g_return_val_if_fail(!error || *error == NULL, result); /* if error is passed, its contents must
                                                            be NULL */
  const int max_dotters = 5;
  GList *selectedSeqs = blxWindowGetSelectedSeqsByType(blxWindow, BLXSEQUENCE_MATCH);

  const int numSeqsSelected = g_list_length(selectedSeqs);
  
  if (numSeqsSelected < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }
  else
    {
      gboolean ok = TRUE;

      if (numSeqsSelected > max_dotters)
        {
          char *msg = g_strdup_printf("You are about to open %d dotters - continue?", numSeqsSelected);

          if (runConfirmationBox(blxWindow, "Run multiple dotters", msg) != GTK_RESPONSE_ACCEPT)
            ok = FALSE;
            
          g_free(msg);
        }

      if (ok)
        {
          /* Call dotter on each selected sequence */
          result = TRUE;
          GList *item = selectedSeqs;

          for ( ; item; item = item->next)
            {
              const BlxSequence *selectedSeq = (const BlxSequence*)(item->data);
              result &= callDotterOnSelectedSeq(blxWindow, hspsOnly, startup_sleep, refType, selectedSeq, error);
            }
        }
    }

  return result;
}


/* Call dotter on the manually-entered adhoc sequence. Returns true if dotter 
 * was called; false if we quit trying. */
gboolean callDotterOnAdhocSeq(DotterDialogData *dialogData, GError **error)
{
  gboolean result = FALSE;
  g_return_val_if_fail(!error || *error == NULL, result); /* if error is passed it must be NULL */
  
  GtkWidget *blxWindow = dialogData->blxWindow;
  BlxContext *bc = blxWindowGetContext(blxWindow);
  const gboolean transcript = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->transcriptButton));
  
  /* We will display the active strand as the main strand in dotter */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);

  /* Get the sequence DNA (and name, if in fasta format) from the text entry */
  char *dotterSName = NULL;
  char *dotterSSeq = NULL;
  textViewGetSeqDetails(dialogData->adhocSeqText, &dotterSSeq, &dotterSName);
  
  if (!dotterSSeq || *dotterSSeq == 0)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "Please enter a sequence into the entry box.\n");
      return FALSE;
    }

  if (!dotterSName)
    {
      /* If we were given a raw sequence we won't have a name: just set a dummy name */
      dotterSName = g_strdup("Unknown sequence");
    }

  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT, dotterZoom = 0;
  GError *rangeError = NULL;
  
  gboolean ok = getDotterRange(blxWindow, dialogData->matchType, dialogData->refType, &dotterStart, &dotterEnd, &dotterZoom, &rangeError);

  if (!ok)
    {
      prefixError(rangeError, "Error calculating dotter range. ");
      g_propagate_error(error, rangeError);
      return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const int frame = 1;
  IntRange dotterRange;
  dotterRange.set(dotterStart, dotterEnd);
  GError *seqError = NULL;
  const char *refSeqName = NULL;

  char *refSeqSegment = getDotterRefSeq(blxWindow, &dotterRange, frame, transcript, &refSeqName, &seqError);
  
  if (!refSeqSegment)
    {
      g_propagate_error(error, seqError);
      return FALSE;
    }
  else
    {
      /* If there was an error set but the sequence was still returned then it's a non-critical warning */
      reportAndClearIfError(&seqError, G_LOG_LEVEL_WARNING);
    }

  IntRange sRange = {1, (int)strlen(dotterSSeq)};
  
  const int offset = dotterRange.min() - 1;
  const BlxStrand refSeqStrand = blxWindowGetActiveStrand(blxWindow);
  
  if (transcript)
    g_message("Calling dotter on %s vs %s\n", refSeqName, dotterSName);
  else
    g_message("Calling dotter on %s [%d,%d] vs %s\n", refSeqName, dotterStart, dotterEnd, dotterSName);
  
  g_debug("reference sequence: name =  %s, offset = %d\n"
          "    match sequence: name =  %s, offset = %d\n", 
          refSeqName, offset, dotterSName, 0);
  
  const gboolean revHozScale = (refSeqStrand == BLXSTRAND_REVERSE);
  const gboolean revVertScale = FALSE; /* don't rev match seq scale, because it would show in dotter with -ve coords, but blixem always shows +ve coords */
  const gboolean clipRange = !transcript; /* don't clip the range if using a transcript range */

  /* If dottering vs a specific transcript, get it now */
  const BlxSequence *transcriptSeq = NULL;
  if (dialogData->refType == BLXDOTTER_REF_TRANSCRIPT)
    transcriptSeq = blxWindowGetSelectedTranscript(blxWindow);

  result = callDotterExternal(blxWindow, bc, dotterZoom, FALSE, dialogData->sleep,
                              refSeqName, &dotterRange, refSeqSegment, refSeqStrand, revHozScale,
                              dotterSName, &sRange, dotterSSeq, qStrand, revVertScale,
                              clipRange, &bc->refSeqRange, transcriptSeq, error);

  /* dotter takes ownership of dotterSSeq but not dotterSName, so free it */
  g_free(dotterSName);
  
  return result;
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
static gboolean callDotterOnSelf(DotterDialogData *dialogData, GError **error)
{
  GtkWidget *blxWindow = dialogData->blxWindow;
  BlxContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the auto range, if requested */
  int dotterStart = UNSET_INT;
  int dotterEnd = UNSET_INT;
  int dotterZoom = 0;
  
  GError *tmpError = NULL;
  if (!getDotterRange(blxWindow, dialogData->matchType, dialogData->refType, &dotterStart, &dotterEnd, &dotterZoom, &tmpError))
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);
  const int frame = 1;
  IntRange qRange;
  qRange.set(dotterStart, dotterEnd);
  const char *refSeqName = NULL;
    
  const gboolean transcript = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->transcriptButton));
  char *refSeqSegment = getDotterRefSeq(blxWindow, &qRange, frame, transcript, &refSeqName, &tmpError);

  if (!refSeqSegment)
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
  char *dotterSSeq = g_strdup(refSeqSegment);

  const gboolean revScale = (qStrand == BLXSTRAND_REVERSE);
  
  if (transcript)
    g_message("Calling dotter on %s vs itself\n", refSeqName);
  else
    g_message("Calling dotter on %s [%d,%d] vs itself\n", refSeqName, dotterStart, dotterEnd);

  /* If dottering vs a specific transcript, get it now */
  const BlxSequence *transcriptSeq = NULL;
  if (dialogData->refType == BLXDOTTER_REF_TRANSCRIPT)
    transcriptSeq = blxWindowGetSelectedTranscript(blxWindow);

  callDotterExternal(blxWindow, bc, dotterZoom, FALSE, dialogData->sleep,
                     refSeqName, &qRange, refSeqSegment, qStrand, revScale,
                     refSeqName, &qRange, dotterSSeq, qStrand, revScale,
                     FALSE, &bc->refSeqRange, transcriptSeq, error);

  return TRUE;
}

