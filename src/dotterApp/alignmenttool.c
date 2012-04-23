/*  File: alignmenttool.c
 *  Author: Gemma Barson, 2010-09-02
 *  Copyright (c) 2010 - 2012 Genome Research Ltd
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
 * Description: Creates the alignment tool window for a Dotter window.
 *              The alignment tool shows the sequence data for the currently-
 *              selected coordinates in the main Dotter window that it is
 *              associated with. The two windows work together so that when the
 *              selected coords are changed on one window they are updated on
 *              the other. The alignment tool is owned by the Dotter window and
 *              will be destroyed when the Dotter window is closed.
 *----------------------------------------------------------------------------
 */

#include <dotterApp/dotter_.h>
#include <seqtoolsUtils/utilities.h>
#include <gtk/gtk.h>
#include <math.h>
#include <string.h>
#include <gdk/gdkkeysyms.h>


static int atob[]       /* OLD (starting at 1) ASCII-to-binary translation table  (Inherited from blast) */
= {
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,24,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,
NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,

NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};



#define DEFAULT_ALIGNMENT_LENGTH                    125   /* default number of sequence chars to display in the alignment tool */
#define MAX_ALIGNMENT_LENGTH                        501   /* max number of sequence chars to display in the alignment tool */
#define MIN_ALIGNMENT_LENGTH                        0     /* min number of sequence chars to display in the alignment tool */
#define SELECTED_COORD_MARKER_HEIGHT                12    /* the height of the marker that indicates the currently-selected coord */


/* Function pointer for expose functions */
typedef gboolean (*ExposeFunc)(GtkWidget *widget, GdkEventExpose *event, gpointer data);


typedef struct _SequenceProperties
{
  const char *seqName;
  const char *sequence;
  BlxSeqType seqType;               /* whether this sequence is in nucleotide or peptide coords */
  BlxStrand strand;
  int frame;
  IntRange *displayRange;           /* the displayed range; pointer to refDisplayRange or matchDisplayRange */
  IntRange *fullRange;              /* the full range of the sequence; pointer to refSeqRange or matchSeqRange */
  gboolean scaleReversed;           /* whether the scale for this sequence is shown reversed (i.e. high-to-low rather than low-to-high) */
  GSList *compSeqs;                 /* list of other sequence widgets that this sequence will be compared against */
  gboolean isMatchSeq;              /* true if this is the match seq, false if it's the ref seq */
} SequenceProperties;


typedef struct _AlignmentToolProperties
{
  int alignmentLen;                 /* the number of coords wide the alignment too displays */
  IntRange refDisplayRange;         /* the current ref seq range displayed in the tool */
  IntRange matchDisplayRange;       /* the current match seq range displayed in the tool */

  gboolean dragging;                /* true when a drag is in progress */
  int dragStart;                    /* set to the start coord where the user initiates a drag */

  GtkWidget *selectionWidget;       /* text in this widget should be highlighted (null if none) */
  IntRange selectionRange;
  
  DotterWindowContext *dotterWinCtx;
  GtkActionGroup *actionGroup;
} AlignmentToolProperties;



/* Local function declarations */
static void                        onCloseMenu(GtkAction *action, gpointer data);
static void                        onPrintMenu(GtkAction *action, gpointer data);
static void                        onSetLengthMenu(GtkAction *action, gpointer data);
static void                        onCopyHCoordMenu(GtkAction *action, gpointer data);
static void                        onCopyVCoordMenu(GtkAction *action, gpointer data);
static void                        onCopySelnMenu(GtkAction *action, gpointer data);
static void                        onCopySelnCoordsMenu(GtkAction *action, gpointer data);
static void                        onClearSelnMenu(GtkAction *action, gpointer data);
static void                        drawSequence(GdkDrawable *drawable, GtkWidget *widget, GtkWidget *alignmentTool);
static void                        drawSequenceHeader(GtkWidget *widget, GtkWidget *alignmentTool, GdkDrawable *drawable, const gboolean horizontal);
static int                         getSequenceOffset(SequenceProperties *properties, DotterContext *dc);
static int                         getCoordAtPos(const int x, GtkWidget *sequenceWidget, GtkWidget *alignmentTool);
static void                        connectSequenceSignals(GtkWidget *widget, GtkWidget *alignmentTool);
static void                        sequenceInitiateDragging(GtkWidget *sequenceWidget, GtkWidget *alignmentTool, const int x);
static void                        sequenceFinishDragging(GtkWidget *sequenceWidget, GtkWidget *alignmentTool, const int x);
static void                        highlightSequenceBase(SequenceProperties *seq1, AlignmentToolProperties *atProperties, DotterWindowContext *dwc, const int displayIdx, const int seq1Idx, const int seq1Start, const gboolean highlight, GdkGC *gc, GdkDrawable *drawable);
static void                        selectVisibleSequence(GtkWidget *sequenceWidget, GtkWidget *alignmentTool);
static int                         getDisplayStart(SequenceProperties *properties, DotterContext *dc);
static char*                       getSequenceBetweenCoords(GtkWidget *sequenceWidget, const int startCoord, const int endCoord, DotterWindowContext *dwc);
static void                        clearSequenceSelection(GtkWidget *alignmentTool);

/* Menu builders - standard menu entries */
static const GtkActionEntry alignmentToolMenuEntries[] = {
{ "Close",          NULL, "_Close tool",              "<control>W", "Close the alignment tool",             G_CALLBACK(onCloseMenu)},
{ "Print",          NULL, "_Print...",                "<control>P", "Print the alignment tool window",      G_CALLBACK(onPrintMenu)},
{ "SetLength",      NULL, "_Set alignment length",    "<control>L", "Set the length of the alignment tool", G_CALLBACK(onSetLengthMenu)},
{ "CopyHCoord",     NULL, "Copy _horizontal coord",   NULL,         "Copy the current horizontal sequence coord to the clipboard", G_CALLBACK(onCopyHCoordMenu)},
{ "CopyVCoord",     NULL, "Copy _vertical coord",     NULL,         "Copy the current vertical sequence coord to the clipboard", G_CALLBACK(onCopyVCoordMenu)},
{ "CopySeln",       NULL, "Copy selectio_n",          "<control>C", "Copy the current selection to the clipboard", G_CALLBACK(onCopySelnMenu)},
{ "ClearSeln",      NULL, "C_lear current selection",  "Escape","Clear the current selection", G_CALLBACK(onClearSelnMenu)},
{ "CopySelnCoords", NULL, "Copy selection coor_ds",   "<shift><control>C","Copy the start/end coords of the current selection to the clipboard", G_CALLBACK(onCopySelnCoordsMenu)}
};


/* This defines the layout of the menu */
static const char alignmentToolMenuDescription[] =
"<ui>"
"  <popup name='MainMenu' accelerators='true'>"
"      <menuitem action='CopyHCoord'/>"
"      <menuitem action='CopyVCoord'/>"
"      <menuitem action='CopySeln'/>"
"      <menuitem action='CopySelnCoords'/>"
"      <menuitem action='ClearSeln'/>"
"      <separator/>"
"      <menuitem action='Close'/>"
"      <menuitem action='Print'/>"
"      <separator/>"
"      <menuitem action='SetLength'/>"
"  </popup>"
"</ui>";




/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static SequenceProperties* sequenceGetProperties(GtkWidget *widget)
{
  return widget ? (SequenceProperties*)(g_object_get_data(G_OBJECT(widget), "SequenceProperties")) : NULL;
}

static void onDestroySequenceWidget(GtkWidget *widget)
{
  SequenceProperties *properties = sequenceGetProperties(widget);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "SequenceProperties", NULL);
    }
}

static void sequenceCreateProperties(GtkWidget *widget, 
                                     const char *seqName,
                                     const char *sequence,
                                     const BlxSeqType seqType,
                                     const BlxStrand strand,
                                     const int frame,
                                     IntRange *displayRange,
                                     IntRange *fullRange,
                                     gboolean scaleReversed,
                                     GSList *compSeqs,
                                     gboolean isMatchSeq)
{
  if (widget)
    {
      SequenceProperties *properties = g_malloc(sizeof *properties);
      
      properties->seqName = seqName;
      properties->sequence = sequence;
      properties->seqType = seqType;
      properties->strand = strand;
      properties->frame = frame;
      properties->displayRange = displayRange;
      properties->fullRange = fullRange;
      properties->scaleReversed = scaleReversed;
      properties->compSeqs = compSeqs;
      properties->isMatchSeq = isMatchSeq;
      
      g_object_set_data(G_OBJECT(widget), "SequenceProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroySequenceWidget), NULL); 
    }
}


static AlignmentToolProperties* alignmentToolGetProperties(GtkWidget *widget)
{
  return widget ? (AlignmentToolProperties*)(g_object_get_data(G_OBJECT(widget), "AlignmentToolProperties")) : NULL;
}

static void onDestroyAlignmentTool(GtkWidget *widget)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(widget);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "AlignmentToolProperties", NULL);
    }
}

static void alignmentToolCreateProperties(GtkWidget *widget, DotterWindowContext *dotterWinCtx)
{
  if (widget)
    {
      AlignmentToolProperties *properties = g_malloc(sizeof *properties);
    
      properties->dotterWinCtx = dotterWinCtx;
      properties->actionGroup = NULL;
      properties->alignmentLen = DEFAULT_ALIGNMENT_LENGTH;
      properties->refDisplayRange.min = 0;
      properties->refDisplayRange.max = 20;
      properties->matchDisplayRange.min = 0;
      properties->matchDisplayRange.max = 20; /* to do: put real values in here for the ranges */

      properties->dragging = FALSE;
      properties->dragStart = 0;
      properties->selectionWidget = NULL;
      properties->selectionRange.min = UNSET_INT;
      properties->selectionRange.max = UNSET_INT;
      
      g_object_set_data(G_OBJECT(widget), "AlignmentToolProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyAlignmentTool), NULL); 
    }
}


/***********************************************************
 *                       Events                            *
 ***********************************************************/

/* Called when the alignment tool window changes size */
static void onSizeAllocateAlignmentTool(GtkWidget *alignmentTool, GtkAllocation *allocation, gpointer data)
{
  callFuncOnAllChildWidgets(alignmentTool, widgetClearCachedDrawable);
  gtk_widget_queue_draw(alignmentTool);
}


/* Expose function for a widget containing a section of sequence. SequenceProperties should
 * be set on any widget that this is to be called for. */
static gboolean onExposeSequence(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkWindow *drawable = widgetGetDrawable(widget);
  
  if (!drawable)
    {
      drawable = createBlankPixmap(widget);
      drawSequence(drawable, widget, alignmentTool);
    }
  
  if (drawable)
    {
      GdkGC *gc = gdk_gc_new(widget->window);
      gdk_draw_drawable(widget->window, gc, drawable, 0, 0, 0, 0, -1, -1);
      g_object_unref(gc);
    }
  
  return TRUE;
}


/* Expose function for a widget containing a reference sequence header. */
static gboolean onExposeRefSequenceHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkDrawable *drawable = widgetGetDrawable(widget);

  if (!drawable)
    {
      drawable = createBlankPixmap(widget);
      drawSequenceHeader(widget, alignmentTool, drawable, TRUE);
    }
  
  if (drawable)
    {
      GdkGC *gc = gdk_gc_new(widget->window);
      gdk_draw_drawable(widget->window, gc, drawable, 0, 0, 0, 0, -1, -1);
      g_object_unref(gc);
    }
  
  return TRUE;
}


/* Expose function for a widget containing a match sequence header. */
static gboolean onExposeMatchSequenceHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkDrawable *drawable = widgetGetDrawable(widget);
  
  if (!drawable)
    {
      drawable = createBlankPixmap(widget);
      drawSequenceHeader(widget, alignmentTool, drawable, FALSE);
    }
  
  if (drawable)
    {
      GdkGC *gc = gdk_gc_new(widget->window);
      gdk_draw_drawable(widget->window, gc, drawable, 0, 0, 0, 0, -1, -1);
      g_object_unref(gc);
    }
  
  return TRUE;
}


/* This should be called when the range displayed by the alignment tool has changed */
static void onAlignmentToolRangeChanged(GtkWidget *alignmentTool)
{
  /* Clear any cached drawables on this or any child widgets */
  callFuncOnAllChildWidgets(alignmentTool, widgetClearCachedDrawable);
  gtk_widget_queue_draw(alignmentTool);
}


/* Create the menu */
static GtkWidget* createAlignmentToolMenu(GtkWidget *window, GtkActionGroup **actionGroup_out)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  gtk_action_group_add_actions (action_group, alignmentToolMenuEntries, G_N_ELEMENTS (alignmentToolMenuEntries), window);
  enableMenuAction(action_group, "CopySeln", FALSE);
  enableMenuAction(action_group, "CopySelnCoords", FALSE);
  enableMenuAction(action_group, "ClearSeln", FALSE);

  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);
  
  GError *error = NULL;
  const char *menuDescription = alignmentToolMenuDescription;
  
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, menuDescription, -1, &error))
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }

  if (actionGroup_out)
    *actionGroup_out = action_group;
  
  return gtk_ui_manager_get_widget (ui_manager, "/MainMenu");
}


/* General mouse button handler */
static gboolean onButtonPressAlignmentTool(GtkWidget *alignmentTool, GdkEventButton *event, gpointer data)
{
  gboolean handled = TRUE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      handled = TRUE;
    }
  
  
  return handled;
}


/* Mouse button handler when clicking on a sequence */
static gboolean onButtonPressSequence(GtkWidget *sequenceWidget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      GtkWidget *alignmentTool = GTK_WIDGET(data);
      sequenceInitiateDragging(sequenceWidget, alignmentTool, event->x);
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* left double-click */
    {
      GtkWidget *alignmentTool = GTK_WIDGET(data);
      selectVisibleSequence(sequenceWidget, alignmentTool);
      handled = TRUE;
    }

  return handled;
}


/* Mouse button handler when releasing from a click on a sequence widget */
static gboolean onButtonReleaseSequence(GtkWidget *sequenceWidget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->button == 1) /* left click */
    {
      GtkWidget *alignmentTool = GTK_WIDGET(data);
      sequenceFinishDragging(sequenceWidget, alignmentTool, event->x);
      handled = TRUE;
    }
    
  return handled;
}


static gboolean onMouseMoveSequence(GtkWidget *sequenceWidget, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->state & GDK_BUTTON1_MASK) /* left drag */
    {
      /* Update the current selection range */
      GtkWidget *alignmentTool = GTK_WIDGET(data);
      
      AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
      const int newCoord = getCoordAtPos(event->x, sequenceWidget, alignmentTool);
      
      intrangeSetValues(&atProperties->selectionRange, atProperties->dragStart, newCoord);
      widgetClearCachedDrawable(sequenceWidget, NULL);
      gtk_widget_queue_draw(sequenceWidget);
      
      handled = TRUE;
    }
  
  return handled;
}


/***********************************************************
 *                       Get/set functions                 *
 ***********************************************************/

/* Set the alignment range centre based on the currently-selected match/ref seq coord. */
void updateAlignmentRange(GtkWidget *alignmentTool, DotterWindowContext *dwc)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  DotterContext *dc = dwc->dotterCtx;

  /* Re-create the range, centred on the set coordinate and with the alignment tool's alignment 
   * length. Note that the length is the number of display chars but the reference sequence is
   * in nucleotide coords so needs converting. (The match sequence is always in display coords so
   * will be correct whether we're displaying nucleotides or peptides.) */
  int len = properties->alignmentLen * getResFactor(dc, TRUE);
  int offset = ceil((double)properties->alignmentLen / 2.0) * getResFactor(dc, TRUE);
  properties->refDisplayRange.min = dwc->refCoord - offset;
  properties->refDisplayRange.max = properties->refDisplayRange.min + len;

  len = properties->alignmentLen * getResFactor(dc, FALSE);
  offset = ceil((double)properties->alignmentLen / 2.0) * getResFactor(dc, FALSE);
  properties->matchDisplayRange.min = dwc->matchCoord - offset;
  properties->matchDisplayRange.max = properties->matchDisplayRange.min + len;
  
  onAlignmentToolRangeChanged(alignmentTool);
}


/* Set the alignment length */
static void setAlignmentLength(GtkWidget *alignmentTool, const int length, GError **error)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);

  if (length > MIN_ALIGNMENT_LENGTH && length <= MAX_ALIGNMENT_LENGTH)
    {
      /* The display code assumes we always have an odd length; if we
       * are given an even length, round it up */
      if (length % 2 == 0)
        properties->alignmentLen = length + 1;
      else
        properties->alignmentLen = length;
    }
  else
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_INVALID_LENGTH, "Length %d is not in valid range %d -> %d.\n",
                  length, MIN_ALIGNMENT_LENGTH, MAX_ALIGNMENT_LENGTH);
    }

  updateAlignmentRange(alignmentTool, properties->dotterWinCtx);
}

/***********************************************************
 *                          Menus                          *
 ***********************************************************/

static void onCloseMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  
  if (alignmentTool)
    {
      gtk_widget_hide(alignmentTool);
    }
}


/* Called when the user selects the print menu option */
static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  blxPrintWidget(alignmentTool, NULL, GTK_WINDOW(alignmentTool), &dwc->printSettings, &dwc->pageSetup, NULL, TRUE, PRINT_FIT_BOTH);
}


/* Callback from the SetLength menu */
static gboolean onAlignmentLenChanged(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  
  const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
  const int newVal = convertStringToInt(text);
  
  GError *error = NULL;
  setAlignmentLength(alignmentTool, newVal, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      result = FALSE;
    }
  
  return result;
}


/* Set-length menu item: pops up a dialog asking the user to enter the alignment length */
static void onSetLengthMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  
  if (alignmentTool)
    {
      GtkWidget *parent = gtk_widget_get_toplevel(alignmentTool);
      
      GtkWidget *dialog = gtk_dialog_new_with_buttons("Dotter - Set alignment length", 
                                                      GTK_WINDOW(parent), 
                                                      GTK_DIALOG_DESTROY_WITH_PARENT,
                                                      GTK_STOCK_CANCEL,
                                                      GTK_RESPONSE_REJECT,
                                                      GTK_STOCK_OK,
                                                      GTK_RESPONSE_ACCEPT,
                                                      NULL);
      
      gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
      
      GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;
      AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
      
      GtkWidget *entry = gtk_entry_new();
      gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

      char *displayText = convertIntToString(properties->alignmentLen);
      gtk_entry_set_text(GTK_ENTRY(entry), displayText);
      g_free(displayText);
      
      widgetSetCallbackData(entry, onAlignmentLenChanged, alignmentTool);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
      
      gtk_box_pack_start(GTK_BOX(contentArea), gtk_label_new("Enter the new length:"), FALSE, FALSE, 0);
      gtk_box_pack_start(GTK_BOX(contentArea), entry, TRUE, TRUE, 0);

      gtk_widget_show_all(dialog);
    }
}


/* Utility to copy an integer value as a string to the primary clipboard */
static void copyIntToPrimaryClipboard(const int val)
{
  char *displayText = convertIntToString(val);
  setPrimaryClipboardText(displayText);
  g_free(displayText); 
}

/* Callback called when the user selects the 'copy horizontal coord' menu option */
static void onCopyHCoordMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  
  copyIntToDefaultClipboard(properties->dotterWinCtx->refCoord);
  copyIntToPrimaryClipboard(properties->dotterWinCtx->refCoord);
}

/* Callback called when the user selects the 'copy vertical coord' menu option */
static void onCopyVCoordMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);

  copyIntToDefaultClipboard(properties->dotterWinCtx->matchCoord);
  copyIntToPrimaryClipboard(properties->dotterWinCtx->matchCoord);
}

/* Callback called when the user selects the 'copy selection' menu option */
static void onCopySelnMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);

  if (atProperties->selectionWidget)
    {
      /* copy the selection to the clipboard */
      char *text = getSequenceBetweenCoords(atProperties->selectionWidget,
                                            atProperties->selectionRange.min,
                                            atProperties->selectionRange.max,
                                            atProperties->dotterWinCtx);

      setDefaultClipboardText(text);
      setPrimaryClipboardText(text);
      
      g_free(text);
    }
}

/* Callback called when the user selects the 'copy selection coords' menu option */
static void onCopySelnCoordsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);

  if (atProperties->selectionWidget)
    {
      SequenceProperties *properties = sequenceGetProperties(atProperties->selectionWidget);
      
      const int start = (properties->strand  == BLXSTRAND_REVERSE ? atProperties->selectionRange.max : atProperties->selectionRange.min);
      const int end = (properties->strand  == BLXSTRAND_REVERSE ? atProperties->selectionRange.min : atProperties->selectionRange.max);
      
      char *text = g_strdup_printf("%d, %d", start, end);

      setDefaultClipboardText(text);
      setPrimaryClipboardText(text);
      
      g_free(text);
    }
}

/* Callback called when the user selects the 'clear selection' menu option */
static void onClearSelnMenu(GtkAction *action, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  clearSequenceSelection(alignmentTool);
}


/***********************************************************
 *                       Initialisation                    *
 ***********************************************************/

/* Create a widget that will display the reference sequence for the given
 * frame/strand. Add it to the given row in the given table and increment
 * the row number. Also add it to the given refSeqList */
static void createRefSeqWidget(GtkWidget *alignmentTool,
                               AlignmentToolProperties *properties,
                               const int frame,
                               const BlxStrand strand,
                               const char strandChar,
                               GSList *matchSeqList,
                               GtkTable *table,
                               const int xpad,
                               const int ypad,
                               int *row,
                               GSList **refSeqList)
{
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* Create a label containing the name */
  char *text = blxprintf("%s (%c%d):", dc->refSeqName, strandChar, frame);
  GtkWidget *label = createLabel(text, 0.0, 0.0, TRUE, TRUE);
  gtk_widget_modify_font(label, dc->fontDesc);
  gtk_misc_set_alignment(GTK_MISC(label), 0, 0);
  g_free(text);
  gtk_table_attach(table, label, 1, 2, *row, *row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  /* Find which sequence is applicable (fwd or rev strand or translated peptide seq) */
  char *sequence = dc->refSeq;
  
  if (dc->displaySeqType == BLXSEQ_DNA && strand == BLXSTRAND_REVERSE)
    sequence = dc->refSeqRev;
  else if (dc->blastMode == BLXMODE_BLASTX)
    sequence = dc->peptideSeqs[frame - 1];
  
  /* Create the widget that will render the sequence data. */
  GtkWidget *refSeqWidget = gtk_drawing_area_new();
  gtk_widget_set_size_request(refSeqWidget, -1, roundNearest(dc->charHeight));
  *refSeqList = g_slist_append(*refSeqList, refSeqWidget);
  
  sequenceCreateProperties(refSeqWidget, dc->refSeqName, sequence, dc->refSeqType, strand, 
                           frame, &properties->refDisplayRange, &dc->refSeqFullRange,
                           dc->hozScaleRev, matchSeqList, FALSE);

  gtk_table_attach(table, refSeqWidget, 2, 3, *row, *row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_widget_add_events(refSeqWidget, GDK_EXPOSURE_MASK);

  connectSequenceSignals(refSeqWidget, alignmentTool);
  
  *row += 1;
}


/* Create a drawing-area widget and add it to the given table in the
 * given row. Increments the row number when done. the given expose
 * function will draw the widget. */
                          
static void createDrawingAreaWidget(DotterContext *dc,
                                    GCallback exposeFunc,
                                    gpointer exposeFuncData,
                                    const int height,
                                    GtkTable *table,
                                    const int xpad,
                                    const int ypad,
                                    int *row)
{
  GtkWidget *widget = gtk_drawing_area_new();
  gtk_widget_set_size_request(widget, -1, height);
  gtk_table_attach(table, widget, 2, 3, *row, *row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  g_signal_connect(G_OBJECT(widget), "expose-event", exposeFunc, exposeFuncData);

  *row += 1;
}



/* Create the section of the alignment tool for the given strand */
static GtkWidget* createAlignmentToolSection(BlxStrand strand,
                                             GtkWidget *alignmentTool,
                                             AlignmentToolProperties *properties)
{
  /* We'll pack everything into a table */
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  const int numRows = dc->numFrames + 3;  /* 2 rows for the headers, 1 for the sseq and 1 for each ref seq frame */
  const int numCols = 2;              /* 2 columns: label and sequence */
  const int xpad = 2;
  const int ypad = 0;
  const char strandChar = (strand == BLXSTRAND_FORWARD ? '+' : '-');
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  int row = 0;

  /* Create the widget that will render the match sequence. Need to do this first so we can pass
   * it as data to the ref sequence widgets. */
  GtkWidget *matchSeqWidget = gtk_drawing_area_new();
  gtk_widget_set_size_request(matchSeqWidget, -1, roundNearest(dc->charHeight));
  
  /* We need to create a list of match seqs (just one entry) and ref seqs to pass as data to the
   * ref seq and match seq respectively. */
  GSList *matchSeqList = NULL;
  matchSeqList = g_slist_append(matchSeqList, matchSeqWidget);
  
  GSList *refSeqList = NULL;
  
  /* Create a header line for this strand of the reference sequence. This widget will display
   * the currently-selected ref seq coord. */
  createDrawingAreaWidget(dc, G_CALLBACK(onExposeRefSequenceHeader), alignmentTool,
                          roundNearest(dc->charHeight) + SELECTED_COORD_MARKER_HEIGHT,
                          table, xpad, ypad, &row);
    
  /* Add a line for each frame of the reference sequence */
  int frame = 1;
  for (frame = 1; frame <= dc->numFrames; ++frame)
    {
      createRefSeqWidget(alignmentTool, properties, frame, strand, strandChar, matchSeqList, table, xpad, ypad, &row, &refSeqList);
    }
  
  /* Add a label for the match sequence */
  char *text = blxprintf("%s:", dc->matchSeqName);
  GtkWidget *label = createLabel(text, 0.0, 0.0, TRUE, TRUE);
  gtk_widget_modify_font(label, dc->fontDesc);
  gtk_misc_set_alignment(GTK_MISC(label), 0, 0);
  g_free(text);
  gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  /* Now add the data to the match-sequence widget and add it to the bottom row of the table */
  char *matchSequence = dc->matchSeqStrand == BLXSTRAND_REVERSE ? dc->matchSeqRev : dc->matchSeq;
  
  sequenceCreateProperties(matchSeqWidget, dc->matchSeqName, matchSequence, dc->matchSeqType, dc->matchSeqStrand,
                           1, &properties->matchDisplayRange, &dc->matchSeqFullRange, dc->vertScaleRev, refSeqList, TRUE);
  
  gtk_table_attach(table, matchSeqWidget, 2, 3, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_widget_add_events(matchSeqWidget, GDK_EXPOSURE_MASK);

  connectSequenceSignals(matchSeqWidget, alignmentTool);
    
  ++row;
  
  /* Create the bottom label. This will display the currently-selected match seq coord.  */
  createDrawingAreaWidget(dc, G_CALLBACK(onExposeMatchSequenceHeader), alignmentTool,
                          roundNearest(dc->charHeight) + SELECTED_COORD_MARKER_HEIGHT + roundNearest(dc->charHeight), /* extra charheight for spacing */
                          table, xpad, ypad, &row);
  
  return GTK_WIDGET(table);
}


GtkWidget* createAlignmentTool(DotterWindowContext *dotterWinCtx)
{
  DEBUG_ENTER("createAlignmentTool");

  GtkWidget *alignmentTool = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(alignmentTool), "Dotter - Alignment Tool");

  gtk_window_set_default_size(GTK_WINDOW(alignmentTool), 1160, -1);
  
  /* We'll put everything in a vbox, inside a scrolled window, inside a frame */  
  GtkWidget *frame = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
  gtk_container_add(GTK_CONTAINER(alignmentTool), frame);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
    
  /* Remember the headers so we can store them in the properties. We'll need to update them
   * when the range updates since they contain the centre coord of the current display range. */
  alignmentToolCreateProperties(alignmentTool, dotterWinCtx);
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* Put the forward ref seq strand on top, unless the display is reversed */
  BlxStrand qStrand = dc->hozScaleRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
  GtkWidget *section1 = createAlignmentToolSection(qStrand, alignmentTool, properties);
  gtk_box_pack_start(GTK_BOX(vbox), section1, FALSE, FALSE, 0);

  /* Only show the other reference sequence strand if we're displaying sequences as nucleotides */
  if (dc->displaySeqType == BLXSEQ_DNA)
    {
      qStrand = (qStrand == BLXSTRAND_FORWARD ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
      GtkWidget *section2 = createAlignmentToolSection(qStrand, alignmentTool, properties);
      gtk_box_pack_start(GTK_BOX(vbox), section2, FALSE, FALSE, 0);
    }

  /* Create the right-click menu */
  GtkWidget *menu = createAlignmentToolMenu(alignmentTool, &properties->actionGroup);
  
  gtk_widget_add_events(alignmentTool, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(alignmentTool, GDK_KEY_PRESS_MASK);

  g_signal_connect(G_OBJECT(alignmentTool), "button-press-event", G_CALLBACK(onButtonPressAlignmentTool), menu);
  g_signal_connect(G_OBJECT(alignmentTool), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
  g_signal_connect(G_OBJECT(alignmentTool), "size-allocate", G_CALLBACK(onSizeAllocateAlignmentTool), NULL);
  gtk_widget_show_all(alignmentTool);
  
  onAlignmentToolRangeChanged(alignmentTool);
  
  DEBUG_EXIT("createAlignmentTool returning ");
  return alignmentTool;
}


/***********************************************************
 *                           Drawing                       *
 ***********************************************************/

/* Draw the sequence data for the given sequence-widget */
static void drawSequence(GdkDrawable *drawable, GtkWidget *widget, GtkWidget *alignmentTool)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
  DotterWindowContext *dwc = atProperties->dotterWinCtx;
  DotterContext *dc = atProperties->dotterWinCtx->dotterCtx;

  GdkGC *gc = gdk_gc_new(drawable);

  /* Get the sequence info for this widget */
  SequenceProperties *seq1 = sequenceGetProperties(widget);
  const int seq1Len = strlen(seq1->sequence);
  
  /* If text on this widget is selected, then we'll highlight any bases within
   * the selection range. */
  const gboolean highlight = (atProperties->selectionWidget == widget);

  /* Get the sequence coord at the start of the display (leftmost edge) */
  int seq1Start = getDisplayStart(seq1, dc);
  const int seq1Offset = getSequenceOffset(seq1, dc);
  
  /* Loop through each display coord (0-based from left edge of display) */
  int displayIdx = 0;
  char seq1Text[atProperties->alignmentLen + 1];

  for ( ; displayIdx <= atProperties->alignmentLen; ++displayIdx)
    {
      seq1Text[displayIdx] = ' ';
    
      /* Get the 0-based index into sequence 1 and extract the character at this index */
      const int seq1Idx = displayIdx - seq1Offset;
      
      if (seq1Idx >= 0 && seq1Idx < seq1Len)
        {
          seq1Text[displayIdx] = seq1->sequence[seq1Idx];
          highlightSequenceBase(seq1, atProperties, dwc, displayIdx, seq1Idx, seq1Start, highlight, gc, drawable);
        }
    }
  
  /* terminate the display text string */
  seq1Text[displayIdx] = '\0';
  
  /* Draw the sequence text. */
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, seq1Text);
  pango_layout_set_font_description(layout, dc->fontDesc);
  
  if (layout)
    {
      gtk_paint_layout(widget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, widget, NULL, 0, 0, layout);
      g_object_unref(layout);
    }
  
  g_object_unref(gc);
}


/* Utility called by drawSequenceHeader to draw a marker line at the given coords */
static void drawSequenceHeaderMarker(GdkDrawable *drawable, const int x, const int y, const gdouble charWidth)
{
  /* Draw a marker line below where the text will go */
  const int x1 = x + roundNearest((charWidth / 2.0));
  
  GdkGC *gc = gdk_gc_new(drawable);
  gdk_draw_line(drawable, gc, x1, y, x1, y +  + SELECTED_COORD_MARKER_HEIGHT);
  g_object_unref(gc);
}


/* Utility called by drawSequenceHeader to draw the header text at the given coords */
static void drawSequenceHeaderText(GtkWidget *widget, 
                                   GdkDrawable *drawable, 
                                   const int x, 
                                   const int y, 
                                   const int coordIn, 
                                   DotterContext *dc,
                                   const gboolean horizontal)
{
  int coord = getDisplayCoord(coordIn, dc, horizontal);
  char *displayText = convertIntToString(coord);
  
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, displayText);
  pango_layout_set_font_description(layout, dc->fontDesc);
  
  if (layout)
    {
      /* Offset the text so that the middle of the text is lined up with the coord of interest */
      const int offset = ceil((((gdouble)numDigitsInInt(coord) / 2.0) - 1) * dc->charWidth);
    
      gtk_paint_layout(widget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, widget, NULL, x - offset, y, layout);
      g_object_unref(layout);
    }
  
  g_free(displayText);
  
}


/* Draw a sequence header that displays the centre coord of the given display range. If
 * markerFirst is true the marker is drawn above the text, otherwise below. */
static void drawSequenceHeader(GtkWidget *widget, 
                               GtkWidget *alignmentTool,
                               GdkDrawable *drawable, 
                               const gboolean horizontal)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = dwc->dotterCtx;
  
  /* Find the coord to display */
  const int coord = horizontal ? dwc->refCoord : dwc->matchCoord;

  const IntRange const *displayRange = horizontal ? &properties->refDisplayRange : &properties->matchDisplayRange;
  
  /* Find the position to display at. Find the position of the char at this coord */
  const int displayIdx = convertToDisplayIdx(coord - displayRange->min, horizontal, dc, 1, NULL) - 1;
  int x = (int)((gdouble)displayIdx * dc->charWidth);
  int y = 0;

  if (horizontal)
    {
      /* For the ref sequence, draw the marker below the label. */
      drawSequenceHeaderText(widget, drawable, x, y, coord, dc, horizontal);
      y += roundNearest(dc->charHeight);
      drawSequenceHeaderMarker(drawable, x, y, dc->charWidth);
    }
  else
    {
      /* For the match sequence, draw the marker above the label */
      drawSequenceHeaderMarker(drawable, x, y, dc->charWidth);
      y += SELECTED_COORD_MARKER_HEIGHT;
      drawSequenceHeaderText(widget, drawable, x, y, coord, dc, horizontal);
    }
}



/***********************************************************
 *                      Utility functions                  *
 ***********************************************************/

/* Return the coord that the alignment tool display starts at (adjusted for strand direction
 * and converted from dna to peptide coords if applicable). */
static int getDisplayStart(SequenceProperties *properties, DotterContext *dc)
{
  const gboolean forward = (properties->strand == BLXSTRAND_FORWARD);
  int displayStart = forward ? properties->displayRange->min : properties->displayRange->max;

  displayStart = convertToDisplayIdx(displayStart, !properties->isMatchSeq, dc, properties->frame, NULL);

  return displayStart;
}


/* Return the coord that the alignment tool display ends at (adjusted for strand direction
 * and converted from dna to peptide coords if applicable). */
static int getDisplayEnd(SequenceProperties *properties, DotterContext *dc)
{
  const gboolean forward = (properties->strand == BLXSTRAND_FORWARD);
  int displayEnd = forward ? properties->displayRange->max : properties->displayRange->min;

  displayEnd = convertToDisplayIdx(displayEnd, !properties->isMatchSeq, dc, properties->frame, NULL);

  return displayEnd;
}


/* Calculate how many bases into the alignment tool's display range the given sequence starts.
 * Also return the coord that the alignment tool display starts at (adjusted for strand direction
 * and converted from dna to peptide coords if applicable). */
static int getSequenceOffset(SequenceProperties *properties, DotterContext *dc)
{
  const gboolean forward = (properties->strand == BLXSTRAND_FORWARD);

  /* Get the start coord of the sequence data */
  int seqStart = forward ? properties->fullRange->min : properties->fullRange->max;
  
  /* Offset this coord so that we have the first coord that starts our reading frame */
  int reqdFrame = properties->frame;
  int frame = UNSET_INT;
  convertToDisplayIdx(seqStart, !properties->isMatchSeq, dc, 1, &frame);
  
  int offset = reqdFrame - frame;
  if (offset < 0)
    {
      offset += dc->numFrames;
    }

  seqStart = forward ? seqStart + offset : seqStart - offset;
  
  /* Convert to display coords. */
  seqStart = convertToDisplayIdx(seqStart, !properties->isMatchSeq, dc, reqdFrame, NULL);
  
  /* Calculate the offset. Note that this may be positive or negative because it
   * can be in either direction. */
  const int displayStart = getDisplayStart(properties, dc);
  int result = forward ? seqStart - displayStart - 1 : displayStart - seqStart;

  return result;
}


/* Get the sequence coordinate at the given x position within the given sequence widget */
static int getCoordAtPos(const int x, GtkWidget *sequenceWidget, GtkWidget *alignmentTool)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
  SequenceProperties *properties = sequenceGetProperties(sequenceWidget);
  DotterContext *dc = atProperties->dotterWinCtx->dotterCtx;

  /* Get the 0-based character index from the leftmost edge of the widget */
  const int displayIdx = (int)((double)x / (double)dc->charWidth);

  /* Get the coord of the first (leftmost) character in the display */
  int displayStart = getDisplayStart(properties, dc);

  /* Apply the offset to get the 0-based coord into the sequence, and add the range
   * min to get the real coord */
  int coord = 0;
  if (properties->strand == BLXSTRAND_REVERSE)
    coord = displayStart - displayIdx;
  else
    coord = displayStart + displayIdx + 1;

  /* Bounds-limit the result */
  boundsLimitValue(&coord, properties->fullRange);
  
  return coord;
}


/* Connect signals for a sequence widget */
static void connectSequenceSignals(GtkWidget *widget, GtkWidget *alignmentTool)
{
  gtk_widget_add_events(widget, GDK_EXPOSURE_MASK);
  gtk_widget_add_events(widget, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(widget, GDK_BUTTON_RELEASE_MASK);
    gtk_widget_add_events(widget, GDK_BUTTON_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(widget), "expose-event", G_CALLBACK(onExposeSequence), alignmentTool);
  g_signal_connect(G_OBJECT(widget), "button-press-event", G_CALLBACK(onButtonPressSequence), alignmentTool);
  g_signal_connect(G_OBJECT(widget), "button-release-event", G_CALLBACK(onButtonReleaseSequence), alignmentTool);
  g_signal_connect(G_OBJECT(widget), "motion-notify-event", G_CALLBACK(onMouseMoveSequence), alignmentTool);
}


/* Get the bit of sequence from the given sequence widget between
 * the given coords. Returns a new string which should be freed by
 * the caller using g_free. */
static char *getSequenceBetweenCoords(GtkWidget *sequenceWidget,
                                      const int startCoord,
                                      const int endCoord,
                                      DotterWindowContext *dwc)
{
  char *result = NULL;
  
  SequenceProperties *properties = sequenceGetProperties(sequenceWidget);
  
  int startIdx = startCoord - properties->fullRange->min;
  const int numChars = endCoord - startCoord + 1;

  //  if (properties->strand == BLXSTRAND_REVERSE)
  //startIdx = properties->fullRange->max - endCoord;
    
  if (properties->strand == BLXSTRAND_REVERSE)
    {
      /* Sequence is stored rev-comp'd but supplied coords are
       * on the forward strand, so invert the coords */
      startIdx = properties->fullRange->max - endCoord;
    }

  result = g_strndup(properties->sequence + startIdx, numChars + 1);
  result[numChars] = '\0';

  return result;
}


/* Set the widget that contains the current sequence selection. */
static void setSelectionWidget(AlignmentToolProperties *atProperties, GtkWidget *selectionWidget)
{
  atProperties->selectionWidget = selectionWidget;
  enableMenuAction(atProperties->actionGroup, "CopySelnCoords", selectionWidget != NULL);
  enableMenuAction(atProperties->actionGroup, "CopySeln", selectionWidget != NULL);
  enableMenuAction(atProperties->actionGroup, "ClearSeln", selectionWidget != NULL);
}


/* This clears the current selection */
static void clearSequenceSelection(GtkWidget *alignmentTool)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);

  if (atProperties->selectionWidget)
    {
      /* clear current selection widget and redraw it */
      GtkWidget *widget = atProperties->selectionWidget;
      setSelectionWidget(atProperties, NULL);
      widgetClearCachedDrawable(widget, NULL);
      gtk_widget_queue_draw(widget);
    }
}


/* Initiate dragging on a sequence widget */
static void sequenceInitiateDragging(GtkWidget *sequenceWidget, GtkWidget *alignmentTool, const int x)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);

  /* flag that we're dragging */
  atProperties->dragging = TRUE;
  atProperties->dragStart = getCoordAtPos(x, sequenceWidget, alignmentTool);

  /* flag that we should highlight the selected sequence (clear any current selection first) */
  clearSequenceSelection(alignmentTool);
  
  atProperties->selectionRange.min = atProperties->dragStart;
  atProperties->selectionRange.max = atProperties->dragStart;

  setSelectionWidget(atProperties, sequenceWidget);
}


/* Finish dragging on a sequence widget */
static void sequenceFinishDragging(GtkWidget *sequenceWidget, GtkWidget *alignmentTool, const int x)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
  DotterWindowContext *dwc = atProperties->dotterWinCtx;

  if (atProperties->dragging)
    {
      /* cancel the dragging and highlighting flags */
      atProperties->dragging = FALSE;
  
      /* get the range of coords that the user dragged over */
      int minCoord = atProperties->dragStart;
      int maxCoord = getCoordAtPos(x, sequenceWidget, alignmentTool);

      if (minCoord > maxCoord)
        {
          int tmp = minCoord;
          minCoord = maxCoord;
          maxCoord = tmp;
        }
      
      /* Get the sequence between these coords and place it on the clipboard */
      char *result = getSequenceBetweenCoords(sequenceWidget, minCoord, maxCoord, dwc);
      setPrimaryClipboardText(result);
      
      g_free(result);
    }
}


/* Select all of the visible text in a sequence widget */
static void selectVisibleSequence(GtkWidget *sequenceWidget, GtkWidget *alignmentTool)
{
  SequenceProperties *properties = sequenceGetProperties(sequenceWidget);
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
  DotterContext *dc = atProperties->dotterWinCtx->dotterCtx;
  
  /* cancel any dragging operation */
  atProperties->dragging = FALSE;

  /* flag that we should highlight the selected sequence (clear any current
   * selection first) and set the coords of the highlighted bit */
  clearSequenceSelection(alignmentTool);

  int start = getDisplayStart(properties, dc);
  int end = getDisplayEnd(properties, dc);

  if (properties->strand == BLXSTRAND_FORWARD)
    {
      ++start;
      ++end;
    }
  
  boundsLimitValue(&start, properties->fullRange);
  boundsLimitValue(&end, properties->fullRange);
  
  intrangeSetValues(&atProperties->selectionRange, start, end);
  setSelectionWidget(atProperties, sequenceWidget);

  /* copy the selection to the primary clipboard */
  char *result = getSequenceBetweenCoords(sequenceWidget,
                                          atProperties->selectionRange.min,
                                          atProperties->selectionRange.max,
                                          atProperties->dotterWinCtx);
  
  setPrimaryClipboardText(result);
  g_free(result);
}


/* Get the colour that the given base in the given sequence should
 * be highlighted according to how well it matches the other sequences
 * at the same display index. */
static DotterColorId getBaseHighlightColor(SequenceProperties *seq1,
                                           DotterContext *dc,
                                           const int displayIdx,
                                           const int seq1Idx,
                                           const int seq1DisplayStart)
{
  /* Compare the base in the current sequence to the base at the
   * same display position in all other sequences and find whether
   * any match. Exact matches take precedence over conserved matches.
   * If none match, just use the background colour (i.e. no highlighting). */
  DotterColorId colorId = DOTCOLOR_BACKGROUND;
 
  GSList *item = seq1->compSeqs;
 
  for ( ; item; item = item->next)
    {
      /* Get info about the other sequence */
      GtkWidget *seq2Widget = GTK_WIDGET(item->data);
      SequenceProperties *seq2 = sequenceGetProperties(seq2Widget);
      const int seq2Len = strlen(seq2->sequence);

      /* Get the sequence coord at the start of the display (leftmost edge) */
      const int seq2Offset = getSequenceOffset(seq2, dc);
    
      /* Get the zero-based index into the sequence and compare the bases to determine the highlight color */
      const int seq2Idx = displayIdx - seq2Offset;
      
      if (seq2Idx >= 0 && seq2Idx < seq2Len)
        {
          if (seq1->sequence[seq1Idx] == seq2->sequence[seq2Idx])
            {
              colorId = DOTCOLOR_MATCH;
              break;
            }
          else if (dc->blastMode != BLXMODE_BLASTN && 
                   dc->matrix[atob[(unsigned int)seq1->sequence[seq1Idx]] - 1 ][atob[(unsigned int)seq2->sequence[seq2Idx]] - 1 ] > 0)
            {
              colorId = DOTCOLOR_CONS;
            }
        }
    }

  return colorId;
}


static void highlightSequenceBase(SequenceProperties *seq1,
                                  AlignmentToolProperties *atProperties,
                                  DotterWindowContext *dwc,
                                  const int displayIdx,
                                  const int seq1Idx,
                                  const int seq1Start,
                                  const gboolean highlight,
                                  GdkGC *gc,
                                  GdkDrawable *drawable)
{
  DotterContext *dc = dwc->dotterCtx;
  
  DotterColorId colorId = getBaseHighlightColor(seq1, dc, displayIdx, seq1Idx, seq1Start);
  const int seq1Coord = (seq1->strand == BLXSTRAND_REVERSE ? seq1Start - displayIdx : seq1Start + displayIdx + 1);
  const gboolean selected = highlight && valueWithinRange(seq1Coord, &atProperties->selectionRange);
  
  /* We don't need to bother drawing the background if it's the standard
   * (non-selected) background color */
  if (colorId != DOTCOLOR_BACKGROUND || selected)
    {
      GdkColor *color = getGdkColor(colorId, dc->defaultColors, selected, dwc->usePrintColors);
      gdk_gc_set_foreground(gc, color);

      const int x = (int)((gdouble)displayIdx * dc->charWidth);
      gdk_draw_rectangle(drawable, gc, TRUE, x, 0, ceil(dc->charWidth), roundNearest(dc->charHeight));
    }
}

