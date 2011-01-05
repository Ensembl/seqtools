/*  File: alignmenttool.c
 *  Author: Gemma Barson, 2010-09-02
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


static int atob[]	/* OLD (starting at 1) ASCII-to-binary translation table  (Inherited from blast) */
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
#define	SELECTED_COORD_MARKER_HEIGHT		    12	  /* the height of the marker that indicates the currently-selected coord */


typedef struct _SequenceProperties
{
  const char *seqName;
  const char *sequence;
  BlxSeqType seqType;               /* whether this sequence is in nucleotide or peptide coords */
  BlxStrand strand;
  int frame;
  int startCoordOffset;             /* how much to offset the first coord in the seq by to give the first base in the first codon in frame 1 */
  IntRange *displayRange;	    /* the displayed range; pointer to refDisplayRange or matchDisplayRange */
  IntRange *fullRange;              /* the full range of the sequence; pointer to refSeqRange or matchSeqRange */
  gboolean scaleReversed;	    /* whether the scale for this sequence is shown reversed (i.e. high-to-low rather than low-to-high) */
  GSList *compSeqs;                 /* list of other sequence widgets that this sequence will be compared against */
  gboolean isMatchSeq;		    /* true if this is the match seq, false if it's the ref seq */
} SequenceProperties;


typedef struct _AlignmentToolProperties
{
  int alignmentLen;                 /* the number of coords wide the alignment too displays */
  IntRange refDisplayRange;	    /* the current ref seq range displayed in the tool */
  IntRange matchDisplayRange;	    /* the current match seq range displayed in the tool */
  
  DotterWindowContext *dotterWinCtx;
  
  GtkWidget *refSeqHeader1;
  GtkWidget *matchSeqHeader1;
  GtkWidget *refSeqHeader2;
  GtkWidget *matchSeqHeader2;
} AlignmentToolProperties;



/* Local function declarations */
static void                       onCloseMenu(GtkAction *action, gpointer data);
static void                       onPrintMenu(GtkAction *action, gpointer data);
static void                       onPrintColorsMenu(GtkAction *action, gpointer data);
static void                       onSetLengthMenu(GtkAction *action, gpointer data);
//static const gchar*               getSequence(GtkWidget *widget, GtkWidget *alignmentTool);
static void                       drawSequence(GdkDrawable *drawable, GtkWidget *widget, GtkWidget *alignmentTool);
static void			  drawSequenceHeader(GtkWidget *widget, GdkDrawable *drawable, DotterContext *dc, const IntRange const *displayRange, const gboolean horizontal);
static int                        convertToDisplayIdx(const int dnaIdx, const gboolean horizontal, DotterContext *dc);


/* Menu builders */
static const GtkActionEntry alignmentToolMenuEntries[] = {
{ "Close",        NULL, "_Close",                   NULL,	"Close the alignment tool",             G_CALLBACK(onCloseMenu)},
{ "Print",        NULL, "_Print",                   NULL,	"Print",                                G_CALLBACK(onPrintMenu)},
{ "PrintColors",  NULL, "Print _colors",            NULL,	"Toggle colors for printing",           G_CALLBACK(onPrintColorsMenu)},
{ "SetLength",    NULL, "_Set alignment length",    NULL,	"Set the length of the alignment tool", G_CALLBACK(onSetLengthMenu)}
};

/* This defines the layout of the menu */
static const char alignmentToolMenuDescription[] =
"<ui>"
"  <popup name='MainMenu'>"
"      <menuitem action='Close'/>"
"      <menuitem action='Print'/>"
"      <menuitem action='PrintColors'/>"
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
                                     const int startCoordOffset,
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
      properties->startCoordOffset = startCoordOffset;
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

static void alignmentToolCreateProperties(GtkWidget *widget,
					  DotterWindowContext *dotterWinCtx,
					  GtkWidget *refSeqHeader1,
					  GtkWidget *matchSeqHeader1,
					  GtkWidget *refSeqHeader2,
					  GtkWidget *matchSeqHeader2)
{
  if (widget)
    {
      AlignmentToolProperties *properties = g_malloc(sizeof *properties);
    
      properties->dotterWinCtx = dotterWinCtx;
      properties->alignmentLen = DEFAULT_ALIGNMENT_LENGTH;
      properties->refDisplayRange.min = 0;
      properties->refDisplayRange.max = 20;
      properties->matchDisplayRange.min = 0;
      properties->matchDisplayRange.max = 20; /* to do: put real values in here for the ranges */
      properties->refSeqHeader1 = refSeqHeader1;
      properties->matchSeqHeader1 = matchSeqHeader1;
      properties->refSeqHeader2 = refSeqHeader2;
      properties->matchSeqHeader2 = matchSeqHeader2;

    g_object_set_data(G_OBJECT(widget), "AlignmentToolProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyAlignmentTool), NULL); 
    }
}


/***********************************************************
 *                       Events				   *
 ***********************************************************/

/* Expose function for a widget containing a section of sequence. SequenceProperties should
 * be set on any widget that this is to be called for. */
static gboolean onExposeSequence(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkWindow *drawable = widget->window;

  if (drawable)
    {
      drawSequence(drawable, widget, alignmentTool);
    }
  
  return TRUE;
}


/* Expose function for a widget containing a reference sequence header. */
static gboolean onExposeRefSequenceHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkWindow *drawable = widget->window;
  
  if (drawable)
    {
      AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
      drawSequenceHeader(widget, drawable, properties->dotterWinCtx->dotterCtx, &properties->refDisplayRange, TRUE);
    }
  
  return TRUE;
}


/* Expose function for a widget containing a match sequence header. */
static gboolean onExposeMatchSequenceHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *alignmentTool = GTK_WIDGET(data);
  GdkWindow *drawable = widget->window;
  
  if (drawable)
    {
      AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
      drawSequenceHeader(widget, drawable, properties->dotterWinCtx->dotterCtx, &properties->matchDisplayRange, FALSE);
    }
  
  return TRUE;
}


/* This should be called when the range displayed by the alignment tool has changed */
static void onAlignmentToolRangeChanged(GtkWidget *alignmentTool)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);

  if (properties)
    {
      char *refCoordText = blxprintf("%d", getRangeCentre(&properties->refDisplayRange));
      char *matchCoordText = blxprintf("%d", getRangeCentre(&properties->matchDisplayRange));
    
      if (properties->refSeqHeader1)
        gtk_widget_queue_draw(properties->refSeqHeader1);
      
      if (properties->refSeqHeader2)
        gtk_widget_queue_draw(properties->refSeqHeader2);
      
      if (properties->matchSeqHeader1)
        gtk_widget_queue_draw(properties->matchSeqHeader1);

      if (properties->matchSeqHeader2)
        gtk_widget_queue_draw(properties->matchSeqHeader2);

      g_free(refCoordText);
      g_free(matchCoordText);
    }
}


/* Create the menu */
static GtkWidget* createAlignmentToolMenu(GtkWidget *window)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  gtk_action_group_add_actions (action_group, alignmentToolMenuEntries, G_N_ELEMENTS (alignmentToolMenuEntries), window);
  
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
  
  return gtk_ui_manager_get_widget (ui_manager, "/MainMenu");
}


/* Mouse button handler */
static gboolean onButtonPressAlignmentTool(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
    }
  
  return TRUE;
}


/***********************************************************
 *                       Get/set functions                 *
 ***********************************************************/

/* Set the alignment range centre based on the currently-selected match/ref seq coord. */
void updateAlignmentRange(GtkWidget *alignmentTool, DotterWindowContext *dwc)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);

  /* Re-create the range, centred on the set coordinate and with the alignment tool's alignment 
   * length. Note that the length is the number of display chars but the reference sequence is
   * in nucleotide coords so needs converting. (The match sequence is always in display coords so
   * will be correct whether we're displaying nucleotides or peptides.) */
  centreRangeOnCoord(&properties->refDisplayRange, dwc->refCoord, properties->alignmentLen * getResFactor(dwc->dotterCtx, TRUE));
  centreRangeOnCoord(&properties->matchDisplayRange, dwc->matchCoord, properties->alignmentLen * getResFactor(dwc->dotterCtx, FALSE));
  
//  boundsLimitRange(&properties->refDisplayRange, &dc->refSeqRange);
//  boundsLimitRange(&properties->matchDisplayRange, &dc->matchSeqRange);

  onAlignmentToolRangeChanged(alignmentTool);
  gtk_widget_queue_draw(alignmentTool);
}


/* Set the alignment length */
static void setAlignmentLength(GtkWidget *alignmentTool, const int length, GError **error)
{
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  
  if (length > MIN_ALIGNMENT_LENGTH && length <= MAX_ALIGNMENT_LENGTH)
    {
      /* Input length is in display coords: convert to nucleotide coords */
      properties->alignmentLen = length;
    }
  else
    {
      g_set_error(error, DOTTER_ERROR, DOTTER_ERROR_INVALID_LENGTH, "Length %d is not in valid range %d -> %d.\n",
                  length, MIN_ALIGNMENT_LENGTH, MAX_ALIGNMENT_LENGTH);
    }

  updateAlignmentRange(alignmentTool, properties->dotterWinCtx);
  
  gtk_widget_queue_draw(alignmentTool);
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


static void onPrintMenu(GtkAction *action, gpointer data)
{
}


static void onPrintColorsMenu(GtkAction *action, gpointer data)
{
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


/***********************************************************
 *                       Initialisation                    *
 ***********************************************************/

/* Create the section of the alignment tool for the given strand */
static GtkWidget* createAlignmentToolSection(BlxStrand strand,
					     DotterWindowContext *dwc,
					     GtkWidget **refSeqHeader,
					     GtkWidget **matchSeqHeader,
					     GtkWidget *alignmentTool,
					     AlignmentToolProperties *properties)
{
  /* We'll pack everything into a table */
  DotterContext *dc = dwc->dotterCtx;
  const int numRows = dc->numFrames + 3;  /* 2 rows for the headers, 1 for the sseq and 1 for each ref seq frame */
  const int numCols = 2;	      /* 2 columns: label and sequence */
  const int xpad = 2;
  const int ypad = 0;
  const char strandChar = (strand == BLXSTRAND_FORWARD ? '+' : '-');
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  int row = 0;

  /* Calculate the reading frame of the first and last coord in the reference sequence. */
  const int resFactor = getResFactor(dc, TRUE);
  const int reqdFrame = (dc->hozScaleRev ? resFactor : 1);

  const int startCoord = (dc->hozScaleRev ? dc->refSeqFullRange.max - (resFactor - 1) : dc->refSeqFullRange.min);
  int startCoordFrame = startCoord % resFactor;
  if (startCoordFrame < 1)
    {
      startCoordFrame += resFactor;
    }

  int startCoordOffset = reqdFrame - startCoordFrame;
  if (startCoordOffset < 0)
    {
      startCoordOffset += resFactor;
    }
  else if (startCoordOffset > resFactor - 1)
    {
      startCoordOffset -= resFactor;
    }
  
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
  *refSeqHeader = gtk_drawing_area_new();
  gtk_widget_set_size_request(*refSeqHeader, -1, roundNearest(dc->charHeight) + SELECTED_COORD_MARKER_HEIGHT);
  gtk_table_attach(table, *refSeqHeader, 2, 3, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  g_signal_connect(G_OBJECT(*refSeqHeader), "expose-event", G_CALLBACK(onExposeRefSequenceHeader), alignmentTool);
  ++row;
  
  /* Add a line for each frame of the reference sequence */
  int frame = 1;
  for (frame = 1; frame <= dc->numFrames; ++frame)
    {
      /* Create a label containing the name */
      char *text = blxprintf("%s (%c%d):", dc->refSeqName, strandChar, frame);
      GtkWidget *label = gtk_label_new(text);
      gtk_widget_modify_font(label, dc->fontDesc);
      gtk_misc_set_alignment(GTK_MISC(label), 0, 0);
      g_free(text);
      gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
      
      /* Find which sequence is applicable (fwd or rev strand or translated peptide seq) */
      char *sequence = dc->refSeq;
      
      if (strand == BLXSTRAND_REVERSE && dc->displaySeqType == BLXSEQ_DNA)
        sequence = dc->refSeqRev;
      else if (dc->blastMode == BLXMODE_BLASTX)
        sequence = dc->peptideSeqs[frame - 1];
      
      /* Create the widget that will render the sequence data. */
      GtkWidget *refSeqWidget = gtk_drawing_area_new();
      gtk_widget_set_size_request(refSeqWidget, -1, roundNearest(dc->charHeight));
      refSeqList = g_slist_append(refSeqList, refSeqWidget);

      sequenceCreateProperties(refSeqWidget, dc->refSeqName, sequence, dc->refSeqType, strand, 
                               frame, startCoordOffset, &properties->refDisplayRange, &dc->refSeqFullRange,
                               dc->hozScaleRev, matchSeqList, FALSE);
      
      gtk_table_attach(table, refSeqWidget, 2, 3, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_widget_add_events(refSeqWidget, GDK_EXPOSURE_MASK);
      g_signal_connect(G_OBJECT(refSeqWidget), "expose-event", G_CALLBACK(onExposeSequence), alignmentTool);
      
      ++row;
    }
  
  /* Add a label for the match sequence */
  char *text = blxprintf("%s:", dc->matchSeqName);
  GtkWidget *label = gtk_label_new(text);
  gtk_widget_modify_font(label, dc->fontDesc);
  gtk_misc_set_alignment(GTK_MISC(label), 0, 0);
  g_free(text);
  gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  /* Now add the data to the match-sequence widget and add it to the bottom row of the table */
  char *matchSequence = dc->matchSeqStrand == BLXSTRAND_REVERSE ? dc->matchSeqRev : dc->matchSeq;
  
  sequenceCreateProperties(matchSeqWidget, dc->matchSeqName, matchSequence, dc->matchSeqType, dc->matchSeqStrand,
                           1, 0, &properties->matchDisplayRange, &dc->matchSeqFullRange, dc->vertScaleRev, refSeqList, TRUE);
  
  gtk_table_attach(table, matchSeqWidget, 2, 3, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  gtk_widget_add_events(matchSeqWidget, GDK_EXPOSURE_MASK);
  g_signal_connect(G_OBJECT(matchSeqWidget), "expose-event", G_CALLBACK(onExposeSequence), alignmentTool);
  
  ++row;
  
  /* Create the bottom header */
  *matchSeqHeader = gtk_drawing_area_new();
  gtk_widget_set_size_request(*matchSeqHeader, -1, roundNearest(dc->charHeight) + SELECTED_COORD_MARKER_HEIGHT + roundNearest(dc->charHeight)); /* add an extra charheight just for spacing */
  gtk_table_attach(table, *matchSeqHeader, 2, 3, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
  g_signal_connect(G_OBJECT(*matchSeqHeader), "expose-event", G_CALLBACK(onExposeMatchSequenceHeader), alignmentTool);
  ++row;
  
  return GTK_WIDGET(table);
}


GtkWidget* createAlignmentTool(DotterWindowContext *dotterWinCtx)
{
  DEBUG_ENTER("createAlignmentTool");

  GtkWidget *alignmentTool = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(alignmentTool), "Dotter - Alignment Tool");

//  const int height = getHeightFromBlastMode(watsonOnly, blastMode);
  gtk_window_set_default_size(GTK_WINDOW(alignmentTool), 1160, -1);
  
  /* We'll put everything in a vbox, inside a scrolled window, inside a frame */  
  GtkWidget *frame = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
  gtk_container_add(GTK_CONTAINER(alignmentTool), frame);
  
//  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
//  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
//  gtk_container_add(GTK_CONTAINER(frame), scrollWin);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
//  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), vbox);
    
  /* Remember the headers so we can store them in the properties. We'll need to update them
   * when the range updates since they contain the centre coord of the current display range. */
  alignmentToolCreateProperties(alignmentTool, dotterWinCtx, NULL, NULL, NULL, NULL);
  AlignmentToolProperties *properties = alignmentToolGetProperties(alignmentTool);
  DotterContext *dc = properties->dotterWinCtx->dotterCtx;
  
  /* Put the forward ref seq strand on top, unless the display is reversed */
  BlxStrand qStrand = dc->hozScaleRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD; //dc->refSeqStrand;
  GtkWidget *section1 = createAlignmentToolSection(qStrand, dotterWinCtx, &properties->refSeqHeader1, &properties->matchSeqHeader1, alignmentTool, properties);
  gtk_box_pack_start(GTK_BOX(vbox), section1, FALSE, FALSE, 0);

  /* Only show the other reference sequence strand if we're displaying sequences as nucleotides */
  if (dc->displaySeqType == BLXSEQ_DNA)
    {
      qStrand = (qStrand == BLXSTRAND_FORWARD ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
      GtkWidget *section2 = createAlignmentToolSection(qStrand, dotterWinCtx, &properties->refSeqHeader2, &properties->matchSeqHeader2, alignmentTool, properties);
      gtk_box_pack_start(GTK_BOX(vbox), section2, FALSE, FALSE, 0);
    }

  /* Create the right-click menu */
  GtkWidget *menu = createAlignmentToolMenu(alignmentTool);

  gtk_widget_add_events(alignmentTool, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(alignmentTool), "button-press-event", G_CALLBACK(onButtonPressAlignmentTool), menu);
  g_signal_connect(G_OBJECT(alignmentTool), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
  gtk_widget_show_all(alignmentTool);
  
  onAlignmentToolRangeChanged(alignmentTool);
  
  DEBUG_EXIT("createAlignmentTool returning ");
  return alignmentTool;
}


/***********************************************************
 *                       Internal routines                 *
 ***********************************************************/

///* Get the visible segment of the sequence in the given sequence-widget */
//static const gchar* getSequence(GtkWidget *widget, GtkWidget *alignmentTool)
//{
//  SequenceProperties *properties = sequenceGetProperties(widget);
//  /* Find the segment of the ref seq to display. Offset the start by the given amount (to make sure we
//   * start at the first codon in frame 1) and by the frame number so we start at the beginning of
//   * the first codon in this particular frame. */
//  const int plusMin = properties->scaleReversed ? -1 : 1;
//  const int offset = properties->startCoordOffset + properties->frame - 1;
//  const int coord1 = properties->displayRange->min + plusMin * offset;
//  const int coord2 = properties->displayRange->max + plusMin * offset;
//  
//  IntRange qRange;
//  intrangeSetValues(&qRange, coord1, coord2);
//
//  /* If we're out of the valid sequence range, shift by a whole number of codons until we're in range */
//  const int resFactor = getResFactor(dc, !properties->isMatchSeq);
//  
//  while (qRange.min < properties->fullRange->min)
//    qRange.min += resFactor;
//  
//  while (qRange.max > properties->fullRange->max)
//    qRange.max -= resFactor;
//  
//  GError *error = NULL;
//
//  gchar *segmentToDisplay = getSequenceSegment(properties->sequence,
//                                               &qRange,
//                                               properties->strand, 
//                                               properties->seqType,                     /* input seq type */
//                                               dc->displaySeqType,                      /* result seq type */
//                                               properties->frame, 
//                                               dc->numFrames,
//                                               properties->fullRange,
//                                               dc->blastMode,
//                                               dc->geneticCode,
//                                               properties->scaleReversed,               /* whether the display is reversed */
//                                               properties->strand == BLXSTRAND_REVERSE, /* reverse the sequence if we want the reverse strand */
//                                               TRUE,                                    /* always complement if it's the reverse strand */
//                                               &error);
//
//  if (!segmentToDisplay)
//    {
//      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
//    }
//  
//  return segmentToDisplay;
//}


/* Draw the sequence data for the given sequence-widget */
static void drawSequence(GdkDrawable *drawable, GtkWidget *widget, GtkWidget *alignmentTool)
{
  AlignmentToolProperties *atProperties = alignmentToolGetProperties(alignmentTool);
  DotterContext *dc = atProperties->dotterWinCtx->dotterCtx;

  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *matchColor = getGdkColor(DOTCOLOR_MATCH, dc->defaultColors, FALSE, dc->usePrintColors);
  GdkColor *consColor = getGdkColor(DOTCOLOR_CONS, dc->defaultColors, FALSE, dc->usePrintColors);

  /* Get the sequence info for this widget */
  SequenceProperties *seq1Properties = sequenceGetProperties(widget);
  const gchar *seq1 = seq1Properties->sequence;
  const int seq1Len = strlen(seq1);
  const gboolean seq1Fwd = (seq1Properties->strand == BLXSTRAND_FORWARD);
  
  /* Get the position of the first valid coord in the display range, which might not be at
   * the start of the display if the sequence doesn't extend the full extent. */
  int seq1Start = seq1Fwd ? seq1Properties->fullRange->min : seq1Properties->fullRange->max;
  int seq1DisplayStart = seq1Fwd ? max(seq1Start, seq1Properties->displayRange->min) : min(seq1Start, seq1Properties->displayRange->max);
  int seq1Offset = seq1Fwd ? seq1DisplayStart - seq1Properties->displayRange->min : seq1Properties->displayRange->max - seq1DisplayStart;
  
  seq1Start = convertToDisplayIdx(seq1Start, !seq1Properties->isMatchSeq, dc);
  seq1DisplayStart = convertToDisplayIdx(seq1DisplayStart, !seq1Properties->isMatchSeq, dc);
  seq1Offset = convertToDisplayIdx(seq1Offset, !seq1Properties->isMatchSeq, dc);
  
  /* Loop through each display coord */
  int displayIdx = 0;
  char seq1Text[seq1Len];
  int y = 0;
  int x = 0;

  for ( ; displayIdx < atProperties->alignmentLen; ++displayIdx)
    {
      x = (int)((gdouble)displayIdx * dc->charWidth);
      
      /* Loop through the sequences to compare against and color the background of this seq according
       * to how well it matches. Exact matches take precedence over conserved matches. */
      gboolean match = FALSE;
      gboolean consMatch = FALSE;

      GSList *item = seq1Properties->compSeqs;
      
      for ( ; item; item = item->next)
        {
          /* Get info about the other sequence */
          GtkWidget *seq2Widget = GTK_WIDGET(item->data);
          SequenceProperties *seq2Properties = sequenceGetProperties(seq2Widget);
	  const gchar *seq2 = seq2Properties->sequence;
	  const int seq2Len = strlen(seq2);
	  const gboolean seq2Fwd = (seq2Properties->strand == BLXSTRAND_FORWARD);

          /* Get the position of the first coord in this sequences that is within the display range */
          int seq2Start = seq2Fwd ? seq2Properties->fullRange->min : seq2Properties->fullRange->max;
	  int seq2DisplayStart = seq2Fwd ? seq2Properties->displayRange->min : seq2Properties->displayRange->max;
	  boundsLimitValue(&seq2DisplayStart, seq2Properties->fullRange);
          int seq2Offset = seq2Fwd ? seq2DisplayStart - seq2Properties->displayRange->min : seq2Properties->displayRange->max - seq2DisplayStart;

          seq2Start = convertToDisplayIdx(seq2Start, !seq2Properties->isMatchSeq, dc);
  	  seq2DisplayStart = convertToDisplayIdx(seq2DisplayStart, !seq2Properties->isMatchSeq, dc);
	  seq2Offset = convertToDisplayIdx(seq2Offset, !seq2Properties->isMatchSeq, dc);

	  /* Get the real indexes into the sequences and compare the bases to determine the highlight color */
          const int seq1Idx = seq1Fwd ? seq1DisplayStart + (displayIdx - seq1Offset) - seq1Start : seq1Start - (seq1DisplayStart - (displayIdx - seq1Offset + 1));
          const int seq2Idx = seq2Fwd ? seq2DisplayStart + (displayIdx - seq2Offset) - seq2Start : seq2Start - (seq2DisplayStart - (displayIdx - seq2Offset + 1));
          
          seq1Text[displayIdx] = ' ';
          
	  if (seq1Idx >= 0 && seq1Idx < seq1Len)
	    {
              seq1Text[displayIdx] = seq1[seq1Idx];
              
              if (seq2Idx >= 0 && seq2Idx < seq2Len)
                {
                  if (seq1[seq1Idx] == seq2[seq2Idx])
                    {
                      match = TRUE;
                      break;
                    }
                  else if (dc->blastMode != BLXMODE_BLASTN && 
                           dc->matrix[atob[(unsigned int)seq1[seq1Idx]] - 1 ][atob[(unsigned int)seq2[seq2Idx]] - 1 ] > 0)
                    {
                      consMatch = TRUE;
                    }
                }
            }
        }
      
      if (match)
        {
          gdk_gc_set_foreground(gc, matchColor);
          gdk_draw_rectangle(drawable, gc, TRUE, x, y, ceil(dc->charWidth), roundNearest(dc->charHeight));
        }
      else if (consMatch)
        {
          gdk_gc_set_foreground(gc, consColor);
          gdk_draw_rectangle(drawable, gc, TRUE, x, y, ceil(dc->charWidth), roundNearest(dc->charHeight));
        }
    }  
  
  seq1Text[displayIdx] = '\0';
  
  /* Draw the sequence text. */
  x = 0;//seq1Offset * dc->charWidth;
  y = 0;
  
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, seq1Text);
  pango_layout_set_font_description(layout, dc->fontDesc);
  
  if (layout)
    {
      gtk_paint_layout(widget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, widget, NULL, x, y, layout);
      g_object_unref(layout);
    }
}


/* Utility called by drawSequenceHeader to draw a marker line at the given coords */
static void drawSequenceHeaderMarker(GdkDrawable *drawable, const int x, const int y, const gdouble charWidth)
{
  /* Draw a marker line below where the text will go */
  const int x1 = x + roundNearest((charWidth / 2.0));
  
  GdkGC *gc = gdk_gc_new(drawable);
  gdk_draw_line(drawable, gc, x1, y, x1, y +  + SELECTED_COORD_MARKER_HEIGHT);
}


/* Utility called by drawSequenceHeader to draw the header text at the given coords */
static void drawSequenceHeaderText(GtkWidget *widget, 
				   GdkDrawable *drawable, 
				   const int x, 
				   const int y, 
				   const int coord, 
				   const gdouble charWidth,
				   PangoFontDescription *fontDesc)
{
  char *displayText = convertIntToString(coord);
  
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, displayText);
  pango_layout_set_font_description(layout, fontDesc);
  
  if (layout)
    {
    /* Offset the text so that the middle of the text is lined up with the coord of interest */
    const int offset = ceil((((gdouble)numDigitsInInt(coord) / 2.0) - 1) * charWidth);
    
    gtk_paint_layout(widget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, widget, NULL, x - offset, y, layout);
    g_object_unref(layout);
    }
  
  g_free(displayText);
  
}


/* Convert a dna index to display (dna or peptide) index, if applicable. If horizontal is true
 * we have the horizontal (reference) sequence, otherwise the vertical (match) sequence. */
static int convertToDisplayIdx(const int dnaIdx, const gboolean horizontal, DotterContext *dc)
{
  int result = dnaIdx;
  
  /* Match seq coords are always in display coords already, so only do anything if this is the
   * ref seq. Also, we only need to convert if it's peptide-nucleotide match. */
  if (horizontal && dc->blastMode == BLXMODE_BLASTX)
    {
      result /= dc->numFrames;
    }
  
  return result;
}


/* Draw a sequence header that displays the centre coord of the given display range. If
 * markerFirst is true the marker is drawn above the text, otherwise below. */
static void drawSequenceHeader(GtkWidget *widget, 
			       GdkDrawable *drawable, 
                               DotterContext *dc,
			       const IntRange const *displayRange, 
			       const gboolean horizontal)
{
  /* Find the coord to display */
  const int coord = getRangeCentre(displayRange);
  
  /* Find the position to display at. Find the position of the char at this coord */
  int x = (int)((gdouble)convertToDisplayIdx(coord - displayRange->min, horizontal, dc) * dc->charWidth);
  int y = 0;

  if (horizontal)
    {
      /* For the ref sequence, draw the marker below the label. */
      drawSequenceHeaderText(widget, drawable, x, y, coord, dc->charWidth, dc->fontDesc);
      y += roundNearest(dc->charHeight);
      drawSequenceHeaderMarker(drawable, x, y, dc->charWidth);
    }
  else
    {
      /* For the match sequence, draw the marker above the label */
      drawSequenceHeaderMarker(drawable, x, y, dc->charWidth);
      y += SELECTED_COORD_MARKER_HEIGHT;
      drawSequenceHeaderText(widget, drawable, x, y, coord, dc->charWidth, dc->fontDesc);
    }
}

