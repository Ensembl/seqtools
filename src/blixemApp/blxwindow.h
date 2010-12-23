/*  File: blxWindow.h
 *  Author: Gemma Barson, 2009-11-24
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
 * Description: Creates the main Blixem window. Also creates a "context", 
 *              which contains all of the variables associated with a Blixem
 *              session. 
 *
 *              The context could live somewhere else, but was just included
 *              here because there is one context for each Blixem window.
 *              Ideally the context would be created earlier (perhaps in the
 *              main function) but so far that has not been practical.
 *----------------------------------------------------------------------------
 */

#ifndef _blxwindow_included_
#define _blxwindow_included_

#include <gtk/gtk.h>
#include <blixemApp/blixem_.h>
#include <seqtoolsUtils/utilities.h>

/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the BlxViewContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    BLXDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */
    
    BLXDIALOG_HELP,                 /* The Help dialog */
    BLXDIALOG_SETTINGS,             /* The Settings dialog */
    BLXDIALOG_FIND,                 /* The Find dialog */
    BLXDIALOG_GROUPS,               /* The Groups dialog */
    BLXDIALOG_VIEW,                 /* The View dialog */
    BLXDIALOG_DOTTER,               /* The Dotter dialog */
    
    BLXDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } BlxDialogId;



/* A Blixem View context, containing all status information required to draw the blixem view */
typedef struct _BlxViewContext
{
  GtkWidget *statusBar;		    /* The Blixem window's status bar */

  char *refSeq;			    /* The reference sequence (always forward strand, always DNA sequence) */
  const char *refSeqName;	    /* The name of the reference sequence */
  IntRange refSeqRange;		    /* The range of the reference sequence */
  IntRange fullDisplayRange;	    /* The range of the displayed sequence */
  
  BlxBlastMode blastMode;	    /* The type of blast matching that was used */
  BlxSeqType seqType;		    /* The type of sequence, e.g. DNA or peptide */
  char* fetchMode;		    /* The fetch method to use */
  char **geneticCode;		    /* The genetic code used to translate DNA <-> peptide */
  int numFrames;		    /* The number of reading frames */

  MSP *mspList;                          /* List of all MSPs. Obsolete - use featureLists array instead */
  GList* featureLists[BLXMSP_NUM_TYPES];  /* Array indexed by the BlxMspType enum. Each array entry contains a GList of all the MSPs of that type. */

  GList *matchSeqs;		    /* List of all match sequences (as BlxSequences). */
  GSList *supportedTypes;           /* List of supported GFF types */
  const char *paddingSeq;	    /* A sequence of padding characters, used if the real sequence could not be found. All padded MSPs
				     * use this same padding sequence - it is constructed to be long enough for the longest required seq. */
  
  gboolean displayRev;		    /* True if the display is reversed (i.e. coords decrease as you read from left to right, rather than increase). */
  const char *net_id;               /* pfetch-socket net id */
  int port;                         /* pfetch-socket port */
  gboolean external;                /* True if Blixem was run externally or false if it was run internally from another program */
  
  GList *selectedSeqs;		    /* A list of sequences that are selected (as BlxSequences) */
  GList *sequenceGroups;	    /* A list of SequenceGroups */
  SequenceGroup *matchSetGroup;	    /* A special group that can be created/deleted quickly from the 'toggle match set' shortcuts */
  
  gboolean autoDotter;		    /* Whether to use automatic dotter params */
  int dotterStart;		    /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterEnd;		    /* End coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterZoom;		    /* Zoom param to call dotter with */
  
  GArray *defaultColors;	    /* Default colors used by Blixem */
  gboolean usePrintColors;	    /* Whether to use print colors (i.e. black and white) */
  
  gboolean flags[BLXFLAG_NUM_FLAGS];              /* Array of all the flags the user can toggle. Indexed by the BlxFlags enum. */
  GtkWidget *dialogList[BLXDIALOG_NUM_DIALOGS];   /* Array of all the persistent dialogs in the application */
  GSList *spawnedProcesses;			  /* List of processes spawned by Blixem */
} BlxViewContext;


/* Public function declarations */
BlxViewContext*		  blxWindowGetContext(GtkWidget *widget);
gboolean		  blxWindowGetDisplayRev(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetBigPicture(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetDetailView(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetMainMenu(GtkWidget *blxWindow);
BlxBlastMode		  blxWindowGetBlastMode(GtkWidget *blxWindow);
const char*		  blxWindowGetFetchMode(GtkWidget *blxWindow);
IntRange*		  blxWindowGetFullRange(GtkWidget *blxWindow);
IntRange*		  blxWindowGetRefSeqRange(GtkWidget *blxWindow);
const char*		  blxWindowGetRefSeqName(GtkWidget *blxWindow);
BlxSeqType		  blxWindowGetSeqType(GtkWidget *blxWindow);
char**			  blxWindowGetGeneticCode(GtkWidget *blxWindow);
char*			  blxWindowGetRefSeq(GtkWidget *blxWindow);
int			  blxWindowGetNumFrames(GtkWidget *blxWindow);
int			  blxWindowGetDotterStart(GtkWidget *blxWindow);
int			  blxWindowGetDotterEnd(GtkWidget *blxWindow);
int			  blxWindowGetDotterZoom(GtkWidget *blxWindow);
int			  blxWindowGetAutoDotter(GtkWidget *blxWindow);
MSP*			  blxWindowGetMspList(GtkWidget *blxWindow);
GList*			  blxWindowGetAllMatchSeqs(GtkWidget *blxWindow);
GList*			  blxWindowGetSequenceGroups(GtkWidget *blxWindow);
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
BlxStrand		  blxWindowGetActiveStrand(GtkWidget *blxWindow);
gboolean                  blxWindowGetNegateCoords(GtkWidget *blxWindow);

GList*			  blxWindowGetSelectedSeqs(GtkWidget *blxWindow);
void			  blxWindowSelectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void			  blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList);
void			  blxWindowDeselectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void			  blxWindowDeselectAllSeqs(GtkWidget *blxWindow);
gboolean		  blxWindowIsSeqSelected(GtkWidget *blxWindow, const BlxSequence *seq);
void			  blxWindowSetSeqSelected(GtkWidget *blxWindow, BlxSequence *seq, const gboolean selected);
void			  blxWindowSelectionChanged(GtkWidget *blxWindow);
BlxSequence*		  blxWindowGetLastSelectedSeq(GtkWidget *blxWindow);

gboolean                  blxContextIsSeqSelected(BlxViewContext *bc, const BlxSequence *seq);

int			  sequenceGetGroupOrder(GtkWidget *blxWindow, const BlxSequence *seq);
void			  copySelectionToClipboard(GtkWidget *blxWindow);
void			  findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);

void                      refreshDialog(const BlxDialogId dialogId, GtkWidget *blxWindow);
void			  showHelpDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void			  showSettingsDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void			  showViewPanesDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void			  showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups, const gboolean bringToFront);
void			  showFindDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void			  showAboutDialog(GtkWidget *blxWindow);
void                      showInfoDialog(GtkWidget *blxWindow);

void			  blxWindowRedrawAll(GtkWidget *blxWindow);
  
GtkWidget*		  createBlxWindow(CommandLineOptions *options, 
					  const char *paddingSeq, 
                                          GList* featureLists[], 
					  GList *seqList, 
                                          GSList *supportedTypes,
					  const char *net_id, 
					  int port, 
					  const gboolean External);


#endif /* _blxwindow_included_ */
