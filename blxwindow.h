/*
 *  blxwindow.h
 *  acedb
 *
 *  Created by Gemma Barson on 24/11/2009.
 *
 */

#ifndef _blxwindow_included_
#define _blxwindow_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>


/* Struct to hold all the settings that come from the command line options */
typedef struct _CommandLineOptions
{
  char *refSeq;			  /* the section of reference sequence we're viewing */
  const char const *refSeqName;	  /* the name of the reference sequence */
  const int refSeqOffset;	  /* how much to offset the first ref seq coord by */
  const int startCoord;		  /* which coord to start the initial display range at */
  MSP *mspList;			  /* the list of alignments */
  char **geneticCode;		  /* the genetic code */
  
  BlxStrand activeStrand;	  /* which strand will initially be the active one */
  int bigPictZoom;		  /* initial zoom level for the big picture (as a multiple of the initial detail view range) */
  gboolean bigPictON;	    
  gboolean bigPictRev;	    
  BlxColumnId initSortColumn;	  /* initial column to sort by */
  gboolean sortInverted;	  /* whether initial sort order should be inverted */
  gboolean gappedHsp;		  /* whether this is a gapped hsp */
  gboolean hiliteSins;	    
  gboolean dotterFirst;		  /* open dotter when blixem starts */
  gboolean startNextMatch;	  /* start at the coord of the next match from the default start coord */
  gboolean parseFullEmblInfo;     /* parse the full EMBL files on startup to populate additional info like tissue-type */
  BlxBlastMode blastMode;	  /* the blast match mode */
  BlxSeqType seqType;		  /* the type of sequence i.e. DNA or peptide */
  int numFrames;		  /* the number of reading frames */
  char *fetchMode;		  /* the default method for fetching sequences */
} CommandLineOptions;


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

  MSP *mspList;			    /* List of all MSPs. */
  GList *matchSeqs;		    /* List of all match sequences (as BlxSequences). */
  GSList *supportedTypes;           /* List of supported GFF types */
  gboolean gappedHsp;		    
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
gboolean		  blxWindowGetGappedHsp(GtkWidget *blxWindow);
MSP*			  blxWindowGetMspList(GtkWidget *blxWindow);
GList*			  blxWindowGetAllMatchSeqs(GtkWidget *blxWindow);
GList*			  blxWindowGetSequenceGroups(GtkWidget *blxWindow);
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
BlxStrand		  blxWindowGetActiveStrand(GtkWidget *blxWindow);

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
					  GList *seqList, 
                                          GSList *supportedTypes,
					  const char *net_id, 
					  int port, 
					  const gboolean External);


#endif /* _blxwindow_included_ */
