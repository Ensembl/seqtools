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
  
  Strand activeStrand;	    /* which strand will initially be the active one */
  int bigPictZoom;	    /* initial zoom level for the big picture (as a multiple of the initial detail view range) */
  gboolean bigPictON;	    
  gboolean bigPictRev;	    
  SortByType initSortMode;  /* initial field to sort by */
  gboolean sortInverted;    /* whether initial sort order should be inverted */
  gboolean gappedHsp;	    /* whether this is a gapped hsp */
  gboolean hiliteSins;	    
  gboolean dotterFirst;	    /* open dotter when blixem starts */
  gboolean startNextMatch;  /* start at the coord of the next match from the default start coord */
  BlxBlastMode blastMode;   /* the blast match mode */
  BlxSeqType seqType;	    /* the type of sequence i.e. DNA or peptide */
  int numReadingFrames;	    /* the number of reading frames */
  char *fetchMode;    /* the default method for fetching sequences */
} CommandLineOptions;



typedef struct _BlxWindowProperties
  {
    GtkWidget *bigPicture;
    GtkWidget *detailView;
    GtkWidget *mainmenu;
    
    char *refSeq;		    /* The reference sequence (always forward strand, always DNA sequence) */
    const char *refSeqName;	    /* The name of the reference sequence */
    IntRange refSeqRange;	    /* The range of the reference sequence */
    IntRange fullDisplayRange;	    /* The range of the displayed sequence */

    MSP *mspList;		    /* List of all MSPs. */
    BlxBlastMode blastMode;	    /* The type of blast matching that was used */
    char* fetchMode;		    /* The fetch method to use */
    BlxSeqType seqType;		    /* The type of sequence, e.g. DNA or peptide */
    char **geneticCode;		    /* The genetic code used to translate DNA <-> peptide */
    int numReadingFrames;	    /* The number of reading frames */
    gboolean gappedHsp;		    
    const char *paddingSeq;	    /* A sequence of padding characters, used if the real sequence could not be found. All padded MSPs
				     * use this same padding sequence - it is constructed to be long enough for the longest required seq. */
    
    /* DYNAMIC PROPERTIES (these can be changed by the user): */
    
    gboolean strandsToggled;	    /* If true, the reverse strand becomes the 'main' or 'top' strand */
    GList *selectedSeqs;	    /* A list of sequence names that are selected */
    GList *sequenceGroups;	    /* A list of SequenceGroups */
    SequenceGroup *matchSetGroup;   /* A special group that can be created/deleted quickly from the 'toggle match set' shortcuts */

    gboolean autoDotterParams;	    /* Whether to use automatic dotter params */
    int dotterStart;		    /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterEnd;		    /* End coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterZoom;		    /* Zoom param to call dotter with */

    GtkPrintSettings *printSettings;  /* Used so that we can re-use the same print settings as a previous print */
    int lastYEnd;		    /* Keeps track of where the last item ended so we can draw the next one flush to it */
    int lastYStart;		    /* Where the last item started (for drawing multiple items at same y pos) */
    int lastYCoord;		    /* Y coord of last item (so we can check if current item should be at same Y pos) */
  } BlxWindowProperties;


/* Public function declarations */
BlxWindowProperties*	  blxWindowGetProperties(GtkWidget *widget);
gboolean		  blxWindowGetStrandsToggled(GtkWidget *blxWindow);
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
int			  blxWindowGetNumReadingFrames(GtkWidget *blxWindow);
int			  blxWindowGetDotterStart(GtkWidget *blxWindow);
int			  blxWindowGetDotterEnd(GtkWidget *blxWindow);
int			  blxWindowGetDotterZoom(GtkWidget *blxWindow);
int			  blxWindowGetAutoDotter(GtkWidget *blxWindow);
gboolean		  blxWindowGetGappedHsp(GtkWidget *blxWindow);
MSP*			  blxWindowGetMspList(GtkWidget *blxWindow);
GList*			  blxWindowGetSequenceMsps(GtkWidget *blxWindow, const char *seqName);
GList*			  blxWindowGetSequenceGroups(GtkWidget *blxWindow);
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const char *seqName);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
Strand			  blxWindowGetActiveStrand(GtkWidget *blxWindow);

GList*			  blxWindowGetSelectedSeqs(GtkWidget *blxWindow);
void			  blxWindowSelectSeq(GtkWidget *blxWindow, char *seqName, const gboolean updateTrees);
void			  blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqNameList);
void			  blxWindowDeselectSeq(GtkWidget *blxWindow, char *seqName, const gboolean updateTrees);
void			  blxWindowDeselectAllSeqs(GtkWidget *blxWindow, const gboolean updateTrees);
gboolean		  blxWindowIsSeqSelected(GtkWidget *blxWindow, const char *msp);
void			  blxWindowSelectionChanged(GtkWidget *blxWindow, const gboolean updateTrees);

int			  sequenceGetGroupOrder(GtkWidget *blxWindow, const char *seqName);
void			  copySelectionToClipboard(GtkWidget *blxWindow);
void			  findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);

void			  displayHelp(GtkWidget *blxWindow);
void			  showSettingsDialog(GtkWidget *blxWindow);
void			  showViewPanesDialog(GtkWidget *blxWindow);
void			  showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups);
void			  showFindDialog(GtkWidget *blxWindow);

void			  blxWindowRedrawAll(GtkWidget *blxWindow);

gchar*			  getSequenceSegment(GtkWidget *blxWindow, 
					     const char const *dnaSequence,
					     const IntRange const *dnaSequenceRange,
					     const int coord1, 
					     const int coord2,
					     const Strand strand,
					     const BlxSeqType inputCoordType,
					     const int frame,
					     const int numFrames,
					     const gboolean rightToLeft,
					     const gboolean reverseResult,
					     const gboolean allowComplement,
					     const gboolean translateResult);
  
GtkWidget*		  createBlxWindow(CommandLineOptions *options, const char *paddingSeq);


#endif /* _blxwindow_included_ */
