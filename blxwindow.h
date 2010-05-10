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
  int numFrames;	    /* the number of reading frames */
  char *fetchMode;    /* the default method for fetching sequences */
} CommandLineOptions;


/* A Blixem View context, containing all status information required to draw the blixem view */
typedef struct _BlxViewContext
{
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
  GList *matchSeqs;		    /* List of all match sequences (as SequenceStructs). */
  gboolean gappedHsp;		    
  const char *paddingSeq;	    /* A sequence of padding characters, used if the real sequence could not be found. All padded MSPs
				     * use this same padding sequence - it is constructed to be long enough for the longest required seq. */
  
  gboolean displayRev;		    /* True if the display is reversed (i.e. coords increase from right to left). Set when strands are toggled. */
  GList *selectedSeqs;		    /* A list of sequence names that are selected */
  GList *sequenceGroups;	    /* A list of SequenceGroups */
  SequenceGroup *matchSetGroup;	    /* A special group that can be created/deleted quickly from the 'toggle match set' shortcuts */
  
  gboolean autoDotter;		    /* Whether to use automatic dotter params */
  int dotterStart;		    /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterEnd;		    /* End coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterZoom;		    /* Zoom param to call dotter with */
  
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
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const SequenceStruct *seqToFind);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
Strand			  blxWindowGetActiveStrand(GtkWidget *blxWindow);

GList*			  blxWindowGetSelectedSeqs(GtkWidget *blxWindow);
void			  blxWindowSelectSeq(GtkWidget *blxWindow, SequenceStruct *seq, const gboolean updateTrees);
void			  blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList);
void			  blxWindowDeselectSeq(GtkWidget *blxWindow, SequenceStruct *seq, const gboolean updateTrees);
void			  blxWindowDeselectAllSeqs(GtkWidget *blxWindow, const gboolean updateTrees);
gboolean		  blxWindowIsSeqSelected(GtkWidget *blxWindow, const SequenceStruct *seq);
void			  blxWindowSelectionChanged(GtkWidget *blxWindow, const gboolean updateTrees);

int			  sequenceGetGroupOrder(GtkWidget *blxWindow, const SequenceStruct *seq);
void			  copySelectionToClipboard(GtkWidget *blxWindow);
void			  findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);

void			  showHelpDialog(GtkWidget *blxWindow);
void			  showSettingsDialog(GtkWidget *blxWindow);
void			  showViewPanesDialog(GtkWidget *blxWindow);
void			  showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups);
void			  showFindDialog(GtkWidget *blxWindow);
void			  showAboutDialog(GtkWidget *blxWindow);

void			  blxWindowRedrawAll(GtkWidget *blxWindow);

gchar*			  getSequenceSegment(BlxViewContext *bc,
					     const char const *dnaSequence,
					     const int coord1, 
					     const int coord2,
					     const Strand strand,
					     const BlxSeqType inputCoordType,
					     const int frame,
					     const gboolean displayRev,
					     const gboolean reverseResult,
					     const gboolean allowComplement,
					     const gboolean translateResult);
  
GtkWidget*		  createBlxWindow(CommandLineOptions *options, const char *paddingSeq);


#endif /* _blxwindow_included_ */
