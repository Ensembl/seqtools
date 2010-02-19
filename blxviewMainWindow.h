/*
 *  blxviewMainWindow.h
 *  acedb
 *
 *  Created by Gemma Barson on 24/11/2009.
 *
 */

#ifndef _blxview_main_window_included_
#define _blxview_main_window_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>

typedef struct _MainWindowProperties
  {
    GtkWidget *bigPicture;
    GtkWidget *detailView;
    GtkWidget *mainmenu;
    
    char *refSeq;		    /* The reference sequence (always forward strand, always DNA sequence) */
    const char *refSeqName;	    /* The name of the reference sequence */
    char *displaySeq;		    /* The displayed sequence (same as ref seq or converted to peptide sequence) */
    IntRange refSeqRange;	    /* The range of the reference sequence */
    IntRange fullDisplayRange;	    /* The range of the displayed sequence */
    const gboolean gappedHsp;	    
    
    MSP *mspList;		    /* Linked list of match sequences */
    char **geneticCode;		    /* The genetic code used to translate DNA <-> peptide */
    BlxBlastMode blastMode;	    /* The type of blast matching that was used */
    BlxSeqType seqType;		    /* The type of sequence, e.g. DNA or peptide */
    int numReadingFrames;	    /* The number of reading frames */

    gboolean strandsToggled;	    /* If true, the reverse strand becomes the 'main' or 'top' strand */
    GList *selectedSeqs;	    /* A list of sequence names that are selected */
    GList *sequenceGroups;	    /* A list of SequenceGroups */
    
    gboolean autoDotterParams;	    /* Whether to use automatic dotter params */
    int dotterStart;		    /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterEnd;		    /* End coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterZoom;		    /* Zoom param to call dotter with */

    GdkDrawable *drawable;	    /* A bitmap where we'll draw the contents we want to print */
    GtkPrintSettings *printSettings;  /* Used so that we can re-use the same print settings as a previous print */
    int lastYEnd;		    /* Keeps track of where the last item ended so we can draw the next one flush to it */
    int lastYStart;		    /* Where the last item started (for drawing multiple items at same y pos) */
    int lastYCoord;		    /* Y coord of last item (so we can check if current item should be at same Y pos) */    
  } MainWindowProperties;


/* Public function declarations */
MainWindowProperties*	  mainWindowGetProperties(GtkWidget *widget);
gboolean		  mainWindowGetStrandsToggled(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetBigPicture(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetDetailView(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetMainMenu(GtkWidget *mainWindow);
BlxBlastMode		  mainWindowGetBlastMode(GtkWidget *mainWindow);
IntRange*		  mainWindowGetFullRange(GtkWidget *mainWindow);
IntRange*		  mainWindowGetRefSeqRange(GtkWidget *mainWindow);
const char*		  mainWindowGetRefSeqName(GtkWidget *mainWindow);
BlxSeqType		  mainWindowGetSeqType(GtkWidget *mainWindow);
char**			  mainWindowGetGeneticCode(GtkWidget *mainWindow);
char*			  mainWindowGetRefSeq(GtkWidget *mainWindow);
char*			  mainWindowGetDisplaySeq(GtkWidget *mainWindow);
int			  mainWindowGetNumReadingFrames(GtkWidget *mainWindow);
int			  mainWindowGetDotterStart(GtkWidget *mainWindow);
int			  mainWindowGetDotterEnd(GtkWidget *mainWindow);
int			  mainWindowGetAutoDotter(GtkWidget *mainWindow);
gboolean		  mainWindowGetGappedHsp(GtkWidget *mainWindow);
MSP*			  mainWindowGetMspList(GtkWidget *mainWindow);
GList*			  mainWindowGetSequenceMsps(GtkWidget *mainWindow, const char *seqName);
GList*			  mainWindowGetSequenceGroups(GtkWidget *mainWindow);
SequenceGroup*		  mainWindowGetSequenceGroup(GtkWidget *mainWindow, const char *seqName);

GList*			  mainWindowGetSelectedSeqs(GtkWidget *mainWindow);
void			  mainWindowSelectSeq(GtkWidget *mainWindow, char *seqName, const gboolean updateTrees);
void			  mainWindowDeselectSeq(GtkWidget *mainWindow, char *seqName, const gboolean updateTrees);
void			  mainWindowDeselectAllSeqs(GtkWidget *mainWindow, const gboolean updateTrees);
gboolean		  mainWindowIsSeqSelected(GtkWidget *mainWindow, const char *msp);
void			  mainWindowSelectionChanged(GtkWidget *mainWindow, const gboolean updateTrees);

int			  sequenceGetGroupOrder(GtkWidget *mainWindow, const char *seqName);

void			  displayHelp(GtkWidget *mainWindow);

void			  mainWindowRedrawAll(GtkWidget *mainWindow);

gchar*			  getSequenceSegment(GtkWidget *mainWindow, 
					     const char const *sequence,
					     const IntRange const *sequenceRange,
					     const int coord1, 
					     const int coord2,
					     const Strand strand,
					     const BlxSeqType seqType,
					     const int frame,
					     const int numReadingFrames,
					     const gboolean reverse,
					     const gboolean translate);
  
GtkWidget*		  createMainWindow(char *refSeq, 
					   const char const *refSeqName,
					   MSP *mspList, 
					   BlxBlastMode blastMode,
					   BlxSeqType seqType, 
					   int numReadingFrames,
					   char **geneticCode,
					   const int refSeqOffset,
					   const int startCoord,
					   const SortByType sortByType,
					   const gboolean sortInverted,
					   const gboolean gappedHsp);


#endif /* _blxview_main_window_included_ */