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
    
    char *refSeq;		    /* The reference sequence (always forward strand, always DNA sequence) */
    char *displaySeq;		    /* The displayed sequence (same as ref seq or converted to peptide sequence) */
    IntRange refSeqRange;	    /* The range of the reference sequence */
    IntRange fullDisplayRange;	    /* The range of the displayed sequence */
    
    MSP *mspList;		    /* Linked list of match sequences */
    char **geneticCode;		    /* The genetic code used to translate DNA <-> peptide */
    BlxBlastMode blastMode;	    /* The type of blast matching that was used */
    BlxSeqType seqType;		    /* The type of sequence, e.g. DNA or peptide */
    int numReadingFrames;	    /* The number of reading frames */

    gboolean strandsToggled;	    /* If true, the reverse strand becomes the 'main' or 'top' strand */
    GList *selectedMsps;	    /* List of MSPs that are selected */
  } MainWindowProperties;


/* Public function declarations */
void blxHelp(void);

MainWindowProperties*	  mainWindowGetProperties(GtkWidget *widget);
gboolean		  mainWindowGetStrandsToggled(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetBigPicture(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetDetailView(GtkWidget *mainWindow);
BlxBlastMode		  mainWindowGetBlastMode(GtkWidget *mainWindow);
IntRange*		  mainWindowGetFullRange(GtkWidget *mainWindow);
IntRange*		  mainWindowGetRefSeqRange(GtkWidget *mainWindow);
BlxSeqType		  mainWindowGetSeqType(GtkWidget *mainWindow);
char**			  mainWindowGetGeneticCode(GtkWidget *mainWindow);
char*			  mainWindowGetRefSeq(GtkWidget *mainWindow);
char*			  mainWindowGetDisplaySeq(GtkWidget *mainWindow);
int			  mainWindowGetNumReadingFrames(GtkWidget *mainWindow);
GList*			  mainWindowGetSelectedMsps(GtkWidget *mainWindow);

void			  mainWindowSelectMsp(GtkWidget *mainWindow, MSP *msp, const gboolean updateTrees);
void			  mainWindowDeselectMsp(GtkWidget *mainWindow, MSP *msp, const gboolean updateTrees);
void			  mainWindowDeselectAllMsps(GtkWidget *mainWindow, const gboolean updateTrees);
gboolean		  mainWindowIsMspSelected(GtkWidget *mainWindow, MSP *msp);

gchar*			  getSequenceSegment(GtkWidget *mainWindow, 
					   const char const *sequence,
					   const IntRange const *sequenceRange,
					   const int coord1, 
					   const int coord2,
					   const Strand strand,
					   const BlxSeqType seqType,
					   const int frame,
					   const int numReadingFrames,
					   const gboolean reverse);
  
GtkWidget*		  createMainWindow(char *refSeq, 
					   const char const *refSeqName,
					   MSP *mspList, 
					   BlxBlastMode blastMode,
					   BlxSeqType seqType, 
					   int numReadingFrames,
					   char **geneticCode,
					   const int refSeqOffset);


#endif /* _blxview_main_window_included_ */