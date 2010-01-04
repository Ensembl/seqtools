/*
 *  detailview.h
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#ifndef _detail_view_included_
#define _detail_view_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>


typedef struct _DetailViewProperties
  {
    GtkWidget *mainWindow;	/* The main window that this view belongs to */
    GtkWidget *feedbackBox;	/* A text box that feeds back info to the user about the currently selected items */
    GtkAdjustment *adjustment;  /* The scroll adjustment control for the detail view */
    
    GList *fwdStrandTrees;	/* A list of all the trees that show the forward strand of the ref seq */
    GList *revStrandTrees;	/* A list of all the trees that show the reverse strand of the ref seq */
    
    char *refSeq;		/* The reference sequence (forward strand) */
    BlxSeqType seqType;		/* The match type, i.e. dna or peptide */
    int numReadingFrames;	/* The number of reading frames */
        
    IntRange displayRange;	/* The currently-displayed range of bases in the reference sequence */
    int selectedBaseIdx;	/* The currently-selected base in the reference sequence */
    
    const char *fontFamily;	/* The fixed-width font that will be used to display the alignments */
  } DetailViewProperties;


/* Public function declarations */
char*			detailViewGetRefSeq(GtkWidget *detailView);
char*			detailViewGetRefSeqReverseStrand(GtkWidget *detailView);
int			detailViewGetNumReadingFrames(GtkWidget *detailView);
IntRange*		detailViewGetDisplayRange(GtkWidget *detailView);
int			detailViewGetSelectedBaseIdx(GtkWidget *detailView);
int			detailViewGetOldSelectedBaseIdx(GtkWidget *detailView);
GtkAdjustment*		detailViewGetAdjustment(GtkWidget *detailView);
GList*			detailViewGetFwdStrandTrees(GtkWidget *detailView);
GList*			detailViewGetRevStrandTrees(GtkWidget *detailView);
GtkWidget*		detailViewGetFrameTree(GtkWidget *detailView, gboolean forward, int frame);
gboolean		detailViewGetStrandsToggled(GtkWidget *detailView);
BlxBlastMode		detailViewGetBlastMode(GtkWidget *detailView);
const char*		detailViewGetFontFamily(GtkWidget *detailView);

DetailViewProperties*	detailViewGetProperties(GtkWidget *widget);

void			setDetailViewScrollPos(GtkAdjustment *adjustment, int value);
void			scrollDetailViewLeft1(GtkWidget *detailView);
void			scrollDetailViewRight1(GtkWidget *detailView);
void			scrollDetailViewLeftStep(GtkWidget *detailView);
void			scrollDetailViewRightStep(GtkWidget *detailView);
void			scrollDetailViewLeftPage(GtkWidget *detailView);
void			scrollDetailViewRightPage(GtkWidget *detailView);

GtkWidget*		createDetailView(GtkWidget *mainWindow,
					 GtkWidget *panedWidget,
					 GtkAdjustment *adjustment, 
					 GtkWidget *fwdStrandGrid, 
					 GtkWidget *revStrandGrid,
					 MSP *mspList,
					 char *refSeq,
					 BlxBlastMode mode,
					 BlxSeqType seqType,
					 int numReadingFrames,
					 IntRange *refSeqRange);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *mainWindow);


#endif /* _detail_view_included_ */