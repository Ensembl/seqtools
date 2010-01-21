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


#define VERTICAL_SEPARATOR_HEIGHT	2
#define NAME_COLUMN_DEFAULT_WIDTH	100
#define SCORE_COLUMN_DEFAULT_WIDTH	30
#define ID_COLUMN_DEFAULT_WIDTH		30
#define START_COLUMN_DEFAULT_WIDTH	50
#define SEQ_COLUMN_DEFAULT_WIDTH	40
#define END_COLUMN_DEFAULT_WIDTH	50


typedef struct _DetailViewProperties
  {
    GtkWidget *mainWindow;	  /* The main window that this view belongs to */
    GtkCellRenderer *renderer;	  /* The cell renderer that renders the sequences */
    GtkAdjustment *adjustment;	  /* The scroll adjustment control for the detail view */
    GtkWidget *feedbackBox;	  /* A text box that feeds back info to the user about the currently selected items */
    GtkWidget *header;		  /* The header for the detail view trees */
    
    GList *fwdStrandTrees;	  /* A list of all the trees that show the forward strand of the ref seq */
    GList *revStrandTrees;	  /* A list of all the trees that show the reverse strand of the ref seq */
    
    BlxSeqType seqType;		  /* The match type, i.e. dna or peptide */
    int numReadingFrames;	  /* The number of reading frames */
    int verticalSeparator;	  /* The vertical distance between the tree rows */
        
    IntRange displayRange;	  /* The currently-displayed range of bases in the reference sequence */
    int selectedBaseIdx;	  /* The currently-selected base in the reference sequence */
    PangoFontDescription *fontDesc; /* The fixed-width font that will be used to display the alignments */
    
    /* Display colours */
    GdkColor refSeqColour;
    GdkColor refSeqSelectedColour;
    GdkColor matchColour;
    GdkColor matchSelectedColour;
    GdkColor mismatchColour;
    GdkColor mismatchSelectedColour;
    GdkColor exonColour;
    GdkColor exonSelectedColour;
    GdkColor gapColour;
    GdkColor gapSelectedColour;
  } DetailViewProperties;


/* Public function declarations */
char*			detailViewGetRefSeq(GtkWidget *detailView);
char*			detailViewGetDisplaySeq(GtkWidget *detailView);
int			detailViewGetNumReadingFrames(GtkWidget *detailView);
IntRange*		detailViewGetDisplayRange(GtkWidget *detailView);
IntRange*		detailViewGetFullRange(GtkWidget *detailView);
int			detailViewGetSelectedBaseIdx(GtkWidget *detailView);
int			detailViewGetOldSelectedBaseIdx(GtkWidget *detailView);
GtkAdjustment*		detailViewGetAdjustment(GtkWidget *detailView);
GList*			detailViewGetFwdStrandTrees(GtkWidget *detailView);
GList*			detailViewGetRevStrandTrees(GtkWidget *detailView);
GtkWidget*		detailViewGetFrameTree(GtkWidget *detailView, const Strand strand, const int frame);
GList*			detailViewGetStrandTrees(GtkWidget *detailView, const Strand strand);
gboolean		detailViewGetStrandsToggled(GtkWidget *detailView);
BlxBlastMode		detailViewGetBlastMode(GtkWidget *detailView);
PangoFontDescription*	detailViewGetFontDesc(GtkWidget *detailView);
GtkCellRenderer*	detailViewGetRenderer(GtkWidget *detailView);
int			detailViewGetVerticalSeparator(GtkWidget *detailView);
BlxSeqType		detailViewGetSeqType(GtkWidget *detailView);
char**			detailViewGetGeneticCode(GtkWidget *detailView);
IntRange*		detailViewGetRefSeqRange(GtkWidget *detailView);
GtkWidget*	        detailViewGetMainWindow(GtkWidget *detailView);

DetailViewProperties*	detailViewGetProperties(GtkWidget *widget);

GdkColor*		detailViewGetRefSeqColour(GtkWidget *tree);
GdkColor*		detailViewGetRefSeqSelectedColour(GtkWidget *tree);
GdkColor*		detailViewGetMatchColour(GtkWidget *tree);
GdkColor*		detailViewGetMatchSelectedColour(GtkWidget *tree);
GdkColor*		detailViewGetMismatchColour(GtkWidget *tree);
GdkColor*		detailViewGetMismatchSelectedColour(GtkWidget *tree);
GdkColor*		detailViewGetExonColour(GtkWidget *tree);
GdkColor*		detailViewGetExonSelectedColour(GtkWidget *tree);
GdkColor*		detailViewGetGapColour(GtkWidget *tree);
GdkColor*		detailViewGetGapSelectedColour(GtkWidget *tree);

void			setDetailViewScrollPos(GtkAdjustment *adjustment, int value);
void			scrollDetailViewLeft1(GtkWidget *detailView);
void			scrollDetailViewRight1(GtkWidget *detailView);
void			scrollDetailViewLeftStep(GtkWidget *detailView);
void			scrollDetailViewRightStep(GtkWidget *detailView);
void			scrollDetailViewLeftPage(GtkWidget *detailView);
void			scrollDetailViewRightPage(GtkWidget *detailView);

void			updateFeedbackBox(GtkWidget *detailView);

void			detailViewAddMspData(GtkWidget *detailView, MSP *mspList);

GtkWidget*		createDetailView(GtkWidget *mainWindow,
					 GtkWidget *panedWidget,
					 GtkAdjustment *adjustment, 
					 GtkWidget *fwdStrandGrid, 
					 GtkWidget *revStrandGrid,
					 MSP *mspList,
					 BlxBlastMode mode,
					 BlxSeqType seqType,
					 int numReadingFrames,
					 IntRange *refSeqRange);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *mainWindow);


#endif /* _detail_view_included_ */