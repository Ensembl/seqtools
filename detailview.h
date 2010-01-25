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

/* Define the columns. Specify a default width, the display text for the
 * column header, and also a name for the property that will be set in the
 * cell renderer. (The latter cannot contain special characters.) */
#define NAME_COLUMN_DEFAULT_WIDTH	100
#define SCORE_COLUMN_DEFAULT_WIDTH	30
#define ID_COLUMN_DEFAULT_WIDTH		30
#define START_COLUMN_DEFAULT_WIDTH	50
#define SEQ_COLUMN_DEFAULT_WIDTH	40
#define END_COLUMN_DEFAULT_WIDTH	70

#define NAME_COLUMN_HEADER_TEXT		"Name"
#define SCORE_COLUMN_HEADER_TEXT	"Score"
#define ID_COLUMN_HEADER_TEXT		"%Id"
#define START_COLUMN_HEADER_TEXT	"Start"
#define SEQ_COLUMN_HEADER_TEXT		"Sequence"
#define END_COLUMN_HEADER_TEXT		"End"

#define NAME_COLUMN_PROPERTY_NAME	"name"
#define SCORE_COLUMN_PROPERTY_NAME	"score"
#define ID_COLUMN_PROPERTY_NAME		"id"
#define START_COLUMN_PROPERTY_NAME	"start"
#define SEQ_COLUMN_PROPERTY_NAME	"sequence"
#define END_COLUMN_PROPERTY_NAME	"end"


/* This enum declares identifiers for each column in the detail view */
typedef enum
  {
    S_NAME_COL,
    SCORE_COL,
    ID_COL,
    START_COL,
    MSP_COL,
    END_COL,
    
    N_COLUMNS
  } ColumnId;


/* This struct describes a column in the detail view. Multiple widgets (i.e. headers
 * and tree columns) in the detail view must all have columns that share the same
 * properties (namely the column width). */
typedef struct _DetailViewColumnInfo
  {
    ColumnId columnId;		/* the column identifier */
    GtkWidget *headerWidget;	/* the header widget for this column (in the detail-view header) */
    GtkCallback refreshFunc;	/* the function that will be called on the header widget when columns are refreshed */
    char *title;		/* the default column title */
    char *propertyName;		/* the property name (used to set the data for the SequenceCellRenderer) */
    int width;			/* the column width */
  } DetailViewColumnInfo;



/* Essential info required by the the detail view */
typedef struct _DetailViewProperties
  {
    GtkWidget *mainWindow;	  /* The main window that this view belongs to */
    GtkCellRenderer *renderer;	  /* The cell renderer that renders the sequences */
    GtkAdjustment *adjustment;	  /* The scroll adjustment control for the detail view */
    GtkWidget *feedbackBox;	  /* A text box that feeds back info to the user about the currently selected items */
    GList *columnList;		  /* A list of details about all the columns in the detail view */
    
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

int			getDetailViewColumnWidth(GtkWidget *detailView, const ColumnId columnId);

GdkColor*		detailViewGetRefSeqColour(GtkWidget *detailView);
GdkColor*		detailViewGetRefSeqSelectedColour(GtkWidget *detailView);
GdkColor*		detailViewGetMatchColour(GtkWidget *detailView);
GdkColor*		detailViewGetMatchSelectedColour(GtkWidget *detailView);
GdkColor*		detailViewGetMismatchColour(GtkWidget *detailView);
GdkColor*		detailViewGetMismatchSelectedColour(GtkWidget *detailView);
GdkColor*		detailViewGetExonColour(GtkWidget *detailView);
GdkColor*		detailViewGetExonSelectedColour(GtkWidget *detailView);
GdkColor*		detailViewGetGapColour(GtkWidget *detailView);
GdkColor*		detailViewGetGapSelectedColour(GtkWidget *detailView);

void			setDetailViewScrollPos(GtkAdjustment *adjustment, int value);
void			scrollDetailViewLeft1(GtkWidget *detailView);
void			scrollDetailViewRight1(GtkWidget *detailView);
void			scrollDetailViewLeftStep(GtkWidget *detailView);
void			scrollDetailViewRightStep(GtkWidget *detailView);
void			scrollDetailViewLeftPage(GtkWidget *detailView);
void			scrollDetailViewRightPage(GtkWidget *detailView);

void			detailViewSetSelectedBaseIdx(GtkWidget *detailView, const int selectedBaseIdx);
void			updateFeedbackBox(GtkWidget *detailView);

void			detailViewAddMspData(GtkWidget *detailView, MSP *mspList);
void			updateDetailViewFontDesc(GtkWidget *detailView);

GtkWidget*		createDetailView(GtkWidget *mainWindow,
					 GtkWidget *panedWidget,
					 GtkAdjustment *adjustment, 
					 GtkWidget *fwdStrandGrid, 
					 GtkWidget *revStrandGrid,
					 MSP *mspList,
					 BlxBlastMode mode,
					 BlxSeqType seqType,
					 int numReadingFrames,
					 const char const *refSeqName);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *mainWindow);


#endif /* _detail_view_included_ */