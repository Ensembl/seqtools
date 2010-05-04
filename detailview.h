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
#include <SeqTools/utilities.h>


#define	DEFAULT_REF_SEQ_BG_COLOUR	GDK_YELLOW
#define HEADER_CONTAINER_NAME		"header container"
#define SNP_TRACK_HEADER_NAME		"SNP track header"
#define DNA_TRACK_HEADER_NAME		"DNA track header"

/* Define the columns. Specify a default width, the display text for the
 * column header, and also a name for the property that will be set in the
 * cell renderer. (The latter cannot contain special characters.) */
#define NAME_COLUMN_DEFAULT_WIDTH	120
#define SCORE_COLUMN_DEFAULT_WIDTH	40
#define ID_COLUMN_DEFAULT_WIDTH		40
#define START_COLUMN_DEFAULT_WIDTH	50
#define SEQ_COLUMN_DEFAULT_WIDTH	40
#define END_COLUMN_DEFAULT_WIDTH	80

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
    SEQUENCE_COL,
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
    GtkWidget *blxWindow;	  /* The main blixem window that this view belongs to */
    GtkCellRenderer *renderer;	  /* The cell renderer that renders the sequences */
    GtkAdjustment *adjustment;	  /* The scroll adjustment control for the detail view */

    GtkWidget *header;		  /* Contains all the widgets in the detail view header */
    GtkWidget *feedbackBox;	  /* A text box that feeds back info to the user about the currently selected items */
    GList *columnList;		  /* A list of details about all the columns in the detail view */
    GHashTable *seqTable;	  /* Hash table linking a sequence name to a SubjectSequence struct. */
    
    GList *fwdStrandTrees;	  /* A list of all the trees that show the forward strand of the ref seq */
    GList *revStrandTrees;	  /* A list of all the trees that show the reverse strand of the ref seq */
    
    BlxSeqType seqType;		  /* The match type, i.e. dna or peptide */
    int numFrames;	  /* The number of reading frames */
    int cellXPadding;		  /* The x padding between the tree cell background area and their drawing area */
    int cellYPadding;		  /* The y padding between the tree cell background area and their drawing area */
        
    IntRange displayRange;	  /* The currently-displayed range of bases in the reference sequence */
    int selectedBaseIdx;	  /* The currently-selected index in the display range */
    int selectedFrame;		  /* The reading frame to display selected bases for */
    int selectedBaseNum;	  /* The currently-selected base within the selected reading frame */
    int selectedDnaBaseIdx;	  /* The currently-selected index in terms of the DNA sequence */
    PangoFontDescription *fontDesc; /* The fixed-width font that will be used to display the alignments */

    gboolean sortInverted;	  /* Whether the sort operations operate in the reverse direction to their default */
    gboolean highlightDiffs;	  /* Whether the 'highlight differences' option is enabled */
    gboolean showSnpTrack;	  /* Whether the 'Show SNP track' option is enabled */

    /* Cached font sizes, needed often for calculations. */
    int charHeight;
    int charWidth;
    
    /* Display properties */
    GdkColor refSeqColour;	      /* background colour for reference sequence base */
    GdkColor refSeqColourSelected;    /* background colour for reference sequence base (when base selected) */
    GdkColor matchColour;	      /* background colour for base that matches */
    GdkColor matchColourSelected;     /* background colour for base that matches (when base selected) */
    GdkColor mismatchColour;	      /* background colour for base that does not match */
    GdkColor mismatchColourSelected;  /* background colour for base that does not match (when base selected) */
    GdkColor consColour;	      /* background colour for peptide that matches similar type */
    GdkColor consColourSelected;      /* background colour for peptide that matches similar type (when base selected) */
    GdkColor exonColour;	      /* background colour for exon base */
    GdkColor exonColourSelected;      /* background colour for exon base (when base selected) */
    GdkColor insertionColour;	      /* background colour for insertion marker in match sequence */
    GdkColor insertionColourSelected; /* background colour for insertion marker in match sequence (when position selected) */
    GdkColor exonBoundaryColourStart; /* line colour for exon boundaries (marking the start of an exon) */
    GdkColor exonBoundaryColourEnd;   /* line colour for exon boundaries (marking the end of an exon) */
    GdkColor highlightTripletColour;  /* For codon triplets, highlight all the bases in the selected triplet in this colour */
    GdkColor highlightDnaBaseColour;  /* For codon triplets, highlight the specific selected DNA base in this colour */
    GdkColor highlightedRowColour;    /* Background colour for rows that have been grouped for highlighting */
    GdkColor snpColour;		      /* Background colour for SNPs or bases that have SNPs associated with them */
    GdkColor snpColourSelected;	      /* Background colour for SNPs when they are selected */
    
    int exonBoundaryLineWidth;	      /* line width for exon boundaries */
    GdkLineStyle exonBoundaryLineStyleStart; /* line style for exon boundaries (marking the start of an exon) */
    GdkLineStyle exonBoundaryLineStyleEnd;   /* line style for exon boundaries (marking the end of the exon) */
  } DetailViewProperties;


/* Public function declarations */
int			detailViewGetNumFrames(GtkWidget *detailView);
IntRange*		detailViewGetDisplayRange(GtkWidget *detailView);
int			detailViewGetSelectedBaseIdx(GtkWidget *detailView);
int			detailViewGetOldSelectedBaseIdx(GtkWidget *detailView);
GtkAdjustment*		detailViewGetAdjustment(GtkWidget *detailView);
GtkWidget*		detailViewGetTree(GtkWidget *detailView, const Strand strand, const int frame);
GtkWidget*		detailViewGetTreeContainer(GtkWidget *detailView, const Strand strand, const int frame);
gboolean		detailViewGetDisplayRev(GtkWidget *detailView);
PangoFontDescription*	detailViewGetFontDesc(GtkWidget *detailView);
GtkCellRenderer*	detailViewGetRenderer(GtkWidget *detailView);
int			detailViewGetCellXPadding(GtkWidget *detailView);
int			detailViewGetCellYPadding(GtkWidget *detailView);
BlxSeqType		detailViewGetSeqType(GtkWidget *detailView);
IntRange*		detailViewGetRefSeqRange(GtkWidget *detailView);
GtkWidget*	        detailViewGetBlxWindow(GtkWidget *detailView);
int			detailViewGetCharWidth(GtkWidget *detailView);
int			detailViewGetCharHeight(GtkWidget *detailView);
GList*			detailViewGetSequenceMsps(GtkWidget *detailView, const char *seqName);
GHashTable*		detailViewGetSeqTable(GtkWidget *detailView);
gboolean		detailViewGetMatchesSquashed(GtkWidget *detailView);
gboolean		detailViewGetSortInverted(GtkWidget *detailView);
gboolean		detailViewGetHighlightDiffs(GtkWidget *detailView);
gboolean		detailViewGetShowSnpTrack(GtkWidget *detailView);
GList*			detailViewGetColumnList(GtkWidget *detailView);
DetailViewColumnInfo*	detailViewGetColumnInfo(GtkWidget *detailView, const ColumnId columnId);

DetailViewProperties*	detailViewGetProperties(GtkWidget *widget);

int			detailViewGetColumnWidth(GtkWidget *detailView, const ColumnId columnId);

int			getBaseIndexAtColCoords(const int x, const int y, const int charWidth, const IntRange const *displayRange);

GdkColor*		detailViewGetSnpColour(GtkWidget *detailView, const gboolean selected);

void			prevMatch(GtkWidget *detailView, GList *seqNameList);
void			nextMatch(GtkWidget *detailView, GList *seqNameList);
void			firstMatch(GtkWidget *detailView, GList *seqNameList);
void			lastMatch(GtkWidget *detailView, GList *seqNameList);
void			goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType);
void			setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void			setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void			scrollDetailViewLeft1(GtkWidget *detailView);
void			scrollDetailViewRight1(GtkWidget *detailView);
void			scrollDetailViewLeftStep(GtkWidget *detailView);
void			scrollDetailViewRightStep(GtkWidget *detailView);
void			scrollDetailViewLeftPage(GtkWidget *detailView);
void			scrollDetailViewRightPage(GtkWidget *detailView);

void			detailViewSortByType(GtkWidget *detailView, const SortByType sortByType);

void			zoomDetailView(GtkWidget *detailView, const gboolean zoomIn);
void			detailViewSetSelectedBaseIdx(GtkWidget *detailView, const int selectedBaseIdx, const int frame, const int baseNum, const gboolean allowScroll, const gboolean scrollMinimum);
void			updateFeedbackBox(GtkWidget *detailView);
void			toggleStrand(GtkWidget *detailView);

void			detailViewAddMspData(GtkWidget *detailView, MSP *mspList);
void			updateDetailViewFontDesc(GtkWidget *detailView);
void			updateSeqColumnSize(GtkWidget *detailView);
void			resizeDetailViewHeaders(GtkWidget *detailView);
void			refreshDetailViewHeaders(GtkWidget *detailView);

void			detailViewSquashMatches(GtkWidget *detailView, const gboolean squash);
void			detailViewSetSortInverted(GtkWidget *detailView, const gboolean invert);
void			detailViewSetHighlightDiffs(GtkWidget *detailView, const gboolean highlightDiffs);
void			detailViewSetShowSnpTrack(GtkWidget *detailView, const gboolean showSnpTrack);
void			detailViewToggleSnpTrack(GtkWidget *detailView);

GtkWidget*		createSnpTrackHeader(GtkBox *parent, GtkWidget *detailView, const int strand);
void			refreshTextHeader(GtkWidget *widget, gpointer data);
gboolean		onExposeDnaTrack(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);

void			selectClickedSnp(GtkWidget *snpTrack,
					 GtkWidget *colHeader,
					 GtkWidget *detailView, 
					 const int xIn, 
					 const int yIn, 
					 const gboolean reCentre,
					 const gboolean expandSnps,
					 const int clickedBase);

void			seqColHeaderSetRow(GtkWidget *header, const int frame);
int			seqColHeaderGetRow(GtkWidget *header);
int			seqColHeaderGetBase(GtkWidget *header, const int frame, const int numFrames);

GtkWidget*		createDetailView(GtkWidget *blxWindow,
					 GtkWidget *container,
					 GtkAdjustment *adjustment, 
					 GtkWidget *fwdStrandGrid, 
					 GtkWidget *revStrandGrid,
					 MSP *mspList,
					 BlxBlastMode mode,
					 BlxSeqType seqType,
					 int numFrames,
					 const char const *refSeqName,
					 const int startCoord,
					 const gboolean sortInverted,
					 const SortByType sortByType);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *blxWindow);


#endif /* _detail_view_included_ */
