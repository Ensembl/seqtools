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
#include <SeqTools/blxwindow.h>


#define HEADER_CONTAINER_NAME		"header container"
#define SNP_TRACK_HEADER_NAME		"SNP track header"
#define DNA_TRACK_HEADER_NAME		"DNA track header"


/* This struct describes a column in the detail view. Multiple widgets (i.e. headers
 * and tree columns) in the detail view must all have columns that share the same
 * properties (namely the column width). */
typedef struct _DetailViewColumnInfo
  {
    BlxColumnId columnId;		/* the column identifier */
    GtkWidget *headerWidget;	/* the header widget for this column (in the detail-view header) */
    GtkCallback refreshFunc;	/* the function that will be called on the header widget when columns are refreshed */
    char *title;		/* the default column title */
    char *propertyName;		/* the property name (used to set the data for the SequenceCellRenderer) */
    char *sortName;             /* the name to display in the sort-by drop-down box (NULL if the view is not sortable on this column) */
    
    int width;			/* the column width */
    gboolean dataLoaded;        /* whether the data for this column has been loaded from the EMBL file (or tried to be loaded, if it doesn't exist) */
  } DetailViewColumnInfo;


/* This struct contains info about canonical splice sites */
typedef struct _BlxSpliceSite
  {
    char donorSite[3];                  /* The bases expected at the donor splice site */
    char acceptorSite[3];               /* The bases expected at the acceptor splice site */

    char donorSiteRev[3];               /* Same as donorSite but reversed */
    char acceptorSiteRev[3];            /* Same as acceptorSite but reversed */
    
    gboolean bothReqd;                  /* Whether both donor and acceptor sites must match in order to be considered canonical */
  } BlxSpliceSite;


/* Essential info required by the the detail view */
typedef struct _DetailViewProperties
  {
    GtkWidget *blxWindow;	  /* The main blixem window that this view belongs to */
    GtkCellRenderer *renderer;	  /* The cell renderer that renders the sequences */
    GtkAdjustment *adjustment;	  /* The scroll adjustment control for the detail view */

    GtkWidget *header;		  /* Contains all the widgets in the detail view header */
    GtkWidget *feedbackBox;	  /* A text box that feeds back info to the user about the currently selected items */
    GtkWidget *statusBar;	  /* A status bar that feeds back info to the user about the currently moused-over items */
    GList *columnList;		  /* A list of details about all the columns in the detail view */
    
    GList *fwdStrandTrees;	  /* A list of all the trees that show the forward strand of the ref seq */
    GList *revStrandTrees;	  /* A list of all the trees that show the reverse strand of the ref seq */
    
    int cellXPadding;		  /* The x padding between the tree cell background area and their drawing area */
    int cellYPadding;		  /* The y padding between the tree cell background area and their drawing area */
        
    IntRange displayRange;	  /* The currently-displayed range of bases in the reference sequence */
    int selectedBaseIdx;	  /* The currently-selected index in the display range */
    int selectedFrame;		  /* The reading frame to display selected bases for */
    int selectedBaseNum;	  /* The currently-selected base within the selected reading frame */
    int selectedDnaBaseIdx;	  /* The currently-selected index in terms of the DNA sequence */
    BlxStrand selectedStrand;	  /* BlxStrand of the tree that the last-selected  */
    PangoFontDescription *fontDesc; /* The fixed-width font that will be used to display the alignments */

    int snpConnectorHeight;	  /* The height of the connector between the SNP track and the DNA base track */
    int numUnalignedBases;        /* If the display-unaligned-sequence option is on, this specifies how many additional bases to show at each end of the alignment */

    /* Cached font sizes, needed often for calculations. */
    int charHeight;
    int charWidth;
        
    int exonBoundaryLineWidth;		     /* line width for exon boundaries */
    GdkLineStyle exonBoundaryLineStyleStart; /* line style for exon boundaries (marking the start of an exon) */
    GdkLineStyle exonBoundaryLineStyleEnd;   /* line style for exon boundaries (marking the end of the exon) */
    
    GSList *spliceSites;           /* List of splice sites that can be found and highlighted by Blixem */
  } DetailViewProperties;


/* Public function declarations */
int			detailViewGetNumFrames(GtkWidget *detailView);
IntRange*		detailViewGetDisplayRange(GtkWidget *detailView);
int			detailViewGetSelectedBaseIdx(GtkWidget *detailView);
int			detailViewGetOldSelectedBaseIdx(GtkWidget *detailView);
GtkAdjustment*		detailViewGetAdjustment(GtkWidget *detailView);
GtkWidget*		detailViewGetTree(GtkWidget *detailView, const BlxStrand strand, const int frame);
GtkWidget*		detailViewGetTreeContainer(GtkWidget *detailView, const BlxStrand strand, const int frame);
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
int                     detailViewGetNumUnalignedBases(GtkWidget *detailView);
GList*			detailViewGetColumnList(GtkWidget *detailView);
DetailViewColumnInfo*	detailViewGetColumnInfo(GtkWidget *detailView, const BlxColumnId columnId);
int			detailViewGetActiveFrame(GtkWidget *detailView);
BlxStrand		detailViewGetSelectedStrand(GtkWidget *detailView);
void			detailViewSetSelectedStrand(GtkWidget *detailView, BlxStrand strand);

DetailViewProperties*	detailViewGetProperties(GtkWidget *widget);

int			detailViewGetColumnWidth(GtkWidget *detailView, const BlxColumnId columnId);
void                    detailViewGetColumnXCoords(GtkWidget *detailView, const BlxColumnId columnId, IntRange *xRange);

int			getBaseIndexAtColCoords(const int x, const int y, const int charWidth, const IntRange const *displayRange);

void			prevMatch(GtkWidget *detailView, GList *seqList);
void			nextMatch(GtkWidget *detailView, GList *seqList);
void			firstMatch(GtkWidget *detailView, GList *seqList);
void			lastMatch(GtkWidget *detailView, GList *seqList);
void			goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType);
void			setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void			setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void			scrollDetailViewLeft1(GtkWidget *detailView);
void			scrollDetailViewRight1(GtkWidget *detailView);
void			scrollDetailViewLeftStep(GtkWidget *detailView);
void			scrollDetailViewRightStep(GtkWidget *detailView);
void			scrollDetailViewLeftPage(GtkWidget *detailView);
void			scrollDetailViewRightPage(GtkWidget *detailView);

void			detailViewSetSortColumn(GtkWidget *detailView, const BlxColumnId sortColumn);

void			zoomDetailView(GtkWidget *detailView, const gboolean zoomIn);
void			detailViewSetSelectedBaseIdx(GtkWidget *detailView, const int selectedBaseIdx, const int frame, const int baseNum, const gboolean allowScroll, const gboolean scrollMinimum);
void			updateFeedbackBox(GtkWidget *detailView);
void			toggleStrand(GtkWidget *detailView);

void			detailViewAddMspData(GtkWidget *detailView, MSP *mspList);

void			updateDetailViewFontDesc(GtkWidget *detailView);
void			updateSeqColumnSize(GtkWidget *detailView);
void			resizeDetailViewHeaders(GtkWidget *detailView);
void			refreshDetailViewHeaders(GtkWidget *detailView);
void			detailViewRedrawAll(GtkWidget *detailView);

void			detailViewUpdateSquashMatches(GtkWidget *detailView, const gboolean squash);
void			detailViewUpdateSortInverted(GtkWidget *detailView, const gboolean invert);
void			detailViewUpdateShowSnpTrack(GtkWidget *detailView, const gboolean showSnpTrack);
void                    detailViewUpdateLimitUnalignedBases(GtkWidget *detailView, const gboolean limitUnalignedBases);

void                    detailViewSetNumUnalignedBases(GtkWidget *detailView, const int numBases);
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

GHashTable*             getRefSeqBasesToHighlight(GtkWidget *detailView, const IntRange const *displayRange, const BlxSeqType seqType, const BlxStrand strand);

void                    drawColumnSeparatorLine(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const BlxViewContext *bc);
gboolean                onExposeGenericHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);

void			drawHeaderChar(BlxViewContext *bc,
				       DetailViewProperties *properties,
				       const int dnaIdx,
				       const char baseChar,
				       const BlxStrand strand, 
                                       const int frame,
				       const BlxSeqType seqType,
				       const gboolean displayIdxSelected, 
				       const gboolean dnaIdxSelected, 
				       const gboolean showBackground,
				       const gboolean showSnps,
				       const gboolean showCodons,
                                       const BlxColorId defaultBgColor,
				       GdkDrawable *drawable,
				       GdkGC *gc,
				       const int x,
				       const int y,
                                       GHashTable *intronBases);

GtkWidget*		createDetailView(GtkWidget *blxWindow,
                                         GtkContainer *parent,
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
					 const BlxColumnId sortColumn,
                                         const gboolean optionalDataLoaded);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *blxWindow);


#endif /* _detail_view_included_ */
