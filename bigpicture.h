/*
 *  bigpicture.h
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#ifndef _big_picture_included_
#define _big_picture_included_

#include <SeqTools/bigpicturegrid.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/blixem_.h>



typedef struct _BigPictureProperties
  {
    GtkWidget *blxWindow;	/* The main blixem window that this grid belongs to */
    GtkWidget *header;		/* The grid header */
    GtkWidget *fwdStrandGrid;	/* The grid that displays the forward ref seq strand */
    GtkWidget *revStrandGrid;	/* The grid that displays the reverse ref seq strand */
    GtkWidget *fwdExonView;	/* The section showing the exons for the forward ref seq strand */
    GtkWidget *revExonView;	/* The section showing the exons for the reverse ref seq strand */
    
    GSList *roundValues;	/* List of "nice" values to round to, for the display values in the grid header */
    int initialZoom;		/* Multiple to multiply the detail view display range by to get the initial big picture display range */
    
    IntRange displayRange;	/* The currently-displayed range in the big picture */
    
    int numHCells;		/* The number of cells in the grid horizontally */
    int basesPerCell;		/* The number of bases show per horizontal cell */
    int roundTo;		/* The number of bases to round grid lines to the nearest multiple of */

    int numVCells;		/* The number of cells in the grid vertically */
    int idPerCell;		/* The percent ID to show per vertical cell */
    IntRange percentIdRange;	/* The max and min %ID values displayed */

    int previewBoxCentre;	/* The base that the preview box is centered on (or UNSET_INT if no preview box) */
    
    int leftBorderChars;	/* The number of characters in the left border of the big picture grids */
    int charWidth;		/* The width of the characters in the big picture grids */
    int charHeight;		/* The height of the characters in the big picture grids */
    
    int highlightBoxMinWidth;
    int previewBoxLineWidth;
    int highlightBoxYPad;       /* Vertical padding between the highlight box and the grid */
  } BigPictureProperties;

typedef struct _GridHeaderProperties
  {
    GtkWidget *bigPicture;  /* The big picture view that this header belongs to */
    GtkWidget *refButton;   /* A reference button, so we can query properties like its height */
    
    GdkRectangle headerRect; /* The actual drawing area where we'll draw the labels */
    int numHeaderLines;	    /* The number of lines of text in the header labels */
    int markerHeight;	    /* The height of the marker lines between the grid and the header labels */
    int headerYPad;	    /* Y padding around the header buttons */
  } GridHeaderProperties;


/* Public function declarations */
void			      setGdkColorBlack(GdkColor *color);
void			      setGdkColorBlue(GdkColor *color);
void			      setGdkColorYellow(GdkColor *color);
void			      setGdkColorCyan(GdkColor *color);

BigPictureProperties*	      bigPictureGetProperties(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetBlxWindow(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveExonView(GtkWidget *bigPicture);
gboolean		      bigPictureGetDisplayRev(GtkWidget *bigPicture);
IntRange*		      bigPictureGetDisplayRange(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetGridHeader(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetDetailView(GtkWidget *bigPicture);
BlxSeqType		      bigPictureGetSeqType(GtkWidget *bigPicture);
int			      bigPictureGetNumFrames(GtkWidget *bigPicture);
int			      bigPictureGetIdPerCell(GtkWidget *bigPicture);
IntRange*		      bigPictureGetPercentIdRange(GtkWidget *bigPicture);
int			      bigPictureGetNumVCells(GtkWidget *bigPicture);
BlxViewContext*		      bigPictureGetContext(GtkWidget *bigPicture);

void			      bigPictureSetIdPerCell(GtkWidget *bigPicture, const int idPerCell);
void			      bigPictureSetMaxPercentId(GtkWidget *bigPicture, const int newValue);
void			      bigPictureSetMinPercentId(GtkWidget *bigPicture, const int newValue);
void                          bigPictureSetPreviewBoxCentre(GtkWidget *bigPicture, int previewBoxCentre);

void			      calculateGridHeaderBorders(GtkWidget *header);
void			      refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean resizeHighlightBox);
void			      calculateNumVCells(GtkWidget *bigPicture);
void			      bigPictureRedrawAll(GtkWidget *bigPicture);

gint                          convertWidgetPosToBaseIdx(const gint widgetPos, const GdkRectangle const *displayRect, const IntRange const *displayRange);

void                          drawPreviewBox(GtkWidget *bigPicture, GdkDrawable *drawable, GdkGC *gc, GdkRectangle *displayRect, GdkRectangle *highlightRect);
void                          showPreviewBox(GtkWidget *bigPicture, const int x);
void                          acceptAndClearPreviewBox(GtkWidget *bigPicture, const int xCentre, GdkRectangle *displayRect, GdkRectangle *highlightRect);

void			      zoomBigPicture(GtkWidget *bigPicture, const gboolean zoomIn);
void			      zoomWholeBigPicture(GtkWidget *bigPicture);

gdouble			      pixelsPerBase(const gint displayWidth, 
					    const IntRange const *displayRange);

gint			      convertBaseIdxToGridPos(const gint baseIdx, 
						      const GdkRectangle const *gridRect, 
						      const IntRange const *displayRange);

int			      getLeftCoordFromCentre(const int centreCoord, 
						     const int width, 
						     const GdkRectangle *outerRect);

int			      getRightCoordFromCentre(const int centreCoord, 
						      const int width, 
						      const GdkRectangle *outerRect);

void			      refreshGridOrder(GtkWidget *bigPicture);

GtkWidget*		      createBigPicture(GtkWidget *blxWindow,
					       GtkWidget *container,
					       GtkWidget **fwdStrandGrid, 
					       GtkWidget **revStrandGrid,
					       const int bigPictZoom,
					       const int lowestId);


#endif /* _big_picture_included_ */
