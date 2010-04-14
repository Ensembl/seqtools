/*
 *  bigpicturegrid.h
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */


#ifndef _big_picture_grid_included_
#define _big_picture_grid_included_


#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>


#define BIG_PICTURE_GRID_NAME		"BigPictureGrid"
#define DEFAULT_GRID_PERCENT_ID_MAX	100
#define DEFAULT_HIGHLIGHT_BOX_LINE_WIDTH 2


typedef struct _GridProperties
  {
    GtkWidget *bigPicture;   /* The big picture that this grid belongs to */
    
    Strand strand;	     /* Whether this grid shows the fwd or rev strand of the ref sequence. */
    
    int numVCells;	     /* The number of cells in the grid vertically */
    IntRange percentIdRange; /* The currently-displayed range of values along the vertical grid axis */
    int mspLineHeight;	     /* The height of the msp lines */
    
    int gridYPadding;	     /* The y padding around the grid */
    int cellYPadding;	     /* The y padding of the grid cells around the text height */
    int highlightBoxYPad;    /* Vertical padding between the highlight box and the grid */
    
    gulong exposeHandlerId;  /* The handler ID for the expose event */
    gboolean ignoreSizeAlloc; /* Flag to indicate that we should ignore size allocation events */
    
    GdkRectangle gridRect;
    GdkRectangle displayRect;
    GdkRectangle highlightRect;
  } GridProperties;


/* Public function declarations */
GridProperties*	    gridGetProperties(GtkWidget *widget);
Strand		    gridGetStrand(GtkWidget *grid);
GtkWidget*	    gridGetBigPicture(GtkWidget *grid);

void		    calculateGridBorders(GtkWidget *grid);
void		    calculateHighlightBoxBorders(GtkWidget *grid);

void		    callFuncOnAllBigPictureGrids(GtkWidget *widget, 
						 gpointer data);

gint		    convertValueToGridPos(GtkWidget *grid, 
					  const gint value);

void		    showPreviewBox(GtkWidget *grid, const int x);
void		    acceptAndClearPreviewBox(GtkWidget *grid, const int x);

GtkWidget*	    createBigPictureGrid(GtkWidget *bigPicture, 
					 Strand strand);

#endif /* _big_picture_grid_included_ */
