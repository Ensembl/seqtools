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


typedef struct _GridProperties
  {
    GtkWidget *bigPicture;   /* The big picture that this grid belongs to */
    
    BlxStrand strand;	     /* Whether this grid shows the fwd or rev strand of the ref sequence. */
    
    int mspLineHeight;	     /* The height of the msp lines */
    
    int gridYPadding;	     /* The y padding around the grid */
    int cellYPadding;	     /* The y padding of the grid cells around the text height */
    
    gulong exposeHandlerId;  /* The handler ID for the expose event */
    gboolean ignoreSizeAlloc; /* Flag to indicate that we should ignore size allocation events */
    
    GdkRectangle gridRect;
    GdkRectangle displayRect;
    GdkRectangle highlightRect;
  } GridProperties;


/* Public function declarations */
GridProperties*	    gridGetProperties(GtkWidget *widget);
BlxStrand	    gridGetStrand(GtkWidget *grid);
GtkWidget*	    gridGetBigPicture(GtkWidget *grid);

void		    calculateGridBorders(GtkWidget *grid);
void		    calculateHighlightBoxBorders(GtkWidget *grid);

void		    callFuncOnAllBigPictureGrids(GtkWidget *widget, 
						 gpointer data);

gint		    convertValueToGridPos(GtkWidget *grid, 
					  const gint value);

GtkWidget*	    createBigPictureGrid(GtkWidget *bigPicture, 
					 BlxStrand strand);

#endif /* _big_picture_grid_included_ */
