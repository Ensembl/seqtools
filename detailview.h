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
#include <SeqTools/blxview.h>


#define FIXED_WIDTH_FONT		"Courier"

typedef struct _DetailViewProperties
  {
    GtkWidget *refTree;         /* Reference tree so that we can get generic info about the trees when it is the same for all trees */
    GtkAdjustment *adjustment;  /* The scroll adjustment control for the detail view */
    
    char *refSeq;		/* The reference sequence */
    int numReadingFrames;	/* The number of reading frames in the detail view e.g. 1 for DNA matches, 3 for peptides */
    
    IntRange displayRange;	/* The currently-displayed range of bases in the reference sequence */
    int selectedBaseIdx;	/* The currently-selected base in the reference sequence */
  } DetailViewProperties;


/* Public function declarations */
char*			detailViewGetRefSeq(GtkWidget *detailView);
int			detailViewGetNumReadingFrames(GtkWidget *detailView);
IntRange*		detailViewGetDisplayRange(GtkWidget *detailView);
int			detailViewGetSelectedBaseIdx(GtkWidget *detailView);

DetailViewProperties*	detailViewGetProperties(GtkWidget *widget);

void			setDetailViewScrollPos(GtkAdjustment *adjustment, 
					       int value);

GtkWidget*		createDetailView(GtkWidget *container, 
					 GtkAdjustment *adjustment, 
					 const MSP const *mspList,
					 GtkWidget *fwdStrandGrid, 
					 GtkWidget *revStrandGrid,
					 char *refSeq,
					 int numReadingFrames,
					 IntRange *refSeqRange);

GtkWidget*		createDetailViewScrollBar(GtkAdjustment *adjustment, 
						  GtkWidget *mainWindow);

#endif /* _detail_view_included_ */