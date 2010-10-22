/*
 *  seqtoolsExonView.h
 *  seqtools
 *
 *  Created by Gemma Barson on 24/12/2009.
 *
 */

#ifndef _seqtools_exon_view_included_
#define _seqtools_exon_view_included_

#include <gtk/gtk.h>
#include <SeqTools/blxmsp.h>
#include <SeqTools/dotter_.h>

#define SEQTOOLS_EXON_VIEW_NAME		"SeqtoolsExonView"
#define DEFAULT_EXON_HEIGHT			 10
#define DEFAULT_EXON_HEIGHT_BUMPED		 7
#define	DEFAULT_EXON_YPAD			 10
#define	DEFAULT_EXON_YPAD_BUMPED		 4


/* Public function declarations */
GtkWidget *createDotterExonView(GtkWidget *parent, 
				GtkCallback refreshFunc,
				const BlxStrand strand, 
				DotterContext *dc,
				DotterWindowContext *dwc,
				const int width,
				const IntRange const *qRange,
                                const gboolean drawCrosshair,
                                GtkWidget **exonViewOut);

//gboolean	exonViewGetBumped(GtkWidget *exonView);
//void		exonViewSetBumped(GtkWidget *exonView, const gboolean bumped);
//void		exonViewToggleBumped(GtkWidget *exonView);
//
//void            callFuncOnAllChildExonViews(GtkWidget *widget, gpointer data);
void		calculateDotterExonViewHeight(GtkWidget *exonView);
void		calculateDotterExonViewBorders(GtkWidget *exonView, const int width);
void            exonViewSetShowCrosshair(GtkWidget *exonView, const gboolean showCrosshair);

#endif /* _seqtools_exon_view_included_ */