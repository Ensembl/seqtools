/*
 *  exonview.h
 *  blixem
 *
 *  Created by Gemma Barson on 24/12/2009.
 *
 */

#ifndef _exon_view_included_
#define _exon_view_included_

#include <gtk/gtk.h>
#include <blixemApp/blixem_.h>

#define BIG_PICTURE_EXON_VIEW_NAME		"BigPictureExonView"


/* Public function declarations */
GtkWidget*	createExonView(GtkWidget *bigPicture, const BlxStrand strand);

gboolean	exonViewGetExpanded(GtkWidget *exonView);
void		exonViewSetExpanded(GtkWidget *exonView, const gboolean expanded);
void		exonViewToggleExpanded(GtkWidget *exonView);

void            callFuncOnAllBigPictureExonViews(GtkWidget *widget, gpointer data);
void		calculateExonViewHeight(GtkWidget *exonView);
void            calculateExonViewHighlightBoxBorders(GtkWidget *exonView);

#endif /* _exon_view_included_ */
