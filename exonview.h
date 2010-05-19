/*
 *  exonview.h
 *  blixem
 *
 *  Created by Gemma Barson on 24/12/2009.
 *
 */

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>

/* Public function declarations */
GtkWidget*	createExonView(GtkWidget *bigPicture, const BlxStrand strand);

gboolean	exonViewGetExpanded(GtkWidget *exonView);
void		exonViewSetExpanded(GtkWidget *exonView, const gboolean expanded);
void		exonViewToggleExpanded(GtkWidget *exonView);

void		calculateExonViewHeight(GtkWidget *exonView);
