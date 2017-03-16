/*  File: seqtoolsExonView.h
 *  Author: Gemma Barson, 2009-12-24
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: Widget to draw the exons and introns in Dotter.
 *
 *              seqtoolsExonView.c/.h were copied from exonview.c/.h. They
 *              contain much of the same code so the aim was to make these
 *              files generic enough for both Blixem and Dotter to use and to
 *              put them in the seqtoolsUtils library. Unfortunately due to
 *              time constraints this didn't happen, so these files still 
 *              contain some Dotter- and Blixem- specific stuff: in Dotter
 *              we need to draw the crosshair over the exon view; in Blixem
 *              we need to draw the highlight box to indicate the current
 *              detail-view range.
 *----------------------------------------------------------------------------
 */

#ifndef _seqtools_exon_view_included_
#define _seqtools_exon_view_included_

#include <gtk/gtk.h>
#include <seqtoolsUtils/blxmsp.hpp>
#include <dotterApp/dotter_.hpp>

#define SEQTOOLS_EXON_VIEW_NAME		"SeqtoolsExonView"
#define DEFAULT_EXON_HEIGHT			 10
#define DEFAULT_EXON_HEIGHT_BUMPED		 7
#define	DEFAULT_EXON_YPAD			 10
#define	DEFAULT_EXON_YPAD_BUMPED		 4


/* Public function declarations */
GtkWidget *createDotterExonView(GtkWidget *parent, 
				GtkCallback refreshFunc,
                                const gboolean horizontal,
				const BlxStrand strand, 
				DotterWindowContext *dwc,
				const int width,
                                const int height,
				const IntRange* const qRange,
                                const gboolean drawCrosshair,
                                GtkWidget **exonViewOut);

//gboolean	exonViewGetBumped(GtkWidget *exonView);
//void		exonViewSetBumped(GtkWidget *exonView, const gboolean bumped);
void		exonViewToggleBumped(GtkWidget *exonView);
//
//void            callFuncOnAllChildExonViews(GtkWidget *widget, gpointer data);
void		calculateDotterExonViewHeight(GtkWidget *exonView);
void		calculateDotterExonViewBorders(GtkWidget *exonView, const int width, const int height);
void            exonViewSetShowCrosshair(GtkWidget *exonView, const gboolean showCrosshair);
void            exonViewPrepareForPrinting(GtkWidget *exonView);

#endif /* _seqtools_exon_view_included_ */
