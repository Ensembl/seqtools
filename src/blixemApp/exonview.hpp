/*  File: exonview.h
 *  Author: Gemma Barson, 2009-12-24
 *  Copyright [2018] EMBL-European Bioinformatics Institute
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
 * Description: A section of the big picture showing exons, introns and basic
 *              (box-shaped) features for a particular strand of the reference
 *              sequence. The naming is all based around exons but has been
 *              expanded to draw any box-shape feature.
 *
 *              Pending: It would be good to consolidate the Blixem and Dotter
 *              exonview stuff - see the comment in seqtoolsExonView.h
 *----------------------------------------------------------------------------
 */

#ifndef _exon_view_included_
#define _exon_view_included_

#include <gtk/gtk.h>
#include <blixemApp/blixem_.hpp>

#define BIG_PICTURE_EXON_VIEW_NAME		"BigPictureExonView"


/* Public function declarations */
GtkWidget*	createExonView(GtkWidget *bigPicture, const BlxStrand strand);

void            exonViewPrepareForPrinting(GtkWidget *exonView);

gboolean	exonViewGetExpanded(GtkWidget *exonView);
void		exonViewSetExpanded(GtkWidget *exonView, const gboolean expanded);
void		exonViewToggleExpanded(GtkWidget *exonView);

void            callFuncOnAllBigPictureExonViews(GtkWidget *widget, gpointer data);
void		calculateExonViewHeight(GtkWidget *exonView);
void            calculateExonViewHighlightBoxBorders(GtkWidget *exonView);

#endif /* _exon_view_included_ */
