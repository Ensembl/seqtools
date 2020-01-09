/*  File: belvuAlignment.h
 *  Author: Gemma Barson, 2011-04-12
 *  Copyright [2018-2020] EMBL-European Bioinformatics Institute
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
 * Description: Draws the alignment section in the main Belvu window
 *----------------------------------------------------------------------------
 */

#ifndef _belvualignment_h_included_
#define _belvualignment_h_included_

#include <belvuApp/belvu_.hpp>
#include <gtk/gtk.h>

GtkWidget*              createBelvuAlignment(BelvuContext *bc, const char *title, const int wrapWidth);
void                    belvuAlignmentRedrawAll(GtkWidget *belvuAlignment);
void                    belvuAlignmentRefreshAll(GtkWidget *belvuAlignment);
void                    updateOnAlignmentLenChanged(GtkWidget *belvuAlignment);

void                    onBelvuAlignmentFontSizeChanged(GtkWidget *belvuAlignment);

void                    updateOnVScrollSizeChaged(GtkWidget *belvuAlignment);
void                    centerHighlighted(BelvuContext *bc, GtkWidget *belvuAlignment);

void                    removeSelectedSequence(BelvuContext *bc, GtkWidget *belvuAlignment);
void                    removeGappySeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removePartialSeqs(BelvuContext *bc, GtkWidget *belvuAlignment);
void                    removeRedundantSeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removeOutliers(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removeByScore(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);

void                    vScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start);
void                    hScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start);
void                    vScrollPageUpDown(GtkWidget *belvuAlignment, const gboolean up);
void                    hScrollPageLeftRight(GtkWidget *belvuAlignment, const gboolean left);
void                    vScrollUpDown(GtkWidget *belvuAlignment, const gboolean up, const int numRows);
void                    hScrollLeftRight(GtkWidget *belvuAlignment, const gboolean left, const int numChars);

int                     belvuAlignmentGetWidth(GtkWidget *belvuAlignment);

void                    calculateDrawingSizes(GtkWidget *belvuAlignment);
void                    updateHeaderColumnsSize(GtkWidget *belvuAlignment);


#endif /* _belvualignment_h_included_ */
