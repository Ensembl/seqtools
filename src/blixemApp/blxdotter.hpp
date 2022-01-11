/*  File: blxdotter.h
 *  Author: Gemma Barson, 2010-02-03
 *  Copyright [2018-2022] EMBL-European Bioinformatics Institute
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
 * Description: Provides functionality to call Dotter as an external
 *              application.
 *
 *              Determines which sequences to call Dotter with based on the
 *              current selection in Blixem. Provides a dialog allowing the
 *              user to set certain parameters and optionally calculates a
 *              sensible coordinate range to call Dotter with.
 *
 *              Generally Dotter will be called with one alignment sequence vs
 *              a relevant portion of the reference sequence. It can also be
 *              called on the reference sequence vs itself. This is useful if
 *              the reference sequence has been constructed from several
 *              sequences concatenated together and one wishes to identify
 *              alignments between those sequences.
 *----------------------------------------------------------------------------
 */

#ifndef _blx_dotter_h_included_
#define _blx_dotter_h_included_


#include <gtk/gtk.h>
#include <blixemApp/blixem_.hpp>

void			showDotterDialog(GtkWidget *blxWindow, const gboolean bringToFront);
gboolean		callDotterOnSelectedSeqs(GtkWidget *blxWindow, const gboolean hspsOnly, const gboolean sleep, const DotterRefType refType, GError **error);


#endif /* _blx_dotter_h_included_ */
