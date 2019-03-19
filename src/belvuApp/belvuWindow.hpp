/*  File: belvuWindow.h
 *  Author: Gemma Barson, 2011-04-11
 *  Copyright [2018-2019] EMBL-European Bioinformatics Institute
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
 * Description: The main Belvu window
 *----------------------------------------------------------------------------
 */

#ifndef _belvuwindow_h_included_
#define _belvuwindow_h_included_

#include <belvuApp/belvu_.hpp>
#include <gtk/gtk.h>

gboolean              createBelvuWindow(BelvuContext *bc, BlxMessageData *msgData);
void                  showAnnotationWindow(BelvuContext *bc);

GtkUIManager*         createUiManager(GtkWidget *window, BelvuContext *bc, GtkActionGroup **actionGroupOut);
GtkWidget*            createBelvuMenu(GtkWidget *window, const char *path, GtkUIManager *ui_manager);

GtkActionGroup*       belvuWindowGetActionGroup(GtkWidget *belvuWindow);

gboolean              onButtonPressBelvu(GtkWidget *window, GdkEventButton *event, gpointer data);

void                  onRowSelectionChanged(BelvuContext *bc);
void                  onColSelectionChanged(BelvuContext *bc);
void                  onTreeOrderChanged(BelvuContext *bc);

void                  incrementFontSize(BelvuContext *bc);
void                  decrementFontSize(BelvuContext *bc);

#endif /* _belvuwindow_h_included_ */
