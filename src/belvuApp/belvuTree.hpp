/*  File: belvuTree.h
 *  Author: Gemma Barson, 2011-05-06
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
 * Description: Takes care of drawing the tree windows for belvu
 *----------------------------------------------------------------------------
 */


#ifndef _belvutree_h_included_
#define _belvutree_h_included_


#include <belvuApp/belvu_.hpp>
#include <gtk/gtk.h>


#define BELVU_TREE_WINDOW_NAME                  "BelvuTreeWindow"
#define DEFAULT_TREE_SCALE_CORR			0.3	/* default scale for methods using distance correction */
#define DEFAULT_TREE_SCALE_NON_CORR		1.0	/* default scale for methods not using distance correction */


GtkWidget*                createAndShowBelvuTree(BelvuContext *bc, const gboolean isMainTree);
GtkWidget*                createBelvuTreeWindow(BelvuContext *bc, Tree *tree, const gboolean isMainTree);
void                      destroyTreeContents(Tree *tree);
void                      destroyTree(Tree **tree);

void                      belvuTreeRemakeTree(GtkWidget *belvuTree);
void                      onBelvuTreeFontSizeChanged(GtkWidget *belvuTree);

GtkActionGroup*           belvuTreeGetActionGroup(GtkWidget *belvuTree);

void                      treeBootstrap(BelvuContext *bc);
void                      belvuTreeRedrawAll(gpointer belvuTree, gpointer data);
BelvuContext*             belvuTreeGetContext(GtkWidget *belvuTree);

void                      showTreeSettingsDialog(GtkWidget *window, BelvuContext *bc);

GtkWidget*                createTreeSettingsDialogContent(BelvuContext *bc, GtkWidget *dialog, const gboolean isMainTree,
                                                          double *treeScale, double *lineWidth,
                                                          gboolean *showBranchLen, gboolean *showOrganism,
                                                          BelvuPickMode *pickMode, BelvuBuildMethod *buildMethod, BelvuDistCorr *distCorr);


void                      treePrintOrthologs(BelvuContext *bc, GtkWidget *treeWindow);

#endif /* _belvutree_h_included_ */
