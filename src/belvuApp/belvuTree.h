/*  File: belvuTree.h
 *  Author: Gemma Barson, 2011-05-06
 *  Copyright (c) 2011 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
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


#include <belvuApp/belvu_.h>
#include <gtk/gtk.h>



void                      createAndShowBelvuTree(BelvuContext *bc);
void                      createBelvuTreeWindow(BelvuContext *bc, TreeNode *treeHead);
void                      showTreeSettingsDialog(GtkWidget *belvuTree);

GtkWidget*                createTreeSettingsDialogContent(BelvuContext *bc, GtkWidget *dialog, 
                                                          double *treeScale, double *lineWidth, 
                                                          gboolean *showBranchLen, gboolean *showOrganism,
                                                          BelvuPickMode *pickMode);

void                      treeBootstrap(BelvuContext *bc);


#endif /* _belvutree_h_included_ */