/*  File: bigpicturegrid.h
 *  Author: Gemma Barson, 2009-11-23
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
 * Description: Creates a grid on the Big Picture section. A grid shows the
 *              alignments on a particular strand of the reference sequence.
 *              The x coord of the grid is the reference sequence coord and the
 *              y coord is the %ID of the alignment. An alignment is shown as
 *              a horizontal line between the start and end coords of the
 *              match, at the y coord corresponding to the relevant %ID.
 *----------------------------------------------------------------------------
 */

#ifndef _big_picture_grid_included_
#define _big_picture_grid_included_


#include <gtk/gtk.h>
#include <blixemApp/blixem_.hpp>
#include <seqtoolsUtils/utilities.hpp>


#define BIG_PICTURE_GRID_NAME		"BigPictureGrid"


class GridProperties
{
public:
  GridProperties(GtkWidget *widget_in,
                 GtkWidget *bigPicture_in,
                 gulong exposeHandlerId_in,
                 BlxStrand strand_in);

  GtkWidget *widget;       /* The grid widget */
  GtkWidget *bigPicture;   /* The big picture that this grid belongs to */

  BlxStrand strand;	     /* Whether this grid shows the fwd or rev strand of the ref sequence. */

  int mspLineHeight;	     /* The height of the msp lines */

  int gridYPadding;	     /* The y padding around the grid */

  gulong exposeHandlerId;  /* The handler ID for the expose event */
  gboolean ignoreSizeAlloc; /* Flag to indicate that we should ignore size allocation events */

  GdkRectangle gridRect;
  GdkRectangle displayRect;
  GdkRectangle highlightRect;
};


/* Public function declarations */
GridProperties*	    gridGetProperties(GtkWidget *widget);
BlxStrand	    gridGetStrand(GtkWidget *grid);
GtkWidget*	    gridGetBigPicture(GtkWidget *grid);

void		    calculateGridBorders(GtkWidget *grid);
void		    calculateGridHighlightBoxBorders(GtkWidget *grid);

void		    callFuncOnAllBigPictureGrids(GtkWidget *widget,
						 gpointer data);

gint		    convertValueToGridPos(GtkWidget *grid,
					  const gdouble value);

void                gridPrepareForPrinting(GtkWidget *grid);

GtkWidget*	    createBigPictureGrid(GtkWidget *bigPicture,
					 BlxStrand strand);

#endif /* _big_picture_grid_included_ */
