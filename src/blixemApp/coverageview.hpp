/*  File: coverageview.h
 *  Author: Gemma Barson, 2011-03-21
 *  Copyright (c) 2011 - 2012 Genome Research Ltd
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
 * Description: This view shows a plot of the number of alignments at each
 *              coordinate in the reference sequence.
 *----------------------------------------------------------------------------
 */
 
#ifndef _coverage_view_h_included_
#define _coverage_view_h_included_

#include "blixemApp/blixem_.hpp"
#include <gtk/gtk.h>

class CoverageViewProperties
{
public:
  CoverageViewProperties(GtkWidget *widget_in, 
                         GtkWidget *blxWindow_in, 
                         BlxContext *bc);

  /* Access */
  GtkWidget* widget();
  double depthPerCell();
  int maxLabeledDepth();

  /* Modify */
  gboolean setDepthPerCell(const double depthPerCell);
  void setDisplayRange(const IntRange *displayRange);

  /* Events */
  gboolean expose(GdkEventExpose *event, gpointer data);
  gboolean buttonPress(GdkEventButton *event, gpointer data);
  gboolean buttonRelease(GdkEventButton *event, gpointer data);
  gboolean mouseMove(GdkEventMotion *event, gpointer data);
  gboolean scroll(GdkEventScroll *event, gpointer data);

  /* Updates */
  void updateDepth();
  void calculateBorders();
  void calculateHighlightBoxBorders();

  /* Drawing */
  void draw(GdkDrawable *drawable);
  void redraw();
  void prepareForPrinting();

private:
  double charWidth() const;
  void drawPlot(GdkDrawable *drawable);
  void recalculate();

  GtkWidget *m_widget;      /* The coverage view widget */
  GtkWidget *m_blxWindow;   /* The main blixem window */
  BlxContext *m_bc;
  const IntRange *m_displayRange; /* The current range to display */

  int m_viewYPadding;	     /* The y padding around the view rect */
  double m_numVCells;	     /* The number of cells to show vertically */
  gdouble m_rangePerCell;    /* The range of depth values shown per grid cell on the plot */
    
  GdkRectangle m_viewRect;   /* The rectangle we draw in */
  GdkRectangle m_displayRect; /* The total display area */
  GdkRectangle m_highlightRect; /* The area that is highlighted (which indicates the detail-view
                                 range) */

  int *m_maxDepth;
};


CoverageViewProperties*     createCoverageView(GtkWidget *blxWindow, BlxContext *bc);


#endif /* _coverage_view_h_included_ */
