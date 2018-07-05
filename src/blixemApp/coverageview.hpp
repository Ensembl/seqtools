/*  File: coverageview.h
 *  Author: Gemma Barson, 2011-03-21
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
 * Description: This view shows a plot of the number of alignments at each
 *              coordinate in the reference sequence.
 *----------------------------------------------------------------------------
 */

#ifndef _coverage_view_h_included_
#define _coverage_view_h_included_

#include "blixemApp/blixem_.hpp"
#include <gtk/gtk.h>


class BlxPanel;



class CoverageViewProperties
{
public:
  CoverageViewProperties(GtkWidget *widget_in,
                         GtkWidget *blxWindow_in,
                         BlxContext *bc);

  /* Access */
  GtkWidget* widget();
  const IntRange *displayRange();
  const IntRange *highlightRange();
  double contentXPos() const;
  double contentWidth() const;
  double depthPerCell();
  int maxLabeledDepth();

  /* Modify */
  gboolean setDepthPerCell(const double depthPerCell);

  /* Events */
  gboolean expose(GdkEventExpose *event, gpointer data);
  gboolean buttonPress(GdkEventButton *event, gpointer data);
  gboolean buttonRelease(GdkEventButton *event, gpointer data);
  gboolean mouseMove(GdkEventMotion *event, gpointer data);
  gboolean scroll(GdkEventScroll *event, gpointer data);

  /* Updates */
  void setPanel(BlxPanel *panel);
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
  BlxContext *m_bc;
  BlxPanel *m_panel;        /* The parent panel that the coverage view is part of */

  GtkWidget *m_blxWindow;   /* The main blixem window */

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
