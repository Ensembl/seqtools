/*  File: blxpanel.hpp
 *  Author: Gemma Barson, 2016-04-07
 *  Copyright [2018] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2016 Genome Research Ltd
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
 * Description: A base class for the big-picture and detail-view panels
 *----------------------------------------------------------------------------
 */

#ifndef _blx_panel_included_
#define _blx_panel_included_


#include <gtk/gtk.h>
#include <seqtoolsUtils/utilities.hpp>


class CoverageViewProperties;
class BlxContext;
class IntRange;



class BlxPanel
{
public:
  // Construct/destruct
  BlxPanel();

  BlxPanel(GtkWidget *widget,
           GtkWidget *blxWindow,
           BlxContext *bc_in,
           CoverageViewProperties *coverageViewP,
           int previewBoxCentre);

  virtual ~BlxPanel();

  // Access
  GtkWidget* widget() const;
  GtkWidget* blxWindow() const;
  virtual double charWidth() const;
  virtual double charHeight() const;
  GtkWidget* coverageView();
  CoverageViewProperties *coverageViewProperties();
  const GList *columnList() const;

  virtual double contentXPos() const = 0;
  virtual double contentWidth() const = 0;
  virtual const IntRange* highlightRange() const { return NULL; };

  // Update
  virtual void refreshDisplayRange(const bool keepCentered);
  virtual void onHighlightBoxMoved(const int displayIdx, const BlxSeqType seqType);

  // Draw
  void startPreviewBox(const int x, const gboolean init, const int offset);
  void finishPreviewBox(const int xCentreIn, GdkRectangle *displayRect, GdkRectangle *highlightRect);
  void drawPreviewBox(GdkDrawable *drawable, GdkRectangle *displayRect, GdkRectangle *highlightRect);


  // Public member variables
  IntRange displayRange;       /* The currently-displayed range of ref-sequence bases in the panel */

private:
  GtkWidget *m_widget;           /* The widget that draws this panel  */
  GtkWidget *m_blxWindow;      /* The main blixem window that this panel belongs to */
  BlxContext *m_bc;
  CoverageViewProperties *m_coverageViewP;

  bool m_previewBoxOn;         /* True if currently displaying the preview box */
  int m_previewBoxCentre;      /* The base that the preview box is centered on */
  int m_previewBoxOffset;      /* Can be used to offset the actual preview box position from previewBoxCentre */
  int m_previewBoxLineWidth;
};



#endif // _blx_panel_included_
