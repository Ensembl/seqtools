/*  File: blxpanel.cpp
 *  Author: Gemma Barson, 2016-04-07
 *  Copyright [2018-2023] EMBL-European Bioinformatics Institute
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
 * Description: See blxpanel.hpp
 *----------------------------------------------------------------------------
 */

#include <blixemApp/blxpanel.hpp>
#include <blixemApp/blxcontext.hpp>
#include <blixemApp/coverageview.hpp>



// Unnamed namespace
namespace
{

#define DEFAULT_PREVIEW_BOX_LINE_WIDTH  1


/* Convert an x coord in the given rectangle to a base index (in nucleotide coords) */
static gint convertRectPosToBaseIdx(const gint x,
                                    const GdkRectangle* const displayRect,
                                    const IntRange* const dnaDispRange,
                                    const gboolean displayRev)
{
  gint result = UNSET_INT;

  gdouble distFromEdge = (gdouble)(x - displayRect->x);
  int basesFromEdge = (int)(distFromEdge / pixelsPerBase(displayRect->width, dnaDispRange));

  if (displayRev)
    {
      result = dnaDispRange->max() - basesFromEdge;
    }
  else
    {
      result = dnaDispRange->min() + basesFromEdge;
    }

  return result;
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the
 * left edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
int getLeftCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int leftCoord = centreCoord - roundNearest((double)width / 2.0);

  if (outerRect)
    {
      if (leftCoord < outerRect->x)
	leftCoord = outerRect->x;
      else
	{
	  int leftCoordMax = outerRect->x + outerRect->width - width;

	  if (leftCoord > leftCoordMax)
	    {
	      leftCoord = leftCoordMax;
	    }
	}
    }

  return leftCoord;
}


} // Unnamed namespace





/* Default constructor */
BlxPanel::BlxPanel() :
  displayRange(0, 0),
  m_widget(NULL),
  m_blxWindow(NULL),
  m_bc(NULL),
  m_coverageViewP(NULL),
  m_previewBoxOn(false),
  m_previewBoxCentre(0),
  m_previewBoxOffset(0),
  m_previewBoxLineWidth(DEFAULT_PREVIEW_BOX_LINE_WIDTH)
{
};

/* Constructor */
BlxPanel::BlxPanel(GtkWidget *widget_in,
                   GtkWidget *blxWindow_in,
                   BlxContext *bc_in,
                   CoverageViewProperties *coverageViewP_in,
                   int previewBoxCentre_in) :
  displayRange(0, 0),
  m_widget(widget_in),
  m_blxWindow(blxWindow_in),
  m_bc(bc_in),
  m_coverageViewP(coverageViewP_in),
  m_previewBoxOn(false),
  m_previewBoxCentre(previewBoxCentre_in),
  m_previewBoxOffset(0),
  m_previewBoxLineWidth(DEFAULT_PREVIEW_BOX_LINE_WIDTH)
{
}


/* Destructor */
BlxPanel::~BlxPanel()
{
}

GtkWidget* BlxPanel::widget() const
{
  return m_widget;
}

/* Return the main blixem window */
GtkWidget* BlxPanel::blxWindow() const
{
  return m_blxWindow;
}

/* Return the coverage view class */
CoverageViewProperties* BlxPanel::coverageViewProperties()
{
  return m_coverageViewP;
}

/* Return the coverage view widget */
GtkWidget* BlxPanel::coverageView()
{
  GtkWidget *result = NULL;

  if (m_coverageViewP)
    result = m_coverageViewP->widget();

  return result;
}

double BlxPanel::charWidth() const
{
  double result = 0.0;

  if (m_bc)
    result = m_bc->charWidth();

  return result;
}

double BlxPanel::charHeight() const
{
  double result = 0.0;

  if (m_bc)
    result = m_bc->charHeight();

  return result;
}

const GList* BlxPanel::columnList() const
{
  const GList *result = NULL;

  if (m_bc)
    result = m_bc->columnList;

  return result;
}


void BlxPanel::refreshDisplayRange(const bool keepCentered)
{
}


/* Show a preview box centred on the given x coord. If init is true, we're initialising
 * a drag; otherwise, we're already dragging */
void BlxPanel::startPreviewBox(const int x, const gboolean init, const int offset)
{
  if (init)
    {
      m_previewBoxOn = TRUE;
      m_previewBoxOffset = offset;

      static GdkCursor *cursor = NULL;
      cursor = gdk_cursor_new(GDK_FLEUR);

      gdk_window_set_cursor(m_blxWindow->window, cursor);
    }

  /* We might get called by a drag operation where the preview box drag has not been
   * initialised; just ignore it */
  if (!m_previewBoxOn)
    return;

  /* Whether dragging or initialising, we need to update the position */
  m_previewBoxCentre = x + m_previewBoxOffset;

  /* Refresh all child widgets, and also the coverage view (which may not
   * be a direct child of the panel widget) */
  gtk_widget_queue_draw(widget());
  gtk_widget_queue_draw(coverageView());
}


/* Scroll the panel so that it is centred on the current preview box position, and clear
 * the preview box. (Do nothing if preview box is not currently shown.)  */
void BlxPanel::finishPreviewBox(const int xCentreIn, GdkRectangle *displayRect, GdkRectangle *highlightRect)
{
  gdk_window_set_cursor(m_blxWindow->window, NULL);

  /* If we're not currently drawing a preview box then nothing to do */
  if (!m_previewBoxOn)
    return;

  /* Apply any offset required to get the real centre coord */
  const int xCentre = xCentreIn + m_previewBoxOffset;

  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&displayRange, m_bc->seqType, m_bc->numFrames, m_bc->displayRev, &m_bc->refSeqRange, &dnaDispRange);

  /* Find the base index where the new scroll range will start. This is the leftmost
   * edge of the preview box if numbers increase in the normal left-to-right direction,
   * or the rightmost edge if the display is reversed. */
  const int x = getLeftCoordFromCentre(xCentre, highlightRect->width, displayRect);
  int baseIdx = convertRectPosToBaseIdx(x, displayRect, &dnaDispRange, m_bc->displayRev);

  /* Subtract 1 if the display is reversed to give the base to the right of x, rather than the base to the left of x */
  if (m_bc->displayRev)
    --baseIdx;

  /* Reset the preview box status. We don't need to un-draw it because we
   * will redraw the whole big picture below. */
  m_previewBoxOn = FALSE;
  m_previewBoxOffset = 0;

  /* Perform any required updates following the change to the highlight box position. The base
   * index is in terms of the nucleotide coords so we need to convert to display coords */
  const int displayIdx = convertDnaIdxToDisplayIdx(baseIdx, m_bc->seqType, 1, m_bc->numFrames, m_bc->displayRev, &m_bc->refSeqRange, NULL);
  onHighlightBoxMoved(displayIdx, m_bc->seqType);

  /* Re-centre the panel */
  refreshDisplayRange(TRUE);

  gtk_widget_queue_draw(widget());
}

void BlxPanel::onHighlightBoxMoved(const int displayIdx, const BlxSeqType seqType)
{
}

/* Draw the preview box on the given drawable within the boundaries of the given displayRect.
 * The boundaries of the preview box are given by highlightRect.
 * Only does anything if the preview box centre is set. */
void BlxPanel::drawPreviewBox(GdkDrawable *drawable,
                              GdkRectangle *displayRect,
                              GdkRectangle *highlightRect)
{
  /* If not currently drawing the preview box then nothing to do */
  if (!m_previewBoxOn)
    {
      return;
    }

  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&displayRange, m_bc->seqType, m_bc->numFrames, m_bc->displayRev, &m_bc->refSeqRange, &dnaDispRange);

  /* Find the x coord for the left edge of the preview box (or the right edge, if
   * the display is right-to-left). */
  int x = getLeftCoordFromCentre(m_previewBoxCentre, highlightRect->width, displayRect);

  /* Convert it to the base index and back again so that we get it rounded to the position of
   * the nearest base. */
  int baseIdx = convertRectPosToBaseIdx(x, displayRect, &dnaDispRange, m_bc->displayRev);
  int xRounded = convertBaseIdxToRectPos(baseIdx, displayRect, &dnaDispRange, TRUE, m_bc->displayRev, TRUE);

  /* The other dimensions of the preview box are the same as the current highlight box. */
  GdkRectangle previewRect = {xRounded, highlightRect->y, highlightRect->width, highlightRect->height};

  GdkColor *previewBoxColor = getGdkColor(BLXCOLOR_PREVIEW_BOX, m_bc->defaultColors, FALSE, m_bc->usePrintColors);

  drawHighlightBox(drawable, &previewRect, m_previewBoxLineWidth, previewBoxColor);
}
