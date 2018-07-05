/*  File: bigpicture.h
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
 * Description: The "big picture" section displays an overview of the entire
 *              portion of reference sequence that was passed to Blixem. The
 *              reference sequence coordinates are shown along the top.
 *
 *              Alignments are shown as lines and are placed in a scaled grid
 *              (one for each strand) against their %ID (see bigpicturegrid.h).
 *
 *              Exons and introns are drawn in standard notation in the
 *              transcript section (see exonview.h).
 *----------------------------------------------------------------------------
 */

#ifndef _big_picture_included_
#define _big_picture_included_

#include <blixemApp/bigpicturegrid.hpp>
#include <blixemApp/blxwindow.hpp>
#include <blixemApp/blixem_.hpp>
#include <blixemApp/blxpanel.hpp>



class BigPictureProperties : public BlxPanel
{
public:
  // Constructors
  BigPictureProperties(GtkWidget *bigPicture_in,
                       GtkWidget *blxWindow_in,
                       BlxContext *bc,
                       CoverageViewProperties *coverageViewP_in,
                       GtkWidget *header_in,
                       GtkWidget *fwdStrandGrid_in,
                       GtkWidget *revStrandGrid_in,
                       GtkWidget *fwdExonView_in,
                       GtkWidget *revExonView_in,
                       int previewBoxCentre_in,
                       const IntRange* const initRange_in,
                       const IntRange* const fullRange_in,
                       const int initialZoom_in,
                       const gdouble lowestId_in);

  ~BigPictureProperties();

  // Access
  double contentXPos() const;
  double contentWidth() const;

  // Query
  const IntRange *highlightRange() const;

  // Update
  void refreshDisplayRange(const bool keepCentered);
  void onHighlightBoxMoved(const int displayIdx, const BlxSeqType seqType);

  // Member variables
  GtkWidget *header;		/* The grid header */
  GtkWidget *fwdStrandGrid;	/* The grid that displays the forward ref seq strand */
  GtkWidget *revStrandGrid;	/* The grid that displays the reverse ref seq strand */
  GtkWidget *fwdExonView;	/* The section showing the exons for the forward ref seq strand */
  GtkWidget *revExonView;	/* The section showing the exons for the reverse ref seq strand */

  GSList *roundValues;	        /* List of "nice" values to round to, for the display values in the grid header */
  int initialZoom;		/* Multiple to multiply the detail view display range by to get the initial big picture display range */

  int numHCells;		/* The number of cells in the grid horizontally */
  int basesPerCell;		/* The number of bases show per horizontal cell */
  int roundTo;	           	/* The number of bases to round grid lines to the nearest multiple of */

  int numVCells;		/* The number of cells in the grid vertically */
  gdouble idPerCell;		/* The percent ID to show per vertical cell */
  DoubleRange percentIdRange;	/* The max and min %ID values displayed */

  int leftBorderChars;          /* The number of characters in the left border of the big picture grids */

};

class GridHeaderProperties
{
public:
  GtkWidget *widget;          /* The grid header widget */
  GtkWidget *bigPicture;      /* The big picture view that this header belongs to */
  GtkWidget *refButton;       /* A reference button, so we can query properties like its height */

  GdkRectangle headerRect;    /* The actual drawing area where we'll draw the labels */
  int numHeaderLines;	      /* The number of lines of text in the header labels */
  int markerHeight;	      /* The height of the marker lines between the grid and the header labels */
  int headerYPad;	      /* Y padding around the header buttons */
};


/* Public function declarations */
void			      setGdkColorBlack(GdkColor *color);
void			      setGdkColorBlue(GdkColor *color);
void			      setGdkColorYellow(GdkColor *color);
void			      setGdkColorCyan(GdkColor *color);

BigPictureProperties*	      bigPictureGetProperties(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetBlxWindow(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveExonView(GtkWidget *bigPicture);
gboolean		      bigPictureGetDisplayRev(GtkWidget *bigPicture);
IntRange*		      bigPictureGetDisplayRange(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetGridHeader(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetDetailView(GtkWidget *bigPicture);
BlxSeqType		      bigPictureGetSeqType(GtkWidget *bigPicture);
int			      bigPictureGetNumFrames(GtkWidget *bigPicture);
gdouble			      bigPictureGetIdPerCell(GtkWidget *bigPicture);
DoubleRange*		      bigPictureGetPercentIdRange(GtkWidget *bigPicture);
int			      bigPictureGetNumVCells(GtkWidget *bigPicture);
BlxContext*		      bigPictureGetContext(GtkWidget *bigPicture);

gboolean                      bigPictureSetIdPerCell(GtkWidget *bigPicture, const gdouble idPerCell);
gboolean		      bigPictureSetMaxPercentId(GtkWidget *bigPicture, const gdouble newValue);
gboolean		      bigPictureSetMinPercentId(GtkWidget *bigPicture, const gdouble newValue);

void			      calculateGridHeaderBorders(GtkWidget *header);
void			      refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean keepCentered);
void                          calculateBigPictureCellSize(GtkWidget *bigPicture, BigPictureProperties *properties);
void			      calculateNumVCells(GtkWidget *bigPicture);
void			      bigPictureRedrawAll(GtkWidget *bigPicture);
void                          bigPicturePrepareForPrinting(GtkWidget *bigPicture);

void                          scrollBigPictureLeftStep(GtkWidget *bigPicture);
void                          scrollBigPictureRightStep(GtkWidget *bigPicture);

gint			      bigPictureGetCellHeight(GtkWidget *bigPicture);

void			      drawVerticalGridLines(GdkRectangle *drawingRect, GdkRectangle *highlightRect,
						    const int yPadding, BlxContext *bc,
						    BigPictureProperties *bpProperties, GdkDrawable *drawable);

void			      drawHorizontalGridLines(GtkWidget *widget, GtkWidget *bigPicture,
						      GdkRectangle *drawingRect, BlxContext *bc,
						      BigPictureProperties *bpProperties, GdkDrawable *drawable,
						      const gint numCells, const gdouble rangePerCell,
						      const gdouble maxVal, const gboolean abbrev, const char *unit);

void			      zoomBigPicture(GtkWidget *bigPicture, const gboolean zoomIn);
void			      zoomWholeBigPicture(GtkWidget *bigPicture);

int			      getRightCoordFromCentre(const int centreCoord,
						      const int width,
						      const GdkRectangle *outerRect);

void			      refreshGridOrder(GtkWidget *bigPicture);

GtkWidget*		      createBigPicture(GtkWidget *blxWindow,
                                               BlxContext *bc,
					       GtkContainer *parent,
					       GtkWidget **fwdStrandGrid,
					       GtkWidget **revStrandGrid,
                                               const IntRange* const initRange,
                                               const IntRange* const fullRange,
					       const int bigPictZoom,
					       const gdouble lowestId);


#endif /* _big_picture_included_ */
