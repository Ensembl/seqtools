/*  File: dotter_.h
 *  Author: esr, 1999-08-26
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
 * Description: Internal header for Dotter package-wide code
 *----------------------------------------------------------------------------
 */

#ifndef _dotter_p_h_included_
#define _dotter_p_h_included_

#include <dotterApp/dotter.hpp>
#include <seqtoolsUtils/version.hpp>
#include <config.h>

#define NR                    23                        /* Not A residue */
#define NA                    24                        /* Not A residue */
#define NUM_READING_FRAMES    3                         /* number of reading frames for protein sequences */

#define DEFAULT_MAJOR_TICK_HEIGHT                   10    /* default height for major tick marks on the scale */
#define DEFAULT_MINOR_TICK_HEIGHT                   5     /* default height for minor tick marks on the scale */
#define DEFAULT_X_PADDING                           0     /* default horizontal padding around the dot plot */
#define DEFAULT_Y_PADDING                           0     /* default vertical padding around the dot plot */
#define SCALE_LINE_WIDTH                            1     /* the line width in pixels of the scale boundary and tick marks */
#define MAX_MATRIX_NAME_LENGTH                      500   /* max length of the matrix name */
#define NUM_COLORS                                  256   /* the number of colors in our greyscale */

extern char *stdcode1[];        /* 1-letter amino acid translation code */


/*            dotter program version and information.                        */
#define DOTTER_TITLE   "Dotter"
#define DOTTER_PREFIX  "Dotter - "
#define DOTTER_PREFIX_ABBREV "D: "
#define DOTTER_DESC    "Dot-matrix plotter for detailed comparision of two sequences."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define DOTTER_VERSION_STRING      SEQTOOLS_VERSION
#define DOTTER_PACKAGE_VERSION     UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, SEQTOOLS_VERSION)
#define DOTTER_TITLE_STRING        UT_MAKE_TITLE_STRING(DOTTER_TITLE, DOTTER_VERSION_STRING)
#define DOTTER_VERSION_COMPILE     DOTTER_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define DOTTER_COPYRIGHT_STRING    UT_MAKE_COPYRIGHT_STRING("2010-2015")
#define DOTTER_WEBSITE_STRING      "http://www.sanger.ac.uk/resources/software/seqtools/"
#define DOTTER_LICENSE_STRING      UT_MAKE_LICENCE_STRING(DOTTER_TITLE)



/* Dotter error domain */
#define DOTTER_ERROR g_quark_from_string("Dotter")

/* Error codes */
typedef enum
  {
    DOTTER_ERROR_INVALID_LENGTH,         /* an invalid alignment length was specified */
    DOTTER_ERROR_INVALID_WIN_SIZE,       /* an invalid sliding window size was specified */
    DOTTER_ERROR_OPENING_FILE,           /* error opening file */
    DOTTER_ERROR_READING_FILE,           /* error reading file */
    DOTTER_ERROR_SAVING_FILE             /* error saving file */
  } DotterError;


typedef enum _DotterHspMode
  {
    DOTTER_HSPS_OFF = 0,            /* HSPs not shown. Make this option 0 so we can use it as a boolean to check if HSPs are shown or not. */
    DOTTER_HSPS_LINE,               /* draw HSPs as solid lines */
    DOTTER_HSPS_FUNC,               /* draw HSPs as a function of their score */
    DOTTER_HSPS_GREYSCALE           /* HSPs are drawn as greyscale, overriding the dot-plot */
  } DotterHspMode;

/* The following are used to define default colors for certain types of features in Dotter.
 * One of several different actual colors from the BlxColor struct may be used depending
 * on state, e.g. we use a different color if "print colors" (i.e. black and
 * white mode) is on. */

typedef enum
{
  DOTCOLOR_MIN,                             /* dummy value so that we don't get a zero ID */

  DOTCOLOR_BACKGROUND,                      /* background color for widgets */
  DOTCOLOR_MATCH,                           /* background color for an exact match */
  DOTCOLOR_CONS,                            /* background color for a conserved match */
  DOTCOLOR_MISMATCH,                        /* background color for a mismatch */
  DOTCOLOR_EXON_FILL,                       /* fill color for an exon in the exon view */
  DOTCOLOR_EXON_LINE,                       /* line color for an exon in the exon view */
  DOTCOLOR_CDS_FILL,                        /* fill color for a CDS in the exon view (coding region) */
  DOTCOLOR_CDS_LINE,                        /* line color for a CDS in the exon view (coding region) */
  DOTCOLOR_UTR_FILL,                        /* fill color for an UTR in the exon view (non-coding/untranslated region) */
  DOTCOLOR_UTR_LINE,                        /* line color for an UTR in the exon view (non-coding/untranslated region) */
  DOTCOLOR_CROSSHAIR,                       /* color of the dotplot crosshair */
  DOTCOLOR_GRID,                            /* color of the dotplot grid */
  DOTCOLOR_BORDER,                          /* color of the dotplot border */
  DOTCOLOR_MARKER_FILL,                     /* fill color of the position markers on the greyramp */
  DOTCOLOR_MARKER_LINE,                     /* line color of the position markers on the greyramp */
  DOTCOLOR_THRESHOLD_MARKER,                /* line color of the threshold marker on the greyramp */
  DOTCOLOR_BREAKLINE,                       /* the color of break-lines between sequences */
  DOTCOLOR_CANONICAL,
  DOTCOLOR_NON_CANONICAL,

  DOTCOLOR_NUM_COLORS
} DotterColorId;


/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the DotterWindowContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    DOTDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */

    DOTDIALOG_HELP,                 /* The Help dialog */
    DOTDIALOG_SETTINGS,             /* The Settings dialog */

    DOTDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } DotterDialogId;


#define CONS_MATRIX_SIZE 24

/* General properties for the dotter session. These settings are common to all dotter windows
 * open for the same process (i.e. if you start a sub-dotter from within dotter it will inherit
 * these properties). */
typedef struct _DotterContext
{
  BlxBlastMode blastMode;                   /* the blast mode used */
  BlxSeqType displaySeqType;                /* whether the display shows sequences as nucleotide or peptide */
  int numFrames;                            /* the number of reading frames */
  char **geneticCode;                       /* code used to translate nucleotide triplets to amino acids */
  gint32 matrix[CONS_MATRIX_SIZE][CONS_MATRIX_SIZE];        /* matrix for determining conserved matches */
  char *matrixName;                         /* matrix name */
  MSP *mspList;                             /* list of all MSPs in dotter */
  GList *seqList;                           /* list of all matches sequences as BlxSequences */
  GSList *windowList;                       /* list of all windows that use this context */

  gboolean watsonOnly;                      /* true if we only display the watson strand */
  gboolean crickOnly;                       /* true if we only display the crick strand */

  PangoFontDescription *fontDesc;           /* fixed-width font to use for alignment */
  gdouble charWidth;                        /* the fixed-width font width */
  gdouble charHeight;                       /* the fixed-width font height */

  char *refSeqName;                         /* the reference sequence name */
  char *refSeq;                             /* the passed-in reference sequence */
  char *refSeqRev;                          /* the reverse-complement of refSeq (or NULL if not applicable) */
  IntRange refSeqFullRange;                 /* the full reference sequence range passed to dotter on startup */
  BlxStrand refSeqStrand;                   /* which strand of the ref seq we were passed */
  BlxSeqType refSeqType;                    /* whether refSeq is in nucleotide or peptide coords */

  char* peptideSeqs[NUM_READING_FRAMES];    /* Array of 3 strings containing the 3 frame translations of the reference sequence */

  char *matchSeqName;                       /* the match sequence name */
  char *matchSeq;                           /* the match sequence */
  char *matchSeqRev;                        /* the revsere-complemented match sequence (or NULL if not applicable) */
  IntRange matchSeqFullRange;               /* the full match sequence range passed to dotter on startup */
  BlxStrand matchSeqStrand;                 /* which strand of the match seq we were passed */
  BlxSeqType matchSeqType;                  /* whether matchSeq is in nucleotide or peptide coords */

  gboolean hozScaleRev;                     /* true if horizontal coords should increase from right-to-left rather than left-to-right */
  gboolean vertScaleRev;                    /* true if vertical coords should increase from bottom-to-top rather than top-to-bottom */
  gboolean negateCoords;                    /* negate displayed coords if the scale is reversed, i.e. so coords still appear to increase left-to-right */
  gboolean displayMirror;                   /* whether to display a mirror image in self comparisons */
  gboolean abbrevTitle;                     /* abbreviate window titles to save space */

  double memoryLimit;                       /* maximum Mb allowed for dotplot */

  int scaleWidth;                           /* width of the dotplot scale */
  int scaleHeight;                          /* height of the dotplot scale */

  GArray *defaultColors;                    /* Default colors used by Dotter */

  BlxMessageData *msgData;                  /* Data to be passed to message handlers */
} DotterContext;


/* Struct containing properties specific to a particular dotter window.  */
typedef struct _DotterWindowContext
  {
    DotterContext *dotterCtx;

    IntRange refSeqRange;                     /* displayed ref seq range currently displayed in this window */
    IntRange matchSeqRange;                   /* displayed match seq range currently displayed in this window */

    int refCoord;                             /* currently-selected ref seq coord */
    int matchCoord;                           /* currently-selected match seq coord */

    gdouble zoomFactor;                       /* we divide by this to scale the dotplot, i.e. <1 zooms in, >1 zooms out */
    gboolean selfComp;                        /* true if we're comparing the same portion ref seq against itself */

    GtkWidget *dialogList[DOTDIALOG_NUM_DIALOGS]; /* list of all persistent dialogs for this window */

    GtkPrintSettings *printSettings;          /* print settings */
    GtkPageSetup *pageSetup;                  /* page setup for printing */
    gboolean usePrintColors;                  /* whether to use print colours */

    GtkUIManager *uiManager;                  /* the ui manager for this dotter window */
    GtkActionGroup *actionGroup;              /* the menu actions for this window */
  } DotterWindowContext;


class DotplotScaleProperties
{
public:
  GtkWidget *widget;                  /* The dotplot scale widget */
  int basesPerMark;                   /* how many bases per major tickmark */
  int basesPerSubmark;                /* how many bases per minor tickmark */
  int numMarks;                       /* how many major tickmarks there are */
  int numSubmarks;                    /* how many minor tickmarks there are */
  int firstSubmarkPos;                /* the position of the first tickmark */
  int firstSubmarkCoord;              /* the coordinate the first tickmark represents */
  int startCoord;                     /* the overall start coord */
  int endCoord;                       /* the overall end coord */
  GArray *labels;
};


class DotplotProperties
{
public:
  GtkWidget *widget;                  /* The dotplot widget */
  DotterWindowContext *dotterWinCtx;  /* pointer to the dotter context for the window that this tool belongs to */
  GdkRectangle plotRect;              /* the space where the dot plot will be */
  GtkWidget *hozExons1;               /* shows main strand exons for horizontal sequence */
  GtkWidget *hozExons2;               /* shows other strand exons for horizontal sequence */
  GtkWidget *vertExons1;              /* shows main strand exons for vertical sequence */
  GtkWidget *vertExons2;              /* shows other strand exons for vertical sequence */

  gulong greyMap[NUM_COLORS];         /* maps weight -> pixel value. fixed mapping in pseudo colour displays
                                         variable mapping in truecolor displays */
  GdkColor greyRamp[NUM_COLORS];      /* 256 grey colors, black->white, only used in true color displays */
  GdkColormap *colorMap;              /* the greyramp colormap */

  int imageWidth;
  int imageHeight;
  GdkImage *image;                    /* the greyramp image */

  double expResScore;
  int pixelFac;
  int slidingWinSize;

  /* Dynamic properties: */
  unsigned char *pixelmap;            /* source data for drawing the dot-plot */
  unsigned char *hspPixmap;           /* source data for drawing the HSP dot-plot */

  gboolean crosshairOn;               /* whether to show the crosshair that marks the position of the currently-selected coord */
  gboolean crosshairCoordsOn;         /* whether to display the crosshair label */
  gboolean crosshairFullscreen;       /* whether to show the crosshair over the whole widget or just within the dot-plot rectangle */

  gboolean pixelmapOn;                /* whether to show the dot-plot pixelmap or not */
  DotterHspMode hspMode;              /* how (and whether) to show high-scoring pairs from Blast */

  gboolean gridlinesOn;               /* whether to show grid lines */
  gboolean breaklinesOn;              /* whether to show break-lines between sequences */
  gboolean hozLabelsOn;               /* whether to show labels for features on the horizontal sequence */
  gboolean vertLabelsOn;              /* whether to show labels for features on the vertical sequence */

  GdkPoint dragStart;                 /* start point for mid-click drag */
  GdkPoint dragEnd;                   /* end point for mid-click drag */

  const char *exportFileName;         /* file to export the plot to if in batch "export to pdf" mode */

  DotplotScaleProperties hozScale;    /* horizontal scale properties */
  DotplotScaleProperties vertScale;   /* vertical scale properties */
};




const char*         dotterGetAppName(void);
const char*         dotterGetTitlePrefix(DotterContext *dc);
const char*         dotterGetCopyrightString(void);
const char*         dotterGetWebSiteString(void);
const char*         dotterGetCommentsString(void);
const char*         dotterGetLicenseString(void);
const char*         dotterGetVersionString(void);

int                 winsizeFromlambdak(int mtx[CONS_MATRIX_SIZE][CONS_MATRIX_SIZE],
                                       int *tob, int abetsize, const char *qseq, const char *sseq,
                                       double *exp_res_score, double *Lambda);

void                argvAdd(int *argc, char ***argv, const char *s);


/* dotter.c */
int                 getResFactor(DotterContext *dc, const gboolean horizontal);
int                 getDisplayCoord(const int coordIn, DotterContext *dc, const gboolean horizontal);
void                copyIntToDefaultClipboard(const int val);
int                 getStartCoord(DotterWindowContext *dwc, const gboolean horizontal);
int                 getEndCoord(DotterWindowContext *dwc, const gboolean horizontal);
int                 getSelectedCoord(DotterWindowContext *dwc, const gboolean horizontal);
void                dotterEnableSelectionMenus(DotterWindowContext *dwc, const gboolean enable);

void                callDotterInternal(DotterContext *dc,
                                       const IntRange* const refSeqRange,
                                       const IntRange* const matchSeqRange,
                                       const gdouble zoomFactor,
                                       const gboolean breaklinesOn);

void                dotterSetToggleMenuStatus(DotterWindowContext *dwc,
                                              const char *menuItem,
                                              const gboolean enable);


/* greyramptool.c */
GtkWidget*          createGreyrampTool(DotterWindowContext *dwc, const int bottomVal, const int topVal, const gboolean swapValues, GtkWidget **greyrampWindow_out);
GtkWidget*          createGreyrampToolMinimised(DotterWindowContext *dwc, const int blackPoint, const int whitePoint);
void                registerGreyrampCallback(GtkWidget *greyramp, GtkWidget *widget, GtkCallback func);
void                updateGreyMap(GtkWidget *greyramp);

/* alignmenttool.c */
GtkWidget*          createAlignmentTool(DotterWindowContext *dotterWinCtx, GtkWidget **alignmentWindow_out);
void                updateAlignmentRange(GtkWidget *alignmentTool, DotterWindowContext *dwc);
void                alignmentToolSetSpliceSitesOn(GtkWidget *alignmentTool, const gboolean spliceSitesOn);
gboolean            alignmentToolGetSpliceSitesOn(GtkWidget *alignmentTool);
void                alignmentToolRedrawAll(GtkWidget *alignmentTool);
void                alignmentToolCopySeln(GtkWidget* alignmentTool);
void                alignmentToolCopySelnCoords(GtkWidget *alignmentTool);
void                alignmentToolClearSequenceSelection(GtkWidget *alignmentTool);

/* dotplot.c */
GtkWidget*          createDotplot(DotterWindowContext *dwc,
                                  const char *loadFileName,
                                  const char *saveFileName,
                                  const DotterSaveFormatType saveFormat,
                                  const char *exportFileName,
                                  const gboolean hspsOn,
                                  const gboolean breaklinesOn,
                                  const char *initWinsize,
                                  const int pixelFacIn,
                                  const int zoomFacIn,
                                  const int qcenter,
                                  const int scenter,
                                  GtkWidget **dotplot);

void                dotplotPrepareForPrinting(GtkWidget *dotplot);
int                 convertToDisplayIdx(const int dnaIdx, const gboolean horizontal, DotterContext *dc, const int frame, int *baseNum);
int                 convertToDnaIdx(const int displayIdx, const gboolean horizontal, DotterContext *dc, const int frame, int baseNum);
int                 getDotplotWidth(GtkWidget *dotplot, DotplotProperties *properties);
int                 getDotplotHeight(GtkWidget *dotplot, DotplotProperties *properties);
DotplotProperties*  dotplotGetProperties(GtkWidget *widget);
void                dotplotToggleBumpExons(GtkWidget *dotplot);

DotterHspMode       dotplotGetHspMode(GtkWidget *dotplot);
int                 dotplotGetSlidingWinSize(GtkWidget *dotplot);
gboolean            dotplotSetSlidingWinSize(GtkWidget *dotplot, const int newValue, GError **error);
void                dotplotSetBreaklinesOn(GtkWidget *dotplot, const gboolean breaklinesOn);
void                dotplotSetHozLabelsOn(GtkWidget *dotplot, const gboolean breaklinesOn);
void                dotplotSetVertLabelsOn(GtkWidget *dotplot, const gboolean breaklinesOn);
int                 dotplotGetImageWidth(GtkWidget *dotplot);
int                 dotplotGetImageHeight(GtkWidget *dotplot);
double              dotplotGetExpResScore(GtkWidget *dotplot);
int                 dotplotGetPixelFac(GtkWidget *dotplot);
void                dotplotUpdateGreymap(GtkWidget *dotplot, gpointer data);
void                dotplotUpdateCrosshair(GtkWidget *dotplot, DotterWindowContext *dwc);
void                dotplotTogglePixelmap(GtkWidget *dotplot);
void                dotplotToggleGrid(GtkWidget *dotplot);
void                toggleCrosshairOn(GtkWidget *dotplot);
void                toggleCrosshairCoordsOn(GtkWidget *dotplot);
void                toggleCrosshairFullscreen(GtkWidget *dotplot);
void                setHspMode(GtkWidget *dotplot, DotterHspMode hspMode);
void                refreshDotplot(GtkWidget *dotplot);
void                redrawDotplot(GtkWidget *dotplot);
void                recalcDotplot(GtkWidget *dotplot);
void                savePlot(GtkWidget *dotplot, DotplotProperties *properties,
                             const char *saveFileName, DotterSaveFormatType saveFormat, GError **error);
void                exportPlot(GtkWidget *dotplot, GtkWindow *window, const char *exportFileName, GError **error);
void                loadPlot(GtkWidget *dotplot, const char *loadFileName, GError **error);

GList*              dotterCreateColumns();

#endif /* _dotter_p_h_included_ */
