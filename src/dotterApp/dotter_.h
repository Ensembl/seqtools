/*  File: dotter_.h
 *  Author: esr
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
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
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Private include file for dotter
 * Exported functions: 
 * HISTORY:
 * Last edited: Sep 15 08:36 2006 (edgrif)
 * Created: Thu Aug 26 17:17:58 1999 (fw)
 * CVS info:   $Id: dotter_.h,v 1.7 2010-11-03 13:39:09 gb10 Exp $
 *-------------------------------------------------------------------
 */

#ifndef _dotter_p_h_included_
#define _dotter_p_h_included_

#include <dotterApp/dotter.h>

#define NR                    23 		        /* Not A residue */
#define NA                    24 		        /* Not A residue */
#define NUM_READING_FRAMES    3                         /* number of reading frames for protein sequences */

#define DEFAULT_MAJOR_TICK_HEIGHT                   10    /* default height for major tick marks on the scale */
#define DEFAULT_MINOR_TICK_HEIGHT                   5     /* default height for minor tick marks on the scale */
#define DEFAULT_X_PADDING                           0     /* default horizontal padding around the dot plot */
#define DEFAULT_Y_PADDING                           0     /* default vertical padding around the dot plot */
#define SCALE_LINE_WIDTH                            1     /* the line width in pixels of the scale boundary and tick marks */
#define MAX_MATRIX_NAME_LENGTH                      500   /* max length of the matrix name */

extern char *stdcode1[];        /* 1-letter amino acid translation code */


/*            dotter program version and information.                        */
#define DOTTER_TITLE   "Dotter program"
#define DOTTER_DESC    "Dot-matrix plotter for detailed comparision of two sequences."

#define DOTTER_VERSION 4
#define DOTTER_RELEASE 0
#define DOTTER_UPDATE  0
#define DOTTER_VERSION_NUMBER	   UT_MAKE_VERSION_NUMBER(DOTTER_VERSION, DOTTER_RELEASE, DOTTER_UPDATE)
#define DOTTER_VERSION_STRING	   UT_MAKE_VERSION_STRING(DOTTER_VERSION, DOTTER_RELEASE, DOTTER_UPDATE)
#define DOTTER_TITLE_STRING	   UT_MAKE_TITLE_STRING(DOTTER_TITLE, DOTTER_VERSION, DOTTER_RELEASE, DOTTER_UPDATE)
#define DOTTER_VERSION_COMPILE	   DOTTER_VERSION_STRING "  " __TIME__ " "__DATE__

#define DOTTER_COPYRIGHT_STRING	   "Copyright (c) 2010: Genome Research Ltd."
#define DOTTER_WEBSITE_STRING	   ""
#define DOTTER_LICENSE_STRING	   "Dotter is distributed under the GNU Public License, see http://www.gnu.org/copyleft/gpl.txt"

#define DOTTER_AUTHOR_LIST	   " Originally written by Erik Sonnhammer and Richard Durbin",\
" Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define DOTTER_AUTHOR_TEXT	   " Originally written by Erik Sonnhammer and Richard Durbin\n" \
" Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define DOTTER_COMMENTS_STRING(TITLE, VERSION, RELEASE, UPDATE)	\
"("DOTTER_TITLE_STRING", "					\
"compiled on - "__DATE__" "__TIME__")\n\n"			\
DOTTER_AUTHOR_TEXT "\n"


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
  DOTCOLOR_MARKER_FILL,                     /* fill color of the position markers on the greyramp */
  DOTCOLOR_MARKER_LINE,                     /* line color of the position markers on the greyramp */
  DOTCOLOR_THRESHOLD_MARKER,                /* line color of the threshold marker on the greyramp */
  
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


/* General properties for the dotter session. These settings are common to all dotter windows
 * open for the same process (i.e. if you start a sub-dotter from within dotter it will inherit
 * these properties). */
typedef struct _DotterContext
{
  BlxBlastMode blastMode;		    /* the blast mode used */
  BlxSeqType displaySeqType;                /* whether the display shows sequences as nucleotide or peptide */
  int numFrames;			    /* the number of reading frames */
  char **geneticCode;                       /* code used to translate nucleotide triplets to amino acids */
  int matrix[24][24];                       /* matrix for determining conserved matches */
  char *matrixName;                         /* matrix name */
  MSP *mspList;                             /* list of all MSPs in dotter */
  GList *seqList;			    /* list of all matches sequences as BlxSequences */
  GSList *windowList;			    /* list of all windows that use this context */

  gboolean watsonOnly;                      /* true if we only display the watson strand */
  gboolean crickOnly;                       /* true if we only display the crick strand */
  
  PangoFontDescription *fontDesc;	    /* fixed-width font to use for alignment */
  gdouble charWidth;                        /* the fixed-width font width */
  gdouble charHeight;                       /* the fixed-width font height */
  
  char *refSeqName;			    /* the reference sequence name */
  char *refSeq;				    /* the reference sequence forward strand */
  char *refSeqRev;                          /* the reference sequence reverse-complemented (NULL if not applicable) */
  IntRange refSeqFullRange;                 /* the full reference sequence range passed to dotter on startup */
  BlxStrand refSeqStrand;                   /* which strand of the ref seq we were passed */
  BlxSeqType refSeqType;                    /* whether refSeq is in nucleotide or peptide coords */

  char* peptideSeqs[NUM_READING_FRAMES];    /* Array of 3 strings containing the 3 frame translations of the reference sequence */
  
  char *matchSeqName;			    /* the match sequence name */
  char *matchSeq;			    /* the match sequence */
  char *matchSeqRev;			    /* the revsere-complemented match sequence (or NULL if not applicable) */
  IntRange matchSeqFullRange;               /* the full match sequence range passed to dotter on startup */
  BlxStrand matchSeqStrand;                 /* which strand of the match seq we were passed */
  BlxSeqType matchSeqType;                  /* whether matchSeq is in nucleotide or peptide coords */
  
  gboolean hozScaleRev;			    /* true if horizontal coords should increase from right-to-left rather than left-to-right */
  gboolean vertScaleRev;		    /* true if vertical coords should increase from bottom-to-top rather than top-to-bottom */
  gboolean selfComp;                        /* true if we're comparing the ref seq against itself */
  gboolean displayMirror;                   /* whether to display a mirror image in self comparisons */
  
  double memoryLimit;                       /* maximum Mb allowed for dotplot */
  
  int scaleWidth;			    /* width of the dotplot scale */
  int scaleHeight;			    /* height of the dotplot scale */

  GArray *defaultColors;		    /* Default colors used by Dotter */
  gboolean usePrintColors;		    /* Whether to use print colors (i.e. black and white) */
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
    
    GtkWidget *dialogList[DOTDIALOG_NUM_DIALOGS]; /* list of all persistent dialogs for this window */
  } DotterWindowContext;



int		    winsizeFromlambdak(int mtx[24][24], int *tob, int abetsize, const char *qseq, const char *sseq, 
				       double *exp_res_score, double *Lambda);

void		    argvAdd(int *argc, char ***argv, char *s);


/* dotter.c */
int		    getResFactor(DotterContext *dc, const gboolean horizontal);
void                callDotterInternal(DotterContext *dc, 
                                       const IntRange const *refSeqRange,
                                       const IntRange const *matchSeqRange,
                                       const gdouble zoomFactor);

/* greyramptool.c */
GtkWidget*	    createGreyrampTool(DotterContext *dc, const int bottomVal, const int topVal, const gboolean swapValues);
void                registerGreyrampCallback(GtkWidget *greyramp, GtkWidget *widget, GtkCallback func);
void                updateGreyMap(GtkWidget *greyramp);

/* alignmenttool.c */
GtkWidget*	    createAlignmentTool(DotterWindowContext *dotterWinCtx);
void                updateAlignmentRange(GtkWidget *alignmentTool, DotterWindowContext *dwc);


/* dotplot.c */
GtkWidget*	    createDotplot(DotterWindowContext *dwc, 
                                  const char *loadFileName,
                                  const char *saveFileName,
                                  const gboolean hspsOn,
                                  const char *initWinsize,
                                  const int pixelFacIn,
                                  const int zoomFacIn,
                                  const int qcenter,
                                  const int scenter,
                                  GtkWidget **dotplot);

DotterHspMode       dotplotGetHspMode(GtkWidget *dotplot);
int                 dotplotGetSlidingWinSize(GtkWidget *dotplot);
void                dotplotSetSlidingWinSize(GtkWidget *dotplot, const int newValue, GError **error);
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
void                redrawDotplot(GtkWidget *dotplot);
void                refreshDotplot(GtkWidget *dotplot);
void                savePlot(GtkWidget *dotplot, const char *saveFileName, GError **error);
void                loadPlot(GtkWidget *dotplot, const char *loadFileName, GError **error);


#endif /* _dotter_p_h_included_ */
