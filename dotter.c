/*  File: dotter.c
 *  Author: Erik Sonnhammer
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
 *      Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *      Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: See comment block below.
 * Exported functions: See dotter.h
 * HISTORY:
 * Last edited: Aug 26 09:08 2009 (edgrif)
 * * May 13 11:29 1999 (edgrif): Add minor fix from Christian Iseli
 * * Mar 17 16:24 1999 (edgrif): Fixed bug which crashed xace when a
 *              negative alignment length was given.
 * Created: Wed Mar 17 16:23:21 1999 (edgrif)
 * CVS info:   $Id: dotter.c,v 1.13 2010-08-31 15:46:31 gb10 Exp $
 *-------------------------------------------------------------------
 */

/*
   DOTTER - Sequence-sequence dotplot in a pixel background image

-------------------------------------------------------------
|  File: dotter.c                                           |
|  Author: Erik Sonnhammer                                  |
|  Copyright (C) E Sonnhammer and R Durbin, 1994            |
-------------------------------------------------------------

   Memory requirements:
   The pixelmap + (alphabetsize+4) x qlen + 2 x slen


   Date   Modification
--------  ---------------------------------------------------
93-09-04  Created
...
94-10-09  Added DNA support, changed crosshair to BoxShift calls
94-10-10  Got rid of *rawdata since we can't store the full n*m matrix
	  Implemented zoomfactor everywhere
	  ScoreVector speedup trick
94-10-12  Reverse-complement for DNA (watson = top/forward strand)
94-10-13  Saving, loading and batch running
94-11-01  Re-calling dotter with middle mouse rectangle.
          Tidy scales even with offset.
	  Zoom factor according to qlen*slen
94-11-11  Introduced header in dotter files [variable(bytes)]:
             version (1), zoom(4), qlen4(4), slen4(4)
          Fixed signed char bug (-1 !> 127) for Sun, Alpha.
          Rewrote inner loop - 3 times faster.
94-11-15  Changed *data to unsigned char.  char* doesn't work for Sun/Sol/OSF.
          Reversed blastx mode
94-12-05  Display of Exons and Introns.
94-12-13  Speedup by better use of scoreVec rows.
95-07-12  [2.0] user-features with start == end drawn as lines.
95-07-13  Finally tracked down malloc problem of graphPixelsRaw().
          This was incorrect in graphRampTool and graphPixles too.
95-08-24  [2.1] Added from command line choosable score matrix file (Both Protein & DNA)
95-08-24  [2.1] Crosshair coordinates in dot-matrix.
95-09-06  [2.2] Karlin-Altschul statistics for estimating best windowsize.
95-09-07  [2.2] Calculation of score -> pixel factor, given MATRIX, and sequences.
95-09-08  [2.2] Dotplot files now with: win, pixelFac, Matrix, Matrixname.
95-11-01  [2.2] Added X options capability.
          Emergency workaround if Karlin/Altschul statistics fail. Limit range to 10-50.
	  Improved zooming in with findCommand().
	  Added Xoptions at zooming.
	  Added MSPlist at zooming.
96-04-16  [2.3] Default usage of 'installed' (private) colormap, by adding "-install"
          to the argument list.  Disable with -i. New graphlib routines by Darren Platt.
96-04-23  [2.3] Added dotterBinary filename to zoom-parameters, since "which" always seems
          to fail from dotter started in a popen, which blocked double zooming.
96-04-24  [2.3] rd changed graphColorSquares to use explicit tints by force majeure.
96-07-18  [2.3] fixed bug (initAlignment called before initCrosshair) that caused crashes.
96-07-24  [2.3] fixed bug in LeftDown that caused box crash when Crosshair wasn't on.
96-07-28  [2.3] Changed to checkmark menus.
                Fixed bug that HSPpixel map didn't get reset before calculation.
96-09-20 (rbrusk): in dotterRedraw(), (re-)initialize  vLineBox and hLineBox for all
			new dotgraph's
97-02-28  [2.4] Fixed bugs in feature drawing routine (incorrect resfac in DNA-protein),
          incorrect parsing of multiple sequence files.
97-03-19  [2.4] Changed findCommand to search path self instead of relying on 'which dotter'.
97-11-19  [2.5] Added featurefile on command line.
          Added series and width drawing of feature segments.
	  Added crosshair-full-window option.
	  Fixed annotation drawing of segment boxes (used to only work on lines).
97-11-19  [2.6] For selfcomparisons, duplication of the mirror image was made default.
                (full selfcomp matrix is now never calculated.  -D has reversed effect).
97-11-19  [2.6] Added selectFeatures Tool to hide/unhide feature series.
98-01-15  [2.7] Changed Feature annotation to msp->desc (no length limit).
                Fixed some bugs in Feature drawing.
		Window sizes automatically to accommodated feature files.
98-07-08  [2.7] Fixed bug in findCommand: run strtok on a copy and not directly on getenv("PATH").
99-07-07  [3.0] Added support for XY curve shapes PARTIAL and INTERPOLATE.
02-03-14  [3.0] Added quotes around sequence names when calling dotterBinary so that names with pipes
                don't cause failures.
02-03-26  [3.1] Added separator line for first sequence in multi-sequence mode.
                Changed header in dotplot from seqname to filename in multi-sequence mode.


          Pending:

	  change filqueryopen to filqueryeditopen (.dotter) when getting new code.

	  Fix revcomp when zooming in.

	  Add to alignment tool the extent and score of the current window.

	  (Raise initial cutoffs when compressed. ?)

	  HSP drawing bugged for reverse matches and for reversed scale
	  (Huge job to fix ... do only if really necessary)

      	  Score matrix for DNA w/ transversions/transitions...? literature?

-------------------------------------------------------------------------------- */

/* CPU usage profiling on SGI:

   cc dotter.c, no optimization
   pixie -o dotter.pixie dotter
   dotter.pixie -wD -b t seq seq
   prof -pixie -h -only calcWindow dotter dotter.Addrs dotter.Counts
*/

#include <ctype.h>
#include <wh/regular.h>
#include <wh/aceio.h>
#include <wh/graph.h>
#include <wh/gex.h>
#include <wh/key.h>
#include <wh/menu.h>
#include <SeqTools/blixem_.h>
#include <SeqTools/dotter_.h>
#include <SeqTools/utilities.h>


/* tint stuff used to be in graph.h, now local - rd 960524
   NB colours as Erik liked them before Jean tinkered!
   could rename #define's more sensibly now
*/

#define TINT_WHITE      0x00
#define TINT_HIGHLIGHT1 0x01	/* highest priority, dna highlighting */
#define TINT_HIGHLIGHT2 0x02	/* highlight friends */
#define TINT_RED        0x04 
#define TINT_LIGHTGRAY  0x08 
#define TINT_MAGENTA    0x10 
#define TINT_CYAN       0x20 
#define TINT_LIGHTGREEN 0x40 
#define TINT_YELLOW     0x80 

static int tints[8] = { LIGHTRED, MIDBLUE, RED, LIGHTGRAY, 
			MAGENTA, CYAN, LIGHTGREEN, YELLOW } ;

#define LeftBorder 70		/* x-pos of y-axis' 0-position */
#define TopBorder 65		/* y-pos of x-axis' 0-position */
#define DEFAULTALIGNLEN 125
#define MAXALIGNLEN 501

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

char *dotterVersion = "3.1",
     *Xoptions=0;


typedef struct
{
  const char *name ;
  int start, end ;
  char strand ;
  MSP *msp_start, *msp_end ;
  float y_pos ;
} GeneDataStruct, *GeneData ;

typedef struct
{
  GList *forward_genes ;
  GList *reverse_genes ;
  char strand  ;
} GeneStrandStruct, *GeneStrand ;



extern void colorMap (void) ;
static void setWindow (void) ;
static void initAlignment(void);
static void keyboard (int key, int modifier_unused);
static void savePlot(void);
static void togglePrintColors(void);
static void setAlnLen(void);
static void drawBlastHSPgray(void);
static void drawBlastHSPline(void);
static void drawBlastHSPlinef(void);
static void dotterRedraw(void);
static void togglePixelmap(void);
static void toggleCrosshair(void);
static void toggleCrosshairPos(void);
static void toggleCrosshairFullscreen(void);
static void toggleGrid(void);
static void clearHSPs(void);
static void initCrosshair(void);
static void callDotterParams(void);
static void readmtx(int mtx[24][24], char *mtxfile);
static void mtxcpy(int mtx[24][24], int BLOSUM62[24][24]);
static void DNAmatrix(int mtx[24][24]);
static void dotterPrint(void);
static void Help(void);
static void dotterDestroy(void) ;
static void loadPlot(char *loadfile) ;
static void initWindow(char *winsize) ;
static void calcWindow(void) ;
static void drawAllFeatures(MSP *msp) ;
static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth) ;
gint compareMSPs(gconstpointer a, gconstpointer b) ;
static void getGenePositionsCB(gpointer data, gpointer user_data) ;
gint compareGenes(gconstpointer a, gconstpointer b) ;
static void setYoffsets(GList *first, float min_offset, float max_offset) ;
static void drawGenesCB(gpointer data, gpointer user_data) ;
static void drawMSPGene(MSP *msp, float y_offset) ;
//static void printMSP(gpointer data, gpointer user_data) ;
//static void printGene(gpointer data, gpointer user_data) ;
static int gArrayGetLen(GArray *array);


#define toggleCrosshairStr    "Crosshair"
#define toggleCrosshairPosStr "Crosshair coordinates"
#define toggleCrosshairFullscreenStr "Crosshair over whole window"
#define toggleGridStr         "Grid"
#define drawBlastHSPgrayStr   "Draw Blast HSPs (gray pixels)"
#define drawBlastHSPlineStr   "Draw Blast HSPs (red lines)"
#define drawBlastHSPlinefStr  "Draw Blast HSPs (colour = f(score))"
#define togglePixelmapStr     "Pixelmap"


static void dotterRampChange(BOOL isDrag);
static MENU dotterMenu ;
static MENUOPT mainMenu[] = {  
  {graphDestroy,     "Quit"},
  {Help,             "Help"},
  {gexRampTool,    "Greyramp Tool"},
  {initAlignment,    "Alignment Tool"},
  {dotterPrint,      "Print"},
  {toggleCrosshair,   toggleCrosshairStr},
  {toggleCrosshairPos,toggleCrosshairPosStr},
  {toggleCrosshairFullscreen, toggleCrosshairFullscreenStr},
  {toggleGrid,        toggleGridStr},
  {selectFeatures,    selectFeaturesStr},
  {menuSpacer,       ""},
  {callDotterParams, "Zoom in with parameter control"},
  {savePlot,         "Save current plot"},
/*    {loadFeaturesPrompt,"Load features from file"},*/
  {setWindow,        "Change size of sliding window"},
  {menuSpacer,       ""},
  {drawBlastHSPgray, drawBlastHSPgrayStr},
  {drawBlastHSPline, drawBlastHSPlineStr},
  {drawBlastHSPlinef,drawBlastHSPlinefStr},
  {clearHSPs,        "Remove HSPs"},
  {togglePixelmap,   togglePixelmapStr},
  {0, 0}
} ;

static MENUOPT alnmenu[] =
{
  {graphDestroy,        "Quit"},
  {graphPrint,          "Print"},
  {togglePrintColors,   "Toggle colours for printing"},
  {setAlnLen,           "Set Alignment length"},
  {0, 0}
};

enum { BLASTNOTHING, BLASTRED, BLASTFUNC, BLASTGREY };

/* Global variables */

static int    MATRIX[24][24],
              qlen,	/* query residues (horizontal sequence) */
              qlen4,	/* query Pixels (pixelmap length) */
              slen,	/* subject residues (vertical sequence) */
              slen4,	/* subject Pixels (pixelmap height) */
              qoffset,  /* Difference between displayed q coord and position in qseq */
              soffset,	/* Difference between displayed s coord and position in sseq  */
              qseqbox, xqseqbox[3], sseqbox, qposbox, sposbox,
              qseqboxCrick, sseqboxCrick, qposboxCrick, sposboxCrick, 
              oldcolor, 
              backgColor=LIGHTGRAY, alnBackgColor=TINT_LIGHTGRAY, 
              vLineBox, hLineBox, CrosshairPosBox,
              win,		/* The length of the sliding window */
              abetsize,		/* The alphabet size (20 or 4) */
              blastp, blastn, blastx,
              zoom,		/* Zoomfactor = 1 ... */
              oldx, oldy,
              selfcomp,
              Display_selfcomp,
              watsonOnly,
              crickOnly,
              resfac,	/* Residue factor. 3 for DNA-Protein, 1 otherwise */
              reversedScale,
              plusmin,
              RightBorder,
              MSPoffset,	/* Difference between real MSP coord and coord stored in MSP */
              CrosshairON = 1,
              CrosshairPosON = 1,
              CrosshairFullscreenON = 1,
              PixelmapON,
              pixelmap_done,
              pixelFac,
              datalen, 
              ALIGNLEN = DEFAULTALIGNLEN,      /* use an odd number please */
              gridOn = 0,
              BlastHSPsOn = 0,
              BlastMode = BLASTNOTHING,
              printMode = 0,
              HSPgaps = 0,
              fsBoxStart,
              fsRightOn = 1,
              fsBottomOn = 1,
              fsAnnRightOn = 1,
              fsAnnBottomOn = 1,
              fsEndLinesOn = 0,
              alignmentInitialized = 0,
              greyRampSwap = 0;

       Graph  dotterGraph=0;
static Graph  alnGraph=0, fsGraph=0;
static UCHAR *data, *HSPpixels=0 ;
static char  *qseq, *sseq, *qname, *sname,
  qpos[10], spos[10],
  qseqDisp[MAXALIGNLEN+1],
  xqseqDisp[3][MAXALIGNLEN+1],
  sseqDisp[MAXALIGNLEN+1],
  qseqDispCrick[MAXALIGNLEN+1],
  sseqDispCrick[MAXALIGNLEN+1],
  qcolors[MAXALIGNLEN+1], xqcolors[3][MAXALIGNLEN+1], scolors[MAXALIGNLEN+1],
  qcolorsCrick[MAXALIGNLEN+1], scolorsCrick[MAXALIGNLEN+1],
  *qrevcompl=0, *pepqseq[3],
  CrosshairPosText[100],
  MATRIX_NAME[81] = "",
  fsPlotHeighttx[10],
  *banner;

char         *dotterBinary=0;

static double rectx, oldrectx, recty, oldrecty, exp_res_score, Lambda;
static float  crossx, crossy, oldLinew;
       float fsPlotHeight=2.0;	/* The height of feature series XY plots */
static int fontw, fonth;		/* fontwidth, fontheight (in pixels) */
static FILE  *saveFil;
static MSP   *MSPlist=0;     /* List of MSPs - the first object contains data */
static MSP   *msp;
static STORE_HANDLE  handle;


int atob_0[]	/* NEW (starting at 0) ASCII-to-binary translation table */
	= {
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,23,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
	14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,
	NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
	14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,

	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
	NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR 
};

int atob[]	/* OLD (starting at 1) ASCII-to-binary translation table  (Inherited from blast) */
	= {
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,24,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,

	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};
char aa_btoa[]	/* binary-to-ASCII translation table */
	= "-ARNDCQEGHILKMFPSTWYVBZX*" ;

/*  BLOSUM62 930809

     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
int BLOSUM62[24][24] = {
  {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4 },
  { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4 },
  { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4 },
  { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
  {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 },
  { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4 },
  { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4 },
  { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4 },
  { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4 },
  { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4 },
  { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4 },
  { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4 },
  { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4 },
  { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4 },
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4 },
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4 },
  { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4 },
  { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4 },
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4 },
  { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
  { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
  { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4 },
  { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1 }
 };



double aafq[20]	/* Amino acid residue frequencies used by S. Altschul */
	= {.081, .057, .045, .054, .015, .039, .061, .068, .022, .057,
	   .093, .056, .025, .040, .049, .068, .058, .013, .032, .067 } ;

#define NN 5

int ntob[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

int ntob_compl[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN, 3,NN, 2,NN,NN,NN, 1,NN,NN,NN,NN,NN,NN, 4,NN,
        NN,NN,NN,NN, 0,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

/* binary-to-ASCII translation table */
char bton[] = "ACGTN*";

GArray *fsArr = NULL;  /* Stores Feature Series - the actual segments are stored
			   as MSPs, using these fields:
			   msp->sframe  = [1..2] sequence
			   msp->qstart = segment start
			   msp->qend   = segment end
			   msp->fs     = Ordinal number of series that this MSP belongs to.
			   msp->fsColor  = color
			   msp->desc   = annotation
			   */





/* 
 *               External routines.
 */



void dotter (char  type,
	     char *opts,
	     const char *queryname,
	     char *queryseq,
	     int   qoff,
	     const char *subjectname,
	      char *subjectseq,
	      int   soff,
	      int   qcenter,
	      int   scenter,
	      char *savefile,
	      char *loadfile,
	      char *mtxfile,
	      float memoryLimit,
	      int   zoomFac,
	      MSP  *MSPs,
	      int   MSPoff,
	      char *winsize,
	      int   pixelFacset)
{
  int  i;

//  messalloccheck();

  /* Reset global statics */
  resfac = PixelmapON = 1;
  blastp = blastn = blastx = selfcomp =
    Display_selfcomp = watsonOnly = crickOnly = reversedScale = pixelmap_done = 0;
  BlastHSPsOn = 0;
  saveFil = 0;

  switch(type) {
  case 'P':  
    blastp = 1; 
    abetsize = 20;  break;
  case 'N':  
    blastn = 1; 
    abetsize = 4;   break;
  case 'X':  
    blastx = 1; 
    resfac = 3; 
    abetsize = 20;  break;
  default: g_error("Invalid sequence type passed to Dotter: %c", type);
  }
  
  /* Option parsing */
  if (opts)
    while (*opts)
      {
	switch (*opts)
	  {
	  case 'D': Display_selfcomp = 1; break;
	  case 'R': reversedScale = 1;    break;
	  case 'H': 
	    BlastHSPsOn = 1;         
	    PixelmapON = 0;             break;
	  case 'W': watsonOnly = 1;       break;
	  case 'C': crickOnly = 1;        break;
	  case 'G': HSPgaps = 1;          break;
	  case 'S': greyRampSwap = 1;     break;
	  case 'L': fsEndLinesOn = 1;     break;
	  }
	opts++;
      }

  plusmin = ( reversedScale ? -1 : 1 );
    
  if (graphActivate(dotterGraph))
    {
      dotterDestroy();

      /* Don't free data or HSPpixels, since every time rawImage 
	 is called with new data, a new XImage struct is added to gSubDev->images.
	 These are later free'd by XDestroyImage in graphSubDevDestroy.

	 However, memory allocated this way never seems to be entirely free'd
	 by graphDestroy, in the sense that it's not always reused. 
	 Maybe an illusion? */
	
      /* free(data);
	 if (HSPpixels) free(HSPpixels);*/
    }

  HSPpixels = 0;
  handle = handleCreate();
  banner = (char *)handleAlloc(0, handle, 1000);

  qname = g_malloc(strlen(queryname)+1); strcpy(qname, queryname);
  qseq = queryseq;
  qoffset = qoff;

  sname = g_malloc(strlen(subjectname)+1); strcpy(sname, subjectname);
  sseq = subjectseq;
  soffset = soff;
  MSPlist = MSPs;
  MSPoffset = MSPoff;
  
  if (!(qlen = strlen(qseq))) g_error("queryseq is empty");
  if (!(slen = strlen(sseq))) g_error("subjectseq is empty");

  for (i = 0; i < qlen; i++) qseq[i] = freeupper(qseq[i]);
  for (i = 0; i < slen; i++) sseq[i] = freeupper(sseq[i]);

  if (!memoryLimit) memoryLimit = 0.5; /* Mb */

  if (!strcmp(qseq, sseq)) selfcomp = 1;

  if (blastn && !watsonOnly) {
    /* Reverse and complement qseq for drawAlignment */
    qrevcompl = handleAlloc(0, handle, qlen+1);
    revComplement(qrevcompl, qseq) ;
    /*for (i = 0; i < qlen; i++) qrevcompl[i] = bton[ntob_compl[qseq[qlen-i-1]]];
      qrevcompl[qlen] = 0;*/
  }
    
  if (blastx) {
    for (i = 0; i < 3; i++)
      pepqseq[i] = blxTranslate(qseq+i, stdcode1);
  }


  /* Get score matrix */
  if (mtxfile) {
    readmtx(MATRIX, mtxfile);
    strncpy(MATRIX_NAME, mtxfile, 80);
  }
  else {
    if (blastn) {
      DNAmatrix(MATRIX);
      strcpy(MATRIX_NAME, "DNA+5/-4");
    }
    else {
      mtxcpy(MATRIX, BLOSUM62);
      strcpy(MATRIX_NAME, "BLOSUM62");
    }
  }


  /* Don't do batch processing if output file can't be opened */
  if (savefile) 
    {
      if (!(saveFil = fopen (savefile, "wb")))
	g_error("Failed to open %s", savefile);
    }
	

  if (loadfile)
    {
      /* Load existing dotplot file */
      loadPlot(loadfile);
    }
  else 
    {
      initWindow(winsize);

      /* Set pixelFac so that exp_res_score is at 1/5 of the range. 
       * This positions exp_res_score at 51.2
       */
      if (pixelFacset)
	pixelFac = pixelFacset;
      else
	pixelFac = 0.2*256/exp_res_score;
	
      if (!zoomFac)
	zoom = (int)sqrt((qlen/resfac/1e6*slen - 1e-6)/memoryLimit) + 1;
      else
	zoom = zoomFac;

      qlen4 = (int)ceil((double)qlen/resfac/zoom);
      if (qlen4 % 4)
	qlen4 += 4-(qlen4 % 4);
      if (qlen/resfac > qlen4*zoom)
	g_error("qlen/resfac > qlen4*zoom (%d > %d (%d*%d))",
	      qlen/resfac, qlen4*zoom, qlen4, zoom);

      slen4 = (int)ceil((double)slen/zoom);
      if (slen4 % 4)
	slen4 += 4-(slen4 % 4);
      if (slen > slen4*zoom)
	g_error("slen > slen4*zoom (%d > %d (%d*%d))", slen, slen4, zoom, slen4*zoom);

      datalen = slen4*qlen4;
      data = (UCHAR *)g_malloc(datalen);
    }

  if (savefile)
    {
      calcWindow();
      savePlot();
    }
  else
    {
      dotterMenu = menuInitialise ("dotter", (MENUSPEC*)mainMenu) ;
      if (!MSPlist)
	menuSuppress (dotterMenu, drawBlastHSPgrayStr);
      menuSetFlags(menuItem(dotterMenu, selectFeaturesStr), MENUFLAG_HIDE);

      if (qcenter < 1 || qcenter > qlen)
	{
	  oldx = qlen/resfac/2;
	  crossx = oldx/zoom + LeftBorder;
	}
      else
	{
	  crossx = LeftBorder + qcenter/resfac/zoom;
	  oldx = qcenter;
	}

      if (scenter < 1 || scenter > slen)
	{
	  oldy = slen/2;
	  crossy = oldy/zoom + TopBorder;
	}
      else
	{
	  crossy = TopBorder + scenter/zoom;
	  oldy = scenter;
	}
	
      RightBorder = LeftBorder-1 + qlen4 - qlen % 4;
	
      if (fsArr) /* feature-series array exists */
	{
	  CrosshairFullscreenON = 1;
	}
	
      if (!loadfile && !BlastHSPsOn)
	{
	  calcWindow();
	}

      dotterRedraw();	/* Note: must be after calcWindow; used to be before in AW to fail early */
    }


  return ;
}



char Seqtype(char *seq)
{
    char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
    char *primenuc    = "ACGTUN";
    char *protonly    = "EFIPQZ";

    /* Simplified version of Sean Eddy's */

    int  pos;
    char c;
    int  po = 0;			/* count of protein-only */
    int  nt = 0;			/* count of t's */
    int  nu = 0;			/* count of u's */
    int  na = 0;			/* count of nucleotides */
    int  aa = 0;			/* count of amino acids */
    int  no = 0;			/* count of others */
  
    /* Look at the first 300 characters
     */
    for (pos = 0; seq[pos] && pos < 300; pos++)
    {
	c = freeupper(seq[pos]);

	if (strchr(protonly, c)) 
	    po++;
	else if (strchr(primenuc, c)) {
	    na++;
	    if (c == 'T') nt++;
	    else if (c == 'U') nu++;
	}
	else if (strchr(aminos, c)) 
	    aa++;
	else if (isalpha(c)) 
	    no++;
    }

    if (po > 0) return 'P';
    else if (na > aa) return 'N';
    else return 'P';
} /* Seqtype */




/* 
 *                Internal routines.
 */


/************************/

/* RMEXP will subtract the expected score to remove background noise
   due to biased composition. Gos's idea - not explored yet.

static void rmExp(void){}
*/


static void menuCheck(MENU menu, int mode, int thismode, char *str)
{
    if (mode == thismode)
	menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
    else
	menuUnsetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);

    menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE);

}
static void setMenuCheckmarks(void)
{
    menuCheck(dotterMenu, 1, CrosshairON, toggleCrosshairStr);
    menuCheck(dotterMenu, 1, CrosshairPosON, toggleCrosshairPosStr);
    menuCheck(dotterMenu, 1, CrosshairFullscreenON, toggleCrosshairFullscreenStr);
    menuCheck(dotterMenu, 1, gridOn, toggleGridStr);

    menuCheck(dotterMenu, BlastMode, BLASTGREY, drawBlastHSPgrayStr);
    menuCheck(dotterMenu, BlastMode, BLASTRED, drawBlastHSPlineStr);
    menuCheck(dotterMenu, BlastMode, BLASTFUNC, drawBlastHSPlinefStr);

    menuCheck(dotterMenu, 1, PixelmapON, togglePixelmapStr);

    graphNewMenu(dotterMenu);
}


static void Help(void)
{
    graphMessage (messprintf("\
Dotter, version %s\n\
Copyright (C) Erik Sonnhammer\n\
and Richard Durbin, 1995\n\
\n\n\
Sliding window length = %d\n\
Pixel values = %d x score/residue\n\
Matrix = %s\n\
Zoom (compression) factor = %d\n\
\n\
LEFT MOUSE BUTTON:\n\
  Position crosshair.\n\
\n\
MIDDLE MOUSE BUTTON:\n\
  Zoom in region.\n\
\n\
RIGHT MOUSE BUTTON:\n\
  Menu.\n\
\n\n\
KEYSTROKES:\n\
  Arrow keys: move up/down/left/right\n\
  < > : move along diagonals\n\
  { } : move along reverse diagonals\n\
\n\n\
RESIDUE COLOURS in alignment tool:\n\
  Cyan      = Identical Residue.\n\
  DarkBlue  = Positive Score.\n\
  No colour = Negative score.\n", dotterVersion, win, pixelFac, MATRIX_NAME, zoom));
}


/* REVERSEBYTES changes order of n bytes at location ptr
*/
void reversebytes(void *ptr, int n)
{ 
  static char copy[256], *cp;
  int  i;

  cp = ptr;
  memcpy(copy, ptr, n);  /* Note: strcpy doesn't work - stops at \0 */
  
  for(i=0; i<n; i++) *cp++ = copy[n-i-1];
}


/* SAVEPLOT writes a 1 byte fileformat first, various params and then the dot-matrix
*/
static void savePlot(void)
{
    static char 
	dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE] ;
    int 
	i, j,
	batch=1,
	MNlen, MNlenSave,
	mtx;
    UCHAR 
	format = 2;


    if (!saveFil) {
	if (!(saveFil = filqueryopen(dirName, fileName, "dotter","w", 
				     "Save as file: (adds .dotter to filename)")))
	    return ;

	batch = 0;
    }

    MNlen = MNlenSave = strlen(MATRIX_NAME);

#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
    reversebytes(&pixelFac, 4);
    reversebytes(&win, 4);
    reversebytes(&MNlenSave, 4);
#endif

    fwrite(&format,    1, 1, saveFil);
    fwrite(&zoom,      1, 4, saveFil);
    fwrite(&qlen4,     1, 4, saveFil);
    fwrite(&slen4,     1, 4, saveFil);
    fwrite(&pixelFac,  1, 4, saveFil); /* New feature of format 2  */
    fwrite(&win,       1, 4, saveFil); /* New feature of format 2  */
    fwrite(&MNlenSave, 1, 4, saveFil); /* New feature of format 2  */
    fwrite(&MATRIX_NAME, 1, MNlen, saveFil); /* New feature of format 2  */
    for (i = 0; i < 24; i++)
	for (j = 0; j < 24; j++) {
	    mtx = MATRIX[i][j];
#ifdef ALPHA
	    reversebytes(&mtx, 4);
#endif
	    fwrite(&mtx, 1, 4, saveFil); /* New feature of format 2  */
	}

#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
    reversebytes(&pixelFac, 4);
    reversebytes(&win, 4);
#endif
    
    fwrite(data, 1, qlen4*slen4, saveFil);

    fclose(saveFil);
    saveFil = 0;
}


static void loadPlot(char *loadfile)
{
    FILE 
	*fil;
    int 
	i, j, n, 
	dotstart = 0,
	MNlen,
	mtx;
    UCHAR 
	format;

    if (!(fil = fopen (loadfile, "rb")))
	g_error("Failed to open %s", loadfile);

    if ((fread(&format, 1, 1, fil)) != 1) g_error("reading file %s", loadfile);
    if ((fread(&zoom,   1, 4, fil)) != 4) g_error("reading file %s", loadfile);
    if ((fread(&qlen4,  1, 4, fil)) != 4) g_error("reading file %s", loadfile);
    if ((fread(&slen4,  1, 4, fil)) != 4) g_error("reading file %s", loadfile);
#ifdef ALPHA
    reversebytes(&zoom, 4);
    reversebytes(&qlen4, 4);
    reversebytes(&slen4, 4);
#endif

    if (format == 1) {
	dotstart = 13;
	/* Don't actually know these variables for sure - guess the most common */
	pixelFac = 50;
	win = 25;
    }
    else if (format == 2) 
    {
	if ((fread(&pixelFac, 1, 4, fil)) != 4) g_error("reading file %s", loadfile);
	if ((fread(&win,      1, 4, fil)) != 4) g_error("reading file %s", loadfile);
	if ((fread(&MNlen, 1, 4, fil)) != 4) g_error("reading file %s", loadfile);
#ifdef ALPHA
	reversebytes(&pixelFac, 4);
	reversebytes(&win, 4);
	reversebytes(&MNlen, 4);
#endif
	if ((fread(&MATRIX_NAME, 1, MNlen, fil)) != MNlen) g_error("reading file %s", loadfile);
	MATRIX_NAME[MNlen] = 0;
	for (i = 0; i < 24; i++)
	    for (j = 0; j < 24; j++) {
		if ((fread(&mtx, 1, 4, fil)) != 4) g_error("reading file %s", loadfile);
#ifdef ALPHA
		reversebytes(&mtx, 4);
#endif
		MATRIX[i][j] = mtx;
	    }
	dotstart = MNlen + 2329;
    }
    else 
	g_error("Unknown dotter file format version: %d", format);

    fseek(fil, 0, SEEK_END);
    n = ftell(fil);

    if (n-dotstart != qlen4*slen4)
	g_error("Wrong number of pixels in %s: %d. Expected %d * %-d = %d\n", 
	      loadfile, n, qlen4, slen4, qlen4*slen4);

    datalen = slen4*qlen4;
    data = (UCHAR *)g_malloc(datalen);
    fseek(fil, dotstart, SEEK_SET);

    if ((n = fread(data, 1, qlen4*slen4, fil)) != qlen4*slen4)
	g_error("Read wrong number of pixels from %s: %d. Expected %d * %-d = %d\n", 
	      loadfile, n, qlen4, slen4, qlen4*slen4);
    
    fclose(fil);

    printf("I read your dotter dotplot file %s.\n", loadfile);
    if (format == 1) 
	printf("It was in the old format 1, so the windowsize and pixel factor had to be guessed.\n");
    fflush(stdout);
}

/*
void loadFeaturesDotter(FILE* fil)
{
    loadFeatures(fil, &MSPlist);

    dotterRedraw();
    if (fsGraph) selectFeatures();
}


static void loadFeaturesPrompt(void)
{
    FILE *fil;
    static char dirName[DIR_BUFFER_SIZE], fileName[FIL_BUFFER_SIZE];

    if (!(fil = filqueryopen(dirName, fileName, "","r", "Load file: "))) return;

    loadFeaturesDotter(fil);
}
*/


/* INITMATRIX calculates the min, max and mean score of the given matrix.
   Archaic and not needed any more.
* /
static void initMatrix (void)
{
    int i, j, x, sum ;

    minScore = maxScore = sum = 0 ;
    for (i = 0 ; i < abetsize ; i++)
	for (j = 0 ; j < abetsize ; j++)
	    { 
		x = MATRIX[i][j];
		if (x < minScore) minScore = x;
		if (x > maxScore) maxScore = x;
		sum += x;
	    }
	
    meanScore = (float)sum/(abetsize*abetsize);

    if (blastn)
	printf("DNA");
    else
	printf("Protein");

    printf(" score matrix: Max score= %d, Min score= %d, Mean score= %.2f\n", 
	   maxScore, minScore, meanScore);
}
*/


/* CALCWINDOW puts the max diagonal for each pixel into *data
 ************************************************************/
static void calcWindow(void)
{
    STORE_HANDLE calcHandle;

    register int 
	q, s, qmax,     /* Loop variables */
	*newsum,	/* New sum pointer */
	*oldsum,	/* Old sum pointer */
	*delrow,	/* Pointer to scoreVec, row to remove */
	*addrow,	/* Pointer to scoreVec, row to add */
	dotpos, dotposq, dotposs;
    
    int 
	**scoreVec = NULL,     /* Array of precalculated scores for qseq residues */
        *sIndex,	/* Array of binary coded sseq residues */
	ql,	        /* query position in local submatrix (of one pixel) */
	sl,	        /* subject position in local submatrix (of one pixel) */
	i, min, sec, frame, pepqseqlen, 
	*zero = NULL, *sum1 = NULL, *sum2 = NULL,
	win2 = win/2, 
	val;

    float 
	dots,	        /* total number of dots (millions) */
	speed;          /* speed in Mdots/seconds */

    UCHAR 
	*dot, dotValue;

    calcHandle = handleCreate();
    
    speed = 17.2;  /* SGI MIPS R10000 (clobber) */
    /* speed = 5.7;  DEC Alpha AXP 3000/700 */
    /* speed = 3.7;  SGI R4400: */
    dots = qlen/1e6*slen;
    if (selfcomp) dots /= 2;
    if (blastn && !(watsonOnly || crickOnly)) dots *= 2;
    if (blastx) dots *= 3;

    min = (int)(dots/speed/60);
    sec = (int)(dots/speed) - min*60;
    printf("%d vs. %d residues => %.2f million dots. ",
	   qlen, slen, dots);

    if (min+sec >= 2) {
	printf("(Takes ");

	if (min)
	    printf("%d:%.2d minutes", min, sec);
	else 
	    printf("%d seconds", sec);
    
	printf(" on an SGI MIPS R10000)");
    }

    printf("\n");
    fflush(stdout);

    /* Initialize lookup tables for faster execution */
    sIndex = (int *) handleAlloc(0, calcHandle, slen * sizeof(int)) ;

    if (blastp)
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 25*sizeof(int*));
	for (i = 0; i < 25; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen*sizeof(int));
	
	for (i = 0; i < 24; i++)
	    for (q = 0; q < qlen; q++)
		scoreVec[i][q] = MATRIX[i][atob[(int)qseq[q]]-1];
	
	/* Non-protein symbols in scorevector */
	for (q = 0; q < qlen; q++) scoreVec[24][q] = MATRIX[23][23];
	
	for (s = 0; s < slen; s++) sIndex[s] = atob[(int)sseq[s]]-1;
    }
    else if (blastn) 
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen * sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 6*sizeof(int*));
	
	for (i = 0; i < 6; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen*sizeof(int));
	/* Fill the scoreVec below depending on which strand */

	    
	for (s = 0 ; s < slen ; s++) sIndex[s] = ntob[(int)sseq[s]];
    }
    
    /* Reset background */
    for (i = 0; i < slen4*qlen4 ; i++) data[i] = 0;


    if (blastx)
    {
	zero = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	sum1 = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	sum2 = (int *)handleAlloc(0, calcHandle, qlen/3*sizeof(int));

	scoreVec = (int**) handleAlloc(0, calcHandle, 25*sizeof(int *));
	for (i = 0; i < 25; i++)
	    scoreVec[i] = (int*) handleAlloc(0, calcHandle, qlen/3*sizeof(int));
	
	/* Non-protein symbols in scorevector */
	for (q = 0; q < qlen/3; q++) scoreVec[24][q] = MATRIX[23][23];
	
	for (s = 0; s < slen; s++) sIndex[s] = atob[(int)sseq[s]]-1;

	for (frame = 0; frame < 3; frame++)
	{
	    pepqseqlen = strlen(pepqseq[frame]);
	    for (i = 0; i < 24; i++)
		for (q = 0; q < pepqseqlen; q++)
		    scoreVec[i][q] = MATRIX[i][atob[(int)pepqseq[frame][q]]-1];

	    for (i = 0; i<pepqseqlen; i++) {
		sum1[i] = 0;
		sum2[i] = 0;
	    }

	    for (s = 0; s < slen; ++s)
	    {   
		if (s & 1) {
		    newsum = sum1 ;
		    oldsum = sum2 ;
		}
		else {
		    newsum = sum2 ;
		    oldsum = sum1 ;
		}

		if (s >= win) delrow = scoreVec[sIndex[s-win]];
		else delrow = zero;

		addrow = scoreVec[sIndex[s]];
		*newsum = *addrow++;
		
		qmax = min(win, pepqseqlen);
		for (q = 1; q < qmax ; ++q)
		    *++newsum = *oldsum++ + *addrow++;
		
		qmax = pepqseqlen;
		for ( ; q < qmax ; ++q) {
		    *++newsum = *oldsum++ + *addrow++ - *delrow++ ;
		    if (*newsum > 0 && s >= win) 
		    {
			dotposq = (q-win2)/zoom;
			dotposs = (s-win2)/zoom;
			    
			/* Only fill half the submatrix */
			ql = q-win2 - dotposq*zoom;
			sl = s-win2 - dotposs*zoom;
			if (sl >= ql)
			{
			    dotpos = qlen4*dotposs + dotposq;

			    if (dotpos < 0 || dotpos >= datalen) {
				g_critical ( "Pixel out of bounds (%d) in blastx: %d\n",
					datalen-1, dotpos);
			    }
			    else {
				val = *newsum * pixelFac / win;
				dotValue = (val > 255 ? 255 : (UCHAR)val);
				dot = &data[dotpos];
				if (dotValue > *dot) *dot = dotValue;
			    }
			}
		    }
		}
	    }
	}
    }

    if (blastp || (blastn && !crickOnly)) 
    {
	if (blastn)
	    for (i = 0; i < 6; i++)
		for (q = 0; q < qlen; q++)
		    scoreVec[i][q] = MATRIX[i][ntob[(int)qseq[q]]];

	for (s = 0; s < slen; ++s)
	{ 
	    if (s & 1) {
		newsum = sum1 ;
		oldsum = sum2 ;
	    }
	    else {
		newsum = sum2 ;
		oldsum = sum1 ;
	    }

	    if (s >= win) delrow = scoreVec[sIndex[s-win]];
	    else delrow = zero;

	    addrow = scoreVec[sIndex[s]];
	    *newsum = *addrow++;

	    qmax = min(win, qlen);
	    for (q = 1; q < qmax ; ++q)
		*++newsum = *oldsum++ + *addrow++;

	    qmax = (selfcomp ? s+1 : qlen);
	    for ( ; q < qmax ; ++q) 
	    {
		*++newsum = *oldsum++ + *addrow++ - *delrow++;
		if (*newsum > 0 && s >= win) 
		{
		    dotposq = (q-win2)/zoom;
		    dotposs = (s-win2)/zoom;
		    
		    /* Only fill half the submatrix */
		    ql = q-win2 - dotposq*zoom;
		    sl = s-win2 - dotposs*zoom;
		    if (sl >= ql)
		    {
			dotpos = qlen4*dotposs + dotposq;

			if (dotpos < 0 || dotpos > datalen-1) {
			    g_critical ( "Pixel out of bounds (%d) in blastp/blastn-forw: %d\n", 
				    datalen-1, dotpos);
			}
			else {
			    /* Keep the max dot value of all diagonals in this pixel */
			    val = *newsum * pixelFac / win;
			    dotValue = (val > 255 ? 255 : (UCHAR)val);
			    dot = &data[dotpos];
			    if (dotValue > *dot) *dot = dotValue;
			}
		    }
		}
	    }
	}
    }


    if (blastn && !watsonOnly) 
    {
	if (blastn)
	    for (i = 0; i < 6; i++)
		for (q = 0; q < qlen; q++)
		    scoreVec[i][q] = MATRIX[i][ntob_compl[(int)qseq[q]]];

	for (i = 0; i<qlen; i++) {
	    sum1[i] = 0;
	    sum2[i] = 0;
	}

	for (s = slen-1; s >= 0; --s)
	{ 
	    if (s & 1) {
		newsum = sum1 ;
		oldsum = sum2 ;
	    }
	    else {
		newsum = sum2 ;
		oldsum = sum1 ;
	    }

	    if (s < slen-win) delrow = scoreVec[sIndex[s+win]];
	    else delrow = zero;

	    addrow = scoreVec[sIndex[s]];
	    *newsum = *addrow++;

	    qmax = min(win, qlen);
	    for (q = 1; q < qmax ; ++q)
		*++newsum = *oldsum++ + *addrow++;

	    qmax = (selfcomp ? s+1 : qlen);
	    for ( ; q < qmax ; ++q) {
		*++newsum = *oldsum++ + *addrow++ - *delrow++ ;
		if (*newsum > 0 && s <= slen-win) 
		{
		    dotposq = (q-win2)/zoom;

		    dotposs = (s+win2)/zoom;

		    /* Only fill half the submatrix */
		    ql = q-win2 - dotposq*zoom;

		    /* Set the origin (0,0) to the bottom left corner of submatrix
		       Ugly but correct. Zoom = pixels/submatrix */
		    sl = zoom-1 - (s+win2 - dotposs*zoom);
		    
		    if (sl >= ql)
		    {
			dotpos = qlen4*dotposs + dotposq;
			
			if (dotpos < 0 || dotpos >= datalen) {
			    g_critical ( "Pixel out of bounds (%d) in blastn-rev: %d\n",
				    datalen-1, dotpos);
			}
			else {
			    val = *newsum * pixelFac / win;
			    dotValue = (val > 255 ? 255 : (UCHAR)val);
			    dot = &data[dotpos];
			    if (dotValue > *dot) *dot = dotValue;
			}
		    }
		}
	    }
	}
    }

    if (selfcomp && Display_selfcomp) {
	/* Copy mirror image */

	int dotposCopy;
	
	for (s = 0; s < slen4; ++s) { 
	    for (q = 0; q < s ; ++q) {
		
		dotpos = qlen4*s + q;
		dotposCopy = qlen4*q + s;

		if (dotpos < 0 || dotpos >= datalen)
		    g_critical ( "Source pixel out of bounds (%d) in mirrorCopy: %d\n",
				datalen-1, dotpos);
		if (dotposCopy < 0 || dotposCopy >= datalen)
		    g_critical ( "Destination pixel out of bounds (%d) in mirrorCopy: %d\n",
				datalen-1, dotposCopy);
		data[dotposCopy] = data[dotpos];
	    }
	}
    }

    handleDestroy(calcHandle);
//    messalloccheck();
    pixelmap_done = 1;
}


void findScaleUnit (float cutoff, float *u, float *sub)
{
  float unit = *u ;
  float subunit = *u ;

  if (cutoff < 0)
    cutoff = -cutoff ;

  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5000001 ;	/* safe rounding */
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }
  subunit /= 10 ;
  if (subunit > *sub)
    *sub = subunit ;
  *u = unit ;
}


static void drawScale(void)
{
  int i, x, y,
    TickStart;			/* First minor tick   */
  float unit=1,		/* Major tick unit */
    subunit=1;		/* First Major tick   */
  
  graphRectangle(LeftBorder-1, TopBorder-1, LeftBorder-1+qlen4, TopBorder-1+slen4);
  
  /* QSEQ */
  /* min 100 pixels between major ticks */
  findScaleUnit(100*zoom*resfac, &unit, &subunit);
  y = TopBorder -1;
  
  graphTextFormat(FIXED_WIDTH);
  
  if (!qoffset) 
    TickStart = 0;
  else 
    TickStart = subunit*ceil(qoffset/subunit);
  
  for (i = TickStart; i < qlen+qoffset; i += subunit)
    {
      if (reversedScale) 
	x = RightBorder - (i-qoffset)/zoom/resfac;
      else 
	x = LeftBorder-1 + ceil((float)(i-qoffset)/zoom/resfac);
      
      graphLine(x, y, x, y-5); /* Minor tick */
      
      /* Major tick */
      if (!(i % (int)unit))
	{
	  graphLine(x, y, x, y-10);
	  graphText(messprintf("%d", i), x-4, y-25);
	}
      
      if (gridOn && x != LeftBorder-1)
	{
	  graphColor(LIGHTRED);
	  graphLine(x, y+1, x, y+slen4-1);
	  graphColor(BLACK);
	}
    }
  
    /* SSEQ */
    /* min 50 pixels between major ticks */
    unit = subunit = 1;
    findScaleUnit(50*zoom, &unit, &subunit);
    x = LeftBorder - 1;

    if (!soffset) 
	TickStart = 0;
    else 
	TickStart = subunit*ceil(soffset/subunit);

    for (i = TickStart; i < slen+soffset; i += subunit) 
    {
	y = TopBorder-1 + (i-soffset)/zoom;

	graphLine(x, y, x-5, y); /* Minor tick */

	/* Major tick */
	if (!(i % (int)unit))
	  {
	    graphLine(x, y, x-10, y);
	    graphText(messprintf("%7d", i), 0, y-6);
	  }

	if (gridOn && y != TopBorder-1)
	  {
	    graphColor(LIGHTRED);
	    graphLine(x+1, y, x+qlen4-1, y);
	    graphColor(BLACK);
	  }
    }

    graphTextFormat(PLAIN_FORMAT);
}


/* Utility to return the Feature Series that the given MSP belongs to */
static FeatureSeries* mspGetFeatureSeries(const MSP *msp)
{
  return msp->fs;
}


/* Returns true if this MSP is part of a Feature Series and that feature series is currently displayed */
static gboolean mspShowFs(const MSP *msp)
{
  FeatureSeries *fs = mspGetFeatureSeries(msp);
  return (fs && fs->on);
}


/* Return the y position of the lower edge of a feature-series MSP, given its height.
   Calculate it if it is not set yet. Updates the max y coord to be the lower edge of
   this MSP. Returns 0 if the MSP is not in a Feature Series */
float mspGetFsBottomEdge(MSP *msp, float *maxy, float height)
{
  float result = 0;
  
  FeatureSeries *fs = mspGetFeatureSeries(msp);
  
  if (fs)
    {
      result = fs->y;

      if (!result)
        {
          *maxy += height;
          fs->y = *maxy;
          result = *maxy;
        }
    }
  
  return result;
}


/* Return the x position of the rightmost edge of a feature-series MSP, given its height.
   Calculae it if it is not yet set. Updates the max x coord to be the rightmost edge of
   this MSP. Returns 0 if the MSP is not in a FeatureSeries */
float mspGetFsRightEdge(MSP *msp, float *maxx, float height)
{
  float result = 0;

  FeatureSeries *fs = mspGetFeatureSeries(msp);

  if (fs)
    {
      result = fs->x;
      
      if (!result)
        {
          *maxx += height;
          fs->x = *maxx;
          result = *maxx;
        }
    }
  
  return result;
}


static float seq2graphX(int pos)
{
    float x;

    x = ceil((float)(pos + MSPoffset - qoffset)/zoom/resfac);

    if (reversedScale)
	x = RightBorder - x;
    else
	x += LeftBorder-1;

    return x;
}

static float seq2graphY(int pos)
{
    float y;

    y = ceil((float)(pos - soffset)/zoom);
    
    y += TopBorder-1;

    return y;
}


static void XdrawSEG(MSP *msp, float offset)
{
    /* Horizontal axis */

    float  
	sx = seq2graphX(mspGetQStart(msp)), 
	ex = seq2graphX(mspGetQEnd(msp));
    
    offset += TopBorder + slen4 -1;
    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);

    if (fsEndLinesOn) {
	graphLine(sx, TopBorder, sx, TopBorder+slen4-2);
	graphLine(ex, TopBorder, ex, TopBorder+slen4-2);
    }

    graphFillRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);
    graphColor(BLACK);
    graphRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);

    graphColor(BLACK);
    if (fsAnnBottomOn && msp->desc) {
	graphText(msp->desc, sx, offset);
    }    
    graphColor(oldcolor); graphLinewidth(oldLinew);
}


static void XdrawSEGxy(MSP *msp, float offset)
{
    int i, inNotFilled=0, descShown=0;
    float  
	x, y, 
	xold=0, yold=0;

    offset += TopBorder + slen4 -1;
    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);

    for (i = 0; i < qlen; i++)
      {
	const int xyVal = g_array_index(msp->xy, int, i);

	if (xyVal == XY_NOT_FILLED)
	  {
	    inNotFilled = 1;
	  }
	else
	  {
	    x = seq2graphX(i);
	    y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
	    
	    if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE))
	      {
	        graphLine(xold, yold, x, y);
		
		if (fsAnnBottomOn && msp->desc && !descShown)
		  {
		    int linecolor = graphColor(BLACK);
		    graphText(msp->desc, xold, offset);
		    graphColor(linecolor);
		    descShown = 1;
		  }
	      }
		
	    xold = x;
	    yold = y;
	    inNotFilled = 0;
	  }
      }
    
    graphColor(oldcolor); graphLinewidth(oldLinew);
}


static void YdrawSEG(MSP *msp, float offset)
{
    /* Vertical axis */

    float  
	sx = seq2graphY(mspGetQStart(msp)), 
	ex = seq2graphY(mspGetQEnd(msp));
    
    offset += LeftBorder + qlen4 -1;
    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);

    if (fsEndLinesOn) {
	graphLine(LeftBorder, sx, LeftBorder+qlen4-2, sx);
	graphLine(LeftBorder, ex, LeftBorder+qlen4-2, ex);
    }

    graphFillRectangle(offset, sx, offset - msp->score/100.0 * fonth, ex);
    graphColor(BLACK);
    graphRectangle    (offset, sx, offset - msp->score/100.0 * fonth, ex);

    graphColor(BLACK);
    if (fsAnnRightOn && msp->desc) {
	graphText(msp->desc, offset, sx);
    }    
    graphColor(oldcolor); graphLinewidth(oldLinew);
}


static void YdrawSEGxy(MSP *msp, float offset)
{
    int i, inNotFilled=0, descShown=0;
    float  
	x, y, 
	xold=0, yold=0;

    offset += LeftBorder + qlen4 -1;
    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);

    for (i = 0; i < qlen; i++)
      {
	const int xyVal = g_array_index(msp->xy, int, i);
	
	if (xyVal == XY_NOT_FILLED)
	  {
	    inNotFilled = 1;
	  }
	else
	  {
	    x = seq2graphY(i);
	    y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
	    
	    if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE)) 
	      {
		graphLine(yold, xold, y, x);
		
		if (fsAnnRightOn && msp->desc && !descShown) 
		  {
		    int linecolor = graphColor(BLACK);
		    graphText(msp->desc, offset, xold);
		    graphColor(linecolor);
		    descShown = 1;
		  }
	      }

	    xold = x;
	    yold = y;
	    inNotFilled = 0;
	  }
      }
    
    graphColor(oldcolor); graphLinewidth(oldLinew);
}


static int isHorizontalMSP(MSP *msp) {

    if (!msp->qname)
      g_error("No qname set in MSP - I need this to associate it with one of the sequences");

    if (!strcasecmp(msp->qname, qname) || !strcmp(msp->qname, "@1"))
        return 1;
    else
        return 0;
}


static int isVerticalMSP(MSP *msp) {

    if (!msp->qname) 
      g_error("No qname set in MSP - I need this to associate it with one of the sequences");

    if (!strcasecmp(msp->qname, sname) || !strcmp(msp->qname, "@2"))
        return 1;
    else
        return 0;
}



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
static void drawGenes(MSP *msp)
{
  float sx, ex, sy, ey, midx, midy, x, y, 
    height,		/* box height/width */
    boxHeight=10,
    textHeight,
    oldLinew;
  int i;
  float posx=0, posy=0;		/* Next series to be drawn */

  float forward_y, reverse_y ;


  forward_y = TopBorder + slen4 + 10 ;
  reverse_y = forward_y + 30 ;


  if (fsArr)
    {
      /* Loop through each feature-series in the feature-series array and set coords to 0 */
      for (i = 0; i < gArrayGetLen(fsArr); i++)
	{
	  FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
	  fs->y = fs->x = 0;
	}
    }

  textHeight = fonth;
  boxHeight = fonth;

  oldLinew = graphLinewidth(1);


  BlxStrand strand = BLXSTRAND_NONE;
  
  if (selfcomp || !reversedScale) 
    strand = BLXSTRAND_FORWARD;
  else 
    strand = BLXSTRAND_REVERSE;


  for (; msp; msp = msp->next)
    {    
      height = boxHeight;

      if (msp->score < 0.0 || mspHasFs(msp))
	{
	  sx = seq2graphX(mspGetQStart(msp));
	  ex = seq2graphX(mspGetQEnd(msp));

	  if (msp->qStrand != strand)
	    y = reverse_y ;
	  else
	    y = forward_y ;

	  if (msp->score == -1.0) /* EXON */
	    {
	      oldcolor = graphColor(BLUE); 
	      graphRectangle(sx, y, ex, y + height);
	    }
	  else if (msp->score == -2.0) /* INTRON */
	    {
	      oldcolor = graphColor(BLUE); 
	      midx = 0.5 * (sx + ex) ;
	      graphLine (sx, y + height/2, midx, y) ;
	      graphLine (ex, y + height/2, midx, y) ;
	    }
	  else if (mspHasFs(msp)) /* FEATURE SEGMENT - COLOURED BOX OR LINE */
	    {
              if (!mspShowFs(msp))
		continue;

	      /* Adjust height to score */
	      if (msp->score > 0)
		{
		  height = boxHeight * msp->score/100.0;
		}

	      if (fsBottomOn && (isHorizontalMSP(msp) || (selfcomp && isVerticalMSP(msp))))
		{
		  /* HORIZONTAL AXIS (X) */
		    
		  if (msp->type == BLXMSP_XY_PLOT)
		    {
		      XdrawSEGxy(msp, mspGetFsRightEdge(msp, &posx, fonth*(fsPlotHeight+1)));
		    }
		  else if (msp->type == BLXMSP_FS_SEG)
		    {
		      XdrawSEG(msp, mspGetFsRightEdge(msp, &posx, boxHeight+textHeight));
		    }
		}

	      if (fsRightOn &&  (isVerticalMSP(msp) || (selfcomp && isHorizontalMSP(msp))))
		{
		/* VERTICAL AXIS (Y) */

		if (msp->type == BLXMSP_XY_PLOT)
		  {
		    YdrawSEGxy(msp, mspGetFsBottomEdge(msp, &posy, fonth*(fsPlotHeight+1)));
		  }
		else if (msp->type == BLXMSP_FS_SEG)
		  {
		    YdrawSEG(msp, mspGetFsBottomEdge(msp, &posy, boxHeight+textHeight));
		  }
		}
	    }


	  if (selfcomp) /* Draw the boxes on vertical axes too */
	    {
	      sy = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/zoom);
	      ey = ceil((float)(mspGetQEnd(msp)+MSPoffset - qoffset)/zoom);
		
	      sy += TopBorder-1;
	      ey += TopBorder-1;
		
	      x = LeftBorder + qlen4 + 10;
              
	      if (msp->qStrand != strand)
                {
                  x += 20;
                }
		
	      if (msp->score == -1.0) /* EXON */
	        {
		  oldcolor = graphColor(BLUE); 
		  graphRectangle(x, sy, x + height, ey);
		}
	      else if (msp->score == -2.0) /* INTRON */
		{
		  oldcolor = graphColor(BLUE); 
		  midy = 0.5 * (sy + ey) ;
		  graphLine (x + height/2, sy, x, midy) ;
		  graphLine (x + height/2, ey, x, midy) ;
		}
	    }

	  graphColor(oldcolor); 
	  graphLinewidth(oldLinew);
	}
    }

  return ;
}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


static void drawAllFeatures(MSP *msp)
{
  float sx, ex, 
    y, 
    height,		/* box height/width */
    boxHeight=10,
    textHeight,
    oldLinew;
  int i;
  float posx=0, posy=0;		/* Next series to be drawn */

  int top, bottom ;
  float feature_boundary, feature_top, feature_bottom, feature_strand ;
  float forward_y, reverse_y, depth ;
  float old_line_width ;


  /* Set forward/reverse strand gene drawing positions. */
  graphGetBounds(&top, &bottom) ;
  feature_boundary = 3.0 ;				    /* Allows for line thickness etc... */
  feature_top = TopBorder + slen4 ;
  feature_bottom = bottom ;
  feature_strand = feature_top + ((feature_bottom - feature_top + 1) / 2) ;

  forward_y = feature_top + feature_boundary ;
  reverse_y = feature_strand + feature_boundary ;
  depth = (feature_strand - feature_boundary) - (feature_top + feature_boundary) - fonth ;


  /* Draw a strand separator. */
  old_line_width = graphLinewidth(2) ;
  graphLine(LeftBorder - 1, feature_strand, LeftBorder - 1 + qlen4, feature_strand) ;
  graphLinewidth(old_line_width) ;


  /* Now draw the genes. */
  drawGenes(msp, forward_y, reverse_y, depth) ;


  if (fsArr)
    {
      /* Loop through each feature-series in the feature-series array and set coords to 0 */
      for (i = 0; i < gArrayGetLen(fsArr); i++)
	{
	  FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
	  fs->y = fs->x = 0;
	}
    }

  textHeight = fonth;
  boxHeight = fonth;

  oldLinew = graphLinewidth(1);

  BlxStrand strand = BLXSTRAND_NONE;
  
  if (selfcomp || !reversedScale) 
    strand = BLXSTRAND_FORWARD;
  else 
    strand = BLXSTRAND_REVERSE;


  for (; msp; msp = msp->next)
    {    
      height = boxHeight;

      if (mspHasFs(msp))
	{
	  sx = seq2graphX(mspGetQStart(msp));
	  ex = seq2graphX(mspGetQEnd(msp));

	  if (msp->qStrand != strand)
	    y = reverse_y ;
	  else
	    y = forward_y ;

          if (!mspShowFs(msp))
	    continue;

	  /* Adjust height to score */
	  if (msp->score > 0)
	    {
	      height = boxHeight * msp->score / 100.0;
	    }

	  if (fsBottomOn && (isHorizontalMSP(msp) || (selfcomp && isVerticalMSP(msp))))
	    {
	      /* HORIZONTAL AXIS (X) */
		    
	      if (msp->type == BLXMSP_XY_PLOT)
		{
		  XdrawSEGxy(msp, mspGetFsRightEdge(msp, &posx, fonth*(fsPlotHeight+1)));
		}
	      else if (msp->type == BLXMSP_FS_SEG)
		{
		  XdrawSEG(msp, mspGetFsRightEdge(msp, &posx, boxHeight+textHeight));
		}
	    }

	  if (fsRightOn &&  (isVerticalMSP(msp) || (selfcomp && isHorizontalMSP(msp))))
	    {
	      /* VERTICAL AXIS (Y) */

	      if (msp->type == BLXMSP_XY_PLOT)
		{
		  YdrawSEGxy(msp, mspGetFsBottomEdge(msp, &posy, fonth*(fsPlotHeight+1)));
		}
	      else if (msp->type == BLXMSP_FS_SEG)
		{
		  YdrawSEG(msp, mspGetFsBottomEdge(msp, &posy, boxHeight+textHeight));
		}
	    }

	  graphColor(oldcolor); 
	  graphLinewidth(oldLinew);
	}
    }

  return ;
}


static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth)
{
  gboolean bump_genes = FALSE ;				    /* Make this user settable... */

  float height,		/* box height/width */
    sy, ey, midy, x, y, oldLinew ;

  height = fonth ;

  oldLinew = graphLinewidth(1);

  BlxStrand strand = BLXSTRAND_NONE;
  if (selfcomp || !reversedScale) 
    strand = BLXSTRAND_FORWARD;
  else 
    strand = BLXSTRAND_REVERSE;


  if (bump_genes)
    {
      MSP *gene_msp, *tmp ;
      GList *exon_intron_list = NULL ;
      GeneStrandStruct strand_genes = {NULL} ;

      /* MSP's are in any old order and only some are exons/introns, here we make a list of
       * intron/exon MSP's and then sort that list into transcripts for drawing... */

      /* Sort all introns/exons by gene and position within gene. */
      gene_msp = tmp = msp ;
      for ( ; tmp ; tmp = tmp->next)
	{
	  if (tmp->score < 0)
	    exon_intron_list = g_list_append(exon_intron_list, tmp) ;
	}
      exon_intron_list = g_list_sort(exon_intron_list, compareMSPs) ;


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
      g_list_foreach(exon_intron_list, printMSP, NULL) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


      /* Produce two lists, one for each strand, of genes sorted by position from the sorted intron/exon list. */
      strand_genes.strand = '+' ;
      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
      strand_genes.forward_genes = g_list_sort(strand_genes.forward_genes, compareGenes) ;
      strand_genes.forward_genes = g_list_first(strand_genes.forward_genes) ;

      strand_genes.strand = '-' ;
      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
      strand_genes.reverse_genes = g_list_sort(strand_genes.reverse_genes, compareGenes) ;
      strand_genes.reverse_genes = g_list_first(strand_genes.reverse_genes) ;

      /* Set y offsets for forward and reverse stranded genes. */
      setYoffsets(strand_genes.forward_genes, forward_y, forward_y + depth) ;
      setYoffsets(strand_genes.reverse_genes, reverse_y, reverse_y + depth) ;


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
      printf("Forward:\n") ;
      g_list_foreach(strand_genes.forward_genes, printGene, NULL) ;
      printf("Reverse:\n") ;
      g_list_foreach(strand_genes.reverse_genes, printGene, NULL) ;
      printf("\n") ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

      /* Draw the strands. */
      oldcolor = graphColor(BLUE); 
      g_list_foreach(strand_genes.forward_genes, drawGenesCB, &strand_genes) ;

      oldcolor = graphColor(MAGENTA); 
      g_list_foreach(strand_genes.reverse_genes, drawGenesCB, &strand_genes) ;

      graphColor(oldcolor); 

    }
  else
    {
      for (; msp; msp = msp->next)
	{    
	  if (msp->score < 0)
	    {
	      if (msp->qStrand != strand)
		y = reverse_y ;
	      else
		y = forward_y ;

	      oldcolor = graphColor(BLUE); 
	      drawMSPGene(msp, y) ;

	      if (selfcomp) /* Draw the boxes on vertical axes too */
		{
		  sy = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/zoom);
		  ey = ceil((float)(mspGetQEnd(msp)+MSPoffset - qoffset)/zoom);
		
		  sy += TopBorder-1;
		  ey += TopBorder-1;
		
		  x = LeftBorder + qlen4 + 10;
		  if (msp->qStrand != strand) x += 20;
		
		  if (msp->score == -1.0) /* EXON */
		    {
		      oldcolor = graphColor(BLUE); 
		      graphRectangle(x, sy, x + height, ey);
		    }
		  else if (msp->score == -2.0) /* INTRON */
		    {
		      oldcolor = graphColor(BLUE); 
		      midy = 0.5 * (sy + ey) ;
		      graphLine (x + height/2, sy, x, midy) ;
		      graphLine (x + height/2, ey, x, midy) ;
		    }
		}

	      graphColor(oldcolor) ;
	    }
	}
    }

  graphLinewidth(oldLinew) ;

  return ;
}


/* A GCompareFunc() to compare gene MSPs.... */
gint compareMSPs(gconstpointer a, gconstpointer b)
{
  int result = 0 ;
  MSP *msp_a = (MSP *)a, *msp_b = (MSP *)b ;

  if (!(result = strcmp(msp_a->sname, msp_b->sname)))
    {
      if (mspGetSStart(msp_a) < mspGetSStart(msp_b))
	result = -1 ;
      else if (mspGetSStart(msp_a) > mspGetSStart(msp_b))
	result = 1 ;
      else
	{
	  /* I actually don't think this should ever happen as it means either there are
	   * duplicates or worse there are overlapping introns/exons within a gene.... */
	  if (mspGetSEnd(msp_a) < mspGetSEnd(msp_b))
	    result = -1 ;
	  else if (mspGetSEnd(msp_a) > mspGetSEnd(msp_b))
	    result = 1 ;
	  else
	    result = 0 ;
	}
    }

  return result ;
}


/* A GFunc() to record start/end of genes. */
static void getGenePositionsCB(gpointer data, gpointer user_data)
{
  MSP *msp = (MSP *)data ;
  GeneStrand strand_genes = (GeneStrand)user_data ;
  GList *gene_list ;
  GeneData gene ;

  /* We need genes added according to strand... */
  if (msp->qframe[1] == strand_genes->strand)
    {
      const char *curr_name = NULL ;
      int curr_length ;

      if (strand_genes->strand == '+')
	gene_list = strand_genes->forward_genes ;
      else
	gene_list = strand_genes->reverse_genes ;

      if (gene_list)
	{
	  /* Look at last element, this is the last gene we added. */
	  gene_list = g_list_last(gene_list) ;

	  curr_name = (((GeneData)(gene_list->data))->name) ;
	  curr_length = strlen(curr_name) - 1 ;
	}

      /* If there's no gene or we are starting a new gene then just add to the list,
       * otherwise record latest msp end position. */
      if (!gene_list || strncmp(mspGetSName(msp), curr_name, curr_length) != 0)
	{
	  gene = g_new0(GeneDataStruct, 1) ;
	  gene->name = mspGetSName(msp) ;
	  gene->start = mspGetSStart(msp) ;
	  gene->end = mspGetSEnd(msp) ;
	  gene->strand = msp->qframe[1] ;
	  gene->msp_start = msp ;

	  gene_list = g_list_append(gene_list, gene) ;
	}
      else
	{
	  gene = (GeneData)(gene_list->data) ;

	  gene->end = mspGetSEnd(msp) ;
	  gene->msp_end = msp ;
	}

      if (strand_genes->strand == '+')
	strand_genes->forward_genes = gene_list ;
      else
	strand_genes->reverse_genes = gene_list ;
    }

  return ;
}



/* A GCompareFunc() to compare genes for position.... */
gint compareGenes(gconstpointer a, gconstpointer b)
{
  int result = 0 ;
  GeneData gene_a = (GeneData)a, gene_b = (GeneData)b ;

  if (gene_a->strand == '+' && gene_b->strand == '-')
    result = -1 ;
  else if (gene_a->strand == '-' && gene_b->strand == '+')
    result = 1 ;
  else
    {
      if (gene_a->start < gene_b->start)
	result = -1 ;
      else if (gene_a->start > gene_b->start)
	result = 1 ;
      else
	result = 0 ;
    }

  return result ;
}



static void setYoffsets(GList *gene_list, float min_offset, float max_offset)
{
  GList *curr_ptr, *next_ptr ;
  float curr_y = 0.0, bump_incr = 0.5 ;


  /* Go through the gene list comparing, reordering and assigning y offsets to the genes
   * so they can be drawn without overlapping. */
  curr_ptr = gene_list ;
  next_ptr = curr_ptr->next ;
  curr_y = min_offset - bump_incr ;

  /* Go through all genes. */
  while (curr_ptr)
    {
      GList *list_ptr = g_list_first(gene_list) ;

      /* For each gene look at all the preceding genes and to decide its y offset so it is
       * bumped correctly. */
      while (list_ptr)
	{
	  GeneData list_gene = (GeneData)list_ptr->data ;
	  GeneData curr_gene = (GeneData)curr_ptr->data ;

	  if (list_ptr == curr_ptr)
	    {
	      /* This gene overlaps all previous ones and needs a new offset. */
	      curr_y += bump_incr ;
	      if (curr_y > max_offset)
		curr_y = max_offset ;

	      curr_gene->y_pos = curr_y ;
	      break ;
	    }
	  else if (curr_gene->start > list_gene->end)
	    {
	      /* This gene coes not overlap the list gene so give is the same offset. */
	      curr_gene->y_pos = list_gene->y_pos ;

	      if (list_ptr->next != curr_ptr)
		{
		  list_ptr = g_list_remove(list_ptr, curr_gene) ;
		  list_ptr = g_list_insert_before(gene_list, list_ptr->next, curr_gene) ;
		}

	      break ;
	    }
	  else
	    {
	      /* This gene overlaps the list gene so move on to the next one. */
	      list_ptr = g_list_next(list_ptr) ;
	    }
	}

      /* update curr/next until we get to the end of the list... */
      if ((curr_ptr = next_ptr))
	next_ptr = curr_ptr->next ;
    }

  return ;
}


/* A GFunc() to record start/end of genes. */
static void drawGenesCB(gpointer data, gpointer user_data_unused)
{
  GeneData gene = (GeneData)data ;
  MSP *msp ;

  msp = gene->msp_start ;
  do
    {
      drawMSPGene(msp, gene->y_pos) ;

      msp = msp->next ;
    } while (msp && msp != gene->msp_end) ;

  return ;
}


static void drawMSPGene(MSP *msp, float y_offset)
{
  float height,		/* box height/width */
    sx, ex, midx ;

  height = fonth ;

  sx = seq2graphX(mspGetQStart(msp));
  ex = seq2graphX(mspGetQEnd(msp));

  if (msp->score == -1.0) /* EXON */
    {
      graphRectangle(sx, y_offset, ex, y_offset + height);
    }
  else if (msp->score == -2.0) /* INTRON */
    {

      midx = 0.5 * (sx + ex) ;
      graphLine (sx, y_offset + height/2, midx, y_offset) ;
      graphLine (ex, y_offset + height/2, midx, y_offset) ;
    }

  return ;
}


//static void printMSP(gpointer data, gpointer user_data)
//{
//  MSP *msp = (MSP *)data ;
//  
//  printf("%s:\t%d,%d\t->\t%d,%d\n",
//	 msp->sname, msp->sstart, msp->send, msp->qstart, msp->qend) ;
//  
//
//  return ;
//}
//

//static void printGene(gpointer data, gpointer user_data)
//{
//  GeneData gene = (GeneData)data ;
//  
//  printf("%s: '%c' %d -> %d   is at position: %f\n",
//	 gene->name, gene->strand, gene->start, gene->end, gene->y_pos) ;
//  
//
//  return ;
//}




static void graphPixelLine(int strength, int sx, int sy, int ex, int ey)
{
    int i, inc, len, dotpos, x, y;
    UCHAR dotValue;

    if (sx < ex) 
	inc = 1;
    else
	inc = -1;

    len = abs(ex - sx +1);
    
    if (strength < 256)
	dotValue = (UCHAR)strength;
    else 
	dotValue = 255;

    for (i = 0; i < len; i++)
    {
	x = sx + i*inc;
	y = sy + i;

	if (x >= 0 && x < qlen4 && sy >= 0 && sy < slen4) 
	{
	    dotpos = qlen4*(y) + x;

	    if (dotpos < 0 || dotpos > datalen-1) {
		g_critical("Pixel out of bounds (0-%d) in graphPixelLine: %d."
                           "Crash imminent.", datalen-1, dotpos);
	    }
	    else
	    {
		if (dotValue > HSPpixels[dotpos]) HSPpixels[dotpos] = dotValue;
	    }
	}
    }
}


static int score2color(gdouble score)
{
    if (score < 75.0) return DARKRED;
    if (score < 100.0) return MAGENTA;
    return RED;
}



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
/* currently unused, was for setting line widths according to score. */
static int score2width(int score)
{
    if (score < 75) return 1;
    if (score < 100) return 2;
    return 3;
}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



int illegalSubmatrix(int sx, int sy, int shift)
{
    int dotposq, dotposs, ql, sl;

    dotposq = (sx+shift)/zoom;
    dotposs = (sy+shift)/zoom;
    
    ql = (sx+shift) - dotposq*zoom;
    sl = (sy+shift) - dotposs*zoom;

    if (sl >= ql) 
	return 0;		/* legal */
    else 
	return 1;		/* illegal */
}


int illegalSubmatrixRev(int sx, int sy, int shift)
{
    int dotposq, dotposs, ql, sl;

    dotposq = (sx-shift)/zoom;
    dotposs = (sy+shift)/zoom;
    
    ql = (sx-shift) - dotposq*zoom;

    /* Set the origin (0,0) to the bottom left corner of submatrix
       Zoom = pixels/submatrix */
    sl = zoom-1 - (sy - dotposs*zoom);

    if (sl >= ql) 
	return 0;		/* legal */
    else 
	return 1;		/* illegal */
}


/* Notes:
 * msp->[qs][start|end] coordinates go [1..n+1] while real dots go [0..n]
 * In other words, msp-coords [6,6] will make a dot at [5,5].
 */
static void drawBlastHSPs(void)
{
  int i,
	sx, ex, sy, ey, /* Screen coordinates [0..n] */
	matchlen, shift;
    char strand;
    const char *MSPsname;
    MSP *msp;

    strand = reversedScale ? '-' : '+';

    if (BlastMode == BLASTNOTHING) BlastMode = BLASTRED;

    if (BlastMode == BLASTGREY) /* Reset background */
	for (i=0; i < datalen; i++) HSPpixels[i] = 0;

    for (msp = MSPlist; msp;  msp = msp->next)
    {    
      if ((MSPsname = strchr(mspGetSName(msp), ':')))
	MSPsname++;
      else
	MSPsname = mspGetSName(msp);

      if (!strcmp(MSPsname, sname))
	{
	  matchlen = mspGetSRangeLen(msp) / zoom;
	    

/* printf("\n%s: %d,%d - %d,%d\n", 
       msp->sname, msp->qstart, msp->sstart, msp->qend, msp->send);
*/

	    if (mspGetQStart(msp) < mspGetQEnd(msp)) {

		sx = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/resfac) -1;
		sy = mspGetSStart(msp) - soffset -1;

		/* Check if we're in an illegal part of submatrix
		 * - We only want to compress in legal submatrix parts */
		shift = 0;
		while(illegalSubmatrix(sx, sy, shift)) shift++;

		sx = (sx + shift)/zoom -shift;
		sy = (sy + shift)/zoom -shift;
		
		ex = sx + (matchlen-1);
	    }
	    else {
		sx = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/resfac) -1;
		sy = mspGetSStart(msp) - soffset -1;

		/* Check if we're in an illegal part of submatrix
		 * - We only want to compress in legal submatrix parts */
		shift = 0;
		while(illegalSubmatrixRev(sx, sy, shift)) shift++;

		sx = (sx - shift)/zoom + shift;
		sy = (sy + shift)/zoom - shift;
		
		ex = sx - (matchlen-1);

		/* sx = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/resfac);
		sx = sx/zoom - 1;
		ex = sx - (matchlen-1);
		sy = (msp->sstart - soffset)/zoom - 1;*/

	    }

	    if (reversedScale) {
		sx = RightBorder-LeftBorder - sx;
		ex = RightBorder-LeftBorder - ex;
	    }

	    ey = sy + matchlen-1;
	    
	    if (BlastMode == BLASTRED) {
		graphColor(RED);
		graphLinewidth(1 /* score2width(msp->score) */);

		graphLine(LeftBorder + sx, TopBorder + sy,
			  LeftBorder + ex, TopBorder + ey);
		graphColor(BLACK);
	    }
	    else if (BlastMode == BLASTFUNC) {
		graphColor(score2color(msp->score));
		graphLinewidth(1);

		graphLine(LeftBorder + sx, TopBorder + sy,
			  LeftBorder + ex, TopBorder + ey);
		graphColor(BLACK);
	    }
	    else if (BlastMode == BLASTGREY) {
		/* strength = win*msp->score/abs(msp->send - msp->sstart); */
		int strength = (int)msp->score;
		graphPixelLine(strength, sx, sy, ex, ey);
	    }
	    else {
		g_critical("Unknown BlastMode = bug -> contact Sanger");
	    }
	}
    }

}

static void drawBlastHSPgray()
{
    BlastMode = BLASTGREY;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void drawBlastHSPline()
{
    BlastMode = BLASTRED;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void drawBlastHSPlinef(void)
{
    BlastMode = BLASTFUNC;

    BlastHSPsOn = 1;
    dotterRedraw();
}
static void clearHSPs(void)
{
    BlastMode = BLASTNOTHING;

    BlastHSPsOn = 0;
    dotterRedraw();
}


/* x and y real sequence coordinates */
static void setQposSpos(int x, int y)
{
  int position ;

  if (reversedScale)
    position = qlen - x*resfac + qoffset ;
  else
    position = x*resfac + 1 + qoffset ;

  sprintf(qpos, "%d", position) ;


  sprintf(spos, "%d", y+1 + soffset);

  return ;  
}


/* DRAWALIGNMENT draws the alignment of the diagonal at the crosshair
 * Note: x and y are in sequence positions (Real sequence, starting at 0)
 */
static void drawAlignment(int x, int y)
{
    int    i, frame;
    char  *q, *s, *qrc;

    if (x < 0) x = 0;
    else if (x > qlen-1) x = qlen-1;
    if (y < 0) y = 0;
    else if (y > slen-1) y = slen-1;

    oldx = x;
    oldy = y;

    if (!graphActivate(alnGraph))
      return;

    q = qseq + x - ALIGNLEN/2;
    s = sseq + y - ALIGNLEN/2;

    setQposSpos(x, y);

    if (blastp || (blastn && !crickOnly)) 
    {
	for (i = 0; i < ALIGNLEN; i++) 
	{
	    qcolors[i] = scolors[i] = alnBackgColor;

	    qseqDisp[i] = sseqDisp[i] = 0;
	    
	    /* Handle sticky ends */
	    if (q+i < qseq || q+i > qseq+qlen-1) {
		qseqDisp[i] = ' ';
		if (s+i >= sseq && s+i <= sseq+slen-1) { 
		    sseqDisp[i] = s[i];
		}
	    }
	    if (s+i < sseq || s+i > sseq+slen-1) {
		sseqDisp[i] = ' ';
		if (q+i >= qseq && q+i <= qseq+qlen-1) { 
		    qseqDisp[i] = q[i];
		}
	    }
	    
	    if (qseqDisp[i] || sseqDisp[i]) continue;
	    
	    qseqDisp[i] = q[i];
	    sseqDisp[i] = s[i];
	    
	    /* Highlight matching residues */
	    if (blastp) {
		if (q[i] == s[i])
		    qcolors[i] = scolors[i] = TINT_CYAN;
		else if ( MATRIX[atob[(int)q[i]]-1 ][ atob[(int)s[i]]-1 ] > 0)
		    qcolors[i] = scolors[i] = TINT_HIGHLIGHT2;
	    }
	    else if (blastn && q[i] == s[i])
		qcolors[i] = scolors[i] = TINT_CYAN;

	}

	graphBoxDraw(qseqbox, BLACK, backgColor);
	graphBoxDraw(sseqbox, BLACK, backgColor);
	
	graphBoxDraw(qposbox, BLACK, backgColor);
	graphBoxDraw(sposbox, BLACK, backgColor);
    }

    /* Crick strand - qseq reversed and complemented */

    if (blastn && !watsonOnly)
    {

	qrc = qrevcompl + qlen-1 - x - ALIGNLEN/2;

	for (i = 0; i < ALIGNLEN; i++) 
	{
	    qcolorsCrick[i] = scolorsCrick[i] = alnBackgColor; /* i.e. Background */

	    qseqDispCrick[i] = sseqDispCrick[i] = 0;

	    /* Handle sticky ends */
	    if (qrc+i < qrevcompl || qrc+i > qrevcompl+qlen-1) {
		qseqDispCrick[i] = ' ';
		if (s+i >= sseq && s+i <= sseq+slen-1) { 
		    sseqDispCrick[i] = s[i];
		}
	    }
	    if (s+i < sseq || s+i > sseq+slen-1) {
		sseqDispCrick[i] = ' ';
		if (qrc+i >= qrevcompl && qrc+i <= qrevcompl+qlen-1) { 
		    qseqDispCrick[i] = qrc[i];
		}
	    }
	    
	    if (qseqDispCrick[i] || sseqDispCrick[i]) continue;
	    
	    qseqDispCrick[i] = qrc[i];
	    sseqDispCrick[i] = s[i];
	    
	    /* Highlight matching residues */
	    if (qrc[i] == s[i])
		qcolorsCrick[i] = scolorsCrick[i] = TINT_CYAN;
	}

	graphBoxDraw(qseqboxCrick, BLACK, backgColor);
	graphBoxDraw(sseqboxCrick, BLACK, backgColor);
	
	graphBoxDraw(qposboxCrick, BLACK, backgColor);
	graphBoxDraw(sposboxCrick, BLACK, backgColor);

    }
    
    if (blastx)
    {
	for (i = 0; i < ALIGNLEN; i++) scolors[i] = alnBackgColor;

	for (frame = 0; frame < 3; frame++)
	{
	    for (i = 0; i < ALIGNLEN; i++) 
	    {
		q = pepqseq[frame] + x - ALIGNLEN/2;

		xqcolors[frame][i] = alnBackgColor;

		xqseqDisp[frame][i] = sseqDisp[i] = 0;
	    
		/* Handle sticky ends */
		if (q+i < pepqseq[frame] || q+i > pepqseq[frame]+qlen/3-1) {
		    xqseqDisp[frame][i] = ' ';
		    if (s+i >= sseq && s+i <= sseq+slen-1) { 
			sseqDisp[i] = s[i];
		    }
		}
		if (s+i < sseq || s+i > sseq+slen-1) {
		    sseqDisp[i] = ' ';
		    if (q+i >= pepqseq[frame] && q+i <= pepqseq[frame]+qlen/3-1) { 
			xqseqDisp[frame][i] = q[i];
		    }
		}
	    
		if (xqseqDisp[frame][i] || sseqDisp[i]) continue;
	    
		xqseqDisp[frame][i] = q[i];
		sseqDisp[i] = s[i];
	    
		/* Highlight matching residues */
		if (q[i] == s[i])
		    xqcolors[frame][i] = scolors[i] = TINT_CYAN;
		else if ( MATRIX[atob[(int)q[i]]-1 ][ atob[(int)s[i]]-1 ] > 0 ) {
		    xqcolors[frame][i] = TINT_HIGHLIGHT2;
		    if (scolors[i] != TINT_CYAN) scolors[i] = TINT_HIGHLIGHT2;
		}
	    }
	    graphBoxDraw(xqseqbox[frame], BLACK, backgColor);
	}
	graphBoxDraw(qposbox, BLACK, backgColor);

	graphBoxDraw(sseqbox, BLACK, backgColor);
	graphBoxDraw(sposbox, BLACK, backgColor);
    }

    graphActivate(dotterGraph);
}


static void togglePrintColors(void)
{
    backgColor = (backgColor == LIGHTGRAY ?  WHITE : LIGHTGRAY);
    alnBackgColor = (alnBackgColor == TINT_LIGHTGRAY ?  TINT_WHITE : TINT_LIGHTGRAY);
    initAlignment();
}


static void dotterPrint(void)
{
    printMode = 1;
    dotterRedraw();

    graphPrint();

    printMode = 0;
    dotterRedraw();
}


static void setAlnLen(void)
{
  ACEIN len_in;

  if (!(len_in = messPrompt ("Give Alignment length", 
			     messprintf("%d", ALIGNLEN), 
			     "iz", 0)))
    return;

  aceInInt (len_in, &ALIGNLEN);
  aceInDestroy (len_in);

  if (!(ALIGNLEN % 2)) ALIGNLEN++;
  
  if (ALIGNLEN > MAXALIGNLEN) ALIGNLEN = MAXALIGNLEN;
  else if (ALIGNLEN < 0) ALIGNLEN = DEFAULTALIGNLEN ;
  
  initAlignment();

  return;
} /* setAlnLen */


static void fsSelFinish(void)
{
#ifndef DOTTER
    if (blixemWindow) blviewRedraw();
#endif
    if (graphActivate(dotterGraph)) dotterRedraw();
}


static void fsSelAll(void)
{
  int i;
  graphActivate(fsGraph);
  
  for (i = 0; i < gArrayGetLen(fsArr); i++) 
    {
      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
      fs->on = 1;
      graphBoxDraw(fsBoxStart+i, WHITE, BLACK);
    }
    
  fsSelFinish();
}

static void fsSelNone(void)
{
  int i;
  graphActivate(fsGraph);
  
  for (i = 0; i < gArrayGetLen(fsArr); i++) 
    {
      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
      fs->on = 0;
      graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
    }
    
  fsSelFinish();
}

static void fsSelNoCurves(void)
{
  int i;
  graphActivate(fsGraph);
  
  for (i = 0; i < gArrayGetLen(fsArr); i++) 
    {
      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
      
      if (fs->xy) 
	{
	  fs->on = 0;
	  graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
	}
    }
    
  fsSelFinish();
}

static void fsSelNoSegments(void)
{
  int i;
  graphActivate(fsGraph);
  
  for (i = 0; i < gArrayGetLen(fsArr); i++) 
    {
      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
      
      if (!fs->xy) 
	{
	  fs->on = 0;
	  graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
	}
    }
    
  fsSelFinish();
}

static void fsToggleBottom(void)
{
    fsBottomOn = !fsBottomOn;
    selectFeatures();
    fsSelFinish();
}
static void fsToggleRight(void)
{
    fsRightOn = !fsRightOn;
    selectFeatures();
    fsSelFinish();
}
static void fsToggleAnnBottom(void)
{
    fsAnnBottomOn = !fsAnnBottomOn;
    selectFeatures();
    fsSelFinish();
}
static void fsToggleAnnRight(void)
{
    fsAnnRightOn = !fsAnnRightOn;
    selectFeatures();
    fsSelFinish();
}
static void fsToggleEndLines(void)
{
    fsEndLinesOn = !fsEndLinesOn;
    selectFeatures();
    fsSelFinish();
}
static void fsSel(int box, double x_unused, double y_unused, int modifier_unused)
{
    if (box-fsBoxStart < 0 || box-fsBoxStart > gArrayGetLen(fsArr))
	return;
    
    FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, box - fsBoxStart);
    int *on = &fs->on;

    graphActivate(fsGraph);


    if (*on) {
	*on = 0;
	graphBoxDraw(box, BLACK, WHITE);
    }
    else {
	*on = 1;
	graphBoxDraw(box, WHITE, BLACK);
    }

    fsSelFinish();
}
static void setfsPlotHeight(char *cp)
{
    fsPlotHeight = atof(cp);
#ifndef DOTTER
    if (blixemWindow) blviewRedraw();
#endif
    if (graphActivate(dotterGraph)) dotterRedraw();
}
    

void selectFeatures(void)
{
  int i, box;
  float y=1.0;

  if (!graphActivate(fsGraph))
    {
      fsGraph = graphCreate (TEXT_SCROLL, "Feature Series Selection Tool", 0, 0, 0.4, 0.6);
    }
  graphPop();
  graphRegister(PICK, fsSel);

  if (dotterGraph) {

    box = graphButton("Show bottom series", fsToggleBottom, 1, y);
    if (!fsBottomOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);
    y += 1.5;

    box = graphButton("Show right series", fsToggleRight, 1, y);
    if (!fsRightOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);
    y += 1.5;

    box = graphButton("Show bottom annotations", fsToggleAnnBottom, 1, y);
    if (!fsAnnBottomOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);
    y += 1.5;

    box = graphButton("Show right annotations", fsToggleAnnRight, 1, y);
    if (!fsAnnRightOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);
    y += 1.5;

    box = graphButton("Draw lines at segment ends", fsToggleEndLines, 1, y);
    if (!fsEndLinesOn) graphBoxDraw(box, BLACK, WHITE);
    else graphBoxDraw(box, WHITE, BLACK);
    y += 2;
  }

  sprintf(fsPlotHeighttx, "%.1f", fsPlotHeight);
  graphText("Height of curves (XY series):", 1, y);
  graphTextScrollEntry (fsPlotHeighttx, 6, 4, 31, y, setfsPlotHeight);
  y += 2;

  graphLinewidth(0.1);
  graphLine(0, y, 2000, y);
  y += 0.5;

  graphText("Pick to select/unselect series", 1, y);
  y += 1.5;
  graphButton("All", fsSelAll, 1, y);
  graphButton("No curves", fsSelNoCurves, 13, y);
  graphButton("No segments", fsSelNoSegments, 24, y);
  fsBoxStart = 1+graphButton("None", fsSelNone, 6, y);
  y += 2;

  graphTextBounds(50, gArrayGetLen(fsArr) * 1.5 + y + 5);

  for (i = 0; i < gArrayGetLen(fsArr); i++)
    {
      float  margin = 0.1;
      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);

      box = graphBoxStart();
      graphText(fs->name, 1, y);      
      graphRectangle(1 - margin, y - margin, 1 + margin + strlen(fs->name), y + 1 + margin);
      graphBoxEnd();

      if (!fs->on) 
	{
	  graphBoxDraw(box, BLACK, WHITE);
	}
      else 
	{
	  graphBoxDraw(box, WHITE, BLACK);
	}
	
      y += 1.5;
    }

  graphRedraw();

  return ;
}


float fsTotalHeight(MSP *msplist)
{
    int i;
    float maxy = 0;
	
    if (!fsArr || !gArrayGetLen(fsArr))
      {
	return 0.0;
      }

    for (i = 0; i < gArrayGetLen(fsArr); i++) 
      {
	FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
	fs->y = fs->x = 0;
      }
	
    for (msp = msplist; msp; msp = msp->next) 
      {
        if (mspShowFs(msp))
	  {
	    if (msp->type == BLXMSP_XY_PLOT) 
	      {
		mspGetFsBottomEdge(msp, &maxy, fsPlotHeight+1);
	      }
	    else if (msp->type == BLXMSP_FS_SEG) 
	      {
		mspGetFsBottomEdge(msp, &maxy, 1+1);
	      }
	  }
      }
    
    return maxy + 2;
}


static void initAlignment(void)
{
    float x, y;
    int i, frame, height=0;

    {
	/* static int warned=0;
	if (!warned) g_critical("The residue colours of the Aligment tool have been corrupted by Jean."
			     " Send complaints to mieg@kaa.crbm.cnrs-mop.fr");
	warned=1;*/
    }

    if (!CrosshairON)
      {
	g_critical("Turn on the crosshair !");
	return;
      }
    
    graphActivate(dotterGraph);
    graphBoxDim (vLineBox, &x, 0, 0, 0);
    graphBoxDim (hLineBox, 0, &y, 0, 0);

    if (!graphActivate(alnGraph))
    {
      BOOL oldMap = graphSetInstallMap(TRUE);
      alnGraph = graphCreate (TEXT_SCROLL, "Dotter - Alignment Tool", 0, 0, 1.25, 0.25);
      graphMenu(alnmenu);
      graphRegister(KEYBOARD, keyboard);
      graphSetInstallMap(oldMap);
    }

    graphPop();
    graphClear();
    
    if (blastp)
	height = 7;
    else if (blastn) {
	if (watsonOnly) height = 7;
	else height = 14;
    }
    else if (blastx) 
	height = 9;

    graphTextBounds(MAXALIGNLEN, height);
    graphColor(backgColor); 
    graphRectangle(0, 0, 1000, 1000);
    graphColor(BLACK);
	
    for (i = 0; i < MAXALIGNLEN; i++)
      {
	qseqDisp[i] = sseqDisp[i] =
	  xqseqDisp[0][i] = xqseqDisp[1][i] = xqseqDisp[2][i] =
	  qseqDispCrick[i] = sseqDispCrick[i] = 0;
      }

    /* Maybe do this instead * /
    qseqDisp = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[0] = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[1] = handleAlloc(0, handle, MAXALIGNLEN+1);
    xqseqDisp[2] = handleAlloc(0, handle, MAXALIGNLEN+1);
    sseqDisp = handleAlloc(0, handle, MAXALIGNLEN+1);
    qseqDispCrick = handleAlloc(0, handle, MAXALIGNLEN+1);
    sseqDispCrick = handleAlloc(0, handle, MAXALIGNLEN+1);
    */

    if (blastx)
    {
	/* Draw static features */
	graphText(messprintf("%s 1:", qname), 0.5, 3);
	graphText(messprintf("%s 2:", qname), 0.5, 4);
	graphText(messprintf("%s 3:", qname), 0.5, 5);
	graphText(messprintf("%s:", sname), 0.5, 6);
	graphLine(NAMESIZE+4 + (float)ALIGNLEN/2, 2, NAMESIZE+4 + (float)ALIGNLEN/2, 8);
	
	/* alignment boxes */
	for (frame = 0; frame < 3; frame++)
	{
	    xqseqbox[frame] = graphBoxStart();
	    graphColorSquares (xqcolors[frame], NAMESIZE+4, 3+frame, ALIGNLEN, 1, tints);
	    graphTextPtr     (xqseqDisp[frame], NAMESIZE+4, 3+frame, ALIGNLEN);
	    graphBoxEnd ();
	}

	sseqbox = graphBoxStart();
	graphColorSquares (scolors, NAMESIZE+4, 6, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDisp, NAMESIZE+4, 6, ALIGNLEN);
	graphBoxEnd ();

	/* coordinate boxes */
	qposbox = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 1, 10);
	graphBoxEnd ();
	sposbox = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 8, 10);
	graphBoxEnd ();
    }

    if (blastp || (blastn && !crickOnly)) 
    {
	/* Draw static features */
	graphText(messprintf("%s:", qname), 0.5, 3);
	graphText(messprintf("%s:", sname), 0.5, 4);
	graphLine(NAMESIZE+2 + (float)ALIGNLEN/2, 2, NAMESIZE+2 + (float)ALIGNLEN/2, 6);
	
	/* alignment boxes */
	qseqbox = graphBoxStart();
	graphColorSquares (qcolors, NAMESIZE+2, 3, ALIGNLEN, 1, tints);
	graphTextPtr     (qseqDisp, NAMESIZE+2, 3, ALIGNLEN);
	graphBoxEnd ();
	sseqbox = graphBoxStart();
	graphColorSquares (scolors, NAMESIZE+2, 4, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDisp, NAMESIZE+2, 4, ALIGNLEN);
	graphBoxEnd ();
	    
	/* coordinate boxes */
	qposbox = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 1, 10);
	graphBoxEnd ();
	sposbox = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 6, 10);
	graphBoxEnd ();
    }

    if (blastn && !watsonOnly)
    {
	/* Draw static features */
	graphText("RevComp:", 0.5, 10);
	graphText(messprintf("%s:", sname), 0.5, 11);
	graphLine(NAMESIZE+2 + (float)ALIGNLEN/2, 9, NAMESIZE+2 + (float)ALIGNLEN/2, 13);
	
	/* alignment boxes */
	qseqboxCrick = graphBoxStart();
	graphColorSquares (qcolorsCrick, NAMESIZE+2, 10, ALIGNLEN, 1, tints);
	graphTextPtr     (qseqDispCrick, NAMESIZE+2, 10, ALIGNLEN);
	graphBoxEnd ();
	sseqboxCrick = graphBoxStart();
	graphColorSquares (scolorsCrick, NAMESIZE+2, 11, ALIGNLEN, 1, tints);
	graphTextPtr     (sseqDispCrick, NAMESIZE+2, 11, ALIGNLEN);
	graphBoxEnd ();
	
	/* coordinate boxes */
	qposboxCrick = graphBoxStart();
	graphTextPtr (qpos, NAMESIZE + ALIGNLEN/2+2, 8, 10);
	graphBoxEnd ();
	sposboxCrick = graphBoxStart();
	graphTextPtr (spos, NAMESIZE + ALIGNLEN/2+2, 13, 10);
	graphBoxEnd ();
    }    
    
    graphBoxDraw(0, backgColor, backgColor);
    graphRedraw();
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
}


static void toggleCrosshair(void)
{
  if (CrosshairON)
    {
      CrosshairON  = 0;
    }
  else
    {
      CrosshairON  = 1;
    }

  dotterRedraw();

  return ;
}

static void toggleCrosshairPos(void)
{
    CrosshairPosON = (CrosshairPosON ? 0 : 1);
    dotterRedraw();
}

static void toggleCrosshairFullscreen(void)
{
    CrosshairFullscreenON = !CrosshairFullscreenON;
    dotterRedraw();
}


static void drawCrosshair(float x, float y)
{
  if (CrosshairON)
    {
      if (CrosshairFullscreenON)
	{
	  graphBoxShift(hLineBox, 0.0, y);
	  graphBoxShift(vLineBox, x, 0.0);
	}
      else
	{
	  graphBoxShift(hLineBox, LeftBorder, y);
	  graphBoxShift(vLineBox, x, TopBorder);
	}

      if (CrosshairPosON)
	{
	  sprintf(CrosshairPosText, "%s, %s", qpos, spos);
	  graphBoxShift(CrosshairPosBox, x+10, y+10);
	}
    }

  crossx = x;
  crossy = y;

  return ;
}


static void dotterRampChange(BOOL isDrag)
{
  if (!isDrag)
    {
      if(gridOn || fsEndLinesOn || BlastHSPsOn)
	dotterRedraw();
      else
	drawAllFeatures(MSPlist);
    }
 
  if (hLineBox && vLineBox)
    drawCrosshair(crossx, crossy);
}

static void boundaries(double *x, double *y)
{
  if (*x < LeftBorder)
    *x = LeftBorder;
  else if (*x > LeftBorder + qlen4)
    *x = LeftBorder + qlen4 ;

  if (*y < TopBorder)
    *y = TopBorder ;
  else if (*y > TopBorder + slen4)
    *y = TopBorder + slen4 ;

  return ;
}


static void LeftDrag (double x, double y) 
{
    boundaries(&x, &y);
    setQposSpos((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
    drawCrosshair(x, y);
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
}


static void LeftDown (double x, double y) 
{ 
    if (!CrosshairON) return;
    
    boundaries(&x, &y);
    setQposSpos((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);
    drawCrosshair(x, y);
    drawAlignment((int)(x-LeftBorder)*zoom, (int)(y-TopBorder)*zoom);

    graphRegister (LEFT_DRAG, LeftDrag);
}


/* KEYBOARD handles x, y coords in sequence units! (Mouse does it in screen units)
*/
static void keyboard (int key, int modifier_unused)
{
    int x, y;

    switch (key) {
    case UP_KEY:   x = oldx;    y = oldy -1;   break;
    case DOWN_KEY: x = oldx;	y = oldy+1;    break;
    case LEFT_KEY: x = oldx-1;	y = oldy;      break;
    case RIGHT_KEY: x = oldx+1; y = oldy;      break;
    case '>':	x = oldx+1 ;	y = oldy+1 ;   break ;
    case '<':	x = oldx-1 ;	y = oldy-1 ;   break ;
    case '}':	x = oldx+1 ;	y = oldy-1 ;   break ;
    case '{':	x = oldx-1 ;	y = oldy+1 ;   break ;
    default: return;
    }

    if (x < 0) x = 0;
    else if (x > qlen/resfac-1) x = qlen/resfac-1;
    if (y < 0) y = 0;
    else if (y > slen-1) y = slen-1;

    graphActivate(dotterGraph);
    setQposSpos(x, y);
    drawCrosshair((float)x/zoom+LeftBorder, (float)y/zoom+TopBorder);
    drawAlignment(x, y);

    oldx = x;
    oldy = y;
}


/* Find an executable and return its complete pathname.
 */
int findCommand (char *command, char **retp)
{
#if !defined(NO_POPEN)
    static char retstr[1025] ;
    char *path, file[1025], retval;
    int found=0;

    /* Don't use csh - fails if the path is not set in .cshrc * /
    if (access(csh, X_OK)) {
	g_critical("Could not find %s", csh);
	return 0;
    }
    if (!(pipe = (FILE *)popen(messprintf("%s -cf \"which %s\"", csh, command), "r"))) {
	return 0;
    }

    while (!feof(pipe))
	fgets(retval, 1024, pipe);
    retval[1024] = 0;
    pclose(pipe);

    if (cp = strchr(retval, '\n')) *cp = 0;
    if (retp) *retp = retval;

    / * Check if whatever "which" returned is an existing and executable file * /
    if (!access(retval, F_OK) && !access(retval, X_OK))
        return 1;
    else
	return 0;
    */
    
    path = g_malloc(strlen(getenv("PATH"))+1);
    /* Don't free 'path' since it changes later on - never mind, 
       we're only calling it once */

    strcpy(path, getenv("PATH"));
    path = strtok(path, ":");
    while (path) {
	strcpy(file, path);
	strcat(file,"/");
	strcat(file, command);
	if (!access(file, F_OK) && !access(file, X_OK)) {
	    found = 1;
	    break;
	}

	path = strtok(0, ":");
    }

    if (found) {
	strcpy(retstr, file);
	retval = 1;
    }
    else {
	strcpy(retstr, "Can't find executable 'dotter' in path");
	retval = 0;
    }

    if (retp) *retp = retstr;
    return retval;

#endif
}

static void stringProtect(FILE *file, const char *string)
{
  const char *cp;
 
  fputc(' ', file);
  fputc('"', file);
  if (string)
    for(cp = string; *cp; cp++)
      {
	if (*cp == '"' || *cp == '$')
	  fputc('$', file);
	fputc(*cp, file);
      }
  fputc('"', file);
  
}

static void callDotter(int dotterZoom, int xstart, int ystart, int xend, int yend)
{
#if !defined(NO_POPEN)
    FILE *pipe;
    static char *noname = "NoName", *qnam, *snam ;
    MSP *msp;

    /* Need a name for correct number of arguments */
    qnam = (*qname ? qname : noname);
    snam = (*sname ? sname : noname);

    if (xstart < 1)  xstart = 1;
    if (xend > qlen) xend = qlen;
    if (ystart < 1)  ystart = 1;
    if (yend > slen) yend = slen;


    /* Open pipe to new dotterBinary */
    if (!dotterBinary) { 
	printf("Looking for Dotter ...\n");
	if (!findCommand("dotter", &(dotterBinary))) {
	    g_critical("Failed to zoom in - %s.  "
		    "($PATH=%s)", dotterBinary, getenv("PATH"));
	    dotterBinary = 0;
	    return;
	}
    }
    printf("Calling %s with region: %d,%d - %d,%d\n", 
	   dotterBinary, xstart, ystart, xend, yend);
    fflush(stdout);

    pipe = (FILE *)popen(messprintf("/bin/csh -cf \"%s -z %d -q %d -s %d -S '%s' %d '%s' %d %s %s\"", 
				    dotterBinary, 
				    dotterZoom, 
				    xstart-1+qoffset, 
				    ystart-1+soffset, 
				    qnam, 
				    xend-xstart+1, 
				    snam, 
				    yend-ystart+1,
				    dotterBinary,
				    (Xoptions ? Xoptions : "")), 
			 "w");

    fwrite(qseq+xstart-1, 1, xend-xstart+1, pipe);
    fwrite(sseq+ystart-1, 1, yend-ystart+1, pipe);

    /* Pass on features */
    for (msp = MSPlist; msp; msp = msp->next) 
      {
	if (msp->type == BLXMSP_FS_SEG)
	  {
	    fprintf(pipe, "%d %f %d %d %d %d", 
		    msp->type,
		    msp->score, 
		    msp->fsColor, 
		    mspGetQStart(msp) +MSPoffset,
		    mspGetQEnd(msp)   +MSPoffset,
		    msp->fs ? msp->fs->order : 0);
	    stringProtect(pipe, mspGetSName(msp));
	    stringProtect(pipe, msp->sframe);
	    stringProtect(pipe, msp->qname);
	    stringProtect(pipe, msp->qframe);
	    stringProtect(pipe, msp->desc);
	    fputc('\n', pipe);
	  }
      }

    fprintf(pipe, "%c\n", EOF);
    fflush(pipe);
#endif
}


static void callDotterParams(void)
{
  static int dotterZoom=0, xstart=0, ystart=0, xend=0, yend=0;
  ACEIN params_in;

  params_in = messPrompt ("Dotter parameters: zoom (compression) factor, "
			  "xstart, ystart, xend, yend", 
			  messprintf("%d %d %d %d %d", 
				     dotterZoom,
				     xstart, ystart, xend, yend),
			  "iiiiiz", 0);
  if (!params_in)
    return;

  aceInInt(params_in, &dotterZoom);
  aceInInt(params_in, &xstart);
  aceInInt(params_in, &ystart);
  aceInInt(params_in, &xend);
  aceInInt(params_in, &yend);
  aceInDestroy (params_in);

  callDotter(dotterZoom, xstart, ystart, xend, yend);

  return;
} /* callDotterParams */


static void MiddleUp (double x, double y) 
{
  int  xstart, ystart, xend, yend, t ;

  boundaries(&x, &y);

  if (oldrectx)
    {
      graphXorLine(rectx, recty, rectx, oldrecty);
      graphXorLine(rectx, recty, oldrectx, recty);
      graphXorLine(oldrectx, recty, oldrectx, oldrecty);
      graphXorLine(rectx, oldrecty, oldrectx, oldrecty);
    }

  if (x < rectx)
    {
      t = rectx;
      rectx = x;
      x = t;
    }

  if (y < recty)
    {
      t = recty;
      recty = y;
      y = t;
    }

  xstart = (rectx-LeftBorder)*zoom*resfac;
  ystart = (recty-TopBorder)*zoom;

  xend = (x-LeftBorder)*zoom*resfac;
  yend = (y-TopBorder)*zoom;

  /* Ignore small mouse moves as they are likely to be accidental or cancelled clicks... */
  if (xend-xstart > 10 && yend-ystart > 10)
    {
      callDotter(0, xstart, ystart, xend, yend) ;
    }

  return ;
}


static void MiddleDrag (double x, double y) 
{
    boundaries(&x, &y);

    if (oldrectx) {
	graphXorLine(rectx, recty, rectx, oldrecty);
	graphXorLine(rectx, recty, oldrectx, recty);
	graphXorLine(oldrectx, recty, oldrectx, oldrecty);
	graphXorLine(rectx, oldrecty, oldrectx, oldrecty);
    }
    
    graphXorLine(rectx, recty, rectx, y);
    graphXorLine(rectx, recty, x, recty);
    graphXorLine(x, recty, x, y);
    graphXorLine(rectx, y, x, y);

    oldrectx = x;
    oldrecty = y;
}

static void MiddleDown (double x, double y) 
{ 
    boundaries(&x, &y);

    rectx = x;
    recty = y;
    oldrectx = oldrecty = 0;

    graphRegister (MIDDLE_DRAG, MiddleDrag);
    graphRegister (MIDDLE_UP, MiddleUp);
}


static void initWindow(char *winsize)
{
    double 
	exp1, exp2, exp3;
    int
	win1, win2, win3;

    /* Call winsizeFromlambdak even if we don't want to set the window
       size in order to get the other parameters (exp_res_score) */

    if (blastx) {
	win1 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[0], sseq, &exp1, &Lambda);
	win2 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[1], sseq, &exp2, &Lambda);
	win3 = winsizeFromlambdak(MATRIX, atob_0, abetsize, pepqseq[2], sseq, &exp3, &Lambda);
	exp_res_score = (exp1 + exp2 + exp3)/3.0;
	win = (win1 + win2 + win3)/3.0;
    }
    else if (blastn)
	win = winsizeFromlambdak(MATRIX, ntob, abetsize, qseq, sseq, &exp_res_score, &Lambda);
    else
	win = winsizeFromlambdak(MATRIX, atob_0, abetsize, qseq, sseq, &exp_res_score, &Lambda);
    


    if (!winsize || freeupper(*winsize) == 'K') {
	if (win < 3) {
	    g_critical("Karlin/Altschul estimate of window size = %d ignored. Using 10 instead.\n", win);
	    win = 10;
	}
	if (win > 50) {
	    g_critical("Karlin/Altschul estimate of window size = %d ignored. Using 50 instead.\n", win);
	    win = 50;
	}
	return;
    }

    if (!atoi(winsize))
	g_error("Bad window size specification: %s", winsize);
    win = atoi(winsize);
}


static void setWindow(void)
{
  ACEIN size_in;

  if (!(size_in = messPrompt ("Give window size", 
			      messprintf("%d", win), 
			      "iz", 0)))
    return;
  
  aceInInt(size_in, &win) ;
  aceInDestroy (size_in);

  calcWindow();
  dotterRedraw();
    
    /*    
       sprintf(banner, "%s (horizontal) vs %s (vertical).  Window = %d, Pixel values = %d x score/residue, Matrix = %s", 
       qname, sname, win, pixelFac, MATRIX_NAME);
       graphBoxDraw(bannerbox, BLACK, WHITE);
       */  
}


static void dotterDestroy(void)
{
  int i;
  /* Free Dotter stuff */
  handleDestroy(handle);
  
  if (blastx)
    {
      for (i = 0; i < 3; i++)
	{
	  g_free(pepqseq[i]);
	}
    }

  /* Free stuff g_malloc'ed in calling routine (usually blixem or dotterMain) */
  g_free(qseq);
  g_free(sseq);

  /* Don't free MSP's since that will screw blixem up !!! */
  
  if (graphActivate(dotterGraph))
    graphClear();

  g_free(data);
  g_free(HSPpixels);
  
  if (graphActivate(alnGraph))
    {
      graphDestroy() ;
      alignmentInitialized = 0 ;
    }

  if (graphActivate(fsGraph))
    graphDestroy();
  /* if (graphActivate(rampGraph)) graphDestroy();*/

  return ;
}


static void initCrosshair(void)
{
  if (!CrosshairON)
    {
      menuSetFlags(menuItem(dotterMenu, toggleCrosshairPosStr), MENUFLAG_HIDE);
      menuSetFlags(menuItem(dotterMenu, toggleCrosshairFullscreenStr), MENUFLAG_HIDE);
    }
  else
    {
      menuUnsetFlags(menuItem(dotterMenu, toggleCrosshairPosStr), MENUFLAG_HIDE);
      menuUnsetFlags(menuItem(dotterMenu, toggleCrosshairFullscreenStr), MENUFLAG_HIDE);

      /* Set up crosshair boxes ***********************/
      graphColor(BLUE);
	
      if (!CrosshairFullscreenON)
	{
	  hLineBox = graphBoxStart();
	  graphLine(LeftBorder, crossy, LeftBorder-2+qlen4, crossy);
	  graphBoxEnd();
	
	  vLineBox = graphBoxStart();
	  graphLine(crossx, TopBorder, crossx, TopBorder-2+slen4);
	  graphBoxEnd();
	}
      else
	{
	  int nx, ny;

	  graphGetBounds(&nx, &ny) ;

	  hLineBox = graphBoxStart();
	  graphLine(0.0, crossy, (float)nx+1000, crossy);
	  graphBoxEnd();
      
	  vLineBox = graphBoxStart();
	  graphLine(crossx, 0.0, crossx, (float)ny+1000);
	  graphBoxEnd();
	}

      if (CrosshairPosON)
	{
	  setQposSpos((int)(crossx-LeftBorder)*zoom, (int)(crossy-TopBorder)*zoom);
	  sprintf(CrosshairPosText, "%s, %s", qpos, spos);
	  CrosshairPosBox = graphBoxStart();
	  graphTextPtr(CrosshairPosText, crossx+10, crossy+10, 25);
	  graphBoxEnd();
	  graphBoxDraw(CrosshairPosBox, BLUE, TRANSPARENT);
	}
	
      graphColor(BLACK);
    }

  return ;
}


static void togglePixelmap(void)
{
    if (!PixelmapON) 
    {
	if (!pixelmap_done) {
	    calcWindow();
	}
	PixelmapON = 1;
    }
    else 
	PixelmapON = 0;

    dotterRedraw();
}


static void toggleGrid(void)
{
    gridOn = (gridOn ? 0 : 1);
    dotterRedraw();
}


static void initGreyramp(void)
{
/*
    float 
	minScore,		/ * score/residue * /
	maxScore;		/ * score/residue * /

    float
	offset;			/ * Disused idea for compensating for raised 
				   noise levels in compressed images * /
*/
    int 
	min = 40, 
	max = 100;

    
    if (greyRampSwap) {

      gexRampSet(min, max);
      gexRampTool();

    }
    else {
      gexRampSet(max, min);
      gexRampTool();

    }

    return;


    /* The pixelFac is set so that a pixel value of 50 = exp_res_score.
       Therefore, the stuff below is obsolete.  * /
    
    minScore = 0.8 * exp_res_score;
    maxScore = 2.0 * exp_res_score;

    offset = log((double)zoom*zoom)/(Lambda*win);

    printf("exp_res_score=%.2f => minScore=%.2f, maxScore=%.2f; offset=%.2f\n", 
	   exp_res_score, minScore, maxScore, offset);

    min = (maxScore  +offset)*(float)pixelFac;
    max = (minScore  +offset)*(float)pixelFac;

    if (min < 0) min = 0;
    if (min > 255) min = 255;

    if (max < 0) max = 0;
    if (max > 255) max = 255;

    graphGreyRamp(max, min);*/
}


static void dotterRedraw(void)
{
  static int bannerPos;

  if (!graphActivate(dotterGraph))
    {
      float w, h, fw, fh;
      int pw, ph;
      BOOL oldMap = graphSetInstallMap(TRUE);


      /* (rbrusk): should initialize these for all new dotterGraphs...*/
      vLineBox = hLineBox = 0 ;

      graphScreenSize(&w, &h, &fw, &fh, &pw, &ph);

      fonth = ph/fh;
      fontw = pw/fw;

      dotterGraph = graphCreate (PIXEL_SCROLL, 
				 messprintf("Dotter %s vs. %s", qname, sname), 
				 0, 0, 
				 (LeftBorder+qlen4+60 + fsTotalHeight(MSPlist)*fonth + (MSPlist ? 40:0))*w/pw, 
				 (TopBorder +slen4+60 + fsTotalHeight(MSPlist)*fonth + (MSPlist ? 40:0))*h/ph) ;

      graphNewMenu (dotterMenu);

      /* The menu toggle stuff is all wrong - there's been a misunderstanding
	 I've done the minium to make it work, but it needs to be fixed 
	 properly - srk */
      /* Why must setMenuCheckmarks() be here? - esr */
      setMenuCheckmarks();
      graphRegister(MIDDLE_DOWN, MiddleDown) ;
      graphRegister(KEYBOARD, keyboard);
      graphRegister(DESTROY, dotterDestroy ); 
      graphRegister(LEFT_DOWN, LeftDown);
      graphRegister(RAMP_CHANGE, dotterRampChange);
      initGreyramp();

      graphActivate(dotterGraph);
      graphSetInstallMap(oldMap);

    }

  graphClear();
  graphPop();

  if (fsArr)
    {
      menuUnsetFlags(menuItem(dotterMenu, selectFeaturesStr), MENUFLAG_HIDE);
    }

    
  /*    sprintf(banner, "%s (horizontal) vs. %s (vertical).  Window = %d, Pixel values = %d x score/residue, Matrix = %s",
	qname, sname, win, pixelFac, MATRIX_NAME);
	graphTextPtr(banner, 5, 5, strlen(banner)+5);
  */
    
  sprintf(banner, "%s (horizontal) vs. %s (vertical)", qname, sname);
  bannerPos = (strlen(banner)*fontw < qlen ? LeftBorder : 10);
  graphText(banner, bannerPos, 5);
  if (!printMode)
    graphButton("About", Help, bannerPos, 20);
    
  graphPixelBounds(LeftBorder+qlen4+60 + fsTotalHeight(MSPlist)*fonth + (MSPlist ? 40:0), 
		   TopBorder +slen4+60 + fsTotalHeight(MSPlist)*fonth + (MSPlist ? 40:0));

  if (PixelmapON && !(BlastHSPsOn && BlastMode == BLASTGREY))
    {
      graphPixelsRaw ((unsigned char *)data, qlen4, slen4, qlen4, LeftBorder, TopBorder);
    }
	
  if (BlastHSPsOn)
    {
      static int gapwarned = 0;

    if (BlastMode == BLASTGREY)
      {
	if (!HSPpixels)
	  { 
	    int i;
	    HSPpixels = (UCHAR *)g_malloc(datalen);
	    for (i=0; i < datalen; i++)
	      HSPpixels[i] = 0;
	  }
	graphPixelsRaw ((unsigned char *)HSPpixels, qlen4, slen4, qlen4, LeftBorder, TopBorder);
      }

    drawBlastHSPs();
    
    if (HSPgaps && !gapwarned)
      {
	graphRedraw();
	g_critical("Note: gapped HSPs are shown ungapped in Dotter.");
	gapwarned = 1;
      }
    }

  initCrosshair();
  drawScale();
  drawAllFeatures(MSPlist);
  
  graphRedraw();

  if (graphActivate(alnGraph)) 
    graphPop();
  else
    {
      if (!alignmentInitialized && !getenv("ACEDB_PROJECT"))
	{
	  initAlignment() ;
	  alignmentInitialized = 1 ;
	}
    }

  gexRampPop();

  if (graphActivate(fsGraph)) 
    graphPop() ;

  graphActivate(dotterGraph);

  return ;
}


static void readmtx(int MATRIX[24][24], char *mtxfile)
{
    FILE *fil;
    int row, col;
    char line[1025] = "#", *p;
    
    if (!(fil = fopen(mtxfile, "r")) &&
	!(fil = fopen(messprintf("%s/%s", getenv("BLASTMAT"), mtxfile), "r")))
	g_error("Failed to open score matrix file %s - not found in ./ or $BLASTMAT/.", mtxfile);
    
    /* Ignore header ... */
    while (!feof(fil) && *line == '#') 
	fgets(line, 1024, fil);

    /* Read in the pairwise scores */
    for (row = 0; row < 24; row++)
    {
	if (!(fgets(line, 1024, fil)))
	    g_error("Wrong number of rows in matrix file: %d (should be 24).", row);

	p = strtok(line, " \t\n");
	for (col = 0; col < 24; col++) 
	{
	    while (*p == '*' || isalpha((int) *p))
		p = strtok(NULL, " \t\n");
          
	    if (!p) 
              g_error("Error on row %d in matrix file.", row);

	    MATRIX[row][col] = atoi(p);

	    p = strtok(NULL, " \t\n");
	}
    }

    printf("I read your score matrix %s.\n", mtxfile);
    fclose(fil);
}

static void mtxcpy(int dest[24][24], int src[24][24])
{
    int i, j;

    for (i = 0 ; i < 24 ; i++)
	for (j = 0 ; j < 24 ; j++)
	    dest[i][j] = src[i][j];
}


static void DNAmatrix(int mtx[24][24])
{
    int i, j;

    for (i = 0 ; i < 6 ; i++)
	for (j = 0 ; j < 6 ; j++) {
	    if ( i < 4 && j < 4) 
		mtx[i][j] = (i == j ? 5 : -4);
	    else 
		mtx[i][j] = -4;
	}
}


/* Use this routine to add -install in programs that use Dotter */
void argvAdd(int *argc, char ***argv, char *s)
{
    char **v;
    int i;

    v = (char **)malloc(sizeof(char *)*(*argc+2));

    /* Copy argv */
    for (i = 0; i < (*argc); i++)
	v[i] = (*argv)[i];

    /* Add s */
    v[*argc] = (char *)malloc(strlen(s)+1);
    strcpy(v[*argc], s);

    (*argc)++;

    /* Terminator - ANSI C standard */
    v[*argc] = 0;

    /* free(*argv); */   /* Too risky */
    
    *argv = v;
}


/* Utility to get the length of the given GArray. Returns 0 if array is null. */
static int gArrayGetLen(GArray *array)
{
  return (array ? array->len : 0);
}

/**************************** eof ***************************/
