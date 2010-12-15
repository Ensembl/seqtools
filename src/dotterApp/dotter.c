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
 * CVS info:   $Id: dotter.c,v 1.21 2010-11-08 15:52:49 gb10 Exp $
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

#include <seqtoolsUtils/version.h>
#include <seqtoolsUtils/utilities.h>
#include <dotterApp/dotter_.h>
#include <dotterApp/seqtoolsExonView.h>
#include <gdk/gdkkeysyms.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* tint stuff used to be in graph.h, now local - rd 960524
   NB colours as Erik liked them before Jean tinkered!
   could rename #define's more sensibly now
*/

//#define TINT_WHITE      0x00
//#define TINT_HIGHLIGHT1 0x01	/* highest priority, dna highlighting */
//#define TINT_HIGHLIGHT2 0x02	/* highlight friends */
//#define TINT_RED        0x04 
//#define TINT_LIGHTGRAY  0x08 
//#define TINT_MAGENTA    0x10 
//#define TINT_CYAN       0x20 
//#define TINT_LIGHTGREEN 0x40 
//#define TINT_YELLOW     0x80 

//static int tints[8] = { LIGHTRED, MIDBLUE, RED, LIGHTGRAY, 
//			MAGENTA, CYAN, LIGHTGREEN, YELLOW } ;

#define MAX_WINDOW_WIDTH_FRACTION             0.7 /* max init width of dotter window as fraction of screen size */
#define MAX_WINDOW_HEIGHT_FRACTION            0.7 /* max init height of dotter window as fraction of screen size */


//typedef struct
//{
//  const char *name ;
//  int start, end ;
//  char strand ;
//  MSP *msp_start, *msp_end ;
//  float y_pos ;
//} GeneDataStruct, *GeneData ;
//
//typedef struct
//{
//  GList *forward_genes ;
//  GList *reverse_genes ;
//  char strand  ;
//} GeneStrandStruct, *GeneStrand ;
//


/* Struct containing properties for a dotter window */
typedef struct _DotterProperties
{
  GtkWidget *greyrampTool;                  /* the greyramp tool */
  GtkWidget *alignmentTool;		    /* the alignment tool */
  GtkWidget *dotplot;                       /* the dotplot drawing area */
  
  DotterWindowContext *dotterWinCtx;
} DotterProperties;


#define DOTTER_HELP_TEXT "\
<b><big>DOTTER</big></b>\n\
A dot-matrix program with dynamic threshold control suited for genomic DNA and protein sequence analysis\n\
\n\
\n\
<span foreground=\"blue\">\
<b><big>What's new</big></b>\n\
\t*\t<b><i>Dotter re-vamp</i></b>: Dotter has been re-written to facilitate future improvements.  Most of the changes so far are under-the-hood, but you will notice a few cosmetic differences and additional shortcuts.\n\
\t*\t<b><i>Menu bar</i></b>: As well as the right-click menu, there is now a menu-bar at the top of the main Dotter window.\n\
\t*\t<b><i>CDS/UTR regions</i></b>: Exons are now separated into CDS and UTR regions: CDS regions are coloured green and UTR red.\n\
\t*\t<b><i>Close all sub-Dotters</i></b>: You can close all related Dotter windows using the Quit menu option or Ctrl-Q.  This will close all sub-Dotters that were started under the same parent Dotter (by middle-dragging to zoom in to a region).  To just close an individual Dotter window, click on the x in the corner of the window or use your system shortcut for closing a window (e.g. Ctrl-W, or Cmd-W on Macs).  Note that the alignment tool and greyramp tool will be destroyed along with their parent - however, if the parent window is still open then these tools can be re-opened using the relevant menu options or keyboard shortcuts.\n\
\t*\t<b><i>Close all Dotters from Blixem</i></b>: All Dotter windows spawned from a Blixem will be closed when that Blixem is closed.\n\
\t*\t<b><i>Keyboard shortcuts</i></b>: The following keyboard shortcuts have been added to show the Alignment tool, main Dotter window or Greyramp tool, respectively: Ctrl-A, Ctrl-D and Ctrl-G.\n\
\t*\t<b><i>Settings dialog</i></b>: The Parameter Control dialog box has been replaced by a more intuitive Settings dialog.  From here you can change the zoom or edit the display range.\n\
</span>\
\n\
\n\
<b><big>Mouse controls</big></b>\n\
 - Left button: position crosshair.\n\
 - Middle button: drag to zoom in to a region.\n\
\n\
\n\
<b><big>Keyboard shortcuts</big></b>\n\
 - Arrow keys: move crosshair one dot in arrow direction\n\
 - &lt; &gt; : move crosshair along diagonals\n\
 - { } : move along reverse diagonals\n\
 - Ctrl-Q : quit Dotter (including any child/parent Dotters)\n\
 - Ctrl-H : show this Help page\n\
 - Ctrl-A : show the alignment tool\n\
 - Ctrl-D : show the Dotter main window\n\
 - Ctrl-G : show the greyramp tool\n\
\n\
Hold down Shift to move by nucleotides rather than whole peptides.\n\
\n\
\n\
<b><big>Settings</big></b>\n\
 - Zoom: enter a higher value to zoom out. A value of 1 means 100%%, 2 means 50%% etc. A fraction of 1 can be entered in order to zoom in (e.g. 0.5 for a 200%% zoom), but the display will appear stretched.\n\
 - Horizontal range: enter the min and max coords to display on the horizontal scale. Note that this will be limited to the horizontal sequence range that Dotter was started up with.\n\
 - Vertical range: enter the min and max coords to display on the vertical scale. Note that this will be limited to the vertical sequence range that Dotter was started up with.\n\
 - Sliding window size: affects cut-off limit for how dots are drawn\n\
\n\
\n\
<b><big>Residue colours (alignment tool)</big></b>\n\
Cyan      = Identical Residue.\n\
DarkBlue  = Positive Score.\n\
No colour = Negative score.\n\
\n\
\n\
<b><big>Session details</big></b>\n\
Sliding window length = %d\n\
Pixel values = %d x score/residue\n\
Matrix = %s\n\
Zoom (compression) factor = %f\n\
"



static void showSettingsDialog(GtkWidget *dotterWindow);
static void readmtx(int mtx[24][24], char *mtxfile);
static void mtxcpy(int mtx[24][24], int BLOSUM62[24][24]);
static void DNAmatrix(int mtx[24][24]);
//static void drawAllFeatures(MSP *msp) ;
//static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth) ;
//gint compareMSPs(gconstpointer a, gconstpointer b) ;
//static void getGenePositionsCB(gpointer data, gpointer user_data) ;
//gint compareGenes(gconstpointer a, gconstpointer b) ;
//static void setYoffsets(GList *first, float min_offset, float max_offset) ;
//static void drawGenesCB(gpointer data, gpointer user_data) ;
//static void drawMSPGene(MSP *msp, float y_offset) ;
//static int gArrayGetLen(GArray *array);

static void		      showGreyrampTool(GtkWidget *dotterWindow);
static void		      showAlignmentTool(GtkWidget *dotterWindow);
static GtkWidget*	      createDotterWindow(DotterContext *dc, DotterWindowContext *dwc, const DotterHspMode hspMode, GtkWidget *greyrampTool, GtkWidget *dotplotContainer, GtkWidget *dotplot);
static DotterContext*         dotterGetContext(GtkWidget *dotterWindow);
static void                   redrawAll(GtkWidget *dotterWindow, gpointer data);
static void                   refreshAll(GtkWidget *dotterWindow, gpointer data);
static gboolean               onKeyPressDotter(GtkWidget *widget, GdkEventKey *event, gpointer data);
static gboolean               onKeyPressDotterCoords(GtkWidget *widget, GdkEventKey *event, gpointer data);

static void createDotterInstance(DotterContext *dotterCtx,
                                 DotterWindowContext *dotterWinCtx,
                                 const char *loadFileName,
                                 const char *saveFileName,
                                 const gboolean hspsOn,
                                 const char* winsizeIn,
                                 const int pixelFacIn,
                                 const int zoomFacIn,
                                 const int qcenter,
                                 const int scenter,
                                 const gboolean greyramSwap);

static void                       onQuitMenu(GtkAction *action, gpointer data);
static void                       onSavePlotMenu(GtkAction *action, gpointer data);
static void                       onPrintMenu(GtkAction *action, gpointer data);
static void                       onSettingsMenu(GtkAction *action, gpointer data);
static void                       onShowGreyrampMenu(GtkAction *action, gpointer data);
static void                       onShowAlignmentMenu(GtkAction *action, gpointer data);
static void                       onToggleCrosshairMenu(GtkAction *action, gpointer data);
static void                       onToggleCoordsMenu(GtkAction *action, gpointer data);
static void                       onToggleFullscreenMenu(GtkAction *action, gpointer data);
static void                       onTogglePixelmapMenu(GtkAction *action, gpointer data);
static void                       onToggleGridMenu(GtkAction *action, gpointer data);
static void                       onHelpMenu(GtkAction *action, gpointer data);
static void                       onAboutMenu(GtkAction *action, gpointer data);

/* Menu builders: the action entry list lists menu actions for all menus */
static const GtkActionEntry menuEntries[] = {
{ "Quit",             NULL, "_Quit\t\t\tCtrl-Q",                 NULL,  "Quit dotter",                G_CALLBACK(onQuitMenu)},
{ "SavePlot",         NULL, "_Save plot\t\tCtrl-S",            NULL,  "Save plot",                  G_CALLBACK(onSavePlotMenu)},
{ "Print",            NULL, "_Print\t\t\tCtrl-P",                NULL,  "Print",                      G_CALLBACK(onPrintMenu)},
{ "Settings",         NULL, "Settings\t\t\tCtrl-S",              NULL,  "Set dotter parameters",      G_CALLBACK(onSettingsMenu)},
{ "ShowGreyramp",     NULL, "_Greyramp tool\tCtrl-G",        NULL,  "Show the greyramp tool",     G_CALLBACK(onShowGreyrampMenu)},
{ "ShowAlignment",    NULL, "_Alignment tool\tCtrl-A",       NULL,  "Show the alignment tool",    G_CALLBACK(onShowAlignmentMenu)},
{ "Help",             NULL, "_Help\t\t\tCtrl-H",	NULL,  "Dotter Help",                G_CALLBACK(onHelpMenu)},
{ "About",            NULL, "_About",			NULL,  "About Dotter",               G_CALLBACK(onAboutMenu)}
};

/* Toggle-able menu entries are listed here: */
static GtkToggleActionEntry toggleMenuEntries[] = {
{ "TogglePixmap",     NULL, "Pixelmap",              NULL,  "Show the pixelmap",              G_CALLBACK(onTogglePixelmapMenu),   TRUE},
{ "ToggleGrid",       NULL, "Gridlines",             NULL,  "Show grid lines",                G_CALLBACK(onToggleGridMenu),       FALSE},
{ "ToggleCrosshair",  NULL, "Crosshair",             NULL,  "Show the crosshair",             G_CALLBACK(onToggleCrosshairMenu),  TRUE},
{ "ToggleCoords",     NULL, "Crosshair label",       NULL,  "Show the crosshair label",       G_CALLBACK(onToggleCoordsMenu),     TRUE},
{ "ToggleFullscreen", NULL, "Crosshair fullscreen",  NULL,  "Show the crosshair full screen", G_CALLBACK(onToggleFullscreenMenu), TRUE}
};

/* Radio-button menu entries are listed here: */
static const GtkRadioActionEntry radioMenuEntries[] = {
{ "HspsOff",    NULL, "HSPs off",                     NULL,  "Hide Blast HSPs",                                    DOTTER_HSPS_OFF},
{ "HspsGrey",   NULL, "Draw HSPs (greyramp)",         NULL,  "Draw Blast HSPs as greyramp",                        DOTTER_HSPS_GREYSCALE},
{ "HspsLine",   NULL, "Draw HSPs (red lines)",        NULL,  "Draw Blast HSPs as solid red lines",                 DOTTER_HSPS_LINE},
{ "HspsFunc",   NULL, "Draw HSPs (color = f(score))", NULL,  "Draw Blast HSPs with color as a function of score",  DOTTER_HSPS_FUNC}
};


/* Menu descriptions - these define the layouts of the individual menus */
static const char fileMenuDescription[] =
"<ui>"
"  <popup name='File'>"
"      <menuitem action='SavePlot'/>"
"      <separator/>"
"      <menuitem action='Print'/>"
"      <separator/>"
"      <menuitem action='Quit'/>"
"  </popup>"
"</ui>";

static const char editMenuDescription[] =
"<ui>"
"  <popup name='Edit'>"
"      <menuitem action='Settings'/>"
"  </popup>"
"</ui>";

static const char viewMenuDescription[] =
"<ui>"
"  <popup name='View'>"
"      <menuitem action='ShowGreyramp'/>"
"      <menuitem action='ShowAlignment'/>"
"      <separator/>"
"      <menuitem action='ToggleCrosshair'/>"
"      <menuitem action='ToggleCoords'/>"
"      <menuitem action='ToggleFullscreen'/>"
"      <separator/>"
"      <menuitem action='TogglePixmap'/>"
"      <menuitem action='ToggleGrid'/>"
"      <separator/>"
"      <menuitem action='HspsOff'/>"
"      <menuitem action='HspsGrey'/>"
"      <menuitem action='HspsLine'/>"
"      <menuitem action='HspsFunc'/>"
"  </popup>"
"</ui>";

static const char helpMenuDescription[] =
"<ui>"
"  <popup name='Help'>"
"      <menuitem action='Help'/>"
"      <menuitem action='About'/>"
"  </popup>"
"</ui>";

/* The right-click context menu offers all the above on a single menu (to maintain historic behaviour) */
static const char contextMenuDescription[] =
"<ui>"
"  <popup name='Context'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
"      <menuitem action='SavePlot'/>"
"      <menuitem action='Print'/>"
"      <separator/>"
"      <menuitem action='Settings'/>"
"      <separator/>"
"      <menuitem action='ShowGreyramp'/>"
"      <menuitem action='ShowAlignment'/>"
"      <separator/>"
"      <menuitem action='ToggleCrosshair'/>"
"      <menuitem action='ToggleCoords'/>"
"      <menuitem action='ToggleFullscreen'/>"
"      <separator/>"
"      <menuitem action='TogglePixmap'/>"
"      <menuitem action='ToggleGrid'/>"
"      <separator/>"
"      <menuitem action='HspsOff'/>"
"      <menuitem action='HspsGrey'/>"
"      <menuitem action='HspsLine'/>"
"      <menuitem action='HspsFunc'/>"
"  </popup>"
"</ui>";



/* Global variables */
char *dotterVersion = DOTTER_VERSION_COMPILE ;

static int    MATRIX[24][24];
//              MSPoffset,	/* Difference between real MSP coord and coord stored in MSP */
//              HSPgaps = 0.
//              fsBoxStart,
//              fsRightOn = 1,
//              fsBottomOn = 1,
//              fsAnnRightOn = 1,
//              fsAnnBottomOn = 1,
//              fsEndLinesOn = 0;

//static Graph  fsGraph=0;

#define MAXALIGNLEN 501

//static char fsPlotHeighttx[10];

       float fsPlotHeight=2.0;	/* The height of feature series XY plots */
static MSP   *MSPlist=0;     /* List of MSPs - the first object contains data */
//static MSP   *msp;

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


//GArray *fsArr = NULL;  /* Stores Feature Series - t he actual segments are stored
//			   as MSPs, using these fields:
//			   msp->sframe  = [1..2] sequence
//			   msp->qstart = segment start
//			   msp->qend   = segment end
//			   msp->fs     = Ordinal number of series that this MSP belongs to.
//			   msp->fsColor  = color
//			   msp->desc   = annotation
//			   */


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

/* Create the colors that Dotter will use for various specific purposes */
static void createDotterColors(DotterContext *dc)
{
  /* Initialise the array with empty BlxColor structs */
  dc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), DOTCOLOR_NUM_COLORS);
  int i = DOTCOLOR_MIN + 1;
  
  for ( ; i < DOTCOLOR_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(dc->defaultColors, *blxColor);
    }
    
  /* matches */
  createBlxColor(dc->defaultColors, DOTCOLOR_MATCH, "Exact match", "Exact match", "#00ffe5", BLX_LIGHT_GREY, "#00c3b0", NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CONS, "Conserved match", "Conserved match", "#78b4f0", BLX_LIGHT_GREY, "#5c98d5", NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_MISMATCH, "Mismatch", "Mismatch", "#cacaca", BLX_WHITE, "#989898", NULL);
  
  /* exons */
  createBlxColor(dc->defaultColors, DOTCOLOR_EXON_FILL, "Exon fill color", "Exon fill color", BLX_YELLOW, BLX_GREY, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_EXON_LINE, "Exon line color", "Exon outline color", BLX_BLUE, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CDS_FILL, "CDS fill color", "Coding section fill color", BLX_PALE_GREEN, BLX_GREY, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CDS_LINE, "CDS line color", "Coding section outline color", BLX_DARK_GREEN, BLX_GREY, BLX_VERY_DARK_GREEN, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_UTR_FILL, "UTR fill color", "Untranslated region fill color", BLX_LIGHT_RED, BLX_GREY, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_UTR_LINE, "UTR line color", "Untranslated region outline color", BLX_DARK_RED, BLX_GREY, BLX_VERY_DARK_RED, NULL);

  /* dot plot */
  createBlxColor(dc->defaultColors, DOTCOLOR_CROSSHAIR, "Crosshair", "Color of the crosshair on the dot plot", BLX_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_GRID, "Grid", "Line color of the grid on the dot plot", BLX_LIGHT_RED, BLX_LIGHT_GREY, NULL, NULL);

  /* greyramp */
  createBlxColor(dc->defaultColors, DOTCOLOR_THRESHOLD_MARKER, "Greyramp threshold marker color", "Outline color of the threshold marker on the greyramp tool", BLX_RED, BLX_BLACK, BLX_GREEN, BLX_GREY);
  createBlxColor(dc->defaultColors, DOTCOLOR_MARKER_LINE, "Greyramp marker outline color", "Outline color of the triangle markers on the greyramp tool", BLX_BLACK, BLX_BLACK, BLX_GREEN, BLX_GREY);
  createBlxColor(dc->defaultColors, DOTCOLOR_MARKER_FILL, "Greyramp marker fill color", "Fill color of the triangle markers on the greyramp tool", BLX_WHITE, BLX_WHITE, NULL, NULL);
}


/* Get the initial zoom factor, calculating it if the passed-in zoom is 0 */
static gdouble getInitZoomFactor(DotterContext *dc, const gdouble zoomFacIn, const int qLen, const int sLen)
{
  DEBUG_ENTER("getInitZoomFactor");

  gdouble result = zoomFacIn;
  
  if (result <= 0)
    {
      if (dc->memoryLimit)
        {
          result = (int)sqrt((qLen / dc->numFrames / 1e6 * sLen - 1e-6) / dc->memoryLimit) + 1;
        }
      else
        {
          g_error("Cannot calculate zoom; division by 0 (memory limit is 0).\n");
        }
    }
  
  DEBUG_EXIT("getInitZoomFactor returning %f", result);
  return result;
}


static DotterContext* createDotterContext(BlxBlastMode blastMode, 
                                          const gboolean batchMode,
					  const char *refSeqName,
					  char *refSeq,
                                          const BlxStrand refSeqStrand,
					  const char *matchSeqName,
					  char *matchSeq,
                                          const BlxStrand matchSeqStrand,
					  const int refSeqOffset,
					  const int matchSeqOffset,
					  const gboolean hozScaleRev,
					  const gboolean vertScaleRev,
					  const gboolean selfComp,
					  const gboolean displayMirror,
					  const gboolean watsonOnly,
					  const gboolean crickOnly,
                                          MSP *mspList,
					  GList *seqList,
                                          int matrix[24][24],
                                          char *matrixName,
                                          const double memoryLimit)
{
  DEBUG_ENTER("createDotterContext");

  DotterContext *result = g_malloc(sizeof *result);
  
  result->blastMode = blastMode;
  result->displaySeqType = (blastMode == BLXMODE_BLASTN) ? BLXSEQ_DNA : BLXSEQ_PEPTIDE;
  result->numFrames = (blastMode == BLXMODE_BLASTX) ? 3 : 1;
  result->geneticCode = stdcode1;
  mtxcpy(result->matrix, matrix);
  result->matrixName = matrixName;
  result->mspList = mspList;
  result->seqList = seqList;
  result->windowList = NULL;
  
  result->watsonOnly = watsonOnly;
  result->crickOnly = crickOnly;
  
  /* Set the fixed-width font (not applicable in batch mode) */
  if (!batchMode)
    {
      GtkWidget *tmp = gtk_window_new(GTK_WINDOW_TOPLEVEL);
      const char *fontFamily = findFixedWidthFont(tmp);
      result->fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(result->fontDesc, pango_font_description_get_size(tmp->style->font_desc));
      getFontCharSize(tmp, result->fontDesc, &result->charWidth, &result->charHeight);
      gtk_widget_destroy(tmp);
    }
  else
    {
      result->fontDesc = NULL;
    }
  
  result->refSeqName = g_strdup(refSeqName);
  result->refSeq = refSeq; /* take ownership of passed-in seq */
  result->refSeqRev = NULL;
  result->refSeqType = (blastMode == BLXMODE_BLASTP ? BLXSEQ_PEPTIDE : BLXSEQ_DNA);
  
  /* for dna ref sequences, reverse-complement the ref seq */
  if (result->refSeqType == BLXSEQ_DNA && result->refSeq)
    {
      result->refSeqRev = g_malloc(strlen(result->refSeq) + 1);
      revComplement(result->refSeqRev, result->refSeq);
    }
  else if (result->refSeq)
    {
      /* Just reverse it */
      result->refSeqRev = g_strdup(result->refSeq);
      g_strreverse(result->refSeqRev);
    }
  
  result->refSeqStrand = refSeqStrand;

  result->matchSeqName = g_strdup(matchSeqName);
  result->matchSeq = matchSeq; /* take ownership of passed-in seq */
  result->matchSeqRev = NULL;
  result->matchSeqType = (blastMode == BLXMODE_BLASTN ? BLXSEQ_DNA : BLXSEQ_PEPTIDE);
  result->matchSeqStrand = matchSeqStrand;
  
  /* Reverse/comp match seq if applicable */
  if (result->matchSeqType == BLXSEQ_DNA && result->matchSeqStrand == BLXSTRAND_REVERSE && result->matchSeq)
    {
      result->matchSeqRev = g_malloc(strlen(result->matchSeq) + 1);
      revComplement(result->matchSeqRev, result->matchSeq);
    }
  else if (result->matchSeqStrand == BLXSTRAND_REVERSE && result->matchSeq)
    {
      /* Peptide sequence. Just reverse */
      result->matchSeqRev = g_strdup(result->matchSeq);
      g_strreverse(result->matchSeqRev);
    }
  
  if (result->blastMode == BLXMODE_BLASTX) 
    {
      /* Create the 3 frame translations (for the strand we're interested in only). */
      char *refSeqToUse = (result->refSeqStrand == BLXSTRAND_REVERSE ? result->refSeqRev : result->refSeq);
      
      int i = 0;
      for (i = 0; i < result->numFrames; i++)
        {
          result->peptideSeqs[i] = blxTranslate(refSeqToUse + i, result->geneticCode);
        }
    }
  
  result->refSeqFullRange.min = refSeqOffset + 1;
  result->refSeqFullRange.max = refSeqOffset + strlen(refSeq);
  result->matchSeqFullRange.min = matchSeqOffset + 1;
  result->matchSeqFullRange.max = matchSeqOffset + strlen(matchSeq);
  
  result->hozScaleRev = hozScaleRev;
  result->vertScaleRev = vertScaleRev;

  result->selfComp = selfComp;
  result->displayMirror = displayMirror;
  
  result->memoryLimit = memoryLimit;
  
  result->defaultColors = NULL;
  result->usePrintColors = FALSE;  
  
  if (!batchMode)
    { 
      createDotterColors(result);
    }
  
  /* Calculate the height and width of the horizontal and vertical scales */
  const int leftBorderChars = max(numDigitsInInt(result->matchSeqFullRange.min), numDigitsInInt(result->matchSeqFullRange.max));
  result->scaleWidth = DEFAULT_MAJOR_TICK_HEIGHT + (roundNearest)((gdouble)leftBorderChars * result->charWidth) + SCALE_LINE_WIDTH;
  result->scaleHeight = DEFAULT_MAJOR_TICK_HEIGHT + roundNearest(result->charHeight) + SCALE_LINE_WIDTH;
  
  DEBUG_EXIT("createDotterContext returning");
  return result;
}


static void destroyDotterContext(DotterContext *dc)
{
  DEBUG_ENTER("destroyDotterContext");

  if (dc)
    {
    if (dc->blastMode == BLXMODE_BLASTX)
      {
	int i = 1;
	for ( ; i < dc->numFrames; i++)
	  {
	    g_free(dc->peptideSeqs[i]);
	    dc->peptideSeqs[i] = NULL;
	  }
      }
    
    /* Free stuff g_malloc'ed in calling routine (usually dotterMain) */
    g_free(dc->refSeq);
    dc->refSeq = NULL;
    
    g_free(dc->refSeqName);
    dc->refSeqName = NULL;
    
    g_free(dc->matchSeq);
    dc->matchSeq = NULL;
    
    g_free(dc->matchSeqName);
    dc->matchSeqName = NULL;

    /* To do: free MSP's */
    
      if (dc->matrixName)
        {
          g_free(dc->matrixName);
          dc->matrixName = NULL;
        }
    }  
  
  DEBUG_EXIT("destroyDotterContext returning ");
}

static void destroyDotterWindowContext(DotterWindowContext *dwc)
{
  DEBUG_ENTER("destroyDotterWindowContext");
  /* nothing to do */
  
  DEBUG_EXIT("destroyDotterWindowContext returning ");
}

static DotterProperties* dotterGetProperties(GtkWidget *dotterWindow)
{
  return dotterWindow ? (DotterProperties*)(g_object_get_data(G_OBJECT(dotterWindow), "DotterProperties")) : NULL;
}

static void onDestroyDotterWindow(GtkWidget *dotterWindow)
{
  DEBUG_ENTER("onDestroyDotterWindow");

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties)
    {
      if (properties->greyrampTool)
        {
          gtk_widget_destroy(properties->greyrampTool);
          properties->greyrampTool = NULL;
        }

      if (properties->alignmentTool)
        {
          gtk_widget_destroy(properties->alignmentTool);
          properties->alignmentTool = NULL;
        }
    
      if (properties->dotterWinCtx)
	{
	  if (properties->dotterWinCtx->dotterCtx && properties->dotterWinCtx->dotterCtx->windowList)
	    {
	      properties->dotterWinCtx->dotterCtx->windowList = g_slist_remove(properties->dotterWinCtx->dotterCtx->windowList, dotterWindow);

	      if (g_slist_length(properties->dotterWinCtx->dotterCtx->windowList) < 1)
		{
		  /* It's the last window in the context, so destroy the context */
		  g_slist_free(properties->dotterWinCtx->dotterCtx->windowList);
		  properties->dotterWinCtx->dotterCtx->windowList = NULL;
		  destroyDotterContext(properties->dotterWinCtx->dotterCtx);
		}
	    }
	
	  destroyDotterWindowContext(properties->dotterWinCtx);
	  properties->dotterWinCtx = NULL;
	}

      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(dotterWindow), "DotterProperties", NULL);
    }
  
  DEBUG_EXIT("onDestroyDotterWindow returning ");
}


/* Close all windows associated with the given dotter context */
static void dotterContextCloseAllWindows(DotterContext *dc)
{
  /* Copy the list because pointers are removed from the original
   * list when we destroy windows. */
  GSList *winList = g_slist_copy(dc->windowList);
  GSList *winItem = winList;
  
  for ( ; winItem; winItem = winItem->next)
    {
      GtkWidget *window = GTK_WIDGET(winItem->data);
      gtk_widget_destroy(window);
    }
  
  g_slist_free(winList);
}


static DotterWindowContext* createDotterWindowContext(DotterContext *dotterCtx,
                                                      const IntRange const *refSeqRange,
                                                      const IntRange const *matchSeqRange,
                                                      const gdouble zoomFacIn)
{
  DEBUG_ENTER("createDotterWindowContext");

  DotterWindowContext *result = g_malloc(sizeof *result);
  
  result->dotterCtx = dotterCtx;
  
  result->refSeqRange.min = refSeqRange->min;
  result->refSeqRange.max = refSeqRange->max;
  result->matchSeqRange.min = matchSeqRange->min;
  result->matchSeqRange.max = matchSeqRange->max;

  result->refCoord = UNSET_INT;
  result->matchCoord = UNSET_INT;
  
  result->zoomFactor = getInitZoomFactor(dotterCtx, zoomFacIn, getRangeLength(refSeqRange), getRangeLength(matchSeqRange));
  
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < DOTDIALOG_NUM_DIALOGS; ++dialogId)
    {
      result->dialogList[dialogId] = NULL;
    }
  
  DEBUG_EXIT("createDotterWindowContext returning");
  return result;
}

/* properties specific to a particular dotter window */
static void dotterCreateProperties(GtkWidget *dotterWindow, 
				   GtkWidget *greyrampTool, 
				   GtkWidget *alignmentTool,
                                   GtkWidget *dotplot,
                                   DotterWindowContext *dotterWinCtx)
{
  DEBUG_ENTER("dotterCreateProperties");

  if (dotterWindow)
    {
      DotterProperties *properties = g_malloc(sizeof *properties);

      properties->greyrampTool = greyrampTool;
      properties->alignmentTool = alignmentTool;
      properties->dotplot = dotplot;
      properties->dotterWinCtx = dotterWinCtx;
      
      g_object_set_data(G_OBJECT(dotterWindow), "DotterProperties", properties);
      g_signal_connect(G_OBJECT(dotterWindow), "destroy", G_CALLBACK(onDestroyDotterWindow), NULL); 
    }
  
  DEBUG_EXIT("dotterCreateProperties returning ");
}


/* Utility to extract the dotter context from the dotter window properties */
static DotterContext* dotterGetContext(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  return properties ? properties->dotterWinCtx->dotterCtx : NULL;
}


/* Perform required updates following a change to the currently-selected coords */
static void updateOnSelectedCoordsChanged(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  /* Make sure the coords are in range */
  boundsLimitValue(&properties->dotterWinCtx->refCoord, &properties->dotterWinCtx->refSeqRange);
  boundsLimitValue(&properties->dotterWinCtx->matchCoord, &properties->dotterWinCtx->matchSeqRange);
  
  /* Update the alignment view and dotplot */
  updateAlignmentRange(properties->alignmentTool, properties->dotterWinCtx);
  
  /* Refresh all widgets */
  refreshAll(dotterWindow, NULL);
}


/* Set the initial currently-selected ref seq and match seq coords (i.e. the coords 
 * where the crosshair is centred) */
static void setInitSelectedCoords(GtkWidget *dotterWindow, const int refCoord, const int matchCoord)
{
  DEBUG_ENTER("setInitSelectedCoords");

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  if (valueWithinRange(refCoord, &dwc->refSeqRange))
    dwc->refCoord = refCoord;
  else
    dwc->refCoord = getRangeCentre(&dwc->refSeqRange);

  if (valueWithinRange(matchCoord, &dwc->matchSeqRange))
    dwc->matchCoord = matchCoord;
  else
    dwc->matchCoord = getRangeCentre(&dwc->matchSeqRange);
  
  updateOnSelectedCoordsChanged(dotterWindow);
  
  DEBUG_EXIT("setInitSelectedCoords returning ");
}


/***********************************************************
 *                    External routines                    *
 ***********************************************************/



void dotter (const BlxBlastMode blastMode,
	     DotterOptions *options,
	     const char *queryname,
	     char *queryseq,
	     int   qoff,
	     const BlxStrand refSeqStrand,
	     const char *subjectname,
	      char *subjectseq,
	      int   soff,
	     const BlxStrand matchSeqStrand,
	      int   qcenter,
	      int   scenter,
	      char *saveFileName,
	      char *loadFileName,
	      char *mtxfile,
	      double memoryLimit,
	      int   zoomFacIn,
	      MSP  *MSPs,
	      GList *seqList,
	      int   MSPoff,
	      char *winsizeIn,
	      int   pixelFacIn)
{
  DEBUG_ENTER("dotter(mode=%d, qname=%s, qoff=%d, qstrand=%d, sname=%s, soff=%d, sstrand=%d)",
          blastMode, queryname, qoff, refSeqStrand, subjectname, soff, matchSeqStrand);
  
  gboolean selfComp = FALSE;
  MSPlist = MSPs;
  
  const int qlen = strlen(queryseq);
  const int slen = strlen(subjectseq);
  
  if (qlen < 1) g_error("queryseq is empty");
  if (slen < 1) g_error("subjectseq is empty");

  int i = 0;
  for (i = 0; i < qlen; i++) queryseq[i] = toupper(queryseq[i]);
  for (i = 0; i < slen; i++) subjectseq[i] = toupper(subjectseq[i]);

  if (!memoryLimit) 
    {
      memoryLimit = 0.5; /* Mb */
    }

  if (!strcmp(queryseq, subjectseq)) 
    selfComp = TRUE;

  /* Get score matrix */
  char *matrixName = g_malloc((MAX_MATRIX_NAME_LENGTH + 1) * sizeof(char));
  
  if (mtxfile)	
    {
      readmtx(MATRIX, mtxfile);
      strncpy(matrixName, mtxfile, MAX_MATRIX_NAME_LENGTH);
    }
  else if (blastMode == BLXMODE_BLASTN) 
    {
      DNAmatrix(MATRIX);
      strcpy(matrixName, "DNA+5/-4");
    }
  else
    {
      mtxcpy(MATRIX, BLOSUM62);
      strcpy(matrixName, "BLOSUM62");
    }

  gboolean batchMode = FALSE;
  
  if (saveFileName) 
    {
      /* If a save file was given, that implies we're in batch mode (i.e. don't display any windows) */
      batchMode = TRUE;
      
      /* Don't do batch processing if output file can't be opened */
      if (!fopen (saveFileName, "wb"))
	g_error("Failed to open %s", saveFileName);
    }
  
  /* Create the main dotter context (shared properties for all dotter windows in this process) */
  DotterContext *dotterCtx = createDotterContext(blastMode, batchMode, queryname, queryseq, refSeqStrand, 
                                                 subjectname, subjectseq, matchSeqStrand,
                                                 qoff, soff, options->hozScaleRev, options->vertScaleRev, 
                                                 selfComp, options->mirrorImage, options->watsonOnly, options->crickOnly,
                                                 MSPlist, seqList, MATRIX, matrixName, memoryLimit);

  /* Create a context specific to the initial dotter window */
  DotterWindowContext *dotterWinCtx = createDotterWindowContext(dotterCtx, &dotterCtx->refSeqFullRange, &dotterCtx->matchSeqFullRange, zoomFacIn);

  /* Create the widgets */
  createDotterInstance(dotterCtx, 
                       dotterWinCtx,
                       loadFileName,
                       saveFileName,
                       options->hspsOnly,
                       winsizeIn,
                       pixelFacIn,
                       zoomFacIn,
                       qcenter,
                       scenter,
                       options->swapGreyramp);

  DEBUG_EXIT("dotter returning");
  return ;
}


/* Create all the widgets for a dotter instance. Uses the existing dotter context. Multiple
 * instances (i.e. multiple dotter windows) can exist that share the same main context but display
 * a different range of coords etc,. This creates the widgets and shows them. */
static void createDotterInstance(DotterContext *dotterCtx,
                                 DotterWindowContext *dotterWinCtx,
                                 const char *loadFileName,
                                 const char *saveFileName,
                                 const gboolean hspsOn,
                                 const char* winsizeIn,
                                 const int pixelFacIn,
                                 const int zoomFacIn,
                                 const int qcenter,
                                 const int scenter,
                                 const gboolean greyrampSwap)
{
  DEBUG_ENTER("createDotterInstance");

  GtkWidget *dotplot = NULL;
  GtkWidget *dotplotWidget = createDotplot(dotterWinCtx, 
                                           loadFileName,
                                           saveFileName,
                                           hspsOn,
                                           winsizeIn,
                                           pixelFacIn,
                                           zoomFacIn,
                                           qcenter,
                                           scenter,
                                           &dotplot);
  
  /* Only create the graphical elements if there is a graphical dotplot widget */
  if (dotplotWidget)
    {
      GtkWidget *greyrampTool = createGreyrampTool(dotterCtx, 40, 100, greyrampSwap);
      registerGreyrampCallback(greyrampTool, dotplot, dotplotUpdateGreymap);
      
      GtkWidget *alignmentTool = createAlignmentTool(dotterWinCtx);
  
      const DotterHspMode hspMode = dotplotGetHspMode(dotplot);
      GtkWidget *dotterWindow = createDotterWindow(dotterCtx, dotterWinCtx, hspMode, greyrampTool, dotplotWidget, dotplot);

      /* Set the handlers for the alignment and greyramp tools. Connect them here so we can pass
       * the main window as data. */
      gtk_widget_add_events(alignmentTool, GDK_KEY_PRESS_MASK);
      gtk_widget_add_events(greyrampTool, GDK_KEY_PRESS_MASK);
      g_signal_connect(G_OBJECT(alignmentTool), "key-press-event", G_CALLBACK(onKeyPressDotterCoords), dotterWindow);
      g_signal_connect(G_OBJECT(greyrampTool), "key-press-event", G_CALLBACK(onKeyPressDotter), dotterWindow);

      /* Keep track of all the windows we create, so that we can destroy the DotterContext when
       * the last one is closed */
      dotterCtx->windowList = g_slist_append(dotterCtx->windowList, dotterWindow);
      
      dotterCreateProperties(dotterWindow, greyrampTool, alignmentTool, dotplot, dotterWinCtx);
      DotterProperties *properties = dotterGetProperties(dotterWindow);
      
      setInitSelectedCoords(dotterWindow, qcenter, scenter);
      
      updateGreyMap(properties->greyrampTool);
    }
  
  DEBUG_EXIT("createDotterInstance returning ");
}


/* Open another dotter window, internal to the existing process (i.e. using the same sequences
 * etc. but just displaying a different range). */
void callDotterInternal(DotterContext *dc, 
                        const IntRange const *refSeqRange,
                        const IntRange const *matchSeqRange,
                        const gdouble zoomFactor)
{
  DotterWindowContext *dwc = createDotterWindowContext(dc, refSeqRange, matchSeqRange, zoomFactor);
  createDotterInstance(dc, dwc, NULL, NULL, FALSE, NULL, 0, 0, 0, 0, FALSE);
}



/* 
 *                Internal routines.
 */


/************************/

/* RMEXP will subtract the expected score to remove background noise
   due to biased composition. Gos's idea - not explored yet.

static void rmExp(void){}
*/


///* Utility to return the Feature Series that the given MSP belongs to */
//static FeatureSeries* mspGetFeatureSeries(const MSP *msp)
//{
//  return msp->fs;
//}


/* Returns true if this MSP is part of a Feature Series and that feature series is currently displayed */
//static gboolean mspShowFs(const MSP *msp)
//{
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//  return (fs && fs->on);
//}


/* Return the y position of the lower edge of a feature-series MSP, given its height.
   Calculate it if it is not set yet. Updates the max y coord to be the lower edge of
   this MSP. Returns 0 if the MSP is not in a Feature Series */
//float mspGetFsBottomEdge(MSP *msp, float *maxy, float height)
//{
//  float result = 0;
//  
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//  
//  if (fs)
//    {
//      result = fs->y;
//
//      if (!result)
//        {
//          *maxy += height;
//          fs->y = *maxy;
//          result = *maxy;
//        }
//    }
//  
//  return result;
//}


/* Return the x position of the rightmost edge of a feature-series MSP, given its height.
   Calculae it if it is not yet set. Updates the max x coord to be the rightmost edge of
   this MSP. Returns 0 if the MSP is not in a FeatureSeries */
//float mspGetFsRightEdge(MSP *msp, float *maxx, float height)
//{
//  float result = 0;
//
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//
//  if (fs)
//    {
//      result = fs->x;
//      
//      if (!result)
//        {
//          *maxx += height;
//          fs->x = *maxx;
//          result = *maxx;
//        }
//    }
//  
//  return result;
//}


//static float seq2graphX(int pos)
//{
//    float x;
//
//    x = ceil((float)(pos + MSPoffset - qoffset)/zoom/resfac);
//
//    if (reversedScale)
//	x = RightBorder - x;
//    else
//	x += LeftBorder-1;
//
//    return x;
//}
//
//static float seq2graphY(int pos)
//{
//    float y;
//
//    y = ceil((float)(pos - soffset)/zoom);
//    
//    y += TopBorder-1;
//
//    return y;
//}


//static void XdrawSEG(MSP *msp, float offset)
//{
//    /* Horizontal axis */
//
//    float  
//	sx = seq2graphX(mspGetQStart(msp)), 
//	ex = seq2graphX(mspGetQEnd(msp));
//    
//    offset += TopBorder + slen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    if (fsEndLinesOn) {
//	graphLine(sx, TopBorder, sx, TopBorder+slen4-2);
//	graphLine(ex, TopBorder, ex, TopBorder+slen4-2);
//    }
//
//    graphFillRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);
//    graphColor(BLACK);
//    graphRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);
//
//    graphColor(BLACK);
//    if (fsAnnBottomOn && msp->desc) {
//	graphText(msp->desc, sx, offset);
//    }    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


//static void XdrawSEGxy(MSP *msp, float offset)
//{
//    int i, inNotFilled=0, descShown=0;
//    float  
//	x, y, 
//	xold=0, yold=0;
//
//    offset += TopBorder + slen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    for (i = 0; i < qlen; i++)
//      {
//	const int xyVal = g_array_index(msp->xy, int, i);
//
//	if (xyVal == XY_NOT_FILLED)
//	  {
//	    inNotFilled = 1;
//	  }
//	else
//	  {
//	    x = seq2graphX(i);
//	    y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
//	    
//	    if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE))
//	      {
//	        graphLine(xold, yold, x, y);
//		
//		if (fsAnnBottomOn && msp->desc && !descShown)
//		  {
//		    int linecolor = graphColor(BLACK);
//		    graphText(msp->desc, xold, offset);
//		    graphColor(linecolor);
//		    descShown = 1;
//		  }
//	      }
//		
//	    xold = x;
//	    yold = y;
//	    inNotFilled = 0;
//	  }
//      }
//    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}
//
//
//static void YdrawSEG(MSP *msp, float offset)
//{
//    /* Vertical axis */
//
//    float  
//	sx = seq2graphY(mspGetQStart(msp)), 
//	ex = seq2graphY(mspGetQEnd(msp));
//    
//    offset += LeftBorder + qlen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    if (fsEndLinesOn) {
//	graphLine(LeftBorder, sx, LeftBorder+qlen4-2, sx);
//	graphLine(LeftBorder, ex, LeftBorder+qlen4-2, ex);
//    }
//
//    graphFillRectangle(offset, sx, offset - msp->score/100.0 * fonth, ex);
//    graphColor(BLACK);
//    graphRectangle    (offset, sx, offset - msp->score/100.0 * fonth, ex);
//
//    graphColor(BLACK);
//    if (fsAnnRightOn && msp->desc) {
//	graphText(msp->desc, offset, sx);
//    }    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}
//
//
//static void YdrawSEGxy(MSP *msp, float offset)
//{
//    int i, inNotFilled=0, descShown=0;
//    float  
//	x, y, 
//	xold=0, yold=0;
//
//    offset += LeftBorder + qlen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    for (i = 0; i < qlen; i++)
//      {
//	const int xyVal = g_array_index(msp->xy, int, i);
//	
//	if (xyVal == XY_NOT_FILLED)
//	  {
//	    inNotFilled = 1;
//	  }
//	else
//	  {
//	    x = seq2graphY(i);
//	    y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
//	    
//	    if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE)) 
//	      {
//		graphLine(yold, xold, y, x);
//		
//		if (fsAnnRightOn && msp->desc && !descShown) 
//		  {
//		    int linecolor = graphColor(BLACK);
//		    graphText(msp->desc, offset, xold);
//		    graphColor(linecolor);
//		    descShown = 1;
//		  }
//	      }
//
//	    xold = x;
//	    yold = y;
//	    inNotFilled = 0;
//	  }
//      }
//    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


//static int isHorizontalMSP(MSP *msp) {
//
//    if (!msp->qname)
//      g_error("No qname set in MSP - I need this to associate it with one of the sequences");
//
//    if (!strcasecmp(msp->qname, qname) || !strcmp(msp->qname, "@1"))
//        return 1;
//    else
//        return 0;
//}
//
//
//static int isVerticalMSP(MSP *msp) {
//
//    if (!msp->qname) 
//      g_error("No qname set in MSP - I need this to associate it with one of the sequences");
//
//    if (!strcasecmp(msp->qname, sname) || !strcmp(msp->qname, "@2"))
//        return 1;
//    else
//        return 0;
//}


//static void drawAllFeatures(MSP *msp)
//{
//  float sx, ex, 
//    y, 
//    height,		/* box height/width */
//    boxHeight=10,
//    textHeight,
//    oldLinew;
//  int i;
//  float posx=0, posy=0;		/* Next series to be drawn */
//
//  int top, bottom ;
//  float feature_boundary, feature_top, feature_bottom, feature_strand ;
//  float forward_y, reverse_y, depth ;
//  float old_line_width ;
//
//
//  /* Set forward/reverse strand gene drawing positions. */
//  graphGetBounds(&top, &bottom) ;
//  feature_boundary = 3.0 ;				    /* Allows for line thickness etc... */
//  feature_top = TopBorder + slen4 ;
//  feature_bottom = bottom ;
//  feature_strand = feature_top + ((feature_bottom - feature_top + 1) / 2) ;
//
//  forward_y = feature_top + feature_boundary ;
//  reverse_y = feature_strand + feature_boundary ;
//  depth = (feature_strand - feature_boundary) - (feature_top + feature_boundary) - fonth ;
//
//
//  /* Draw a strand separator. */
//  old_line_width = graphLinewidth(2) ;
//  graphLine(LeftBorder - 1, feature_strand, LeftBorder - 1 + qlen4, feature_strand) ;
//  graphLinewidth(old_line_width) ;
//
//
//  /* Now draw the genes. */
//  drawGenes(msp, forward_y, reverse_y, depth) ;
//
//
//  if (fsArr)
//    {
//      /* Loop through each feature-series in the feature-series array and set coords to 0 */
//      for (i = 0; i < gArrayGetLen(fsArr); i++)
//	{
//	  FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//	  fs->y = fs->x = 0;
//	}
//    }
//
//  textHeight = fonth;
//  boxHeight = fonth;
//
//  oldLinew = graphLinewidth(1);
//
//  BlxStrand strand = BLXSTRAND_NONE;
//  
//  if (selfcomp || !reversedScale) 
//    strand = BLXSTRAND_FORWARD;
//  else 
//    strand = BLXSTRAND_REVERSE;
//
//
//  for (; msp; msp = msp->next)
//    {    
//      height = boxHeight;
//
//      if (mspHasFs(msp))
//	{
//	  sx = seq2graphX(mspGetQStart(msp));
//	  ex = seq2graphX(mspGetQEnd(msp));
//
//	  if (msp->qStrand != strand)
//	    y = reverse_y ;
//	  else
//	    y = forward_y ;
//
//          if (!mspShowFs(msp))
//	    continue;
//
//	  /* Adjust height to score */
//	  if (msp->score > 0)
//	    {
//	      height = boxHeight * msp->score / 100.0;
//	    }
//
//	  if (fsBottomOn && (isHorizontalMSP(msp) || (selfcomp && isVerticalMSP(msp))))
//	    {
//	      /* HORIZONTAL AXIS (X) */
//		    
//	      if (msp->type == BLXMSP_XY_PLOT)
//		{
//		  XdrawSEGxy(msp, mspGetFsRightEdge(msp, &posx, fonth*(fsPlotHeight+1)));
//		}
//	      else if (msp->type == BLXMSP_FS_SEG)
//		{
//		  XdrawSEG(msp, mspGetFsRightEdge(msp, &posx, boxHeight+textHeight));
//		}
//	    }
//
//	  if (fsRightOn &&  (isVerticalMSP(msp) || (selfcomp && isHorizontalMSP(msp))))
//	    {
//	      /* VERTICAL AXIS (Y) */
//
//	      if (msp->type == BLXMSP_XY_PLOT)
//		{
//		  YdrawSEGxy(msp, mspGetFsBottomEdge(msp, &posy, fonth*(fsPlotHeight+1)));
//		}
//	      else if (msp->type == BLXMSP_FS_SEG)
//		{
//		  YdrawSEG(msp, mspGetFsBottomEdge(msp, &posy, boxHeight+textHeight));
//		}
//	    }
//
//	  graphColor(oldcolor); 
//	  graphLinewidth(oldLinew);
//	}
//    }
//
//  return ;
//}


//static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth)
//{
//  gboolean bump_genes = FALSE ;				    /* Make this user settable... */
//
//  float height,		/* box height/width */
//    sy, ey, midy, x, y, oldLinew ;
//
//  height = fonth ;
//
//  oldLinew = graphLinewidth(1);
//
//  BlxStrand strand = BLXSTRAND_NONE;
//  if (selfcomp || !reversedScale) 
//    strand = BLXSTRAND_FORWARD;
//  else 
//    strand = BLXSTRAND_REVERSE;
//
//
//  if (bump_genes)
//    {
//      MSP *gene_msp, *tmp ;
//      GList *exon_intron_list = NULL ;
//      GeneStrandStruct strand_genes = {NULL} ;
//
//      /* MSP's are in any old order and only some are exons/introns, here we make a list of
//       * intron/exon MSP's and then sort that list into transcripts for drawing... */
//
//      /* Sort all introns/exons by gene and position within gene. */
//      gene_msp = tmp = msp ;
//      for ( ; tmp ; tmp = tmp->next)
//	{
//	  if (tmp->score < 0)
//	    exon_intron_list = g_list_append(exon_intron_list, tmp) ;
//	}
//      exon_intron_list = g_list_sort(exon_intron_list, compareMSPs) ;
//
//
//#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
//      g_list_foreach(exon_intron_list, printMSP, NULL) ;
//#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
//
//
//      /* Produce two lists, one for each strand, of genes sorted by position from the sorted intron/exon list. */
//      strand_genes.strand = '+' ;
//      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
//      strand_genes.forward_genes = g_list_sort(strand_genes.forward_genes, compareGenes) ;
//      strand_genes.forward_genes = g_list_first(strand_genes.forward_genes) ;
//
//      strand_genes.strand = '-' ;
//      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
//      strand_genes.reverse_genes = g_list_sort(strand_genes.reverse_genes, compareGenes) ;
//      strand_genes.reverse_genes = g_list_first(strand_genes.reverse_genes) ;
//
//      /* Set y offsets for forward and reverse stranded genes. */
//      setYoffsets(strand_genes.forward_genes, forward_y, forward_y + depth) ;
//      setYoffsets(strand_genes.reverse_genes, reverse_y, reverse_y + depth) ;
//
//
//#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
//      printf("Forward:\n") ;
//      g_list_foreach(strand_genes.forward_genes, printGene, NULL) ;
//      printf("Reverse:\n") ;
//      g_list_foreach(strand_genes.reverse_genes, printGene, NULL) ;
//      printf("\n") ;
//#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
//
//      /* Draw the strands. */
//      oldcolor = graphColor(BLUE); 
//      g_list_foreach(strand_genes.forward_genes, drawGenesCB, &strand_genes) ;
//
//      oldcolor = graphColor(MAGENTA); 
//      g_list_foreach(strand_genes.reverse_genes, drawGenesCB, &strand_genes) ;
//
//      graphColor(oldcolor); 
//
//    }
//  else
//    {
//      for (; msp; msp = msp->next)
//	{    
//	  if (msp->score < 0)
//	    {
//	      if (msp->qStrand != strand)
//		y = reverse_y ;
//	      else
//		y = forward_y ;
//
//	      oldcolor = graphColor(BLUE); 
//	      drawMSPGene(msp, y) ;
//
//	      if (selfcomp) /* Draw the boxes on vertical axes too */
//		{
//		  sy = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/zoom);
//		  ey = ceil((float)(mspGetQEnd(msp)+MSPoffset - qoffset)/zoom);
//		
//		  sy += TopBorder-1;
//		  ey += TopBorder-1;
//		
//		  x = LeftBorder + qlen4 + 10;
//		  if (msp->qStrand != strand) x += 20;
//		
//		  if (msp->score == -1.0) /* EXON */
//		    {
//		      oldcolor = graphColor(BLUE); 
//		      graphRectangle(x, sy, x + height, ey);
//		    }
//		  else if (msp->score == -2.0) /* INTRON */
//		    {
//		      oldcolor = graphColor(BLUE); 
//		      midy = 0.5 * (sy + ey) ;
//		      graphLine (x + height/2, sy, x, midy) ;
//		      graphLine (x + height/2, ey, x, midy) ;
//		    }
//		}
//
//	      graphColor(oldcolor) ;
//	    }
//	}
//    }
//
//  graphLinewidth(oldLinew) ;
//
//  return ;
//}


/* A GCompareFunc() to compare gene MSPs.... */
//gint compareMSPs(gconstpointer a, gconstpointer b)
//{
//  int result = 0 ;
//  MSP *msp_a = (MSP *)a, *msp_b = (MSP *)b ;
//
//  if (!(result = strcmp(msp_a->sname, msp_b->sname)))
//    {
//      if (mspGetSStart(msp_a) < mspGetSStart(msp_b))
//	result = -1 ;
//      else if (mspGetSStart(msp_a) > mspGetSStart(msp_b))
//	result = 1 ;
//      else
//	{
//	  /* I actually don't think this should ever happen as it means either there are
//	   * duplicates or worse there are overlapping introns/exons within a gene.... */
//	  if (mspGetSEnd(msp_a) < mspGetSEnd(msp_b))
//	    result = -1 ;
//	  else if (mspGetSEnd(msp_a) > mspGetSEnd(msp_b))
//	    result = 1 ;
//	  else
//	    result = 0 ;
//	}
//    }
//
//  return result ;
//}


///* A GFunc() to record start/end of genes. */
//static void getGenePositionsCB(gpointer data, gpointer user_data)
//{
//  MSP *msp = (MSP *)data ;
//  GeneStrand strand_genes = (GeneStrand)user_data ;
//  GList *gene_list ;
//  GeneData gene ;
//
//  /* We need genes added according to strand... */
//  if (msp->qframe[1] == strand_genes->strand)
//    {
//      const char *curr_name = NULL ;
//      int curr_length ;
//
//      if (strand_genes->strand == '+')
//	gene_list = strand_genes->forward_genes ;
//      else
//	gene_list = strand_genes->reverse_genes ;
//
//      if (gene_list)
//	{
//	  /* Look at last element, this is the last gene we added. */
//	  gene_list = g_list_last(gene_list) ;
//
//	  curr_name = (((GeneData)(gene_list->data))->name) ;
//	  curr_length = strlen(curr_name) - 1 ;
//	}
//
//      /* If there's no gene or we are starting a new gene then just add to the list,
//       * otherwise record latest msp end position. */
//      if (!gene_list || strncmp(mspGetSName(msp), curr_name, curr_length) != 0)
//	{
//	  gene = g_new0(GeneDataStruct, 1) ;
//	  gene->name = mspGetSName(msp) ;
//	  gene->start = mspGetSStart(msp) ;
//	  gene->end = mspGetSEnd(msp) ;
//	  gene->strand = msp->qframe[1] ;
//	  gene->msp_start = msp ;
//
//	  gene_list = g_list_append(gene_list, gene) ;
//	}
//      else
//	{
//	  gene = (GeneData)(gene_list->data) ;
//
//	  gene->end = mspGetSEnd(msp) ;
//	  gene->msp_end = msp ;
//	}
//
//      if (strand_genes->strand == '+')
//	strand_genes->forward_genes = gene_list ;
//      else
//	strand_genes->reverse_genes = gene_list ;
//    }
//
//  return ;
//}
//
//

/* A GCompareFunc() to compare genes for position.... */
//gint compareGenes(gconstpointer a, gconstpointer b)
//{
//  int result = 0 ;
//  GeneData gene_a = (GeneData)a, gene_b = (GeneData)b ;
//
//  if (gene_a->strand == '+' && gene_b->strand == '-')
//    result = -1 ;
//  else if (gene_a->strand == '-' && gene_b->strand == '+')
//    result = 1 ;
//  else
//    {
//      if (gene_a->start < gene_b->start)
//	result = -1 ;
//      else if (gene_a->start > gene_b->start)
//	result = 1 ;
//      else
//	result = 0 ;
//    }
//
//  return result ;
//}



//static void setYoffsets(GList *gene_list, float min_offset, float max_offset)
//{
//  GList *curr_ptr, *next_ptr ;
//  float curr_y = 0.0, bump_incr = 0.5 ;
//
//
//  /* Go through the gene list comparing, reordering and assigning y offsets to the genes
//   * so they can be drawn without overlapping. */
//  curr_ptr = gene_list ;
//  next_ptr = curr_ptr->next ;
//  curr_y = min_offset - bump_incr ;
//
//  /* Go through all genes. */
//  while (curr_ptr)
//    {
//      GList *list_ptr = g_list_first(gene_list) ;
//
//      /* For each gene look at all the preceding genes and to decide its y offset so it is
//       * bumped correctly. */
//      while (list_ptr)
//	{
//	  GeneData list_gene = (GeneData)list_ptr->data ;
//	  GeneData curr_gene = (GeneData)curr_ptr->data ;
//
//	  if (list_ptr == curr_ptr)
//	    {
//	      /* This gene overlaps all previous ones and needs a new offset. */
//	      curr_y += bump_incr ;
//	      if (curr_y > max_offset)
//		curr_y = max_offset ;
//
//	      curr_gene->y_pos = curr_y ;
//	      break ;
//	    }
//	  else if (curr_gene->start > list_gene->end)
//	    {
//	      /* This gene coes not overlap the list gene so give is the same offset. */
//	      curr_gene->y_pos = list_gene->y_pos ;
//
//	      if (list_ptr->next != curr_ptr)
//		{
//		  list_ptr = g_list_remove(list_ptr, curr_gene) ;
//		  list_ptr = g_list_insert_before(gene_list, list_ptr->next, curr_gene) ;
//		}
//
//	      break ;
//	    }
//	  else
//	    {
//	      /* This gene overlaps the list gene so move on to the next one. */
//	      list_ptr = g_list_next(list_ptr) ;
//	    }
//	}
//
//      /* update curr/next until we get to the end of the list... */
//      if ((curr_ptr = next_ptr))
//	next_ptr = curr_ptr->next ;
//    }
//
//  return ;
//}
//

///* A GFunc() to record start/end of genes. */
//static void drawGenesCB(gpointer data, gpointer user_data_unused)
//{
//  GeneData gene = (GeneData)data ;
//  MSP *msp ;
//
//  msp = gene->msp_start ;
//  do
//    {
//      drawMSPGene(msp, gene->y_pos) ;
//
//      msp = msp->next ;
//    } while (msp && msp != gene->msp_end) ;
//
//  return ;
//}
//
//
//static void drawMSPGene(MSP *msp, float y_offset)
//{
//  float height,		/* box height/width */
//    sx, ex, midx ;
//
//  height = fonth ;
//
//  sx = seq2graphX(mspGetQStart(msp));
//  ex = seq2graphX(mspGetQEnd(msp));
//
//  if (msp->score == -1.0) /* EXON */
//    {
//      graphRectangle(sx, y_offset, ex, y_offset + height);
//    }
//  else if (msp->score == -2.0) /* INTRON */
//    {
//
//      midx = 0.5 * (sx + ex) ;
//      graphLine (sx, y_offset + height/2, midx, y_offset) ;
//      graphLine (ex, y_offset + height/2, midx, y_offset) ;
//    }
//
//  return ;
//}
//

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


//static void fsSelAll(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      fs->on = 1;
//      graphBoxDraw(fsBoxStart+i, WHITE, BLACK);
//    }
//}
//
//static void fsSelNone(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      fs->on = 0;
//      graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//    }
//}
//
//static void fsSelNoCurves(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      
//      if (fs->xy) 
//	{
//	  fs->on = 0;
//	  graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//	}
//    }
//}
//
//static void fsSelNoSegments(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      
//      if (!fs->xy) 
//	{
//	  fs->on = 0;
//	  graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//	}
//    }
//}
//
//static void fsToggleBottom(void)
//{
//    fsBottomOn = !fsBottomOn;
//    selectFeatures();
//}
//static void fsToggleRight(void)
//{
//    fsRightOn = !fsRightOn;
//    selectFeatures();
//}
//static void fsToggleAnnBottom(void)
//{
//    fsAnnBottomOn = !fsAnnBottomOn;
//    selectFeatures();
//}
//static void fsToggleAnnRight(void)
//{
//    fsAnnRightOn = !fsAnnRightOn;
//    selectFeatures();
//}
//static void fsToggleEndLines(void)
//{
//    fsEndLinesOn = !fsEndLinesOn;
//    selectFeatures();
//}
//static void fsSel(int box, double x_unused, double y_unused, int modifier_unused)
//{
//    if (box-fsBoxStart < 0 || box-fsBoxStart > gArrayGetLen(fsArr))
//	return;
//    
//    FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, box - fsBoxStart);
//    int *on = &fs->on;
//
//    graphActivate(fsGraph);
//
//
//    if (*on) {
//	*on = 0;
//	graphBoxDraw(box, BLACK, WHITE);
//    }
//    else {
//	*on = 1;
//	graphBoxDraw(box, WHITE, BLACK);
//    }
//
//}
//static void setfsPlotHeight(char *cp)
//{
//    fsPlotHeight = atof(cp);
//}
//    
//
//void selectFeatures(void)
//{
//  int i, box;
//  float y=1.0;
//
//  if (!graphActivate(fsGraph))
//    {
//      fsGraph = graphCreate (TEXT_SCROLL, "Feature Series Selection Tool", 0, 0, 0.4, 0.6);
//    }
//  graphPop();
//  graphRegister(PICK, fsSel);
//
//  if (1 /* dotterWindow */) {
//
//    box = graphButton("Show bottom series", fsToggleBottom, 1, y);
//    if (!fsBottomOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show right series", fsToggleRight, 1, y);
//    if (!fsRightOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show bottom annotations", fsToggleAnnBottom, 1, y);
//    if (!fsAnnBottomOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show right annotations", fsToggleAnnRight, 1, y);
//    if (!fsAnnRightOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Draw lines at segment ends", fsToggleEndLines, 1, y);
//    if (!fsEndLinesOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 2;
//  }
//
//  sprintf(fsPlotHeighttx, "%.1f", fsPlotHeight);
//  graphText("Height of curves (XY series):", 1, y);
//  graphTextScrollEntry (fsPlotHeighttx, 6, 4, 31, y, setfsPlotHeight);
//  y += 2;
//
//  graphLinewidth(0.1);
//  graphLine(0, y, 2000, y);
//  y += 0.5;
//
//  graphText("Pick to select/unselect series", 1, y);
//  y += 1.5;
//  graphButton("All", fsSelAll, 1, y);
//  graphButton("No curves", fsSelNoCurves, 13, y);
//  graphButton("No segments", fsSelNoSegments, 24, y);
//  fsBoxStart = 1+graphButton("None", fsSelNone, 6, y);
//  y += 2;
//
//  graphTextBounds(50, gArrayGetLen(fsArr) * 1.5 + y + 5);
//
//  for (i = 0; i < gArrayGetLen(fsArr); i++)
//    {
//      float  margin = 0.1;
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//
//      box = graphBoxStart();
//      graphText(fs->name, 1, y);      
//      graphRectangle(1 - margin, y - margin, 1 + margin + strlen(fs->name), y + 1 + margin);
//      graphBoxEnd();
//
//      if (!fs->on) 
//	{
//	  graphBoxDraw(box, BLACK, WHITE);
//	}
//      else 
//	{
//	  graphBoxDraw(box, WHITE, BLACK);
//	}
//	
//      y += 1.5;
//    }
//
//  graphRedraw();
//
//  return ;
//}
//
//
//float fsTotalHeight(MSP *msplist)
//{
//    int i;
//    float maxy = 0;
//	
//    if (!fsArr || !gArrayGetLen(fsArr))
//      {
//	return 0.0;
//      }
//
//    for (i = 0; i < gArrayGetLen(fsArr); i++) 
//      {
//	FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//	fs->y = fs->x = 0;
//      }
//	
//    for (msp = msplist; msp; msp = msp->next) 
//      {
//        if (mspShowFs(msp))
//	  {
//	    if (msp->type == BLXMSP_XY_PLOT) 
//	      {
//		mspGetFsBottomEdge(msp, &maxy, fsPlotHeight+1);
//	      }
//	    else if (msp->type == BLXMSP_FS_SEG) 
//	      {
//		mspGetFsBottomEdge(msp, &maxy, 1+1);
//	      }
//	  }
//      }
//    
//    return maxy + 2;
//}


/* Callbacks to be called when the dotter parameters have changed. */
static gboolean onZoomFactorChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = FALSE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  gdouble newValue = g_strtod(text, NULL);
  
  if (newValue <= 0)
    {
      g_critical("Zoom factor must be greater than zero.\n");
    }
  else
    {
      properties->dotterWinCtx->zoomFactor = newValue;
      result = TRUE;
    }
  
  return result;
}

static gboolean onQStartChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->refSeqFullRange))
    g_warning("Limiting reference sequence start to range %d -> %d.\n", dwc->dotterCtx->refSeqFullRange.min, dwc->dotterCtx->refSeqFullRange.max);
  
  boundsLimitValue(&newValue, &dwc->dotterCtx->refSeqFullRange);
  
  properties->dotterWinCtx->refSeqRange.min = newValue;
  
  return TRUE;
}

static gboolean onQEndChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->refSeqFullRange))
    g_warning("Limiting reference sequence end to range %d -> %d.\n", dwc->dotterCtx->refSeqFullRange.min, dwc->dotterCtx->refSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->refSeqFullRange);
  
  properties->dotterWinCtx->refSeqRange.max = newValue;
  
  return TRUE;
}

static gboolean onSStartChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->matchSeqFullRange))
    g_warning("Limiting vertical sequence start to range %d -> %d.\n", dwc->dotterCtx->matchSeqFullRange.min, dwc->dotterCtx->matchSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->matchSeqFullRange);
  
  properties->dotterWinCtx->matchSeqRange.min = newValue;
  
  return TRUE;
}

static gboolean onSEndChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->matchSeqFullRange))
    g_warning("Limiting vertical sequence end to range %d -> %d.\n", dwc->dotterCtx->matchSeqFullRange.min, dwc->dotterCtx->matchSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->matchSeqFullRange);
  
  properties->dotterWinCtx->matchSeqRange.max = newValue;
  
  return TRUE;
}

static gboolean onSlidingWinSizeChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = FALSE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);

  GError *error = NULL;
  dotplotSetSlidingWinSize(properties->dotplot, newValue, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  else
    {
      result = TRUE;
    }
  
  return result;
}


/* Create a text entry box initialised with the given double */
static void createTextEntryFromDouble(GtkWidget *dotterWindow,
                                      GtkTable *table, 
                                      const int row,
                                      const int col,
                                      const int xpad,
                                      const int ypad,
                                      const char *mnemonic,
                                      const double value,
                                      BlxResponseCallback callback)
{
  if (mnemonic)
    {
      GtkWidget *label = gtk_label_new_with_mnemonic(mnemonic);
      gtk_misc_set_alignment(GTK_MISC(label), 1, 0);
      gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
    }
  
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);

  /* Only display decimal places if not a whole number */
  const int numDp = value - (int)value > 0 ? 1 : 0;
  
  char *displayText = convertDoubleToString(value, numDp);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 3);

  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  widgetSetCallbackData(entry, callback, dotterWindow);
}

/* Create a text entry box initialised with the given integer */
static void createTextEntryFromInt(GtkWidget *dotterWindow,
                                   GtkTable *table, 
                                   const int row,
                                   const int col,
                                   const int xpad,
                                   const int ypad,
                                   const char *mnemonic,
                                   const int value,
                                   BlxResponseCallback callback)
{
  if (mnemonic)
    {
      GtkWidget *label = gtk_label_new_with_mnemonic(mnemonic);
      gtk_misc_set_alignment(GTK_MISC(label), 1, 0);
      gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
    }
  
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  
  char *displayText = convertIntToString(value);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 3);

  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  widgetSetCallbackData(entry, callback, dotterWindow);
}


/* Callback when we receive a response for the settings dialog */
static void onResponseSettingsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  GtkWidget *dotterWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Destroy if successful */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      redrawAll(dotterWindow, NULL);
      break;
      
    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      redrawAll(dotterWindow, NULL);
      destroy = FALSE;
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      /* If it's a persistent dialog, just hide it, otherwise destroy it */
      const gboolean isPersistent = GPOINTER_TO_INT(data);
      
      if (isPersistent)
        {
          gtk_widget_hide_all(GTK_WIDGET(dialog));
        }
      else
        {
          gtk_widget_destroy(GTK_WIDGET(dialog));
        }
    }
}


/* Pop up a dialog to allow the user to set the dotter parameters */
static void showSettingsDialog(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  const DotterDialogId dialogId = DOTDIALOG_SETTINGS;
  GtkWidget *dialog = getPersistentDialog(dwc->dialogList, dialogId);
  
  if (!dialog)
    {
      /* Create the dialog */
      dialog = gtk_dialog_new_with_buttons("Settings", 
                                           GTK_WINDOW(dotterWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           NULL);
      
      /* These 2 calls are required to make the dialog persistent... */
      addPersistentDialog(dwc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseSettingsDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Clear contents and re-create */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }

  const int numRows = 4;
  const int numCols = 3;
  const int xpad = 2;
  const int ypad = 2;
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), GTK_WIDGET(table));
  
  createTextEntryFromDouble(dotterWindow, table, 1, 2, xpad, ypad, "_Zoom: ", dwc->zoomFactor, onZoomFactorChanged);
  createTextEntryFromInt(dotterWindow, table, 2, 2, xpad, ypad, "_Horizontal range: ", dwc->refSeqRange.min, onQStartChanged);
  createTextEntryFromInt(dotterWindow, table, 2, 3, xpad, ypad, NULL, dwc->refSeqRange.max, onQEndChanged);
  createTextEntryFromInt(dotterWindow, table, 3, 2, xpad, ypad, "_Vertical range: ", dwc->matchSeqRange.min, onSStartChanged);
  createTextEntryFromInt(dotterWindow, table, 3, 3, xpad, ypad, NULL, dwc->matchSeqRange.max, onSEndChanged);
  createTextEntryFromInt(dotterWindow, table, 4, 2, xpad, ypad, "Sliding _window size: ", dotplotGetSlidingWinSize(properties->dotplot), onSlidingWinSizeChanged);
  

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
  
  return;
}


/* Redraw the main dotter window, the alignment tool and the greyramp tool. Re-calculates borders etc. */
static void redrawAll(GtkWidget *dotterWindow, gpointer data)
{
  gtk_widget_queue_draw(dotterWindow);
  
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  /* Check the range values are the correct way round. */
  sortValues(&dwc->refSeqRange.min, &dwc->refSeqRange.max, TRUE);
  sortValues(&dwc->matchSeqRange.min, &dwc->matchSeqRange.max, TRUE);
  
  if (properties)
    {
      gtk_widget_queue_draw(properties->greyrampTool);
      gtk_widget_queue_draw(properties->alignmentTool);
      redrawDotplot(properties->dotplot);
    }
}


/* Refresh the main dotter window, the alignment tool and the greyramp tool. Clears any cached
 * pixmaps but does not recalculate borders etc. */
static void refreshAll(GtkWidget *dotterWindow, gpointer data)
{
  gtk_widget_queue_draw(dotterWindow);
  
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties)
    {
      gtk_widget_queue_draw(properties->greyrampTool);
      gtk_widget_queue_draw(properties->alignmentTool);
      refreshDotplot(properties->dotplot);
    }
}


static void readmtx(int MATRIX[24][24], char *mtxfile)
{
    FILE *fil;
    int row, col;
    char line[1025] = "#", *p;
    
    char *mtxfileText = blxprintf("%s/%s", getenv("BLASTMAT"), mtxfile);
  
    if (!(fil = fopen(mtxfile, "r")) &&
	!(fil = fopen(mtxfileText, "r")))
      {
        char *msg = blxprintf("Failed to open score matrix file %s - not found in ./ or $BLASTMAT/.");
	g_error(msg, mtxfile);
        g_free(msg);
      }
  
    g_free(mtxfileText);
    mtxfileText = NULL;
    
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

    g_message("I read your score matrix %s.\n", mtxfile);
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


/* Utility to get the length of the given GArray. Returns 0 if array is null. */
//static int gArrayGetLen(GArray *array)
//{
//  return (array ? array->len : 0);
//}


/***********************************************************
 *               Show/hide parts of the view               *
 ***********************************************************/

/* Show the greyramp tool, bringing it to the front. Create it if it doesn't exist */
static void showGreyrampTool(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties->greyrampTool && GTK_IS_WIDGET(properties->greyrampTool))
    {
      gtk_widget_show_all(properties->greyrampTool);
      
      if (GTK_IS_WINDOW(properties->greyrampTool))
	{
	  gtk_window_present(GTK_WINDOW(properties->greyrampTool));
	}
    }
  else
    {
      properties->greyrampTool = createGreyrampTool(properties->dotterWinCtx->dotterCtx, 40, 100, FALSE);
    }
}

/* Show the alignment tool, bringing it to the front. Create it if it doesn't exist */
static void showAlignmentTool(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties->alignmentTool && GTK_IS_WIDGET(properties->alignmentTool))
    {
      gtk_widget_show_all(properties->alignmentTool);
      
      if (GTK_IS_WINDOW(properties->alignmentTool))
	{
	  gtk_window_present(GTK_WINDOW(properties->alignmentTool));
	}
    }
  else
    {
      properties->alignmentTool = createAlignmentTool(properties->dotterWinCtx);
    }
}

/* Bring the main dotter window to the front */
static void showDotterWindow(GtkWidget *dotterWindow)
{
  gtk_widget_show_all(dotterWindow);
  gtk_window_present(GTK_WINDOW(dotterWindow));
}


/***********************************************************
 *                       Help Dialog                       *
 ***********************************************************/

/* Returns a string which is the name of the Dotter application. */
static char *dotterGetAppName(void)
{
  return DOTTER_TITLE ;
}

/* Returns a copyright string for the Dotter application. */
static char *dotterGetCopyrightString(void)
{
  return DOTTER_COPYRIGHT_STRING ;
}

/* Returns the Dotter website URL. */
static char *dotterGetWebSiteString(void)
{
  return DOTTER_WEBSITE_STRING ;
}

/* Returns a comments string for the Dotter application. */
static char *dotterGetCommentsString(void)
{
  return DOTTER_COMMENTS_STRING() ;
}

/* Returns a license string for the dotter application. */
static char *dotterGetLicenseString(void)
{
  return DOTTER_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Dotter code. */
static char *dotterGetVersionString(void)
{
  return DOTTER_VERSION_STRING ;
}

/* Shows the 'About' dialog */
static void showAboutDialog(GtkWidget *parent)
{
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (6)
  const gchar *authors[] = {AUTHOR_LIST, NULL} ;
  
  gtk_show_about_dialog(GTK_WINDOW(parent),
			"authors", authors,
			"comments", dotterGetCommentsString(), 
			"copyright", dotterGetCopyrightString(),
			"license", dotterGetLicenseString(),
			"name", dotterGetAppName(),
			"version", dotterGetVersionString(),
			"website", dotterGetWebSiteString(),
			NULL) ;
#endif
  
  return ;
}


static void onResponseHelpDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      destroy = TRUE;
      break;
      
    case GTK_RESPONSE_HELP:
      showAboutDialog(NULL);
      destroy = FALSE;
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      /* If it's a persistent dialog, just hide it, otherwise destroy it */
      const gboolean isPersistent = GPOINTER_TO_INT(data);
      
      if (isPersistent)
        {
          gtk_widget_hide_all(GTK_WIDGET(dialog));
        }
      else
        {
          gtk_widget_destroy(GTK_WIDGET(dialog));
        }
    }
}


static void showHelpDialog(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  const DotterDialogId dialogId = DOTDIALOG_HELP;
  
  GtkWidget *dialog = getPersistentDialog(dwc->dialogList, dialogId);
  
  if (!dialog)
    {
      /* Create the dialog */
      dialog = gtk_dialog_new_with_buttons("Help", 
                                           GTK_WINDOW(dotterWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_ABOUT,
                                           GTK_RESPONSE_HELP,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These 2 calls are required to make the dialog persistent... */
      addPersistentDialog(dwc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

      g_signal_connect(dialog, "response", G_CALLBACK(onResponseHelpDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Clear contents and re-create */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
      
  GdkScreen *screen = gtk_widget_get_screen(dotterWindow);
  const int width = gdk_screen_get_width(screen) * MAX_WINDOW_WIDTH_FRACTION;
  int height = gdk_screen_get_height(screen) * MAX_WINDOW_HEIGHT_FRACTION;

  DotterContext *dc = dwc->dotterCtx;

  char *messageText = blxprintf(DOTTER_HELP_TEXT, 
                                dotplotGetSlidingWinSize(properties->dotplot), 
                                dotplotGetPixelFac(properties->dotplot), 
                                dc->matrixName, 
                                dwc->zoomFactor);
  
  GtkWidget *child = createScrollableTextView(messageText, TRUE, dotterWindow->style->font_desc, TRUE, &height, NULL);
  
  g_free(messageText);
  
  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);
  
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
}


/***********************************************************
 *                          Events                         *
 ***********************************************************/

static void onQuitMenu(GtkAction *action, gpointer data)
{
  /* Close all associated windows */
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterContext *dc = dotterGetContext(dotterWindow);
  dotterContextCloseAllWindows(dc);
}

static void onSavePlotMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  GError *error = NULL;
  savePlot(properties->dotplot, NULL, NULL, &error);
  
  prefixError(error, "Error saving plot. ");
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
}

static void onPrintMenu(GtkAction *action, gpointer data)
{
/* to do: implement this */
//  GtkWidget *dotterWindow = GTK_WIDGET(data);
}

static void onSettingsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showSettingsDialog(dotterWindow);
}

static void onHelpMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showHelpDialog(dotterWindow);
}

static void onAboutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showAboutDialog(dotterWindow);
}

static void onShowGreyrampMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showGreyrampTool(dotterWindow);
}

static void onShowAlignmentMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showAlignmentTool(dotterWindow);
}

static void onToggleCrosshairMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairOn(properties->dotplot);
}

static void onToggleCoordsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairCoordsOn(properties->dotplot);
}

static void onToggleFullscreenMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairFullscreen(properties->dotplot);
}

static void onTogglePixelmapMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  dotplotTogglePixelmap(properties->dotplot);
}

static void onToggleGridMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  dotplotToggleGrid(properties->dotplot);
}

static void onToggleHspMode(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const DotterHspMode hspMode = gtk_radio_action_get_current_value(current);
  setHspMode(properties->dotplot, hspMode);
}

//static void GHelp(GtkButton *button, gpointer data)
//{
//  GtkWidget *dotterWindow = GTK_WIDGET(data);
//  showHelpDialog(dotterWindow);
//}


/* Mouse button handler */
static gboolean onButtonPressDotter(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      GtkMenu *contextMenu = GTK_MENU(data);
      gtk_menu_popup (contextMenu, NULL, NULL, NULL, NULL, event->button, event->time);
    }
  else if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      /* If the dot-plot was clicked the selected coords will have changed. Perform required updates. */
      updateOnSelectedCoordsChanged(window);
      handled = TRUE;
    }
  
  return handled;
}


/* Mouse move handler */
static gboolean onMouseMoveDotter(GtkWidget *window, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->state & GDK_BUTTON1_MASK)  /* left-drag */
    {
      /* If the dot-plot was clicked the selected coords will have changed. Perform required updates. */
      updateOnSelectedCoordsChanged(window);
      handled = TRUE;
    }

  return handled;
}


/* Get the factor to convert between nucleotide and display (e.g. peptide) coords, or 1 if no
 * conversion is necessary */
int getResFactor(DotterContext *dc, const gboolean horizontal)
{
  return horizontal && dc->blastMode == BLXMODE_BLASTX ? dc->numFrames : 1;
}


/* Move the given sequence coord by the given number of coords (which can be negative to move 
 * in the decreasing direction. 'horizontal' indicates whether it's the horizontal or vertical
 * sequence that we're modifying and 'reverse' indicates whether that sequence's scale is shown reversed. */
static void incrementCoord(GtkWidget *dotterWindow, 
			   DotterContext *dc,
			   int *coord, 
			   const gboolean reverse, 
			   const gboolean horizontal, 
			   const gboolean convertCoords,
			   const int numCoords)
{
  const int incValue = convertCoords ? numCoords * getResFactor(dc, horizontal) : numCoords;

  if (reverse)
    {
      *coord -= incValue;
    }
  else
    {
      *coord += incValue;
    }
  
  updateOnSelectedCoordsChanged(dotterWindow);
}


/* Handle Q key press (Ctrl-Q => close all windows) */
static gboolean onKeyPressQ(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    dotterContextCloseAllWindows(dotterGetContext(dotterWindow));
  
  return ctrlModifier;
}

/* Handle H key press (Ctrl-H => show help dialog) */
static gboolean onKeyPressH(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showHelpDialog(dotterWindow);
  
  return ctrlModifier;
}

/* Handle P key press (Ctrl-P => print) */
static gboolean onKeyPressP(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  /* to do: implement this */
  return ctrlModifier;
}


/* Handle S key press (Ctrl-S => show settings dialog) */
static gboolean onKeyPressS(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showSettingsDialog(dotterWindow);
    
  return ctrlModifier;
}

/* Handle G key press (Ctrl-G => show greyramp tool) */
static gboolean onKeyPressG(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showGreyrampTool(dotterWindow);
  
  return ctrlModifier;
}

/* Handle A key press (Ctrl-A => show alignment tool) */
static gboolean onKeyPressA(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showAlignmentTool(dotterWindow);
  
  return ctrlModifier;
}

/* Handle D key press (Ctrl-D => show dotplot) */
static gboolean onKeyPressD(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showDotterWindow(dotterWindow);
  
  return ctrlModifier;
}

/* Handle up/down key presses */
static gboolean onKeyPressUpDown(GtkWidget *dotterWindow, const gboolean isUp, const gboolean modifier)
{
  /* Increment/decrement the vertical (i.e. match) sequence coord */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const int numCoords = isUp ? -1 : 1;
  
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, numCoords);  
  return TRUE;
}

/* Handle left/right key presses */
static gboolean onKeyPressLeftRight(GtkWidget *dotterWindow, const gboolean isLeft, const gboolean modifier)
{
  /* Increment/decrement the horizontal (i.e. reference) sequence coord */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  const int numCoords = isLeft ? -1 : 1;
  
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, numCoords);
  return TRUE;
}

/* Handle comma/period key presses */
static gboolean onKeyPressCommaPeriod(GtkWidget *dotterWindow, const gboolean isComma, const gboolean modifier)
{
  /* Increment/decrement both sequence coords */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, isComma ? -1 : 1);
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, isComma ? -1 : 1);
  return TRUE;
}

/* Handle left/right square bracket presses */
static gboolean onKeyPressLeftRightBracket(GtkWidget *dotterWindow, const gboolean isLeft, const gboolean modifier)
{
  /* Increment the horizontal and decrement the vertical, or vice versa */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, isLeft ? -1 : 1);
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, isLeft ? 1 : -1);
  return TRUE;
}


/* Key presses applicable to all windows */
gboolean onKeyPressDotter(GtkWidget *widget, GdkEventKey *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  
  const gboolean ctrlModifier = (event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK;	
  
  switch (event->keyval)
    {
      case GDK_Q:   /* fall through */
      case GDK_q:   handled = onKeyPressQ(dotterWindow, ctrlModifier);		    break;

      case GDK_H:   /* fall through */
      case GDK_h:   handled = onKeyPressH(dotterWindow, ctrlModifier);		    break;

      case GDK_P:   /* fall through */
      case GDK_p:   handled = onKeyPressP(dotterWindow, ctrlModifier);		    break;

      case GDK_S:   /* fall through */
      case GDK_s:   handled = onKeyPressS(dotterWindow, ctrlModifier);		    break;

      case GDK_G:   /* fall through */
      case GDK_g:   handled = onKeyPressG(dotterWindow, ctrlModifier);		    break;

      case GDK_A:   /* fall through */
      case GDK_a:   handled = onKeyPressA(dotterWindow, ctrlModifier);		    break;

      case GDK_D:   /* fall through */
      case GDK_d:   handled = onKeyPressD(dotterWindow, ctrlModifier);		    break;

      default: break;
  }
  
  return handled;
}


/* Key press handler for dotter windows that move the selected coords (including main window and
 *  alignment tool.). */
gboolean onKeyPressDotterCoords(GtkWidget *widget, GdkEventKey *event, gpointer data)
{
  gboolean handled = onKeyPressDotter(widget, event, data);
  
  if (!handled)
    {
      GtkWidget *dotterWindow = GTK_WIDGET(data);

      const gboolean shiftModifier = (event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK;
      //  const gboolean altModifier = (event->state & GDK_MOD1_MASK) == GDK_MOD1_MASK;
      
      switch (event->keyval)
	{
	  case GDK_Up:            handled = onKeyPressUpDown(dotterWindow, TRUE, shiftModifier);        break;
	  case GDK_Down:          handled = onKeyPressUpDown(dotterWindow, FALSE, shiftModifier);       break;
	    
	  case GDK_Left:          handled = onKeyPressLeftRight(dotterWindow, TRUE, shiftModifier);     break;
	  case GDK_Right:         handled = onKeyPressLeftRight(dotterWindow, FALSE, shiftModifier);    break;
	    
	  case GDK_comma:         handled = onKeyPressCommaPeriod(dotterWindow, TRUE, shiftModifier);   break;
	  case GDK_less:          handled = onKeyPressCommaPeriod(dotterWindow, TRUE, shiftModifier);   break;
	  case GDK_period:        handled = onKeyPressCommaPeriod(dotterWindow, FALSE, shiftModifier);  break;
	  case GDK_greater:       handled = onKeyPressCommaPeriod(dotterWindow, FALSE, shiftModifier);  break;
	    
	  case GDK_bracketleft:   handled = onKeyPressLeftRightBracket(dotterWindow, TRUE, shiftModifier);	  break;
	  case GDK_braceleft:     handled = onKeyPressLeftRightBracket(dotterWindow, TRUE, shiftModifier);	  break;
	  case GDK_bracketright:  handled = onKeyPressLeftRightBracket(dotterWindow, FALSE, shiftModifier);	  break;
	  case GDK_braceright:    handled = onKeyPressLeftRightBracket(dotterWindow, FALSE, shiftModifier);	  break;
	    
	  default: break;
	}
    }
  
  return handled;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Create the UI manager for the menus */
static GtkUIManager* createUiManager(GtkWidget *window, const DotterHspMode hspMode)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  /* If the HSPs are on initially then we're in hspOnly mode, so we don't show the pixmap; therefore,
   * set the toggled state of the pixelmap menu option (first in the toggleMenuEntries list) to be false
   * if HSPs are on. */
  toggleMenuEntries[0].is_active = (hspMode == 0);  
  
  gtk_action_group_add_actions(action_group, menuEntries, G_N_ELEMENTS (menuEntries), window);
  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS (toggleMenuEntries), window);
  gtk_action_group_add_radio_actions(action_group, radioMenuEntries, G_N_ELEMENTS (radioMenuEntries), hspMode, G_CALLBACK(onToggleHspMode), window);
  
  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);
  
  return ui_manager;
}


/* Create a menu. Optionally add it to the given menu bar, if menuBar is not null, with the given label */
static GtkWidget* createDotterMenu(GtkWidget *window, 
                                   const char *menuDescription, 
                                   const char *menuLabel, 
                                   const char *path, 
                                   GtkMenuBar *menuBar,
                                   GtkUIManager *ui_manager)
{
  GError *error = NULL;
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, menuDescription, -1, &error))
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  GtkWidget *menu = gtk_ui_manager_get_widget (ui_manager, path);
  
  if (menuBar)
    {
      GtkWidget *menuItem = gtk_menu_item_new_with_mnemonic(menuLabel);
      gtk_menu_item_set_submenu(GTK_MENU_ITEM(menuItem), menu);
      gtk_menu_bar_append(menuBar, menuItem);
    }
  
  return menu;
}


/* Create the dotter menu bar. */
static GtkWidget *createDotterMenuBar(GtkWidget *window, const DotterHspMode hspMode, GtkUIManager *uiManager)
{
  GtkWidget *menuBar = gtk_menu_bar_new();
  
  createDotterMenu(window, fileMenuDescription, "_File", "/File", GTK_MENU_BAR(menuBar), uiManager);
  createDotterMenu(window, editMenuDescription, "_Edit", "/Edit", GTK_MENU_BAR(menuBar), uiManager);
  createDotterMenu(window, viewMenuDescription, "_View", "/View", GTK_MENU_BAR(menuBar), uiManager);
  createDotterMenu(window, helpMenuDescription, "_Help", "/Help", GTK_MENU_BAR(menuBar), uiManager);
  
  return menuBar;
}


static GtkWidget* createDotterWindow(DotterContext *dc, 
				     DotterWindowContext *dwc,
                                     const DotterHspMode hspMode, 
                                     GtkWidget *greyrampTool, 
                                     GtkWidget *dotplotContainer, 
                                     GtkWidget *dotplot)
{
  DEBUG_ENTER("createDotterWindow");

  GtkWidget *dotterWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  char *title = blxprintf("Dotter %s vs. %s", dc->refSeqName, dc->matchSeqName);
  gtk_window_set_title(GTK_WINDOW(dotterWindow), title);
  g_free(title);
  
  /* Set the message handlers again, this time passing the window and statusbar, now we know them */
  BlxMessageData *msgData = g_malloc(sizeof *msgData);
  msgData->parent = GTK_WINDOW(dotterWindow);
  msgData->statusBar = NULL;
  
  g_log_set_default_handler(defaultMessageHandler, msgData);
  g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, popupMessageHandler, msgData);
  
  /* Create the menu bar, and a right-click context menu */
  GtkUIManager *uiManager = createUiManager(dotterWindow, hspMode);
  GtkWidget *menuBar = createDotterMenuBar(dotterWindow, hspMode, uiManager);
  GtkWidget *contextMenu = createDotterMenu(dotterWindow, contextMenuDescription, "Context", "/Context", NULL, uiManager);
  
  /* We'll put everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(dotterWindow), GTK_WIDGET(vbox));

  gtk_box_pack_start(GTK_BOX(vbox), menuBar, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), dotplotContainer, TRUE, TRUE, 0);

  
  gtk_widget_add_events(dotterWindow, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(dotterWindow, GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(dotterWindow), "key-press-event", G_CALLBACK(onKeyPressDotterCoords), dotterWindow);
  g_signal_connect(G_OBJECT(dotterWindow), "button-press-event", G_CALLBACK(onButtonPressDotter), contextMenu);
  g_signal_connect(G_OBJECT(dotterWindow), "motion-notify-event", G_CALLBACK(onMouseMoveDotter), NULL);
  
  /* Set the default window size based on the dotplot/exon widget size, up to a max based on screen size */
  GdkScreen *screen = gtk_widget_get_screen(dotterWindow);
  const int maxWidth = gdk_screen_get_width(screen) * MAX_WINDOW_WIDTH_FRACTION;
  const int maxHeight = gdk_screen_get_height(screen) * MAX_WINDOW_HEIGHT_FRACTION;
  
  int width = dotplotGetImageWidth(dotplot) + 100;
  int height = dotplotGetImageHeight(dotplot) + 100;
  height += 2 * (DEFAULT_EXON_HEIGHT + (2 * DEFAULT_EXON_YPAD));
  
  width = min(width, maxWidth);
  height = min(height, maxHeight);
  
  gtk_window_set_default_size(GTK_WINDOW(dotterWindow), width, height);
  
  gtk_widget_show_all(dotterWindow);
  
  DEBUG_EXIT("createDotterWindow returning ");
  return dotterWindow;
}




/**************************** eof ***************************/