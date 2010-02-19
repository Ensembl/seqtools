/*

   BLIXEM - BLast matches In an X-windows Embedded Multiple alignment

 -------------------------------------------------------------
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
 * -------------------------------------------------------------
   |  File: blxview.c                                          |
   |  Author: Erik.Sonnhammer                                  |
   |  Copyright (C) E Sonnhammer 1993-1997                     |
   -------------------------------------------------------------
 *
 * Exported functions: See wh/blxview.h
 * HISTORY:
 * Last edited: Sep 10 16:16 2009 (edgrif)
 * * Jan 10 10:35 2002 (edgrif): Fix up socket code and add various
 *              options for better sequence display.
 * previous mods:
  Date      Modification
--------  ---------------------------------------------------
93-02-20  Created
93-05-25  Dispstart/dispend fix by Richard for seq's < displen.
93-07-24  All boxes of a protein turn black when one is picked
          Sorting by protein name or score added
93-11-16  Added picking of Big Picture HSP's and Reverse Strand Big Picture
93-11-17  Added blastn support
94-01-29  Added Highlight sequences by names matching regexp
94-03-07  Fixed window limits. Added initial settings (BigPict, gotoNext)
94-03-08  Added 'always-both-strands' in blastn mode.
94-03-27  Added Tblastn support (works in seqbl mode).
94-12-01  [2.0] Added dotter calling.
95-02-02  Rewrote auxseq and padseq allocation to be fully dynamic and unlimited.
95-02-06  Added Tblastx support (works in seqbl mode).
95-06-01  [2.1] Added entropy display
95-06-23  [2.1] Added Settings window
95-07-21  2.1 announced--------------------------------------
95-08-01  [2.2] Initial Sorting mode on command line.
95-09-15  [2.2] Added Settings pull-down menu for faster manipulation.
95-09-29  [2.2] Improved WWW browser finding with findCommand() - doesn't get fooled so easily.
95-10-04  [2.2] Added "Print whole alignment"
95-10-27  [2.2] Added acedb-fetching at double clicking.
          BLIXEM_FETCH_ACEDB makes this default.
          Reorganised Settings window to Toggles and Menus rows.
95-11-01  [2.3] Changed command line syntax to "-" for piping.
          Added X options capability (-acefont, -font).
96-01-05  [2.3] Force all tblastn HSP's to be qframe +1, to harmonize with MSPcrunch 1.4
          which gives sframe for tblastn (otherwise the output would be dead).
96-02-08  [2.3] Added option -S "Start display at position #" to stand-alone command line.
96-02-08  [2.3] Added checkmarks to pull-down settings menu.
96-03-05  2.3 Announced.
96-03-08  [2.4] Only look for WWW browser once.
96-05-09  [2.4] Proper grayscale print colours.
96-05-09  [2.4] Enabled piping of query sequence too, for Pepmap and WWW calls.
96-08-20  [2.4] Fixed minor bug in squashed mode and added restoring of previous sorting after squash.
97-05-28  [2.4] Fixed parsing to handle gapped matches.
                Added "Highlight lower case residues" for gapped alignments and
                "Show sequence descriptions" (for MSPcrunch 2.1).
		Added setting the color of matching residues in the Settings window.
97-06-17  [2.4] Fixed "Highlight differences" for gapped alignments ('.' -> '-').
                Changed "Highlight lower case residues" to "Highlight subject insertions" and
		set this automatically for gapped alignments.  Works for both lower case and number
		insert markers.
                Changed blviewRedraw to use strlen to accommodate reverse gapped alignments.
		Simplified (and thereby debugged) selection of Big Picture MSPs to be drawn.
		Made Big Picture ON/OFF control more logical and consistent.
		Added a calcID() step to fix sortById() at startup.
		Added "Hilight Upper/Lower case" - useful general purpose feature.
97-10-14  [2.4] Added Description box when picking sequences.
97-03-07  [3.0] Added code for FS Feature Segment data. (shared code with Dotter for
                control and parsing; the display code is unique to Blixem).
                Added Inverted sorting order.
		Fixed bug in coordinate display in tblastn mode.
99-07-07  [3.0] Added msp->shape attribute. Added support for XY curve shapes PARTIAL and INTERPOLATE.
                Overhauled selective drawing of BigPicture MSPs, now simple enough to be bugfree (hopefully).
01-10-05	Added getsseqsPfetch to fetch all missing sseqs in one go via socket connection to pfetch [RD]

 * Created: Thu Feb 20 10:27:39 1993 (esr)
 * CVS info:   $Id: blxview.c,v 1.17 2010-02-19 16:22:03 gb10 Exp $
 *-------------------------------------------------------------------
 */


/*
		Pending:

		        Do exons and introns like SFS segments (i.e. eliminate
			magic scores.  Requires changes to fmapfeatures.c).

Known bugs:
-----------

revcomp broken when called from acedb.  Slen/send problem?

MSP score codes:
----------------
-1 exon                  -> Big picture + Alignment
-2 intron                -> Big picture + Alignment
-3 Any coloured segment  -> Big picture
-4 stringentSEGcolor     -> Big picture
-5 mediumSEGcolor        -> Big picture
-6 nonglobularSEGcolor   -> Big picture
-10 hidden by hand
*/

#include <wh/aceio.h>
#include <wh/key.h>
#include <wh/menu.h>
#include <SeqTools/dotter.h>
#include <wh/dict.h>

#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <gdk/gdkkeysyms.h>

#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blxviewMainWindow.h>

#ifdef ACEDB
#include <wh/display.h>
#endif


/* This is the only place this is set so that you get the same version/compile whether this is
 * compiled stand alone or as part of xace. */
char *blixemVersion = BLIXEM_VERSION_COMPILE ;



/* get rid of these...ugh.....horrible.....horrific.... */
extern void externalCommand (char *command);
//extern int pickMatch (char *cp, char *tp);


typedef struct _BPMSP
{
  char sname[FULLNAMESIZE+1];
  char *desc;
  int box;
//  Graph graph;
  struct _BPMSP *next;
} BPMSP;



#define max(a,b) (((a) > (b)) ? (a) : (b))
#define BPoffset 4
#define NAMESPACE 12
#define SEQ2BP(i) (float)plusmin*(i-BigPictStart-qoffset)*BPx/BigPictLen + BPoffset
#define MAXALIGNLEN 10000
#define separatorwidth 0.5

#define autoDotterParamsStr "Automatic Dotter parameters"
#define SortByScoreStr      "Sort by score"
#define SortByIdStr         "Sort by identity"
#define SortByNameStr       "Sort by name"
#define SortByPosStr        "Sort by position"
#define SortInvStr          "Inverted sorting order"
#define BigPictToggleStr    "Big Picture"
#define BigPictToggleRevStr "Big Picture Other Strand"
#define toggleIDdotsStr     "Highlight differences"
#define squashMatchesStr    "Squash matches"
#define squashFSStr         "Squash features"
#define entropytoggleStr    "Complexity curves"
#define printColorsStr      "B/W Print colours"
#define toggleColorsStr     "No colours"
#define toggleVerboseStr    "Verbose mode"
#define toggleHiliteSinsStr  "Highlight subject insertions"
#define toggleHiliteUpperStr "Highlight upper case"
#define toggleHiliteLowerStr "Highlight lower case"
#define toggleDESCStr        "Show sequence descriptions"
#define toggleMatchPasteStr  "Paste Match Set\t\t\t\tm"
#define toggleMatchClearStr  "Clear Match Set\t\t\t\tm"

//static void blxDestroy(void) ;
//static void blxPrint(void) ;
//static void wholePrint(void) ;

//static void blxPaste(BlxPasteType paste_type) ;
//static void pasteCB(char *text) ;
BOOL parsePasteText(char *text, BlxPasteData paste_data) ;
void freePasteData(BlxPasteData paste_data) ;
static char *parseMatchLines(char *text) ;
static BOOL parseFeatureLine(char *line,
			     char **feature_name_out,
			     int *feature_start_out, int *feature_end_out, int *feature_length_out) ;

//static BOOL setMatchSet(char **matches) ;
//static void clearMatchSet(void) ;

//static void sortToggleInv(void) ;

//static void squashMatches(void) ;
//static void squashFSdo(void) ;

//static void toggleIDdots (void) ;
//static void toggleVerbose(void) ;
//static void toggleHiliteSins(void) ;
//static void toggleHiliteUpper(void) ;
//static void toggleHiliteLower(void) ;
//static void toggleDESC(void) ;
//static void ToggleStrand(void) ;
//static void printColors (void) ;

//static BOOL gotoMatchPosition(char *match, int q_start, int q_end) ;

#if OLD_BLIXEM
static void keyboard(int key, int modifier) ;
static void blviewPick (int box, double x_unused, double y_unused, int modifier_unused);
static void blviewDestroy(GtkWidget *unused) ;
#endif

//static void toggleColors (void);
static void blviewCreate(char *opts, char *align_types, MSP *msplist, char *refSeq, char *refSeqName, const int qOffset, const int startCoord, const SortByType sortByType, const gboolean sortInverted, const gboolean gappedHsp) ;

#if OLD_BLIXEM
static char *get3rd_base(int start, int end, char *q);
#endif

static BOOL haveAllSequences(const MSP const *msplist, DICT *dict) ;
//static void getsseq(MSP *msp, MSP *msplist) ;
//static char *getSeq(char *seqname, char *fetch_prog) ;

//static BOOL smartDotterRange(char *selected_sequence, const MSP const *msp_list, int blastn,
//			     int strand_sign, int view_start, int view_end,
//			     char **dottersseq_out, int *dotter_start_out, int *dotter_end_out) ;

//static char *abbrevTxt(char *text, int max_len) ;

//static void printMSPs(void) ;

//static void toggleMatchSet(void) ;


/*
 *                 Local globals....sigh....
 */

//static void//     BigPictToggle(void),
//            entropytoggle(void),
//            BigPictToggleRev(void),
//            setHighlight(void),
 //           clrHighlight(void),
  /*             dotterPanel(void), */
//            allocAuxseqs(int len),
//            settingsRedraw(void),
//            menuCheck(MENU menu, int mode, int thismode, char *str),
//	    hidePicked(void),
//            setMenuCheckmarks(void);


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
/* currently unused... */
static void pfetchWindow (MSP *msp);
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



/* GLOBAL VARIABLES...really its hard to beat this for a list of globals,
 * one for the Guinness Book of Records perhaps... */

static BPMSP BPMSPlist;// *bpmsp;
//static int  fetchCount;
GtkWidget *blixemWindow = NULL ;
//static GtkWidget *messageWidget;

#if OLD_BLIXEM
//static Graph frameGraph[3];
//static Graph frameGraphQ[3];
#endif

//static Graph settingsGraph=0;

static int   actstart,        /* Active region coordinates relative to query in BASES */
//             actend,
             dispstart;       /* Displayed region relative to query in BASES */
//             displen=240;     /* Displayed sequence length in BASES */
static char  actframe[16]="(+1)";    /* Active frame */
static int   plusmin = 1 ;       /* +1 for top strand, -1 */
static float //queryy,          /* Screen coords of genome seq */
             //separator_y,	/* y coord of previous panel separator */
             //lastExonx,
             //BPboxheight = 5.7,
             //BPboxstart,
             //BPboxwidth,
             BigPictZoom = 10;
//  oldLinew;


//static char *q;    /* The genomic sequence (query=the translated seq) */
//static int   qlen;
//static int   qoffset;         /* Offset to the 'real' DNA start */
//static char *qname_G = NULL ;
//       MSP  *MSPlist,          /* List of MSP's */
//            *msp;
//static MSP  *pickedMSP = NULL ;	      /* Last picked MSP */
static char 
//	    message[1024],
             HighlightSeq[LONG_NAMESIZE+4] = "",
//             searchSeq[NAMESIZE+4] = "",
//            *cp,
            *padseq = 0;
//            *auxseq = 0,
//            *auxseq2 = 0,
//            *dottersseq = NULL,
//  dotterqname[LONG_NAMESIZE] = "",
//             stringentEntropytx[10],
//             mediumEntropytx[10],
//             nonglobularEntropytx[10],
//  sortModeStr[32] = "Identity" ;


//static BlxPasteType paste_type_G = BLXPASTE_INVALID ;

//static int   lastbox = 0;
//static int   colortoggle = 0;
//static int   backgColor = LIGHTGRAY,
//             IDcolor = CYAN,
//             consColor = MIDBLUE,
//             geneColor = BLUE,
//             hiColor = YELLOW,
//             matchSetColor = CERISE,
//	       gridColor = YELLOW,
//             oldcolor;
static int   BigPictON = 1,
             BigPictRev = 0,	/* Draw other strand in Big picture */
//             BigPictStart, BigPictLen, BPbox, BPx,
//             nx, ny,

             blastx = 1, blastp = 0, blastn = 0, tblastn = 0, tblastx = 0,

             symbfact = 3,
//             i,
             start_nextMatch = 0,
             dotter_first = 0,
//             IDdots = 0,
//             squash = 0,
//             squashFS = 1,
//             verbose = 0,
             HiliteSins = 0,
//             HiliteUpperOn = 0,
//             HiliteLowerOn = 0,
//             DESCon = 0,
//             dotterZoom = 0,
//             dotterStart = 0,
//             dotterEnd = 0,
//             dotterHSPs = 0,
//             auxseqlen = 0,
//             smartDotter = 1,
//             entropyOn = 0,
//             stringentEntropycolor = LIGHTGREEN,
//             stringentEntropybox,
//             mediumEntropycolor = GREEN,
//             mediumEntropybox,
//             nonglobularEntropycolor = DARKGREEN,
//             nonglobularEntropybox,
             alphabetsize,
//             stringentEntropywin = 12,
//             mediumEntropywin = 25,
//             nonglobularEntropywin = 45,
//             printColorsOn,
//             wholePrintOn,
//             oneGraph,
//             settingsButton,
             sortInvOn = 0,
             HSPgaps = 0;

static SortByType sortMode = SORTBYUNSORTED ;

static int oldWidth = 0, oldHeight = 0;

#if OLD_BLIXEM
static double oldx, DNAstep;
#endif

/* A stepping stone to having a blixem view context is this global struct. */
//static BlixemViewStruct blixem_context_G = {FALSE} ;



/*
 *                External routines
 */










/*
 *                Internal routines
 */

/* Temporary function while introducing a blixem view context. */
//static BlixemView getBlxViewContext(void)
//{
//  return &blixem_context_G ;
//}


//static void toggleMatchSet(void)
//{
//  BlixemView blxview = getBlxViewContext() ;
//  BOOL result ;
//
//  if (blxview->match_set)
//    {
//      clearMatchSet() ;
//      blxview->match_set = FALSE ;
//
//      result = menuSetLabel(menuItem(blixemMenu, toggleMatchClearStr), toggleMatchPasteStr) ;
//    }      
//  else
//    {
//      /* We cannot reset the menu yet because blxPaste() is asynchronous so we don't know if it worked. */
//      blxPaste(BLXPASTE_MATCHSET) ;
//    }
//
//
//
//  return ;
//}



//static void toggleVerbose(void)
//{
//  verbose = (verbose ? 0 : 1);
//
//  if (verbose)
//    printMSPs() ;
//
//  blviewRedraw();
//
//  return ;  
//}

//static void toggleHiliteSins(void)
//{
//    HiliteSins = (HiliteSins ? 0 : 1);
//    blviewRedraw();
//}

//static void toggleHiliteUpper(void)
//{
//    HiliteUpperOn = (HiliteUpperOn ? 0 : 1);
//    blviewRedraw();
//}

//static void toggleHiliteLower(void) {
//    HiliteLowerOn = (HiliteLowerOn ? 0 : 1);
//    blviewRedraw();
//}

//static void toggleDESC(void) {
//    DESCon = (DESCon ? 0 : 1);
//    blviewRedraw();
//}

/* AGH...strand probably not good..... */
//static BOOL gotoMatchPosition(char *match, int q_start, int q_end)
//{
//  BOOL result = FALSE ;
//  MSP *msp ;
//
//  for (msp = MSPlist; msp ; msp = msp->next)
//    {
//      if (g_ascii_strcasecmp(msp->sname, match) == 0)
//	{
//	  if (msp->qstart == (q_start - qoffset) && msp->qend == (q_end - qoffset))
//	    {
//	      dispstart = q_start - qoffset ;
//
//	      blviewRedraw() ;
//
//	      result = TRUE ;
//
//	      break ;
//	    }
//	}
//    }
//
//  return result ;
//}


//
//#if OLD_BLIXEM
//static void keyboard(int key, int modifier)
//{
//  switch (key)
//    {
//    case '<':
//    case ',':
//      scrollLeft1();
//      break;
//    case '>':
//    case '.':
//      scrollRight1();
//      break;
//
//    case UP_KEY:
//      blviewPick(lastbox - 1, 0, 0, 0);
//      break;
//    case DOWN_KEY:
//      blviewPick(lastbox + 1, 0, 0, 0);
//      break;
//
//    case GDK_G:
//    case GDK_g: 
//      {
//	blxPaste(BLXPASTE_MOVETO) ;
//	break ;
//      }
//
//    case GDK_M:
//    case GDK_m: 
//      {
//	toggleMatchSet() ;
//	break ;
//      }
//
//    default:
//      break;
//    }
//
//  return ;
//}
//#endif


//static void toggleColors (void)
//{
 //   static int oldback, oldgrid, oldID, oldcons, oldgene, oldhi;
//
////    graphActivate(settingsGraph);
//
//    if (!colortoggle) {
//	oldback = backgColor; backgColor = WHITE;
//	oldgrid = gridColor; gridColor = BLACK;
//	oldID = IDcolor; IDcolor = WHITE;
//	oldcons = consColor; consColor = WHITE;
//	oldgene = geneColor; geneColor = BLACK;
//	oldhi = hiColor; hiColor = WHITE;
//	colortoggle = 1;
//    }
//    else {
//	backgColor = oldback;
//	gridColor= oldgrid;
//	IDcolor = oldID;
//	consColor = oldcons;
//	geneColor = oldgene;
//	hiColor = oldhi;
//	colortoggle = 0;
//    }
//    blviewRedraw();
//}


//static void printColors (void)
//{
//    static int oldback, oldgrid, oldID, oldcons, oldgene, oldhi;
//
////    graphActivate(settingsGraph);
//
//    if (!printColorsOn) {
//	oldback = backgColor; backgColor = WHITE;
//	oldgrid = gridColor; gridColor = LIGHTGRAY;
//	oldID = IDcolor; IDcolor = GRAY;
//	oldcons = consColor; consColor = PALEGRAY;
//	oldgene = geneColor; geneColor = BLACK;
//	oldhi = hiColor; hiColor = LIGHTGRAY;
//	printColorsOn = 1;
//    }
//    else {
//	backgColor = oldback;
//	gridColor= oldgrid;
//	IDcolor = oldID;
//	consColor = oldcons;
//	geneColor = oldgene;
//	hiColor = oldhi;
//	printColorsOn = 0;
//    }
//    blviewRedraw();
//}


//static void toggleIDdots (void)
//{
//    IDdots = !IDdots;
//    blviewRedraw();
//}


//static void calcEntropyArray(Array array, int win)
//{
//    int i, j, *rescount;
//    float pi, sum;
//
//    rescount = (int *)g_malloc(24*sizeof(int));
//
//    for (i=0; i < qlen; i++) {
//	rescount[aa_atob[(unsigned int)q[i]]]++;
//	if (i+1 >= win) {
//	    for (sum = j = 0; j < 24; j++)
//		if (rescount[j]) {
//		    pi = (float)rescount[j]/win;
//		    sum += pi*log(pi);
//		}
//	    arr(array, i-win/2, float) = -sum/log(2);
//	    rescount[aa_atob[(unsigned int)q[i+1-win]]]--;
//	}
//    }
//
//    g_free(rescount);

    /* TEST - delete later * /
    for (i=0; i < qlen; i++)
	printf ("%3d  %c  %f\n", i, q[i], arr(stringentEntropyarr, i, float));
	*/
//}

//static void calcEntropyArrays(BOOL force)
//{
//    /* force:
//       FALSE - only calculate if necessary, i.e. first call.
//       TRUE - force (re)calculation.
//       */
//
//    if (!stringentEntropyarr) {
//	calcEntropyArray(stringentEntropyarr = arrayCreate(qlen, float), stringentEntropywin);
//	calcEntropyArray(mediumEntropyarr = arrayCreate(qlen, float), mediumEntropywin);
//	calcEntropyArray(nonglobularEntropyarr = arrayCreate(qlen, float), nonglobularEntropywin);
//    }
//    else if (force) {
//	calcEntropyArray(stringentEntropyarr = arrayCreate(qlen, float), stringentEntropywin);
//	calcEntropyArray(mediumEntropyarr = arrayCreate(qlen, float), mediumEntropywin);
//	calcEntropyArray(nonglobularEntropyarr = arrayCreate(qlen, float), nonglobularEntropywin);
//    }
//}


//static void entropytoggle (void)
//{
//    entropyOn = !entropyOn;
//    if (entropyOn) BigPictON = 1;
//    calcEntropyArrays(FALSE);
//    blviewRedraw();
//}


/* Find the expression 'query' in the string 'text'
 * Return 1 if found, 0 otherwise
 */
//int strMatch(char *text, char *query)
//{
//    /* Non-ANSI bsd way: * /
//       if (re_exec(text) == 1) return 1;
//       else return 0;
//    */
//
//    /* ANSI way: */
//    return pickMatch(text, query);
//}
//
void highlightProteinboxes(BOOL warpScroll)
{//
//  MSP *msp;
//
//  /* Highlight alignment boxes of current search string seq */
//  if (*searchSeq)
//    for (msp = MSPlist; msp ; msp = msp->next)
//      if (msp->box && strMatch(msp->sname, searchSeq))
//	{
//	  graphActivate(msp->graph);
//	  graphBoxDraw(msp->box, BLACK, RED);
//	}
//
//  /* Highlight alignment boxes of currently selected seq */
//  if (!squash)
//    {
//      for (msp = MSPlist; msp ; msp = msp->next)
//	{
//	  if (msp->box && !strcmp(msp->sname, HighlightSeq))
//	    {
//	      float x1, x2, y1, y2;
//
//	      graphActivate(msp->graph);
//	      graphBoxDraw(msp->box, WHITE, BLACK);
//
//	      if (warpScroll)
//		{
//		  /* Scroll the alignment window so that the currently
//		     selected seq is visible. Not that this only happens
//		     in response to clicking on the big picture
//		     when warpScroll is TRUE. */
//		  graphBoxDim(msp->box, &x1, &y1, &x2, &y2);
//		  graphGoto(x1, y1);
//		}
//	    }
//	  }
//    }
//
//  if (BigPictON)
//    {
//      /* Highlight Big Picture boxes of current search string seq */
//      if (*searchSeq)
//	{
//	  for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
//	    {
//	      if (strMatch(bpmsp->sname, searchSeq))
//		{
//		  graphActivate(bpmsp->graph);
//		  graphBoxDraw(bpmsp->box, RED, BLACK);
//		}
//	    }
//	}
//
//      /* Highlight Big Picture boxes of currently selected seq */
//      for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
//	{
//	  if (!strcmp(bpmsp->sname, HighlightSeq))
//	    {
//	      graphActivate(bpmsp->graph);
//	      graphBoxDraw(bpmsp->box, CYAN, BLACK);
//	    }
//	}
//    }
//
//  return ;
}

//static void hidePicked (void)
//{
//  MSP *msp;
//
//  for (msp = MSPlist; msp ; msp = msp->next)
//    if (!strcmp(msp->sname, HighlightSeq)) {
//      msp->id = msp->score;
//      msp->score = -999;
//    }
//  blviewRedraw () ;
//}


/* Callback for when user clicks on a sequence to retrieve the EMBL entry
 * for that sequence. The method of retrieving the sequence can be changed
 * via environment variables.
 *                                                                           */
//#if OLD_BLIXEM
//static void blviewPick(int box, double x_unused, double y_unused, int modifier_unused)/
//{
  //MSP *msp;
////  Graph origGraph = graphActive();
//
//  if (!box)
//    return ;
//
//  if (box == lastbox)
//    {
//      /* Second click - efetch this seq */
//      for (msp = MSPlist; msp ; msp = msp->next)
//	{
//	  if (msp->box == box && msp->graph == graphActive())
//	    break;
//	}
//
//      if (msp && *msp->sname)
//	{
//	  blxDisplayMSP(msp) ;
//	}
//    }
//  else
//    {
//      /* Reset all highlighted boxes */
//      if (!squash)
//	{
//	  for (msp = MSPlist; msp ; msp = msp->next)
//	    {
//	      if (msp->box && !strcmp(msp->sname, HighlightSeq))
//		{
//		  graphActivate(msp->graph);
//		  graphBoxDraw(msp->box, BLACK, backgColor);
//		}
//	    }
//	}
//
//      if (BigPictON)
//	{
//	  for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
//	    {
//	      if (!strcmp(bpmsp->sname, HighlightSeq))
//		{
//		  graphActivate(bpmsp->graph);
//		  graphBoxDraw(bpmsp->box, BLACK, backgColor);
//		}
//	    }
//	}
//
//      /* Find clicked protein  ***********************/
//
////      for (msp = MSPlist; msp ; msp = msp->next)
////	{
////	  if (msp->box == box && msp->graph == origGraph)
////	    break ;
////	}
//
//      if (msp)
//	{
//	  /* Picked box in alignment */
//	  strcpy(HighlightSeq, msp->sname) ;
//	  pickedMSP = msp ;
//	}
//      else if (BigPictON)
//        {
//	  /* look for picked box in BigPicture */
////	  for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
////	    {
////	      if (bpmsp->box == box && bpmsp->graph == origGraph)
////		break;
////	    }
//
//	  if (bpmsp && *bpmsp->sname)
//	    strcpy(HighlightSeq, bpmsp->sname);
//	}
//
//      if (msp || (bpmsp && *bpmsp->sname))
//	{
//	  char *sname ;
//	  char *desc ;
//	  BOOL highlight_protein ;
//	  
//
//	  /* Put description in message box and in paste buffer. */
//	  if (msp)
//	    {
//	      sname = msp->sname ;
//	      desc = msp->desc ;
//	      highlight_protein = FALSE ;
//	    }
//	  else
//	    {
//	      sname = bpmsp->sname ;
//	      desc = bpmsp->desc ;
//	      highlight_protein = TRUE ;
//	    }
//
//	  strncpy(message, sname, 1023) ;
//	  if (desc)
//	    {
//	      strcat(message, " ") ;
//	      strncat(message, desc, 1023-strlen(message)) ;
//	    }
//
//	  graphPostBuffer(sname) ;
//	  
//	  /* Highlight picked box */
////	  graphActivate(blixemGraph) ;
//	  highlightProteinboxes(highlight_protein) ;
//
//	  gtk_entry_set_text(GTK_ENTRY(messageWidget), message) ;
//	  lastbox = box ;
//	}
//    }
//
//  return ;
//}
//#endif


//static void mspcpy(MSP *dest, MSP *src)
//{
//    dest->type   = src->type;
//    dest->score  = src->score;
//    dest->id     = src->id;
//
//    strcpy(dest->qframe, src->qframe);
//    dest->qstart = src->qstart;
//    dest->qend   = src->qend;
//
//    dest->sname  = src->sname;
//    strcpy(dest->sframe, src->sframe);
//    dest->sstart = src->sstart;
//    dest->send   = src->send;
//    dest->sseq   = src->sseq;
//
//    dest->desc   = src->desc;
//    dest->box    = src->box;
//
//    dest->color  = src->color;
//    dest->shape  = src->shape;
//    dest->fs     = src->fs;
//
//    dest->xy     = src->xy;
//
//    dest->gaps   = src->gaps;
//
//#ifdef ACEDB
//    dest->key    = src->key;
//#endif
//}


//static void sortMSPs(int (*func)())
//{
//  MSP tmpmsp, *msp1, *msp2;
//
//  if (!MSPlist)
//    return;
//
//  for (msp1 = MSPlist ; msp1->next ; msp1 = msp1->next )
//    {
//      for (msp2 = msp1->next ; msp2 ; msp2 = msp2->next )
//	{
//	  if ( (*func)(msp1, msp2)*(sortInvOn ? -1 : 1) > 0 )
//	    {
//	      mspcpy(&tmpmsp, msp2);
//	      mspcpy(msp2, msp1);
//	      mspcpy(msp1, &tmpmsp);
//	    }
//	}
//    }
//
//  return ;
//}


/* aghhh, this all needs rewriting, how opaque can you get sigh... */

#define FSpriority {if (FS(msp1) && !FS(msp2)) return -1; else if (FS(msp2) && !FS(msp1)) return 1;}

//static int possort(MSP *msp1, MSP *msp2) {
//    FSpriority
//    return ( (msp1->qstart > msp2->qstart) ? 1 : -1 );
//}
//
//static int namesort(MSP *msp1, MSP *msp2) {
//    FSpriority
//    if (strcmp(msp1->sname, msp2->sname))
//	return strcmp(msp1->sname, msp2->sname);
//    else
//	return possort(msp1, msp2);
//}
//static int scoresort(MSP *msp1, MSP *msp2) {
//    FSpriority
//    return ( (msp1->score < msp2->score) ? 1 : -1 );
//}
//
//static int idsort(MSP *msp1, MSP *msp2)
//{
//  int result = 0 ;
//
//  FSpriority
//
//  result = ( (msp1->id < msp2->id) ? 1 : -1 ) ;
//
//  return result ;
//}


//static void sortToggleInv(void)
//{
//  sortInvOn = !sortInvOn;
//  switch (sortMode)
//    {
//    case SORTBYNAME : sortByName(); break;
//    case SORTBYSCORE : sortByScore(); break;
//    case SORTBYPOS : sortByPos(); break;
//    case SORTBYID : sortById(); break;
//    default: blviewRedraw(); break;
//    }
//
//  return ;
//}

/*
static void incBack(void) {
    backgColor++;
    if (!(backgColor % BLACK)) backgColor++;
    blviewRedraw();
}
static void decBack(void) {
    backgColor--;
    if (!(backgColor % BLACK)) backgColor--;
    blviewRedraw();
}
static void incGrid(void) {
    gridColor++;
    blviewRedraw();
}
static void decGrid(void) {
    gridColor--;
    blviewRedraw();
}
*/

//static void squashMatches(void)
//{
//    static int oldSortMode;
//
//    if (!squash) {
//      oldSortMode = sortMode;
//      sortByName();
//      squash = 1;
//    }
//    else {
//	switch (oldSortMode) {
//	case SORTBYNAME : sortByName(); break;
//	case SORTBYSCORE : sortByScore(); break;
//	case SORTBYPOS : sortByPos(); break;
//	case SORTBYID : sortById(); break;
//	}
//	squash = 0;
//    }
//
//    blviewRedraw();
//
//    return ;
//}


//static void squashFSdo(void)
//{
//    squashFS = !squashFS;
//    blviewRedraw();
//}


//static void wholePrint(void)
//{
   // int
//	tmp,
//	dispstart_save = dispstart,
//	BigPictON_save = BigPictON;
//
//    static int
//	start=1, end=0;
//    ACEIN zone_in;
//
//    if (!end) end = qlen;
//
//    /* Swap coords if necessary */
//    if ((plusmin < 0 && start < end) || (plusmin > 0 && start > end )) {
//	tmp = start;
//	start = end;
//	end = tmp;
//    }
//
//    /* Apply max limit MAXALIGNLEN */
//    if ((abs(start-end)+1) > MAXALIGNLEN*symbfact) {
//	start = dispstart - plusmin*MAXALIGNLEN*symbfact;
//	if (start < 1) start = 1;
//	if (start > qlen) start = qlen;
//
//	end = start + plusmin*MAXALIGNLEN*symbfact;
//	if (end > qlen) end = qlen;
//	if (end < 1) end = 1;
//    }
//
//    if (!(zone_in = messPrompt("Please state the zone you wish to print",
//			       messprintf("%d %d", start, end),
//			       "iiz", 0)))
//      return;
//
//    aceInInt(zone_in, &start);
//    aceInInt(zone_in, &end);
//    aceInDestroy (zone_in);
//
//    dispstart = start;
//    displen = abs(end-start)+1;
//
//    /* Validation */
//    if (plusmin > 0 && start > end) {
//	messout("Please give a range where from: is less than to:");
//	return;
//    }
//    else if (plusmin < 0 && start < end) {
//	messout("Please give a range where from: is more than to:");
//	return;
//    }
//    if (displen/symbfact > MAXALIGNLEN) {
//	messout("Sorry, can't print more than %d residues.  Anyway, think of the paper!", MAXALIGNLEN*symbfact);
//	return;
//    }
//
//    wholePrintOn = 1;
//    BigPictON = 0;
////    oneGraph = 1;
//
//    blviewRedraw();
//    graphPrint();
//
//    /* Restore */
//    wholePrintOn = 0;
//    dispstart = dispstart_save;
//    displen = dispstart_save;
////    oneGraph = 0;
//    BigPictON = BigPictON_save;
//
//    blviewRedraw();
//}

//static void blxPrint(void)
//{
//  oneGraph = 1;
//  blviewRedraw();
//  graphPrint();
//  oneGraph = 0;
//
//  blviewRedraw();
//}


//static void blxDestroy(void)
//{
//  gtk_widget_destroy(blixemWindow) ;
//
//  return ;
//}



/* 
 *       Set of functions to handle getting data from the clipboard following user actions.
 */

/* Called by menu/keyboard code. */
//static void blxPaste(BlxPasteType paste_type)
//{
//  paste_type_G = paste_type ;				    /* acedb callbacks force us to have a global. */
//
//  graphPasteBuffer(pasteCB) ;				    /* get clipboard data. */
//
//  return ;
//}


/* Called by graph asynchronously (when it gets the "select" signal) with data
 * taken from the windowing systems clipboard. */
//static void pasteCB(char *text)
//{
//  BlixemView blxview = getBlxViewContext() ;
//  BlxPasteDataStruct paste_data = {BLXPASTE_INVALID} ;
//
//  paste_data.type = paste_type_G ;
//
//  if (!parsePasteText(text, &paste_data))
//    {
//      messerror("Could not paste from clipboard, unknown format: \"%s\"", (text ? text : "no data on clipboard")) ;
//    }
//  else
//    {
//      if (paste_type_G && (paste_type_G != paste_data.type))
//	{
//	  messerror("Wrong type of clipboard data, expected %d but got %d", paste_type_G, paste_data.type) ;
//	}
//      else
//	{
//	  switch (paste_data.type)
//	    {
//	    case BLXPASTE_MATCHSET:
//	      {
//		if ((blxview->match_set = setMatchSet(paste_data.data.match_set)))
//		  {
//		    menuSetLabel(menuItem(blixemMenu, toggleMatchPasteStr), toggleMatchClearStr) ;
//		  }
//
//		blviewRedraw() ;
//		break ;
//	      }
//	    case BLXPASTE_MOVETO:
//	      {
//		if (!gotoMatchPosition(paste_data.data.move_to.match,
//				       paste_data.data.move_to.start,
//				       paste_data.data.move_to.end))
//		  messerror("Failed to find/move to match \"%s\" at %d %d",
//			    paste_data.data.move_to.match,
//			    paste_data.data.move_to.start, paste_data.data.move_to.end) ;
//		break ;
//	      }
//	    default:
//	      {
//		printf("debug: data unrecognised !!\n") ;
//		break ;
//	      }
//	    }
//	}
//
//      freePasteData(&paste_data) ;
//    }
//
//  paste_type_G = BLXPASTE_INVALID ;			    /* Be sure to reset. */
//
//  return ;
//}


/* Parses clipboard text looking for one of several formats. */
BOOL parsePasteText(char *text, BlxPasteData paste_data)
{
  BOOL result = FALSE ;

  if (text && *text)
    {
      switch (paste_data->type)
	{
	case BLXPASTE_MATCHSET:
	  {
	    char *match_names ;

	    /* Look for a list of match names. */
	    if ((match_names = parseMatchLines(text))
		&& (match_names = g_strstrip(match_names))
		&& (paste_data->data.match_set = g_strsplit_set(match_names, " ", -1))) /* -1 means do all tokens. */
	      result = TRUE ;

	    break ;
	  }
	case BLXPASTE_MOVETO:
	  {
	    /* Look for a match with a start/end. */
	    char *match_name ;
	    int start = 0, end = 0, length = 0 ;
	  
	    if (parseFeatureLine(text, &match_name, &start, &end, &length))
	      {
		paste_data->data.move_to.match = match_name ;
		paste_data->data.move_to.start = start ;
		paste_data->data.move_to.end = end ;
		
		result = TRUE ;
	      }
	    else
	      result = FALSE ;

	    break ;
	  }
	default:
	  {
	    /* Simply return FALSE. */
	    result = FALSE ;
	    break ;
	  }
	}
    }

  return result ;
}

void freePasteData(BlxPasteData paste_data)
{
  switch (paste_data->type)
    {
    case BLXPASTE_MATCHSET:
      {
	g_strfreev(paste_data->data.match_set) ;

	break ;
      }
    case BLXPASTE_MOVETO:
      {
	g_free(paste_data->data.move_to.match) ;

	break ;
      }
    default:
      break ;
    }

  return ;
}

/* Uses parseLine() to look for match lines in text and then
 * extracts the match names into a space separated list. */
static char *parseMatchLines(char *text)
{
  char *match_names = NULL ;
  char **matches ;

  if ((matches = g_strsplit_set(text, "\n", -1)))		    /* -1 means do all tokens. */
    {
      GString *names ;
      char *match = *matches ;
      gboolean string_free = TRUE ;

      names = g_string_sized_new(100) ;

      while (matches && match)
	{
	  char *match_name ;
	  int start = 0, end = 0, length = 0 ;

	  if (parseFeatureLine(match, &match_name, &start, &end, &length))
	    {
	      g_string_append_printf(names, " %s", match_name) ;

	      g_free(match_name) ;
	    }

	  matches++ ;
	  match = *matches ;
	}

      if (names->len)
	{
	  match_names = names->str ;
	  string_free = FALSE ;
	}

      g_string_free(names, string_free) ;
    }

  return match_names ;
}



/* Parses a line of the form:
 * 
 *                       "\"name\" start end (length)"
 *
 * Returns TRUE if all elements there, FALSE otherwise.
 * 
 * We could verify "name" more but probably not worth it. */
static BOOL parseFeatureLine(char *line,
			     char **feature_name_out,
			     int *feature_start_out, int *feature_end_out, int *feature_length_out)
{
  BOOL result = FALSE ;
  int fields ;
  char sequence_name[1000] = {'\0'} ;
  int start = 0, end = 0, length = 0 ;
  char *format_str = "\"%[^\"]\"%d%d (%d)" ;
	  
  if ((fields = sscanf(line, format_str, &sequence_name[0], &start, &end, &length)) == 4)
    {
      if (sequence_name && (start > 0 && start < end) && length > 0)
	{
	  *feature_name_out = g_strdup(sequence_name) ;
	  *feature_start_out = start ;
	  *feature_end_out = end ;
	  *feature_length_out = length ;

	  result = TRUE ;
	}
    }

  return result ;
}





/* "Match set" is a subset of matches supplied by the user which will be highlighted
 * in the matchSetColor. Can be used to identify a subset for any number of
 * reasons, e.g. the set of matches used so far as evidence to support a transcript. */


/* Takes a list of match names and looks for those names in the msp list,
 * returns TRUE if any of the matches found. Error message reports matches
 * not found. */
//static BOOL setMatchSet(char **matches_in)
//{
//  BOOL result = FALSE ;
 // GString *not_found ;
//  char **matches = matches_in ;
//  char *match = *matches ;
//
//  not_found = g_string_sized_new(100) ;
//
//  while (matches && match)
//    {
//      MSP *msp ;
//      BOOL found = FALSE ;
//
//      for (msp = MSPlist; msp; msp = msp->next)
//	{
//	  if (g_ascii_strcasecmp(msp->sname, match) == 0)
//	    {
//	      result = found = msp->in_match_set = TRUE ;
//	    }
//	}
//
//      if (!found)
//	g_string_append_printf(not_found, " %s ", match) ;
//
//      matches++ ;
//      match = *matches ;
//    }
//
//  if (not_found->len)
//    messerror("Match setting: following matches not found in blixem: %s", not_found->str) ;
//
//  g_string_free(not_found, TRUE) ;
//
//  return result ;
//}


/* Reset to no match set. */
//static void clearMatchSet(void)
//{
//  MSP *msp ;
//
//  for (msp = MSPlist; msp; msp = msp->next)
//    {
//      if (msp->in_match_set)
//	msp->in_match_set = FALSE ;
//    }
//
//  blviewRedraw() ;
//
//  return ;
//}


static void setModeP(void)
{
  blastp = 1;
  blastx = blastn = tblastn = tblastx = 0;
  alphabetsize = 24;
  symbfact = 1;
  BigPictZoom = 10;
}

static void setModeN(void)
{
  blastn = 1;
  blastp = blastx = tblastn = tblastx = 0;
  alphabetsize = 4;
  symbfact = 1;
  BigPictZoom = 30;
}

static void setModeX(void) {
    blastx = 1;
    blastp = blastn = tblastn = tblastx = 0;
    alphabetsize = 4;
    symbfact = 3;
    BigPictZoom = 10;
}
static void setModeT(void) {
    tblastn = 1;
    blastp = blastx = blastn = tblastx = 0;
    alphabetsize = 24;
    symbfact = 1;
}
static void setModeL(void) {
    tblastx = 1;
    blastp = blastx = blastn = tblastn = 0;
    alphabetsize = 24;
    symbfact = 3;
    BigPictZoom = 10;
}


static void blxviewInitGlobals(char *seq, char *seqname, int start, int offset, const MSP *msplist)
{
//  q = seq;
//  qlen = actend = strlen(q) ;
//  qname_G = g_strdup(seqname) ;
  
  dispstart = start;
  actstart=1;
//  qoffset = offset;
//  MSPlist = msplist;
  BPMSPlist.next = 0;
  *HighlightSeq = 0;
  blastx = blastp =  blastn = tblastn = tblastx = 0 ;
  sortMode = SORTBYUNSORTED ;
}


static void blxviewGetOpts(char *opts, char *refSeq)
{
  char *opt = opts;
  while (*opt)
    {
      /* Used options: 	 BGILMNPRSTXZ-+brsinp      */
      
      switch (*opt)
      {
	case 'I':
	  sortInvOn = 1;                         break;
	case 'G':
	  /* Gapped HSPs */
	  HiliteSins = HSPgaps = 1;                 break;
	case 'P': setModeP();                       break;
	case 'N': setModeN();                       break;
	case 'X': setModeX();                       break;
	case 'T': setModeT();                       break;
	case 'L': setModeL();                       break;
	case '-':
	  strcpy(actframe, "(-1)");
	  plusmin = -1;                           break;
	case '+':
	  strcpy(actframe, "(+1)");
	  plusmin = 1;                            break;
	case 'B':
	  BigPictON = 1;                          break;
	case 'b':
	  BigPictON = 0;                          break;
	case 'd':
	  dotter_first = 1;                       break;
	case 'i':
	  sortMode = SORTBYID ;                   break;
	case 'M':
	  start_nextMatch = 1;                    break;
	case 'n':
	  sortMode = SORTBYNAME ;                 break;
	case 'p':
	  sortMode = SORTBYPOS ;                  break;
	case 'R':
	  BigPictRev = 1;                         break;
	case 'r':
	  BigPictRev = 0;                         break;
	case 's':
	  sortMode = SORTBYSCORE ;                break ;
	case 'Z':
	  BigPictZoom = strlen(refSeq);              break;
      }
      
      opt++;
    }
  
  if (blastx + blastn + blastp + tblastn + tblastx == 0)
    {
      printf("\nNo blast type specified. Detected ");
      
      if (Seqtype(refSeq) == 'P')
	{
	  printf("protein sequence. Will try to run Blixem in blastp mode\n");
	  setModeP();
	}
      else
	{
	  printf("nucleotide sequence. Will try to run Blixem in blastn mode\n");
	  setModeN();
	}
    }
}


/* Returns the blast mode, i.e. blastx, blastn etc. */
static BlxBlastMode getBlastMode()
{
  BlxBlastMode mode;
  
  if (blastx)
    mode = BLXMODE_BLASTX;
  else if (tblastx)
    mode = BLXMODE_TBLASTX;
  else if (blastn)
    mode = BLXMODE_BLASTN;
  else if (tblastn)
    mode = BLXMODE_TBLASTN;
  else if (blastp)
    mode = BLXMODE_BLASTP;
  
  return mode;
}


/* Returns the type of sequence we're dealing with */
static BlxSeqType getSeqType()
{
  return blastx || tblastx ? BLXSEQ_PEPTIDE : BLXSEQ_DNA;
}


/* Returns the number of reading frames */
static int getNumReadingFrames()
{
  return getSeqType() == BLXSEQ_PEPTIDE ? 3 : 1;
}


/* Find out if we need to fetch any sequences (they may all be contained in the input
 * files), if we do need to, then fetch them by the appropriate method. */
static gboolean blxviewFetchSequences(PfetchParams *pfetch, gboolean External, MSP *msplist)
{
  gboolean status = TRUE;

  /* First, set the fetch mode */
  if (pfetch)
    {
      /* If pfetch struct then this sets fetch mode to pfetch. */
      
      if (blxConfigSetPFetchSocketPrefs(pfetch->net_id, pfetch->port))
	blxSetFetchMode(BLX_FETCH_PFETCH) ;
    }
  else
    {
      blxFindFetchMode() ;
    }
  
  
  /* See if we have any sequences to fetch */
  DICT *dict = dictCreate(128) ;
  if (!haveAllSequences(msplist, dict))
    {
      if (strcmp(blxGetFetchMode(), BLX_FETCH_PFETCH) == 0)
	{
	  /* Fill msp->sseq fast by pfetch if possible
	   * two ways to use this:
	   *    1) call blixem main with "-P node:port" commandline option
	   *    2) setenv BLIXEM_PFETCH to a dotted quad IP address for the
	   *       pfetch server, e.g. "193.62.206.200" = Plato's IP address
	   *       or to the name of the machine, e.g. "plato"
	   *       and optionally setenv BLIXEM_PORT to the port number for
	   *       the pfetch server. */
	  enum {PFETCH_PORT = 22100} ;			    /* default port to connect on */
	  char *net_id = NULL ;
	  int port = PFETCH_PORT ;
	  
	  if (pfetch)
	    {
	      net_id = pfetch->net_id ;
	      if (!(port = pfetch->port))
		port = PFETCH_PORT ;
	    }
	  else if ((net_id = getenv("BLIXEM_PFETCH")))
	    {
	      char *port_str ;
	      
	      port = 0 ;
	      if ((port_str = getenv("BLIXEM_PORT")))
		port = atoi(port_str) ;
	      
	      if (port <= 0)
		port = PFETCH_PORT ;
	    }
	  else
	    {
	      /* Lastly try for a config file. */
	      
	      blxConfigGetPFetchSocketPrefs(&net_id, &port) ;
	    }
	  
	  if (net_id)
	    status = blxGetSseqsPfetch(msplist, dict, net_id, port, External) ;
	}
#ifdef PFETCH_HTML 
      else if (strcmp(blxGetFetchMode(), BLX_FETCH_PFETCH_HTML) == 0)
	{
	  status = blxGetSseqsPfetchHtml(msplist, dict, getSeqType()) ;
	}
#endif
    }
  messfree(dict) ;
  
  return status;
}


/* blxview() can be called either from other functions in the Blixem
 * program itself or directly by functions in other programs such as
 * xace.
 *
 * Interface
 *      opts:  may contain a number of options that tell blixem to
 *             start up with different defaults
 *    pfetch:  if non-NULL, then we use pfetch instead of efetch for
 *             _all_ sequence fetching (use the node/port info. in the
 *             pfetch struct to locate the pfetch server).
 *
 */
int blxview(char *refSeq, char *refSeqName, int start, int qOffset, MSP *msplist,
            char *opts, PfetchParams *pfetch, char *align_types, BOOL External)
{
  if (blixemWindow)
    gtk_widget_destroy(blixemWindow) ;

  if (!External)
    {
      char *config_file = NULL ;
      GError *error = NULL ;
      GKeyFile *blxGetConfig(void) ;

      /* Set up program configuration. */
      if (!(blxGetConfig()) && !blxInitConfig(config_file, &error))
	{
	  messcrash("Config File Error: %s", error->message) ;
	}
    }


  blxviewInitGlobals(refSeq, refSeqName, start, qOffset, msplist);
  blxviewGetOpts(opts, refSeq);
  gboolean status = blxviewFetchSequences(pfetch, External, msplist);
  

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(opts, align_types, msplist, refSeq, refSeqName, qOffset, start, sortMode, sortInvOn, HSPgaps) ;
    }

  return 0;
}




#if OLD_BLIXEM
//static void addgraph(GtkWidget *vbox, BOOL doFrame, Graph graph)
//{
//  GtkWidget *widget;
//
//  if (doFrame)
//    {
//      widget = gtk_frame_new(NULL);
//      gtk_container_add(GTK_CONTAINER(widget), gexGraph2Widget(graph));
//    }
//  else
//    widget = gexGraph2Widget(graph);
//
//  gtk_container_border_width (GTK_CONTAINER (widget), 0);
//
//  gtk_box_pack_start(GTK_BOX(vbox),
//		     widget,
//		     !doFrame, TRUE, 0);
//
//  graphActivate(graph);
//  graphRegister (PICK, blviewPick) ;
//  graphRegister (MIDDLE_DOWN, MiddleDownQ) ;
//  graphRegister (KEYBOARD, keyboard) ;
//  graphNewMenu(blixemMenu);
//
//  return ;
//}
#endif


/* Initialize the display and the buttons */
static void blviewCreate(char *opts, 
			 char *align_types, 
			 MSP *msplist, 
			 char *refSeq, 
			 char *refSeqName, 
			 const int qOffset,
			 const int startCoord,
			 const SortByType sortByType,
			 const gboolean sortInverted,
			 const gboolean gappedHsp)
{
  if (!blixemWindow)
    {
      blixemWindow = createMainWindow(refSeq, 
				      refSeqName, 
				      msplist, 
				      getBlastMode(), 
				      getSeqType(), 
				      getNumReadingFrames(), 
				      stdcode1, 
				      qOffset, 
				      startCoord,
				      sortByType,
				      sortInverted, 
				      gappedHsp);
      
      if (!oldWidth)
	gtk_window_set_default_size(GTK_WINDOW(blixemWindow),
				    (int)(((float)gdk_screen_width())*0.9),
				    (int)(((float)gdk_screen_height())*0.6));
      else
	gtk_window_set_default_size(GTK_WINDOW(blixemWindow),
				    oldWidth, oldHeight);


      BOOL pep_nuc_align = (*opts == 'X' || *opts == 'N') ;
      gtk_window_set_title(GTK_WINDOW(blixemWindow),
			   messprintf("Blixem %s%s%s:   %s",
				      (pep_nuc_align ? "  (" : ""),
				      (align_types ? align_types : (*opts == 'X' ? "peptide" :
								    (*opts == 'N' ? "nucleotide" : ""))),
				      (pep_nuc_align ? " alignment)" : ""),
				      refSeqName));
    }


  char *nameSeparatorPos = (char *)strrchr(refSeqName, '/');
  if (nameSeparatorPos)
    refSeqName = nameSeparatorPos + 1;

  if (dotter_first && msplist && msplist->sname && (msplist->type == HSP || msplist->type == GSP))
    {
      strcpy(HighlightSeq, msplist->sname);
//      callDotter(msplist);
    }

//  if (start_nextMatch)
//    nextMatch();
//  else
//    blviewRedraw();
}

/***********************************************************/


/* BLVIEWDESTROY frees all the allocated memory
   N.B. This memory was allocated in the calling program (acedb)

   WHAT THIS ROUTINE DOES NOT ADDRESS IS RESETTING THE LARGE NUMBER OF GLOBALS
   IN ANY SENSIBLE WAY......NOT IDEAL TO HAVE GLOBALS, EVEN LESS TO TO NOT RESET
   THEM....GRRRRRRRRRRRRRR............


   Could free auxseq, auxseq2 and padseq too, but they'd have to be remalloc'ed
   next time then.
*/
#if OLD_BLIXEM
static void blviewDestroy(GtkWidget *unused)
{
  MSP *msp, *fmsp;
  BPMSP *tmsp;

  g_free(qname_G) ;

  g_free(q);

  /* Free the allocated sequences and names */
  for (msp = MSPlist; msp; msp = msp->next)
    {
      if (msp->sseq && msp->sseq != padseq)
	{
	  for (fmsp = msp->next; fmsp; fmsp = fmsp->next)
	    {
	      if (fmsp->sseq == msp->sseq)
		fmsp->sseq = 0;
	    }

	/* Bug in fmapfeatures.c causes introns to have stale sseq's */
	if (msp->score >= 0)
	  {
	    g_free(msp->sseq);
	    g_free(msp->qname);
	    g_free(msp->sname);
	    g_free(msp->desc);
	    arrayDestroy(msp->gaps);
	    arrayDestroy(msp->xy);
	  }
	}
    }

  for (msp = MSPlist; msp; )
    {
      fmsp = msp;
      msp = msp->next;
      g_free(fmsp);
    }

  for (bpmsp = BPMSPlist.next; bpmsp; )
    {
      tmsp = bpmsp;
      bpmsp = bpmsp->next;
      g_free(tmsp);
    }

  arrayDestroy(stringentEntropyarr);
  arrayDestroy(mediumEntropyarr);
  arrayDestroy(nonglobularEntropyarr);

  blixemWindow = NULL ;
  pickedMSP = NULL ;


  /* Need to start trying to reset some state.... */
  dotterStart = dotterEnd = 0 ;


  return ;
}
#endif


void drawBigPictMSP(MSP *msp, int BPx, char strand)
{//
//  float  msp_y, msp_sx, msp_ex, midx;
//
//  if (FS(msp)) return;
//
//  msp_sx = max((float)plusmin*(msp->qstart-BigPictStart)*BPx/BigPictLen +BPoffset, 4);
//  msp_ex = max((float)plusmin*(msp->qend-BigPictStart)*BPx/BigPictLen +BPoffset, 4);
//
//  if (msp->score == -1) /* EXON */
//    {
//      oldcolor = graphColor(geneColor);
//      oldLinew = graphLinewidth(.15);
//
//      msp_y = 7.9 + queryy ;
//      if (strand == 'R')
//	msp_y += 1.5 ;
//
//      graphRectangle(msp_sx, msp_y, msp_ex, msp_y + 0.7);
//      graphColor(oldcolor);
//      graphLinewidth(oldLinew);
//    }
//  else if (msp->score == -2) /* INTRON */
//    {
//      oldcolor = graphColor(geneColor);
//      oldLinew = graphLinewidth(.15);
//
//      msp_y = 7.9 + queryy ;
//      if (strand == 'R')
//	msp_y += 1.5 ;
//
//      midx = 0.5 * (msp_sx + msp_ex) ;
//      graphLine (msp_sx, msp_y+0.4, midx, msp_y) ;
//      graphLine (msp_ex, msp_y+0.4, midx, msp_y) ;
//      graphColor(oldcolor);
//      graphLinewidth(oldLinew);
//    }
//  else if (msp->score > 0) /* BLAST MATCH */
//    {
//      int colour ;
//
//      if (!msp->sseq)
//	getsseq(msp);
//      if (!msp->id)
//	calcID(msp);
//
//      msp_y = (float)(140 - msp->id)/20 + queryy ;
//
//      if (strand == 'R')
//	msp_y += 9 ;
//
//      if (!bpmsp->next)
//	bpmsp->next = (BPMSP *)g_malloc(sizeof(BPMSP));
//
//      bpmsp = bpmsp->next;
//      strcpy(bpmsp->sname, msp->sname);
//      bpmsp->desc = msp->desc;
//
//      oldLinew = graphLinewidth(0.1);
//      bpmsp->box = graphBoxStart();
//      bpmsp->graph = graphActive();
//
//      graphFillRectangle(msp_sx, msp_y, msp_ex, msp_y+.2);
//
//      /* graphLine doesn't want to be picked */
//
//      graphBoxEnd();
//
//      if (msp->in_match_set)
//	colour = matchSetColor ;
//      else
//	colour = BLACK ;
//      
//      graphBoxDraw(bpmsp->box, colour, TRANSPARENT);
//      graphLinewidth(oldLinew);
//    }
//
//  return ;
}


/* Function: return true if pos is a value in between or equal to start and end
 */
//static int inRange(pos, start, end)
//{
//    if (start < end) {
//        if (pos >= start && pos <= end) return 1;
//	else return 0;
//    }
//    else {
//        if (pos >= end && pos <= start) return 1;
//	else return 0;
//    }
//}


/* Checks if the msp is supposed to be drawn given its frame and position */
//static void selectBigPictMSP(MSP *msp, int BPx, int BigPictStart, int BigPictLen)
//{
//    int BPend = (actframe[1] == '+' ? BigPictStart+BigPictLen : BigPictStart-BigPictLen);
//
//    /* Check if MSP is in range.  There are three cases:
//       1. MSP spans beginning of BP range
//       2. MSP spans end of BP range
//       3. MSP spans entire BP range (i.e. BigPictStart is within the MSP-range)
//       (4. MSP is included in BP range, falls into 1. and 2.)
//    */
//    if (!inRange(msp->qstart, BigPictStart, BPend) &&
//	!inRange(msp->qend, BigPictStart, BPend) &&
//	!inRange(BigPictStart, msp->qstart, msp->qend))
//        return;
//
//
//    if (actframe[1] == msp->qframe[1])
//	drawBigPictMSP(msp, BPx, 'M');
//
//    if (BigPictRev)
//	if (actframe[1] != msp->qframe[1])
//	    drawBigPictMSP(msp, BPx, 'R');
//
//    return;
//
//#ifdef OLD_CODE_TOO_COMPLICATED_AND_BUGGED__LEFT_FOR_AMUSEMENT
//    if (( ((actframe[1] == '+' && msp->frame[1] == '+')) &&
//	  (msp->qstart <  BigPictStart && msp->qend >  BigPictStart || /* left || span */
//	   msp->qstart >= BigPictStart && msp->qend <= BigPictStart+BigPictLen || /* middle */
//	   msp->qstart <  BigPictStart+BigPictLen && msp->qend > BigPictStart+BigPictLen)) /* right || span */
//	||
//	(actframe[1] == '-' && msp->frame[1] == '-' &&
//	 (msp->qend > BigPictStart && msp->qstart < BigPictStart ||
//	  msp->qend <= BigPictStart && msp->qstart >= BigPictStart-BigPictLen ||
//	  msp->qend > BigPictStart-BigPictLen && msp->qstart < BigPictStart-BigPictLen)))
//	drawBigPictMSP(msp, BPx, 'M');
//    if (BigPictRev)
//	if ( (strchr(actframe, '+') && strchr(msp->frame, '-') &&
//	      (msp->qend <  BigPictStart && msp->qstart < BigPictStart+BigPictLen ||
//	       msp->qend >= BigPictStart && msp->qstart <= BigPictStart+BigPictLen ||
//	       msp->qend <  BigPictStart+BigPictLen && msp->qstart > BigPictStart+BigPictLen))
//	     ||
//	     (strchr(actframe, '-') && strchr(msp->frame, '+') &&
//	      (msp->qend >  BigPictStart && msp->qstart < BigPictStart ||
//	       msp->qend <= BigPictStart && msp->qstart >= BigPictStart-BigPictLen ||
//	       msp->qend > BigPictStart-BigPictLen && msp->qstart < BigPictStart-BigPictLen)))
//	    drawBigPictMSP(msp, BPx, 'R');
//#endif
//}


/* Function: Put feature segment on the screen

   Note: this routine does not support reversed mode.  Worry about that later.
*/
//static void drawSEG(MSP *msp, float offset)
//{
//    float
//	msp_sy,
//	msp_ey = queryy + offset-1,
//	msp_sx,
//	msp_ex;
//
//    if (msp->qstart > BigPictStart+BigPictLen-1 ||
//	msp->qend < BigPictStart)
//        return;
//
//    msp_sx = max(SEQ2BP(msp->qstart), 4);
//    msp_ex = max(SEQ2BP(msp->qend+1), 4);
//
//    msp_sy = msp_ey - (float)msp->score/100;
//
//    oldcolor = graphColor(msp->color); oldLinew = graphLinewidth(.1);
//
//    graphFillRectangle(msp_sx, msp_sy, msp_ex, msp_ey);
//    graphColor(BLACK);
//    graphRectangle(msp_sx, msp_sy, msp_ex, msp_ey);
//    graphText(msp->desc, msp_sx, msp_ey);
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


/* Function: put XY curves on the screen.

   Note: this routine does not support reversed mode.  Worry about that later.
*/
//static void drawSEGxy(MSP *msp, float offset)
//{
//    int i, inNotFilled=0, descShown=0;
//    float
//	msp_y = queryy + offset-1,
//	x, y, xold=0, yold=0;
//
//    oldcolor = graphColor(msp->color); oldLinew = graphLinewidth(.25);
//
//    /* Must go through interpolated data outside the visible area in case the interpolation
//       spans the start or the end of the visible area */
//    if (msp->shape == XY_INTERPOLATE) {
//        for (i = 0; i < BigPictStart; i++) {
//	    if (arr(msp->xy, i, int) != XY_NOT_FILLED) {
//	      xold = SEQ2BP(i);
//	      yold = msp_y - (float)arr(msp->xy, i, int)/100*fsPlotHeight;
//	    }
//	}
//    }
//
//    for (i = BigPictStart; i < BigPictStart+BigPictLen-1; i++) {
//	if (arr(msp->xy, i, int) == XY_NOT_FILLED) {
//	    inNotFilled = 1;
//	}
//	else {
//	    x = SEQ2BP(i);
//	    y = msp_y - (float)arr(msp->xy, i, int)/100*fsPlotHeight;
//	    if (xold && (!inNotFilled || msp->shape == XY_INTERPOLATE)) {
//	        if (x != xold || y != yold) graphLine(xold, yold, x, y);
//		if (!descShown && msp->desc) {
//		      int linecolor = graphColor(BLACK);
//		      graphText(msp->desc, (xold > BPoffset ? xold : BPoffset), msp_y);
//		      graphColor(linecolor);
//		      descShown = 1;
//		  }
//	    }
//	    xold = x;
//	    yold = y;
//	    inNotFilled = 0;
//	}
//    }
//
//    /* Draw interpolated data if it spans the end of the visible area */
//    if (msp->shape == XY_INTERPOLATE && xold) {
//        for (; i < qlen; i++) {
//	    if (arr(msp->xy, i, int) != XY_NOT_FILLED) {
//	        x = SEQ2BP(i);
//		y = msp_y - (float)arr(msp->xy, i, int)/100*fsPlotHeight;
//		graphLine(xold, yold, x, y);
//		break;
//	    }
//	}
//    }
//
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


//static void drawEntropycurve(int start, int end, Array array, int win)
//{
//    int i;
//    float x, y, xold=0, yold=0;
//
//    oldLinew = graphLinewidth(.3);
//
//    for (i = start; i < end; i++) {
//	if (i > win/2 && i < qlen - win/2) {
//	    x = SEQ2BP(i);
//	    y = queryy + 9 - arr(array, i, float)*2;
//	    if (xold) graphLine(xold, yold, x, y);
//	    xold = x;
//	    yold = y;
//
//	    /*if (arr(array, i, float) < min) min = arr(array, i, float);
//	    if (arr(array, i, float) > max) max = arr(array, i, float);*/
//	}
//	/**c = q[i];
//	  graphText(c, SEQ2BP(i), queryy+5);*/
//    }
//    /* printf("min = %f , max = %f\n", min, max); */
//
//    graphLinewidth(oldLinew);
//}


/* Draw separator line a la Mosaic <HR> */
//static void drawSeparator(void) {
//    oldLinew = graphLinewidth(separatorwidth);
//    if (squash)
//	graphColor(RED);
//    else
//	graphColor(DARKGRAY);
//    graphLine(0, queryy-1.0, nx+1, queryy-1.0);
//
//    graphLinewidth(.2);
//    graphColor(WHITE);
//    graphLine(0, queryy-1.1+0.5*separatorwidth, nx+1, queryy-1.1+0.5*separatorwidth);
//    graphColor(BLACK);
//    graphLinewidth(oldLinew);
//}


int frame2graphno(int frame)
{
//  if (blastn)
//    {
//      if (frame == 1)
//	return 0;
//      else
//	return 1;
//    }
//  else
//    {
//      if (abs(frame) == 1)
//	return 0;
//      else if (abs(frame) == 2)
//	return 1;
//      else
//	return 2;
//    }
  return 0;
}


/* Draw the given sequence name at the given coordinates. The name is formatted so that
 * it includes the direction of the strand and does not exceed the given length. */
void drawStrandName(char *name, char *frame, int max_len, int x, int y)
{
//  gchar *displayName = g_strconcat(abbrevTxt(name, max_len - 2), 
//                                  (strchr(frame, '+') ? " +" : " -"),
//                                   NULL);
//  graphText(displayName, x, y);
//  g_free(displayName);
}


/* Redraw the entire blixem window. (Call mainWindowRedrawAll directly where possible,
 * rather than this function, which relies on the global variable 'blixemWindow'). */
void blviewRedraw(void)
{
  if (blixemWindow)
    {
      mainWindowRedrawAll(blixemWindow);
    }
}


/* GET3RD_BASE returns every third base from string q
*/
#if OLD_BLIXEM
static char *get3rd_base(int start, int end, char *q)
{
  static int i, len;
  static char *bases = NULL, *aux, *aux2;

  if (start < 1 || end > qlen)
    {
      messerror ( "Genomic sequence coords out of range: %d - %d\n",
		  start, end);
      
      return NULL;
    }

  if (!bases)
    {
      bases = g_malloc(qlen+1);
      aux = g_malloc(qlen+1);
      aux2 = g_malloc(qlen+1);
    }

  len = abs(end-start) + 1;

  if (start < end)
    {
      for (i=start; i < end; i += 3 )
	bases[(i-start)/3] = q[i-1];
    }
  else if (start > end) /* Reverse and complement it first */
    {
      strncpy(aux2, q+end-1, start-end+1);
      aux2[start-end+1] = 0;

      if (!revComplement(aux, aux2))
	  messcrash ("Cannot reverse & complement") ;

      for (i=0; i < len; i += 3)
	bases[i/3] = aux[i];
    }
  else
    {
      return NULL;
    }
      
  bases[len/3] = '\0';

  return bases;
}
#endif


//static void allocAuxseqs(int len)
//{
//    if (auxseq) {
//	g_free(auxseq);
//	g_free(auxseq2);
//    }
//
//    auxseq = g_malloc(len+1);
//    auxseq2 = g_malloc(len+1);
//    auxseqlen = len;
//}




void blxAssignPadseq(MSP *msp, MSP *msplist)
{
    static int padseqlen=0;
    char *oldpadseq;
    int len = max(msp->sstart, msp->send);
    MSP *hsp;

    if (!padseq) {
	padseq = g_malloc(INITDBSEQLEN+1);
	memset(padseq, '-', INITDBSEQLEN);
	padseqlen = INITDBSEQLEN;
    }

    if (len > padseqlen) {
	oldpadseq = padseq;
	g_free(padseq);

	padseq = g_malloc(len+1);
	memset(padseq, '-', len);
	padseqlen = len;

	/* Change all old padseqs to new */
	for (hsp = msplist; hsp ; hsp = hsp->next)
	    if (hsp->sseq == oldpadseq) hsp->sseq = padseq;
    }

    msp->sseq = padseq;
 }



/* GETQSEQ translates a segment of the query seq (with 'sequence' coords = 1...) */
char *getqseq(int start, int end, const char const *refSeq)
{
//  char *query = NULL ;
//  int refSeqLen = strlen(refSeq);
//
//  if (start < 1 || end < 1 || start > refSeqLen || end > refSeqLen)
//    {
//      messout ( "Requested query sequence %d - %d out of available range: 1 - %d\n",
//		start, end, refSeqLen);
//
//      return NULL;
//    }
//
//  if (blastp || tblastn)
//    {
//      query = g_malloc(end-start+2);
//      strncpy(query, q+start-1, end-start+1);
//      query[end-start+1] = 0;
//
//      return query;
//    }
//
//  if (abs(end-start)+1 > auxseqlen)
//    allocAuxseqs(abs(end-start)+1);
//
//  if (start <= end)
//    {
//      strncpy(auxseq, q+start-1, end-start+1);
//      auxseq[end-start+1] = 0;
//    }
//  else if (start > end) /* Reverse and complement it */
//    {
//      strncpy(auxseq2, refSeq + end - 1, start - end + 1);
//      auxseq2[start-end+1] = 0;
//
//      if (!revComplement(auxseq, auxseq2))
//	messcrash ("Cannot reverse & complement") ;
//    }
//  else
//    {
//      return NULL;
//    }
//
//  if (blastn) 
//    {
//      query =  g_malloc(strlen(auxseq) + 1);
//      strcpy(query, auxseq);
//
//      return query;
//    }
//
//  /* Translate the DNA sequence */
//  if (!(query = blxTranslate(auxseq, stdcode1)))
//    messcrash ("Cannot translate the genomic sequence") ;
//
//  return query;
  return NULL;
}


/* GETSSEQ fetches the database sequence from an external database,
 * currently uses either efetch or pfetch. */
//static void getsseq(MSP *msp, MSP *msplist)
//{
//
//  MSP  *auxmsp;
//  int   len;
//  char *fetch_prog ;
//  char *seq_buf = NULL ;
//
//  if (!*msp->sname)
//    {
//      messout ( "Nameless HSP at %d-%d - skipping Efetch\n",
//		msp->qstart+qoffset, msp->qend+qoffset);
//      blxAssignPadseq(msp, msplist);
//
//      return ;
//    }
//
//  fetch_prog = blxGetFetchProg() ;
//
//  if (verbose)
//    printf("%sing %s\n", fetch_prog, msp->sname);
//  else
//    {
//      if (!fetchCount)
//	{
//	  fetchCount++;
//	  printf("\n%sing external sequences", fetch_prog);
//	}
//      printf(".");
//      fflush(stdout);
//    }
//
//  if (msp->score >= 0)
//    {
//      if ((seq_buf = getSeq(msp->sname, blxGetFetchMode())))
//	{
//	  msp->sseq = g_malloc(strlen(seq_buf)+1);
//
//	  /* Harmonize upper and lower case - (t)blastx always translate
//	   * the query to upper case */
//	  if (isupper(*q) || tblastx || blastx)
//	    {
//	      for (i=0; seq_buf[i]; i++)
//		msp->sseq[i] = freeupper(seq_buf[i]);
//	      msp->sseq[i] = 0;
//	    }
//	  else
//	    {
//	      for (i=0; seq_buf[i]; i++)
//		msp->sseq[i] = freelower(seq_buf[i]);
//
//	      msp->sseq[i] = 0;
//	    }
//
//	  /* Check illegal offsets */
//	  len = strlen(msp->sseq);
//	  for (auxmsp = msplist; auxmsp ; auxmsp = auxmsp->next)
//	    if (!strcmp(auxmsp->sname, msp->sname) && auxmsp->send > len )
//	      {
//		printf("%s HSP with offset beyond sequence (%d > %d) - using pads\n",
//		       msp->sname, auxmsp->send, len);
//		blxAssignPadseq(msp, msplist);
//
//		break;
//	      }
//	}
//      else
//	{
//#if !defined(ACEDB)
//	  messout ( "Unable to %s %s - using pads instead\n", fetch_prog, msp->sname);
//#endif
//	  /* Sequence not in database - fill up with pads */
//	  blxAssignPadseq(msp, msplist);
//	}
//
//      
//
//
//
//
//
//      /* Set sseq for all MSPs of this subject, either a padseq or
//       * the sequence for that subject, taking into account the _strand_.
//       * 
//       * THIS IS NOT TOO EFFICIENT...WE COULD BE HASHING HERE....
//       *  */
//      for (auxmsp = msplist ; auxmsp ; auxmsp = auxmsp->next)
//	{
//	  if (strcmp(auxmsp->sname, msp->sname) == 0)
//	    {
//	      /* Either assign */
//	      if (msp->sseq == padseq)
//		blxAssignPadseq(auxmsp, msplist);
//	      else if (MSPSTRAND(auxmsp->sframe) == MSPSTRAND(msp->sframe))
//		auxmsp->sseq = msp->sseq ;
//
//	      calcID(auxmsp);
//	    }
//	}
//    }
//
//
//  return ;
//}



/********************************************************************************
**                            BIG PICTURE ROUTINES                            ***
********************************************************************************/

//static void BigPictToggle (void) {
//    BigPictON = !BigPictON;
//    blviewRedraw();
//}
//
//static void BigPictToggleRev (void) {
//    BigPictRev = !BigPictRev;
//    blviewRedraw();
//}
//
//
//static void zoomOut(void)
//{
//  BigPictZoom *= 2;
//
//  blviewRedraw();
//}
//
//static void zoomIn(void)
//{
//  BigPictZoom /= (float)2;
//  if (BigPictZoom < 1)
//    BigPictZoom = 1;
//
//  blviewRedraw();
//}
//
//static void zoomWhole(void)
//{
//  BigPictZoom = (float)qlen/displen;
//
//  blviewRedraw();
//}


/* If crosshair-coordinates are screwed up, change here!
 ********************************************************/
#if OLD_BLIXEM
static int x_to_residue(float x)
{
  int retval;

  if (blastx || tblastx)
    retval = dispstart + plusmin*(x - NAMESIZE - 22.3)*3;
  else
    retval = dispstart + plusmin*(x - NAMESIZE - 22);


  if (plusmin > 0) {
    if (retval < dispstart) retval = dispstart;
    if (retval > dispstart+displen-1) retval = dispstart+displen-1;
    return retval;
  }
  else {
    if (retval > dispstart-1) retval = dispstart-1;
    if (retval < dispstart-displen) retval = dispstart-displen;
    return retval +1;
  }
}
#endif


//static void displayResidue(double x)
//{
//  int qpos, spos ;
//  static char queryname[NAMESIZE+1] ;
//
//  if (!*queryname)
//    {
//      if (!*qname_G)
//	strcpy(queryname, "Query");
//      else
//	{
//	  char *abbrev ;
//
//	  abbrev = abbrevTxt(qname_G, NAMESIZE) ;
//	  strncpy(queryname, abbrev, NAMESIZE) ;
//	}
//
//      queryname[NAMESIZE] = 0 ;
//    }
//
//
//  qpos = x_to_residue(x);
//
//  if (!pickedMSP)
//    {
//      sprintf(message, "%d   No subject picked", qpos + qoffset) ;
//    }
//  else
//    {
//      if (blastx || tblastx)
//	{
//	  spos = gapCoord(pickedMSP, qpos, 3);
//	  
//	  if (spos < pickedMSP->sstart)
//	    spos = pickedMSP->sstart;
//	  if (spos > pickedMSP->send)
//	    spos = pickedMSP->send;
//	}
//      else if (blastn)
//	{
//	  spos = gapCoord(pickedMSP, qpos, 1);
//
//	  if (spos < pickedMSP->sstart)
//	    spos = pickedMSP->sstart;
//	  if (spos > pickedMSP->send)
//	    spos = pickedMSP->send;
//	}
//      else if (tblastn)
//	{
//	  if (pickedMSP->sstart < pickedMSP->send)
//	    {
//	      spos = (qpos - pickedMSP->qstart)*3 + pickedMSP->sstart;
//
//	      if (spos < pickedMSP->sstart) spos = pickedMSP->sstart;
//	      if (spos > pickedMSP->send)   spos = pickedMSP->send;
//	    }
//	  else
//	    {
//	      spos = pickedMSP->sstart - (qpos - pickedMSP->qstart)*3;
//	      
//	      if (spos > pickedMSP->sstart) spos = pickedMSP->sstart;
//	      if (spos < pickedMSP->send)   spos = pickedMSP->send;
//	    }
//	}
//      else
//	{
//	  spos = qpos - pickedMSP->qstart + pickedMSP->sstart;
//	  
//	  if (spos < pickedMSP->sstart) spos = pickedMSP->sstart;
//	  if (spos > pickedMSP->send)   spos = pickedMSP->send;
//	}
//
//      sprintf(message, "%d   ", qpos + qoffset) ;
//
//      if (HSPgaps)
//	strcat(message, "Gapped HSP - no coords");
//      else
//	strcat(message, messprintf("%s: %d", pickedMSP->sname, spos));
//
//    }
//
//  gtk_entry_set_text(GTK_ENTRY(messageWidget), message);
//
//  return ;
//}


#if OLD_BLIXEM
static void markDNA(double y)
{
//  Graph old = graphActive();
//
//  if (y < 0)
//    return;
//
//  graphActivate(frameGraphQ[0]);
//
//  graphXorLine (NAMESIZE+22, 1.5 + y,
//		NAMESIZE+22+displen/3, 1.5 + y);
//
//  graphActivate(old);
}
#endif



/* Mouse drag on the upper, "big picture", section of the graph. */
#if OLD_BLIXEM
static void MiddleDownBP (double x, double y)
{
  graphXorBox (BPbox, x - BPboxwidth/2, 1.85);

  graphRegister (MIDDLE_DRAG, MiddleDragBP) ;	/* must redo */
  graphRegister (MIDDLE_UP, MiddleUpBP) ;

  oldx = x;
}

static void MiddleDragBP (double x, double y)
{
  graphXorBox (BPbox, oldx - BPboxwidth/2, 1.85);
  graphXorBox (BPbox, x - BPboxwidth/2, 1.85);

  oldx = x ;
}

static void MiddleUpBP (double x, double y)
{
  graphFitBounds (&nx, 0);
  nx -= 2;						    /* scrollbars */

  dispstart = BigPictStart + (plusmin * ((x-BPoffset)/(nx-BPoffset)*BigPictLen - displen/2)) ;
  blviewRedraw();

  return ;
}



/* Mouse drag on the lower, "query", section of the graph. */
static void MiddleDragQ (double x, double y)
{
//  int i, noframes = blastn ? 2 : 3;
//
//  for(i=0; i<noframes; i++)
//    {
//      graphActivate(frameGraphQ[i]);
//      graphXorLine (oldx, 0, oldx, 1000) ;
//      graphXorLine (x, 0, x, 1000);
//      graphActivate(frameGraph[i]);
//      graphXorLine (oldx, 0, oldx, 1000) ;
//      graphXorLine (x, 0, x, 1000);
//    }
//
//  if (blastx || tblastx)
//    {
//      markDNA(DNAstep);
//      DNAstep =  abs ( (x_to_residue(x) - dispstart) % 3);
//      markDNA(DNAstep);
//    }
//
//  graphActivate(blixemGraph);
//
//  displayResidue(x);
//
//  oldx = x ;
}

static void MiddleDownQ (double x, double y)
{//
//  float nyf;
//  int i, noframes = blastn ? 2 : 3;
//
//  graphFitBounds (&nx, 0);
//  graphWhere(0,0,0, &nyf);
//  ny = nyf - 0.5 ;
//  graphRegister(MIDDLE_DRAG, MiddleDragQ) ;	/* must redo */
//  graphRegister(MIDDLE_UP, MiddleUpQ) ;
//
//  for (i = 0 ; i < noframes ; i++)
//    {
//      graphActivate(frameGraphQ[i]);
//      graphXorLine (x, 0, x, 1000);
//      graphActivate(frameGraph[i]);
//      graphXorLine (x, 0, x, 1000);
//    }
//
//  /* Cleanse messagebox from pick */
//  *message = 0 ;
//  gtk_entry_set_text(GTK_ENTRY(messageWidget), message) ;
//
//  graphActivate(frameGraphQ[0]);
//  if (blastx || tblastx)
//    {
//      DNAstep = (x_to_residue(x) - dispstart) % 3;
//      markDNA(DNAstep);
//    }
//
//  displayResidue(x);
//  oldx = x;
//
//  return ;
}

static void MiddleUpQ (double x, double y)
{

//  dispstart = x_to_residue(x) - plusmin*displen/2;
//
//  blviewRedraw();
}
#endif
	


//static void setHighlight(void)
//{
//    static char dfault[64] = "";
//    ACEIN string_in;
//
//    if ((string_in = messPrompt ("String: (wildcards: * ?)",
//				 dfault, "t", 0)))
//      {
//	/* ANSI way */
//	strncpy(searchSeq, aceInWord(string_in), NAMESIZE+3);
//	searchSeq[NAMESIZE+3] = 0;
//	for (cp = searchSeq; *cp ; cp++) *cp = freeupper(*cp);
//
//	/* Non-ANSI bsd way :
//	   if (!re_comp(searchSeq)) fprintf(stderr, "%s\n", re_comp(searchSeq));*/
//
//	strncpy(dfault, searchSeq, 63);
//	dfault[63] = '\0';
//
//	blviewRedraw();
//
//	aceInDestroy (string_in);
//    }
//}


//static void clrHighlight(void)
//{
//    MSP *msp;
//
//    /* Clear highlighted */
//    *searchSeq = *HighlightSeq = 0;
//    pickedMSP = 0;
//
//    /* Unhide hidden matches */
//    for (msp = MSPlist; msp ; msp = msp->next)
//	if (msp->score == -999) {
//	    msp->score = msp->id;
//	    calcID(msp);
//	}
//
//    blviewRedraw();
//}




/************************  BLIXEM SETTINGS  WINDOW  **************************/

//static void blixemConfColourMenu(KEY key, int box)
//{
//  /* Taken from ? */
//  int *colour;
//
//  if (graphAssFind(assVoid(box+2000), &colour))
//    {
//      *colour = key;
//      graphBoxDraw(box, BLACK, *colour);
//      graphRedraw();
//      blviewRedraw();
//    }
//
//  return ;
//}


//static void blixemConfColour(int *colour, int init, float x, float *y, int len, char *text)
//{
//    int box;
//    if (text)
//	graphText(text, x+len+1, *y);
//
//    box = graphBoxStart();
//    graphRectangle(x, *y, x+len, *y+1);
//    graphBoxEnd();
//    graphBoxFreeMenu(box, blixemConfColourMenu, graphColors);
//    *colour = init;
//    graphBoxDraw(box, BLACK, init);
//    graphAssociate(assVoid(box+2000), colour);
//
//    *y += 1.5;
//}


//static void buttonCheck(char* text, void (*func)(void), float x, float *y, int On)
//{
//  char *but_text ;
//
//  but_text = hprintf(0, "%s %s", On ? "*" : " ", text) ;
//  graphButton(but_text, func, x, *y);
//  g_free(but_text) ;
//
//  /* Could be more fancy like this
//     if (On) graphFillRectangle(x-.5, *y, x-2, *y+1.2);
//     else graphRectangle(x-.5, *y, x-2, *y+1.2);*/
//
//  *y += 1.5;
//
//  return ;
//}


//static void graphButtonDisable(char* text, float x, float *y, int On)
//{
//    int box = graphBoxStart();
//    graphText(messprintf("%s %s", On ? "*" : " ", text), x, *y);
//    graphBoxEnd();
//    graphBoxDraw (box, DARKGRAY, WHITE);
//
//    *y += 1.5;
//}


//static void setStringentEntropywin(char *cp)
//{
//    stringentEntropywin = atoi(cp);
//    calcEntropyArrays(TRUE);
//    blviewRedraw();
//}
//static void setMediumEntropywin(char *cp)
//{
//    mediumEntropywin = atoi(cp);
//    calcEntropyArrays(TRUE);
//    blviewRedraw();
//}
//static void setNonglobularEntropywin(char *cp)
//{
//    nonglobularEntropywin = atoi(cp);
//    calcEntropyArrays(TRUE);
//    blviewRedraw();
//}



//static void settingsPick(int box, double x_unused, double y_unused, int modifier_unused)
//{
//    if (box == stringentEntropybox) graphTextScrollEntry(stringentEntropytx,0,0,0,0,0);
//    else if (box == mediumEntropybox) graphTextScrollEntry(mediumEntropytx,0,0,0,0,0);
//    else if (box == nonglobularEntropybox) graphTextScrollEntry(nonglobularEntropytx,0,0,0,0,0);
//}

//static void settingsRedraw(void)
//{
//  float x1=1, x2=35, y;
//
//  static MENUOPT sortMenu[] =
//    {
//      {sortByScore,      "Score"},
//      {sortById,         "Identity"},
//      {sortByName,       "Name"},
//      {sortByPos,        "Position"},
//      {0, 0}
//    };
//
//
////  if (!graphActivate(settingsGraph))
////    return;
//
//  graphClear();
//
//    /* Background */
//    graphBoxDraw(0, backgColor, backgColor);
//    graphFitBounds (&nx, &ny);
//    graphColor(backgColor); graphRectangle(0, 0, nx+100, ny+100);graphColor(BLACK);
//
//    y  = 1;
//    graphText("Toggles:", x1, y);
//    y += 1.5;
//
//    buttonCheck(BigPictToggleStr, BigPictToggle, x1, &y, BigPictON);
//    if (blastn || blastx || tblastx) {
//	if (BigPictON) {
//	    buttonCheck(BigPictToggleRevStr, BigPictToggleRev, x1, &y, BigPictRev);
//	    menuUnsetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
//	}
//	else {
//	    graphButtonDisable(BigPictToggleRevStr, x1, &y, BigPictRev);
//	    menuSetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
//	}
//    }
//    buttonCheck(entropytoggleStr, entropytoggle, x1, &y, entropyOn);
//
//    buttonCheck(toggleDESCStr, toggleDESC, x1, &y, DESCon);
//    buttonCheck(squashFSStr, squashFSdo, x1, &y, squashFS);
//    buttonCheck(squashMatchesStr, squashMatches, x1, &y, squash);
//    buttonCheck(toggleIDdotsStr, toggleIDdots, x1, &y, IDdots);
//    buttonCheck(printColorsStr, printColors, x1, &y, printColorsOn);
//    buttonCheck(SortInvStr, sortToggleInv, x1, &y, sortInvOn);
//
///* Old disused toggles ...? */
///*    buttonCheck("Display colours", toggleColors, x1, &y, );*/
///*    buttonCheck("Printerphilic colours", printColors, x1, &y, printColorsOn);*/
//
//
//    y = 1;
//    graphText("Menus:", x2, y);
//    y += 1.5;
//
//    blixemConfColour(&backgColor, backgColor, x2, &y, 3, "Background colour");
//    blixemConfColour(&gridColor, gridColor, x2, &y, 3, "Grid colour");
//    blixemConfColour(&IDcolor, IDcolor, x2, &y, 3, "Identical residues");
//    blixemConfColour(&consColor, consColor, x2, &y, 3, "Conserved residues");
//
//
//    graphText("Sort HSPs by ", x2, y);
//    graphBoxMenu(graphButton(messprintf("%s", sortModeStr), settingsRedraw, x2+13, y), sortMenu);
//    y += 1.5;
//
//    graphText("Fetch by ", x2, y);
//    graphBoxMenu(graphButton(messprintf("%s", blxGetFetchMode()), settingsRedraw, x2+9, y), blxPfetchMenu()) ;
//    y += 1.5;
//
//    if (entropyOn) {
//	graphText("Complexity curves:", x2, y);
//	y += 1.5;
//
//	sprintf(stringentEntropytx, "%d", stringentEntropywin);
//	stringentEntropybox = graphTextScrollEntry (stringentEntropytx, 6, 3, x2+19, y, setStringentEntropywin);
//	blixemConfColour(&stringentEntropycolor, stringentEntropycolor, x2, &y, 3, "window length:");
//
//	sprintf(mediumEntropytx, "%d", mediumEntropywin);
//	mediumEntropybox = graphTextScrollEntry (mediumEntropytx, 6, 3, x2+19, y, setMediumEntropywin);
//	blixemConfColour(&mediumEntropycolor, mediumEntropycolor, x2, &y, 3, "window length:");
//
//	sprintf(nonglobularEntropytx, "%d", nonglobularEntropywin);
//	nonglobularEntropybox = graphTextScrollEntry (nonglobularEntropytx, 6, 3, x2+19, y, setNonglobularEntropywin);
//	blixemConfColour(&nonglobularEntropycolor, nonglobularEntropycolor, x2, &y, 3, "window length:");
//    }
//
//    graphRedraw();
//
//
//    return ;
//}


/*
static void dotterPanel(void)
{
    graphCreate(TEXT_FIT, "Dotter Panel", 0, 0, .3, .1);
    graphButton("Dotter", callDotter, 1, 1);
    graphButton("Dotter HSPs only",callDotterHSPs, 1, 2.5);
    graphButton("Dotter query vs. itself",callDotterSelf, 1, 4);
    graphRedraw();
}
*/

/* Menu checkmarks */

//static void menuCheck(MENU menu, int mode, int thismode, char *str)
//{
//  if (mode == thismode)
//    menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
//  else
//    menuUnsetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
//
//  return ;
//}


//static void setMenuCheckmarks(void)
//{
//  menuCheck(settingsMenu, sortMode, SORTBYSCORE, SortByScoreStr);
//  menuCheck(settingsMenu, sortMode, SORTBYID, SortByIdStr);
//  menuCheck(settingsMenu, sortMode, SORTBYNAME, SortByNameStr);
//  menuCheck(settingsMenu, sortMode, SORTBYPOS, SortByPosStr);
//  menuCheck(settingsMenu, 1, sortInvOn, SortInvStr);
//  
//  menuCheck(settingsMenu, 1, BigPictON, BigPictToggleStr);
//  menuCheck(settingsMenu, 1, BigPictRev, BigPictToggleRevStr);
//  menuCheck(settingsMenu, 1, IDdots, toggleIDdotsStr);
//  menuCheck(settingsMenu, 1, squash, squashMatchesStr);
//  menuCheck(settingsMenu, 1, squashFS, squashFSStr);
//  menuCheck(settingsMenu, 1, entropyOn, entropytoggleStr);
//  menuCheck(settingsMenu, 1, printColorsOn, printColorsStr);
//  menuCheck(settingsMenu, 1, colortoggle, toggleColorsStr);
//  menuCheck(settingsMenu, 1, verbose, toggleVerboseStr);
//  menuCheck(settingsMenu, 1, HiliteSins, toggleHiliteSinsStr);
//  menuCheck(settingsMenu, 1, HiliteUpperOn, toggleHiliteUpperStr);
//  menuCheck(settingsMenu, 1, HiliteLowerOn, toggleHiliteLowerStr);
//  menuCheck(settingsMenu, 1, DESCon, toggleDESCStr);
//  graphNewBoxMenu(settingsButton, settingsMenu);
//  graphNewMenu(blixemMenu);

//  return ;
//}


/* Checks the MSP list of sequences for blixem to display to see if it contains sequence
 * data for each sequence (this may happen if acedb, for instance, starts blixem up and
 * acedb already contained all the sequences).
 * Returns TRUE if all the sequences are already there, FALSE otherwise.
 * Optionally returns the names of all the sequences that need to be fetched in dict if
 * one is supplied by caller. */
static BOOL haveAllSequences(const MSP const *msplist, DICT *dict)
{
  BOOL result = TRUE ;

  const MSP *msp = NULL;
  for (msp = msplist ; msp ; msp = msp->next)
    {
      if (!msp->sseq && msp->sname && *msp->sname && msp->score >= 0)
	{
	  result = FALSE ;
	  if (dict)
	    dictAdd (dict, msp->sname, 0) ;
	  else
	    break ;					    /* No dict so no need to carry on. */
	}
    }
  /* RD note: I would expect to have looked at msp->type above */

  return result ;
}



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE

/* currently unused.... */
static void pfetchWindow (MSP *msp)
{
  ACEIN pipe;
  char  cmd[50+1] = "pfetch";
  char  title[50+1];
  char *lineBuf, *textBuf;
  STORE_HANDLE handle = handleCreate();

  sprintf(cmd, "pfetch --client=acedb_%s_%s -F '%s' &",
	  getSystemName(), getLogin(TRUE), msp->sname);
  sprintf(title, "pfetch: %s", msp->sname);

  pipe = aceInCreateFromPipe(cmd, "r", NULL, handle);
  if (pipe)
    {
      textBuf = strnew("", handle);
      aceInSpecial(pipe, "\n\t");

      while ((lineBuf = aceInCard(pipe)))
	  textBuf = g_strconcat(textBuf,"\n", lineBuf, NULL);

      /* call gexTextEditor with the pfetch output & no buttons */
      gexTextEditorNew(title,
		       textBuf,
		       0,
		       NULL,
		       NULL,          /* editorOK */
		       NULL,          /* editorCancel */
		       FALSE          /* set to readonly */
		       );
    }
  handleDestroy(handle);
}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

/* Given a string, returns the string if is less than a fudge factor or returns
 * the string abbreviated in the form "xxxxx...yyyyy". The string is a pointer
 * to a static buffer, caller should make a copy of it if they wish to retain
 * it or alter it. */
//static char *abbrevTxt(char *text, int max_len)
//{
//  char *result ;
//  char *abbrev = "<>" ;
//  static char *abbrev_buf = NULL ;
//  static int buf_len, text_len, head_bytes, tail_bytes ;
//
//  /* First time through allocate the reusable buffer, allocate larger one if required. */
//  if (abbrev_buf == NULL)
//    {
//      buf_len = 50 ;					    /* not many sequence names this long. */
//      abbrev_buf = (char *)g_malloc(buf_len) ;
//    }
//  else if (max_len > (buf_len - 1))
//    {
//      g_free(abbrev_buf) ;
//      buf_len = (max_len * 1.5) ;
//      abbrev_buf = (char *)g_malloc(buf_len) ;
//    }
//
//
//  result = abbrev_buf ;
//
//  text_len = strlen(text) ;
//  if (text_len <= max_len)
//    {
//      result = strcpy(result, text) ;
//    }
//  else
//    {
//      char *tail_ptr ;
//
//      /* don't really need to calculate these each time... */
//      head_bytes = (max_len / 2) - 1 ;
//      tail_bytes = max_len - ((max_len / 2) + 1) ;
//
//      tail_ptr = text + text_len - tail_bytes ;		    /* trailing null fudged in here. */
//
//      result = strncpy(result, text, head_bytes) ;
//
//      strcpy((result + head_bytes), "") ;
//
//      result = strcat(result, abbrev) ;
//
//      result = strcat(result, tail_ptr) ;
//    }
//
//  return (result) ;
//}


/* Print out MSP's, probably for debugging.... */
//static void printMSPs(void)
//{
//  MSP *msp;
//
//  for (msp = MSPlist; msp; msp = msp->next)
//    {
//#ifdef ACEDB
//      if (msp->key)
//	printf("%d %s ", msp->key, name(msp->key)) ;
//      else
//	printf("0 NULL_KEY") ;
//#endif
//      printf("%s %d %d %d %d %d %d :%s:\n",
//	     msp->sname, msp->score, msp->id,
//	     msp->qstart+qoffset, msp->qend+qoffset,
//	     msp->sstart, msp->send, (msp->sseq ? msp->sseq : ""));
//    }
//
//  return ;
//}



/***************** end of file ***********************/
