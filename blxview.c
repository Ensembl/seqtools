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
 * CVS info:   $Id: blxview.c,v 1.21 2010-03-11 11:17:14 gb10 Exp $
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
#include <SeqTools/detailview.h>

#ifdef ACEDB
#include <wh/display.h>
#endif


/* This is the only place this is set so that you get the same version/compile whether this is
 * compiled stand alone or as part of xace. */
char *blixemVersion = BLIXEM_VERSION_COMPILE ;



typedef struct _BPMSP
{
  char sname[FULLNAMESIZE+1];
  char *desc;
  int box;
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
static void blviewCreate(char *opts, char *align_types, MSP *msplist, char *refSeq, char *refSeqName, const int qOffset, const int startCoord, const int bigPictZoom, const SortByType sortByType, const gboolean sortInverted, const gboolean gappedHsp) ;

#if OLD_BLIXEM
static char *get3rd_base(int start, int end, char *q);
#endif

static BOOL haveAllSequences(const MSP const *msplist, DICT *dict) ;
//static void getsseq(MSP *msp, MSP *msplist) ;
//static char *getSeq(char *seqname, char *fetch_prog) ;

//static BOOL smartDotterRange(char *selected_sequence, const MSP const *msp_list, int blastn,
//			     int strand_sign, int view_start, int view_end,
//			     char **dottersseq_out, int *dotter_start_out, int *dotter_end_out) ;


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
GtkWidget *blixemWindow = NULL ;



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
      blviewCreate(opts, align_types, msplist, refSeq, refSeqName, qOffset, start, BigPictZoom, sortMode, sortInvOn, HSPgaps) ;
    }

  return 0;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *opts, 
			 char *align_types, 
			 MSP *msplist, 
			 char *refSeq, 
			 char *refSeqName, 
			 const int qOffset,
			 const int startCoord,
			 const int bigPictZoom,
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
				      bigPictZoom,
				      sortByType,
				      sortInverted, 
				      gappedHsp);

      /* Set the initial window size */
      GdkScreen *screen = gtk_widget_get_screen(blixemWindow);
      const int initWidth = oldWidth ? oldWidth : gdk_screen_get_width(screen) * 0.9;
      const int initHeight = oldHeight ? oldHeight : gdk_screen_get_height(screen) * 0.6;
      gtk_window_set_default_size(GTK_WINDOW(blixemWindow), initWidth, initHeight);


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
    {
      refSeqName = nameSeparatorPos + 1;
    }
  
  if (dotter_first && msplist && msplist->sname && (msplist->type == HSP || msplist->type == GSP))
    {
      strcpy(HighlightSeq, msplist->sname);
//      callDotter(msplist);
    }

  if (start_nextMatch)
    {
      /* Set the start coord to be the start of the next MSP on from the default start coord */
      nextMatch(mainWindowGetDetailView(blixemWindow));
    }
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
//static void blviewDestroy(GtkWidget *unused)
//{
//  MSP *msp, *fmsp;
//  BPMSP *tmsp;
//
//  g_free(qname_G) ;
//
//  g_free(q);
//
//  /* Free the allocated sequences and names */
//  for (msp = MSPlist; msp; msp = msp->next)
//    {
//      if (msp->sseq && msp->sseq != padseq)
//	{
//	  for (fmsp = msp->next; fmsp; fmsp = fmsp->next)
//	    {
//	      if (fmsp->sseq == msp->sseq)
//		fmsp->sseq = 0;
//	    }
//
//	/* Bug in fmapfeatures.c causes introns to have stale sseq's */
//	if (msp->score >= 0)
//	  {
//	    g_free(msp->sseq);
//	    g_free(msp->qname);
//	    g_free(msp->sname);
//	    g_free(msp->desc);
//	    arrayDestroy(msp->gaps);
//	    arrayDestroy(msp->xy);
//	  }
//	}
//    }
//
//  for (msp = MSPlist; msp; )
//    {
//      fmsp = msp;
//      msp = msp->next;
//      g_free(fmsp);
//    }
//
//  for (bpmsp = BPMSPlist.next; bpmsp; )
//    {
//      tmsp = bpmsp;
//      bpmsp = bpmsp->next;
//      g_free(tmsp);
//    }
//
//  arrayDestroy(stringentEntropyarr);
//  arrayDestroy(mediumEntropyarr);
//  arrayDestroy(nonglobularEntropyarr);
//
//  blixemWindow = NULL ;
//  pickedMSP = NULL ;
//
//
//  /* Need to start trying to reset some state.... */
//  dotterStart = dotterEnd = 0 ;
//
//
//  return ;
//}



/* Redraw the entire blixem window. (Call mainWindowRedrawAll directly where possible,
 * rather than this function, which relies on the global variable 'blixemWindow'). */
void blviewRedraw(void)
{
  if (blixemWindow)
    {
      mainWindowRedrawAll(blixemWindow);
    }
}


/* This function is called if an MSP sequence cannot be found. It fills the length
 * of the sequence with padding characters (dashes). If there is already a pad sequence
 * that is long enough it uses that (because any out-of-range characters will be ignored);
 * otherwise it creates a new one of the required length, and makes all other padded MSPs
 * point to the new one too, so that we only have to maintain one. */
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



/* If crosshair-coordinates are screwed up, change here!
 ********************************************************/
//static int x_to_residue(float x)
//{
//  int retval;
//
//  if (blastx || tblastx)
//    retval = dispstart + plusmin*(x - NAMESIZE - 22.3)*3;
//  else
//    retval = dispstart + plusmin*(x - NAMESIZE - 22);
//
//
//  if (plusmin > 0) {
//    if (retval < dispstart) retval = dispstart;
//    if (retval > dispstart+displen-1) retval = dispstart+displen-1;
//    return retval;
//  }
//  else {
//    if (retval > dispstart-1) retval = dispstart-1;
//    if (retval < dispstart-displen) retval = dispstart-displen;
//    return retval +1;
//  }
//}


//static void markDNA(double y)
//{
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
