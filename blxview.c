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
96-05-09  [2.4] Proper grayscale print colors.
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
 * CVS info:   $Id: blxview.c,v 1.40 2010-06-14 11:14:39 gb10 Exp $
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
-3 Any colored segment  -> Big picture
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
#include <SeqTools/blxwindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/blxdotter.h>

#ifdef ACEDB
#include <wh/display.h>
#endif



#define MAXALIGNLEN 10000



/* This is the only place this is set so that you get the same version/compile whether this is
 * compiled stand alone or as part of xace. */
char *blixemVersion = BLIXEM_VERSION_COMPILE ;


static void            blviewCreate(char *opts, char *align_types, const char *paddingSeq, GList *seqList, CommandLineOptions *options) ;
static GList*          getEmptySeqs(GList *inputList);


/* GLOBAL VARIABLES... sigh... */

GtkWidget *blixemWindow = NULL ;
static char *padseq = 0;



/*
 *                Internal routines
 */


//static int possort(MSP *msp1, MSP *msp2) {
//    return ( (msp1->qstart > msp2->qstart) ? 1 : -1 );
//}
//
//static int namesort(MSP *msp1, MSP *msp2) {
//    if (strcmp(msp1->sname, msp2->sname))
//	return strcmp(msp1->sname, msp2->sname);
//    else
//	return possort(msp1, msp2);
//}
//static int scoresort(MSP *msp1, MSP *msp2) {
//    return ( (msp1->score < msp2->score) ? 1 : -1 );
//}
//
//static int idsort(MSP *msp1, MSP *msp2)
//{
//  int result = 0 ;
//
//  result = ( (msp1->id < msp2->id) ? 1 : -1 ) ;
//
//  return result ;
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
//	g_message("Please give a range where from: is less than to:");
//	return;
//    }
//    else if (plusmin < 0 && start < end) {
//	g_message("Please give a range where from: is more than to:");
//	return;
//    }
//    if (displen/symbfact > MAXALIGNLEN) {
//	g_message("Sorry, can't print more than %d residues.  Anyway, think of the paper!", MAXALIGNLEN*symbfact);
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



static void setModeP(CommandLineOptions *options)
{
  options->blastMode = BLXMODE_BLASTP;
  options->bigPictZoom = 10;
}

static void setModeN(CommandLineOptions *options)
{
  options->blastMode = BLXMODE_BLASTN;
  options->bigPictZoom = 30;
}

static void setModeX(CommandLineOptions *options)
{
  options->blastMode = BLXMODE_BLASTX;
  options->bigPictZoom = 10;
}

static void setModeT(CommandLineOptions *options)
{
  options->blastMode = BLXMODE_TBLASTN;
  options->bigPictZoom = 10;
}

static void setModeL(CommandLineOptions *options)
{
  options->blastMode = BLXMODE_TBLASTX;
  options->bigPictZoom = 10;
}


static void blxviewGetOpts(char *opts, char *refSeq, CommandLineOptions *options)
{
  char *opt = opts;
  while (*opt)
    {
      /* Used options: 	 BGILMNPRSTXZ-+brsinp      */
      
      switch (*opt)
      {
	case 'I':
	  options->sortInverted = TRUE;             break;
	case 'G':
	  /* Gapped HSPs */
	  options->hiliteSins = options->gappedHsp = TRUE;    break;
	case 'P': setModeP(options);                break;
	case 'N': setModeN(options);                break;
	case 'X': setModeX(options);                break;
	case 'T': setModeT(options);                break;
	case 'L': setModeL(options);                break;
	case '-':
	  options->activeStrand = BLXSTRAND_REVERSE;   break;
	case '+':
	  options->activeStrand = BLXSTRAND_FORWARD;   break;
	case 'B':
	  options->bigPictON = TRUE;              break; 
	case 'b':
	  options->bigPictON = FALSE;             break;
	case 'd':
	  options->dotterFirst = TRUE;            break;
	case 'i':
	  options->initSortMode = BLXSORT_ID ;      break;
	case 'M':
	  options->startNextMatch = 1;            break;
	case 'n':
	  options->initSortMode = BLXSORT_NAME ;    break;
	case 'p':
	  options->initSortMode = BLXSORT_POS ;     break;
	case 'R':
	  options->bigPictRev = 1;                break; /* to do: this is currently not implemented. not sure what it should do */
	case 'r':
	  options->bigPictRev = 0;                break;
	case 's':
	  options->initSortMode = BLXSORT_SCORE ;   break;
	case 'Z':
	  options->bigPictZoom = strlen(refSeq);  break;
      }
      
      opt++;
    }
  
  if (options->blastMode == BLXMODE_UNSET)
    {
      printf("\nNo blast type specified. Detected ");
      
      if (Seqtype(refSeq) == 'P')
	{
	  printf("protein sequence. Will try to run Blixem in blastp mode\n");
	  setModeP(options);
	}
      else
	{
	  printf("nucleotide sequence. Will try to run Blixem in blastn mode\n");
	  setModeN(options);
	}
    }
  
  /* Determine the sequence type and number of reading frames based on the blast mode */
  if (options->blastMode == BLXMODE_BLASTX || options->blastMode == BLXMODE_TBLASTX)
    {
      options->seqType = BLXSEQ_PEPTIDE;
      options->numFrames = 3;
    }
  else
    {
      options->seqType = BLXSEQ_DNA;
      options->numFrames = 1;
    }
}


/* Find out if we need to fetch any sequences (they may all be contained in the input
 * files), if we do need to, then fetch them by the appropriate method. */
static gboolean blxviewFetchSequences(PfetchParams *pfetch, 
                                      gboolean External, 
                                      GList *seqList, /* list of BlxSequence structs for all required sequences */
                                      CommandLineOptions *options)
{
  gboolean success = TRUE;

  /* First, set the fetch mode */
  if (pfetch)
    {
      /* If pfetch struct then this sets fetch mode to pfetch. */
      
      if (blxConfigSetPFetchSocketPrefs(pfetch->net_id, pfetch->port))
	{
	  options->fetchMode = BLX_FETCH_PFETCH;
	}
    }
  else
    {
      blxFindInitialFetchMode(options->fetchMode);
    }
  
  
  /* See if we have any sequences to fetch (i.e. any that do not have sequence data already) */
  GList *seqsToFetch = getEmptySeqs(seqList);
  
  if (g_list_length(seqsToFetch) > 0)
    {
      if (strcmp(options->fetchMode, BLX_FETCH_PFETCH) == 0)
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
	    {
	      success = blxGetSseqsPfetch(seqsToFetch, net_id, port, External) ;
	    }
	}
#ifdef PFETCH_HTML 
      else if (strcmp(options->fetchMode, BLX_FETCH_PFETCH_HTML) == 0)
	{
	  success = blxGetSseqsPfetchHtml(seqsToFetch, options->seqType) ;
	}
#endif

      g_list_free(seqsToFetch);
    }
  
  return success;
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
int blxview(char *refSeq, char *refSeqName, int start, int qOffset, MSP *msplist, GList *seqList,
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
	  g_error("Config File Error: %s\n", error->message) ;
	}
    }

  char fetchMode[32] = BLX_FETCH_EFETCH;
  const int startCoord = start < 1 ? 1 : start;
  
  CommandLineOptions options = {refSeq, refSeqName, qOffset, startCoord, msplist, stdcode1,
				BLXSTRAND_FORWARD, 10, TRUE, FALSE, BLXSORT_ID, FALSE, FALSE,
				FALSE, FALSE, FALSE, BLXMODE_UNSET, BLXSEQ_INVALID, 1, fetchMode};
  
  blxviewGetOpts(opts, refSeq, &options);
  
  gboolean status = blxviewFetchSequences(pfetch, External, seqList, &options);
  

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(opts, align_types, padseq, seqList, &options) ;
    }

  return 0;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *opts, 
			 char *align_types, 
			 const char *paddingSeq,
                         GList *seqList,
			 CommandLineOptions *options)
{
  if (!blixemWindow)
    {
      /* Create the window */
      blixemWindow = createBlxWindow(options, paddingSeq, seqList);

      BOOL pep_nuc_align = (*opts == 'X' || *opts == 'N') ;
      gtk_window_set_title(GTK_WINDOW(blixemWindow),
			   messprintf("Blixem %s%s%s:   %s",
				      (pep_nuc_align ? "  (" : ""),
				      (align_types ? align_types : (*opts == 'X' ? "peptide" :
								    (*opts == 'N' ? "nucleotide" : ""))),
				      (pep_nuc_align ? " alignment)" : ""),
				      options->refSeqName));
    }

  char *nameSeparatorPos = (char *)strrchr(options->refSeqName, '/');
  if (nameSeparatorPos)
    {
      options->refSeqName = nameSeparatorPos + 1;
    }
  
  if (options->dotterFirst && options->mspList && options->mspList->sname &&
      (options->mspList->type == BLXMSP_HSP || options->mspList->type == BLXMSP_GSP))
    {
      blxWindowSelectSeq(blixemWindow, options->mspList->sSequence);
      
      GError *error = NULL;
      callDotter(blixemWindow, FALSE, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  if (options->startNextMatch)
    {
      /* Set the start coord to be the start of the next MSP on from the default start coord */
      nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
    }
}

/***********************************************************
 *            Create/destroy MSPs and BlxSequences
 ***********************************************************/

/* Add or create a BlxSequence struct, creating the BlxSequence if one does not
 * already exist for the MSP's sequence name. Seperate BlxSequence structs are created
 * for the forward and reverse strands of the same sequence. The passed-in sequence 
 * should always be forwards, and we reverse complement it here if we need the 
 * reverse strand. Returns the new BlxSequence */
BlxSequence* addBlxSequence(MSP *msp, BlxStrand strand, GList **seqList, char *sequence, GError **error)
{
  char *name = mspGetSeqName(msp);

  /* If this is an exon or intron the match strand is not applicable. The exon should 
   * be in the same direction as the ref seq, so use the ref seq strand. */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      strand = msp->qStrand;
    }
  
  /* See if this strand for this sequence already exists. */
  BlxSequence *blxSeq = findBlxSequence(*seqList, name, strand);
  
  if (!blxSeq)
    {
      /* Create a new BlxSequence, and take ownership of the passed in sequence (if any) */
      blxSeq = createEmptyBlxSequence(name);
      *seqList = g_list_prepend(*seqList, blxSeq);
      
      blxSeq->strand = strand;
      
      /* Set the sequence data (if relevant and if passed) */
      blxSeq->sequenceReqd = mspIsBlastMatch(msp) || mspIsSnp(msp);
    }
  
  g_free(name);
  
  /* Add the MSP to the BlxSequence's list */
  blxSeq->mspList = g_list_prepend(blxSeq->mspList, msp);
  
  /* Set a pointer to the BlxSequence from the MSP */
  msp->sSequence = blxSeq;
  
  /* Add the sequence data */
  addBlxSequenceData(blxSeq, sequence, error);
  
  return blxSeq;
}


/* Allocates memory for an MSP and initialise all its fields to the given values.
 * Adds it into the given MSP list and makes the end pointer ('lastMsp') point to the new
 * end of the list. Returns a pointer to the newly-created MSP. Also creates a BlxSequence
 * for this MSP's sequence name (or adds the MSP to the existing one, if it exists already), 
 * and adds that BlxSequence to the given seqList. Takes ownership of 'sequence'. */
MSP* createNewMsp(MSP **lastMsp, 
                  MSP **mspList,
                  GList **seqList,
                  const BlxMspType mspType,
                  const int score,
                  const char *qName,
                  const int qStart,
                  const int qEnd,
                  const BlxStrand qStrand,
                  const int qFrame,
                  const char *sName,
                  const int sStart,
                  const int sEnd,
                  BlxStrand sStrand,
                  char *sequence,
                  char *opts, 
                  GError **error)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  msp->next = NULL;
  msp->type = mspType;
  msp->score = score;
  
  msp->qname = g_strdup(qName);
  
  msp->qFrame = qFrame;
  msp->qStrand = qStrand;
  msp->qstart = qStart;
  msp->qend = qEnd;
  
  msp->sname = g_ascii_strup(sName, -1);
  msp->sstart = sStart;
  msp->send = sEnd;
  
  msp->desc = NULL;
  
  msp->fs = NULL;
  msp->fsColor = 0;
  msp->fsShape = BLXCURVE_BADSHAPE;
  msp->xy = NULL;
  
  msp->gaps = NULL;
  
#ifdef ACEDB
  msp->key = 0;
#endif
  
  /* ID will be calculated later */
  msp->id = 0;
  
  /* For exons and introns, the s strand is not applicable. We always want the exon
   * to be in the same direction as the ref sequence, so set the match seq strand to be 
   * the same as the ref seq strand */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      sStrand = qStrand;
    }
  
  /* Sort coords according to the strand direction */
  sortValues(&msp->qstart, &msp->qend, qStrand == BLXSTRAND_FORWARD);
  sortValues(&msp->sstart, &msp->send, sStrand == qStrand);
  
  /* Dotter still uses the old text versions of strand and frame, so create them here */
  sprintf(msp->qframe, "(%c%d)", getStrandAsChar(qStrand), qFrame);
  sprintf(msp->sframe, "(%c%d)", getStrandAsChar(sStrand), 1);
  
  /* Add/create a BlxSequence for this MSP's sequence name */
  addBlxSequence(msp, sStrand, seqList, sequence, error);
  
  /* Add it to the list */
  if (!*mspList) 
    {
      *mspList = msp; /* first entry in list */
    }
  
  if (*lastMsp)
    {
      (*lastMsp)->next = msp; /* add to end of list */
    }
  
  *lastMsp = msp;
  
  if (error && *error)
    {
      prefixError(*error, "Error creating MSP (ref seq='%s' [%d - %d %s], match seq = '%s' [%d - %d %s]). ",
                  qName, qStart, qEnd, msp->qframe, sName, sStart, sEnd, msp->sframe);
    }
  
  return msp;
}


/* Allocate memory for an MSP and initialise all its fields to relevant 'empty' values.
 * Add it into the given MSP list and make the end pointer ('lastMsp') point to the new
 * end of the list. Returns a pointer to the newly-created MSP */
MSP* createEmptyMsp(MSP **lastMsp, MSP **mspList)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  msp->next = NULL;
  msp->type = BLXMSP_INVALID;
  msp->score = 0;
  msp->id = 0;
  
  msp->qname = NULL;
  
  msp->qframe[0] = '(';
  msp->qframe[1] = '+';
  msp->qframe[2] = '1';
  msp->qframe[3] = ')';
  msp->qframe[4] = '\0';
  
  msp->qstart = 0;
  msp->qend = 0;
  
  msp->sSequence = NULL;
  msp->sname = NULL;
  
  msp->slength = 0;
  msp->sstart = 0;
  msp->send = 0;
  
  msp->sframe[0] = '(';
  msp->sframe[1] = '+';
  msp->sframe[2] = '1';
  msp->sframe[3] = ')';
  msp->sframe[4] = '\0';
  
  msp->desc = NULL;
  
  msp->fs = NULL;
  msp->fsColor = 0;
  msp->fsShape = BLXCURVE_BADSHAPE;
  
  msp->xy = NULL;
  msp->gaps = NULL;
  
#ifdef ACEDB
  msp->key = 0;
#endif
  
  /* Add it to the list */
  if (!*mspList) 
    {
      /* Nothing in the list yet: make this the first entry */
      *mspList = msp;
    }
  
  if (*lastMsp)
    {
      /* Tag it on to the end of the list */
      (*lastMsp)->next = msp;
    }
  
  /* Make the 'lastMsp' pointer point to the new end of the list */
  *lastMsp = msp;
  
  return msp;
}


void destroyMspData(MSP *msp)
{
  g_free(msp->qname);
  msp->qname = NULL;
  
  g_free(msp->sname);
  msp->sname = NULL;
  
  g_free(msp->desc);
  msp->desc = NULL;
  
  g_slist_free(msp->gaps);
  msp->gaps = NULL;
  
  arrayDestroy(msp->xy);
  msp->xy = NULL;
}


/* Destroy all of the MSPs */
void destroyMspList(MSP **mspList)
{
  /* Free the allocated sequences and names */
  MSP *msp = NULL;
  for (msp = *mspList; msp; msp = msp->next)
    {
      destroyMspData(msp);
    }
  
  /* Now free the MSPs themselves. */
  MSP *fmsp = NULL;
  for (msp = *mspList; msp; )
    {
      fmsp = msp;
      msp = msp->next;
      g_free(fmsp);
    }
  
  *mspList = NULL;
  
  return ;
}


/* Destroy all of the BlxSequences */
void destroyBlxSequenceList(GList **seqList)
{
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      destroyBlxSequence(blxSeq);
    }
  
  g_list_free(*seqList);
  *seqList = NULL;
}


/* Redraw the entire blixem window. (Call blxWindowRedrawAll directly where possible,
 * rather than this function, which relies on the global variable 'blixemWindow'). */
void blviewRedraw(void)
{
  if (blixemWindow)
    {
      blxWindowRedrawAll(blixemWindow);
    }
}


void blviewResetGlobals()
{
  blixemWindow = NULL;
}


/* Checks the list of sequences for blixem to display to see if they contain sequence
 * data (this may happen if acedb, for instance, starts blixem up and acedb already contained
 * all the sequences). Returns a list of all those that do not have sequence data. */
static GList* getEmptySeqs(GList *inputList)
{
  GList *resultList = NULL;

  /* Loop through the input list */
  GList *inputItem = inputList;
  
  for ( ; inputItem; inputItem = inputItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(inputItem->data);
      
      /* Only alignments need a sequence */
      if (blxSeq->sequenceReqd && blxSeq->sequence == NULL)
        {
          resultList = g_list_prepend(resultList, blxSeq);
        }
    }

  return resultList;
}


/***************** end of file ***********************/
