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
 * CVS info:   $Id: blxview.c,v 1.63 2010-08-27 12:25:14 gb10 Exp $
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

#include <SeqTools/dotter.h>

#include <math.h>
#include <string.h>
#include <gdk/gdkkeysyms.h>

#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blxwindow.h>
#include <SeqTools/detailview.h>
#include <SeqTools/blxdotter.h>

#ifdef ACEDB
#include <wh/regular.h>
#endif



#define MAXALIGNLEN 10000


/* A struct used to record warning/error messages */
typedef struct _BlxMessage
  {
    char *text;
    time_t timestamp;
    GLogLevelFlags logLevel;
  } BlxMessage;



/* This is the only place this is set so that you get the same version/compile whether this is
 * compiled stand alone or as part of xace. */
char *blixemVersion = BLIXEM_VERSION_COMPILE ;


static void            blviewCreate(char *opts, char *align_types, const char *paddingSeq, GList *seqList, GSList *supportedTypes, CommandLineOptions *options, const char *net_id, int port, const gboolean External) ;
static MSP*            createEmptyMsp(MSP **lastMsp, MSP **mspList);
static void            finaliseBlxSequences(MSP **mspList, GList **seqList, char *opts, const int offset);
static void            processGeneName(BlxSequence *blxSeq);
static void            displayMessageAsPopup(const gchar *message, GLogLevelFlags log_level, gpointer data);
static void            displayMessageAsList(GSList *messageList, const gboolean bringToFront, gpointer data);
static BlxMessage*     createBlxMessage(const char *text, const GLogLevelFlags logLevel);
static void            destroyBlxMessage(BlxMessage **blxMessage);
static GSList**        getMessageList();
void                   destroyMessageList();
static char*           blxMessageGetDisplayText(const BlxMessage *msg, const gboolean incTimestamp);
static void            addMessagesToBuffer(GSList *messageList, GtkTextBuffer *textBuffer, GtkTextTag *normalTag, GtkTextTag *highlightTag);
static gboolean        getUseScrolledMessages();
static void            setUseScrolledMessages(const gboolean newValue);
static void            onSetUseScrolledMessages(GtkWidget *button, const gint responseId, gpointer data);
static void            onSetUsePopupMessages(GtkWidget *button, const gint responseId, gpointer data);
static char*           getDialogIcon(GLogLevelFlags log_level);
static void            printMessageToStatusbar(const gchar *message, gpointer data);


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
      /* Used options: 	 BFGILMNPRSTXZ-+bgrstimnop      */
      
      switch (*opt)
      {
	case 'G':
	  /* Gapped HSPs */
	  options->hiliteSins = options->gappedHsp = TRUE;    break;
          
	case 'P': setModeP(options);                  break;
	case 'N': setModeN(options);                  break;
	case 'X': setModeX(options);                  break;
	case 'T': setModeT(options);                  break;
	case 'L': setModeL(options);                  break;
          
	case '-':
	  options->activeStrand = BLXSTRAND_REVERSE; break;
	case '+':
	  options->activeStrand = BLXSTRAND_FORWARD; break;
	case 'B':
	  options->bigPictON = TRUE;                 break; 
	case 'b':
	  options->bigPictON = FALSE;                break;
	case 'd':
	  options->dotterFirst = TRUE;               break;
        case 'F':
          options->parseFullEmblInfo = TRUE;         break;
	case 'M':
	  options->startNextMatch = 1;               break;
	case 'R':
	  options->bigPictRev = 1;                   break; /* to do: this is currently not implemented. not sure what it should do */
	case 'r':
	  options->bigPictRev = 0;                   break;
	case 'Z':
	  options->bigPictZoom = strlen(refSeq);     break;

          /* Sorting... */
        case 'I':
	  options->sortInverted = TRUE;                   break;
	case 'i':
	  options->initSortColumn = BLXCOL_ID ;           break;
	case 'n':
	  options->initSortColumn = BLXCOL_SEQNAME ;      break;
	case 'p':
	  options->initSortColumn = BLXCOL_START ;        break;
	case 's':
	  options->initSortColumn = BLXCOL_SCORE ;        break;
	case 't':
	  options->initSortColumn = BLXCOL_TISSUE_TYPE ;  break;
	case 'm':
	  options->initSortColumn = BLXCOL_STRAIN ;       break;
	case 'g':
	  options->initSortColumn = BLXCOL_GENE_NAME ;    break;
	case 'o':
	  options->initSortColumn = BLXCOL_ORGANISM ;     break;
      }
      
      opt++;
    }

  if (!options->refSeq)
    {
      g_error("Error: reference sequence not specified.\n");
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
  
  /* For DNA matches, always get the full EMBL file, because it doesn't seem to be any slower
   * (in fact, it seems to be quicker).  For protein matches, only get the full file if requested. */
  if (options->seqType == BLXSEQ_DNA)
    {
      options->parseFullEmblInfo = TRUE;
    }
}


/* Performs the work of fetching the given list of sequence via pfetch */
static gboolean pfetchSequences(GList *seqsToFetch, 
                                GList *seqList, 
                                const char *net_id, 
                                int port,
                                const gboolean parseOptionalData,
                                const gboolean parseSequenceData,
                                const gboolean External,
                                const BlxSeqType seqType,
                                GError **error)
{
  gboolean success = FALSE;

  if (net_id && port != UNSET_INT)
    {
      if (parseSequenceData && !parseOptionalData)
        {
          /* Just parse the minimal EMBL file containing just the fasta sequence */
          success = populateFastaDataPfetch(seqsToFetch, net_id, port, External, seqType, error) ;
        }
      else if (parseOptionalData)
        {
          /* This first call tries fetching the full EMBL entry and parses
           * additional information such as organism and tissue-type */
          success = populateFullDataPfetch(seqsToFetch, net_id, port, External, seqType, error) ;
          
          /* If we were requested to populate the sequence fasta data but there are still 
           * some empty sequences (which is likely in protein matches because pfetch -F 
           * fails for protein variants), try fetching the EMBL entry containing just the 
           * fasta sequence. */
          if (success && parseSequenceData)
            {
              GList *remainingSeqs = getSeqsToPopulate(seqList, TRUE, FALSE);
              
              if (g_list_length(remainingSeqs) > 0)
                {
                  success = populateFastaDataPfetch(remainingSeqs, net_id, port, External, seqType, error) ;
                }
              
              g_list_free(remainingSeqs);
            }
        }
    }
  else
    {
      g_set_error(error, BLX_ERROR, 1, "pfetch is not initialised: net_id = %s, port = %d.\n", net_id, port);
    }
 
  return success;   
}


#ifdef PFETCH_HTML
static gboolean pfetchSequencesHttp(GList *seqsToFetch, 
                                    GList *seqList, 
                                    const gboolean parseOptionalData,
                                    const gboolean parseSequenceData,
                                    const BlxSeqType seqType,
                                    GError **error)
{
  gboolean success = FALSE;

  if (parseSequenceData && !parseOptionalData)
    {
      success = populateSequenceDataHtml(seqsToFetch, seqType, FALSE);
    }
  else if (parseOptionalData)
    {
      /* This first call tries fetching the full EMBL entry and parses
       * additional information such as organism and tissue-type */
      success = populateSequenceDataHtml(seqsToFetch, seqType, TRUE);
      
      /* If we were requested to populate the sequence fasta data but there are still 
       * some empty sequences (which is likely in protein matches because pfetch -F 
       * fails for protein variants), try fetching the EMBL entry containing just the 
       * fasta sequence. */
      if (success && parseSequenceData)
        {
          GList *remainingSeqs = getSeqsToPopulate(seqList, TRUE, FALSE);
          
          if (g_list_length(remainingSeqs) > 0)
            {
              success = populateSequenceDataHtml(remainingSeqs, seqType, FALSE) ;
            }
          
          g_list_free(remainingSeqs);
        }
    }

  return success;
}
#endif


/* Utility to get the parent sequence of the given variant, if it is not already set.
 * returns true if the returned parent is not null. */
static gboolean getParent(BlxSequence *variant, BlxSequence **parent, GList *allSeqs)
{
  if (*parent == NULL)
    {
      *parent = blxSequenceGetVariantParent(variant, allSeqs);
    }
  
  return (*parent != NULL);
}


/* This function checks if the given sequence is missing its optional data (such as organism
 * and gene name) and, if so, looks for the parent sequence and copies the data from there. */
static void populateMissingDataFromParent(BlxSequence *curSeq, GList *seqList)
{
  BlxSequence *parent = NULL;
  
  /* For each optional column, check if the data in the BlxSequence is null */
  if (!curSeq->organism && getParent(curSeq, &parent, seqList) && parent->organism && parent->organism->str)
    {
      curSeq->organism = g_string_new(parent->organism->str);
    }
  
  if (!curSeq->geneName && getParent(curSeq, &parent, seqList) && parent->geneName && parent->geneName->str)
    {
      curSeq->geneName = g_string_new(parent->geneName->str);
    }
  
  if (!curSeq->tissueType && getParent(curSeq, &parent, seqList) && parent->tissueType && parent->tissueType->str)
    {
      curSeq->tissueType = g_string_new(parent->tissueType->str);
    }
  
  if (!curSeq->strain && getParent(curSeq, &parent, seqList) && parent->strain && parent->strain->str)
    {
      curSeq->strain = g_string_new(parent->strain->str);
    }
}


/* This function actually performs the work of fetching the given list of sequences.
 * The first argument is the list of sequences to fetch and the second is the list of all
 * sequences. */
gboolean fetchSequences(GList *seqsToFetch, 
                        GList *seqList,
                        char *fetchMode, 
                        const BlxSeqType seqType,
                        const char *net_id, 
                        int port, 
                        const gboolean parseOptionalData,
                        const gboolean parseSequenceData,
                        const gboolean External,
                        GError **error)
{
  gboolean success = TRUE;
  
  if (g_list_length(seqsToFetch) > 0)
    {
      if (strcmp(fetchMode, BLX_FETCH_PFETCH) == 0)
        {  
          success = pfetchSequences(seqsToFetch, seqList, net_id, port, parseOptionalData, parseSequenceData, External, seqType, error);
	}
#ifdef PFETCH_HTML 
      else if (strcmp(fetchMode, BLX_FETCH_PFETCH_HTML) == 0)
	{
          success = pfetchSequencesHttp(seqsToFetch, seqList, parseOptionalData, parseSequenceData, seqType, error);
	}
#endif
      else
        {
          g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", fetchMode);
        }
    }
    
  
  /* Once we've fetched all the sequences we need to do some post-processing */
  GList *seqItem = seqList;
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      populateMissingDataFromParent(blxSeq, seqList);
      processGeneName(blxSeq);
    }
  
  return success;
}


/* Find out if we need to fetch any sequences (they may all be contained in the input
 * files), if we do need to, then fetch them by the appropriate method. */
static gboolean blxviewFetchSequences(PfetchParams *pfetch, 
                                      gboolean External, 
                                      const gboolean parseFullEmblInfo,
                                      const BlxSeqType seqType,
                                      GList *seqList, /* list of BlxSequence structs for all required sequences */
                                      char **fetchMode,
                                      const char **net_id,
                                      int *port)
{
  gboolean success = TRUE;

  setupFetchMode(pfetch, fetchMode, net_id, port);
  
  /* Fetch any sequences that do not have their sequence data already populated (or
   * optional data too, if requested). */
  GList *seqsToFetch = getSeqsToPopulate(seqList, TRUE, parseFullEmblInfo);

  GError *error = NULL;
  success = fetchSequences(seqsToFetch, seqList, *fetchMode, seqType, *net_id, *port, parseFullEmblInfo, TRUE, External, &error);
  
  g_list_free(seqsToFetch);
  
  if (error && !success)
    {
      prefixError(error, "Error fetching sequences. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  else if (error)
    {
      prefixError(error, "Warning: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
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
gboolean blxview(char *refSeq, 
                 char *refSeqName, 
                 int start, 
                 int qOffset, 
                 MSP *msplist, 
                 GList *seqList,
                 GSList *supportedTypes,
                 char *opts, 
                 PfetchParams *pfetch, 
                 char *align_types, 
                 gboolean External)
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
				BLXSTRAND_FORWARD, 10, TRUE, FALSE, BLXCOL_ID, FALSE, FALSE,
				FALSE, FALSE, FALSE, FALSE, BLXMODE_UNSET, BLXSEQ_INVALID, 1, fetchMode};
  
  blxviewGetOpts(opts, refSeq, &options);
  
  const char *net_id = NULL;
  int port = UNSET_INT;
  gboolean status = blxviewFetchSequences(pfetch, External, options.parseFullEmblInfo, options.seqType, seqList, &options.fetchMode, &net_id, &port);
  
  if (status)
    {
      /* Construct missing data and do any other required processing now we have all the sequence data */
      finaliseBlxSequences(&msplist, &seqList, opts, qOffset);
    }

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(opts, align_types, padseq, seqList, supportedTypes, &options, net_id, port, External) ;
    }
  
  return status;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *opts, 
			 char *align_types, 
			 const char *paddingSeq,
                         GList *seqList,
                         GSList *supportedTypes,
			 CommandLineOptions *options,
                         const char *net_id,
                         int port,
                         const gboolean External)
{
  if (!blixemWindow)
    {
      /* Create the window */
      blixemWindow = createBlxWindow(options, paddingSeq, seqList, supportedTypes, net_id, port, External);

      gboolean pep_nuc_align = (*opts == 'X' || *opts == 'N') ;
      
      char *title = blxprintf("Blixem %s%s%s:   %s",
                              (pep_nuc_align ? "  (" : ""),
                              (align_types ? align_types : (*opts == 'X' ? "peptide" :
                                                            (*opts == 'N' ? "nucleotide" : ""))),
                              (pep_nuc_align ? " alignment)" : ""),
                              options->refSeqName);
      
      gtk_window_set_title(GTK_WINDOW(blixemWindow), title);
      g_free(title);
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

/* Compare the start position in the ref seq of two MSPs. Returns a negative value if a < b; zero
 * if a = b; positive value if a > b. Secondarily sorts by type in the order that types appear in 
 * the BlxMspType enum. */
static gint compareFuncMspPos(gconstpointer a, gconstpointer b)
{
  gint result = 0;

  const MSP const *msp1 = (const MSP const*)a;
  const MSP const *msp2 = (const MSP const*)b;
  
  if (msp1->qRange.min == msp2->qRange.min)
    {
      /* Sort by type. Lower type numbers should appear first. */
      result = msp2->type - msp1->type;
    }
  else 
    {
      result = msp1->qRange.min -  msp2->qRange.min;
    }
  
  return result;
}

/* Add or create a BlxSequence struct, creating the BlxSequence if one does not
 * already exist for the MSP's sequence name. Seperate BlxSequence structs are created
 * for the forward and reverse strands of the same sequence. The passed-in sequence 
 * should always be forwards, and we reverse complement it here if we need the 
 * reverse strand. Returns the new BlxSequence */
BlxSequence* addBlxSequence(const char *name, const char *idTag, BlxStrand strand, GList **seqList, char *sequence, MSP *msp, GError **error)
{
  BlxSequence *blxSeq = NULL;
  
  if (name || idTag)
    {
      /* If this is an exon or intron the match strand is not applicable. The exon should 
       * be in the same direction as the ref seq, so use the ref seq strand. */
      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
        {
          strand = msp->qStrand;
        }
    
      /* See if this strand for this sequence already exists. Horrible hack for backwards compatibility:
       * if the msp is an exon/intron, cut off the old-style 'x' or 'i' postfix from the name, if it has one. */
      char *seqName = msp ? mspGetExonTranscriptName(msp) : g_strdup(name);
      blxSeq = findBlxSequence(*seqList, seqName, idTag, strand);
      
      if (!blxSeq)
        {
          /* Create a new BlxSequence, and take ownership of the passed in sequence (if any) */
          blxSeq = createEmptyBlxSequence(seqName, idTag, NULL);
          *seqList = g_list_prepend(*seqList, blxSeq);
          
          blxSeq->strand = strand;
          
          /* Set whether the sequence data is required by any of this sequence's MSPs */
          blxSeq->sequenceReqd |= mspIsBlastMatch(msp) || mspIsPolyATail(msp) || mspIsSnp(msp);
        }
      
      if (seqName && !blxSeq->fullName)
	{
	  /* It's possible that the BlxSequence was created without a name if we found an
	   * unnamed child exon before we found the parent transcript, so set the name if we have it. */
	  blxSequenceSetName(blxSeq, seqName);
	}
      
      if (msp)
        {
          /* Add the MSP to the BlxSequence's list. Keep it sorted by position. */
          blxSeq->mspList = g_list_insert_sorted(blxSeq->mspList, msp, compareFuncMspPos);
          msp->sSequence = blxSeq;
        }
      
      /* Add the sequence data */
      addBlxSequenceData(blxSeq, sequence, error);
      
      if (msp)
        {
          g_free(seqName);
        }
    }
  else
    {
      g_set_error(error, BLX_ERROR, 1, "Sequence name or parent ID must be set.\n");
    }
  
  return blxSeq;
}


/* Exons and UTRs don't have phase, but we want to display them in the same reading frame
 * as the CDS, if there is one, so if we have a corresponding CDS copy its reading frame data. */
static void copyCdsReadingFrame(MSP *exon, MSP *cds, MSP *utr)
{
  if (cds && exon)
    {
      exon->qFrame = cds->qFrame;
      exon->phase = cds->phase;
    }

  if (cds && utr)
    {
      utr->qFrame = cds->qFrame;
      utr->phase = cds->phase;
    }
}


/* Utility used by constructTranscriptData to create a missing exon/cds/utr given 
 * two others out of the three - i.e. if we have an overlapping exon and cds we can
 * construct the corresponding utr. If created, the new msp is added to the given 
 * BlxSequence and the  MSP list. If a CDS is given and no UTR exists, assume the exon
 * spans the entire CDS (and similarly if a UTR is given but no CDS exists) */
static void createMissingExonCdsUtr(MSP **exon, MSP **cds, MSP **utr, 
                                    BlxSequence *blxSeq, MSP **lastMsp, MSP **mspList, GList **seqList, 
                                    char *opts, GError **error)
{
  BlxMspType newType = BLXMSP_INVALID;
  MSP **ptrToUpdate = NULL;
  int newStart = UNSET_INT;
  int newEnd = UNSET_INT;
  int newPhase = UNSET_INT;

  if (!*exon && (*cds || *utr))
    {
      /* The exon spans both the cds and utr */
      if (*cds && *utr)
        {
          newStart = min((*cds)->qRange.min, (*utr)->qRange.min);
          newEnd = max((*cds)->qRange.max, (*utr)->qRange.max);
          newPhase = (*cds)->phase;
          newType = BLXMSP_EXON;
          ptrToUpdate = exon;
        }
      else if (*cds)
        {
          newStart = (*cds)->qRange.min;
          newEnd = (*cds)->qRange.max;
          newPhase = (*cds)->phase;
          newType = BLXMSP_EXON;
          ptrToUpdate = exon;
        }
      else
        {
          newStart = (*utr)->qRange.min;
          newEnd = (*utr)->qRange.max;
          newPhase = (*utr)->phase;
          newType = BLXMSP_EXON;
          ptrToUpdate = exon;
        }
    }
  else if (!*cds && *exon && *utr)
    {
      /* The cds is the range of the exon that is not in the utr */
      if ((*utr)->qRange.max < (*exon)->qRange.max)
        {
          newStart = (*utr)->qRange.max + 1;
          newEnd = (*exon)->qRange.max;
          newPhase = (*exon)->phase;
          newType = BLXMSP_CDS;
          ptrToUpdate = cds;
        }
      else if ((*exon)->qRange.min < (*utr)->qRange.min)
        {
          newStart = (*exon)->qRange.min;
          newEnd = (*utr)->qRange.min - 1;
          newPhase = (*exon)->phase;
          newType = BLXMSP_CDS;
          ptrToUpdate = cds;
        }
    }
  else if (!*utr && *exon && *cds)
    {
      /* The utr is the range of the exon that is not in the cds */
      if ((*exon)->qRange.min < (*cds)->qRange.min)
        {
          newStart = (*exon)->qRange.min;
          newEnd = (*cds)->qRange.min - 1;
          newPhase = (*cds)->phase;
          newType = BLXMSP_UTR;
          ptrToUpdate = utr;
        }
      else if ((*cds)->qRange.max < (*exon)->qRange.max)
        {
          newStart = (*cds)->qRange.max + 1;
          newEnd = (*exon)->qRange.max;
          newPhase = (*cds)->phase;
          newType = BLXMSP_UTR;
          ptrToUpdate = utr;
        }
    }
  
  if (newType != BLXMSP_INVALID)
    {
      /* Create the new exon/cds/utr */
      DEBUG_OUT("Creating MSP for transcript '%s' of type %d.\n", blxSeq->fullName, newType);
      
      GError *tmpError = NULL;
      
      *ptrToUpdate = createNewMsp(lastMsp, mspList, seqList, newType, NULL, 0,  newPhase, blxSeq->idTag,
                                  NULL, newStart, newEnd, blxSeq->strand, UNSET_INT, blxSeq->fullName,
                                  UNSET_INT, UNSET_INT, blxSeq->strand, NULL, opts, &tmpError);
      
      if (tmpError)
        {
          prefixError(tmpError, "Error constructing missing exon/cds/utr [type='%d']", newType);
          g_propagate_error(error, tmpError);
        }
    }
}


/* Construct any missing transcript data, i.e.
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructTranscriptData(BlxSequence *blxSeq, MSP **lastMsp, MSP **mspList, GList **seqList, char *opts)
{
  GError *tmpError = NULL;

  const MSP *prevMsp = NULL;
  const MSP *prevExon = NULL;
  
  MSP *curExon = NULL;
  MSP *curCds = NULL;
  MSP *curUtr = NULL;

  /* Loop through all MSPs on this sequence (which must be sorted by position on the
   * ref seq - createNewMsp automatically sorts them for us) and create any missing 
   * exon/cds/utr/introns. */
  GList *mspItem = blxSeq->mspList;
  gboolean finished = FALSE;
  
  while (!finished)
    {
      MSP *msp = mspItem ? (MSP*)(mspItem->data) : NULL;
              
      /* Only consider exons and introns */
      if (mspIsExon(msp) || mspIsIntron(msp) || !msp)
        {
          /* See if there was a gap between this exon and the previous one. There's a gap if
           * we have two exons with space between them, or if we're at the first or last exon 
           * and there's a gap to the end of the transcript. */
          gboolean foundGap = FALSE;
          
          if (msp && 
              ((prevMsp && prevMsp->qRange.max < msp->qRange.min) ||
               (!prevMsp && blxSeq->qRange.min < msp->qRange.min) ||
               ((mspItem->next == NULL && blxSeq->qRange.max > msp->qRange.max))))
            {
              foundGap = TRUE;
            }

          if (foundGap || msp == NULL)
            {
              /* We've found a gap between exons, or reached the end. First, see if the current exon/cds or utr
               * is missing and construct it if possible. Also do this if we're at the last MSP. */
              createMissingExonCdsUtr(&curExon, &curCds, &curUtr, blxSeq, lastMsp, mspList, seqList, opts, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
              copyCdsReadingFrame(curExon, curCds, curUtr);

              /* Create an intron in the gap, unless there's already one here */
              if (prevExon && curExon && !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  createNewMsp(lastMsp, mspList, seqList, BLXMSP_INTRON, NULL, UNSET_INT, UNSET_INT, blxSeq->idTag,
                               NULL, prevExon->qRange.max, curExon->qRange.min, blxSeq->strand, UNSET_INT, blxSeq->fullName,
                               UNSET_INT, UNSET_INT, blxSeq->strand, NULL, opts, &tmpError);
                  
                  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
                }

              /* We're moving past the gap, so reset the pointers */
              prevExon = curExon;
              curExon = NULL;
              curCds = NULL;
              curUtr = NULL;
            }
          
          if (msp && msp->type == BLXMSP_EXON)
            curExon = msp;
          else if (msp && msp->type == BLXMSP_CDS)
            curCds = msp;
          else if (msp && msp->type == BLXMSP_UTR)
            curUtr = msp;

          /* Remember the last MSP we saw */
          prevMsp = msp;
          
          /* Proceed to the next MSP. We allow an extra loop with a NULL mspItem at the end, and then finish. */
          if (mspItem)
            mspItem = mspItem->next;
          else
            finished = TRUE;
        } 
      else
        {
          /* Something that's not an exon/intron. Skip this BlxSequence. */
          finished = TRUE;
        }
    }
}


/* Given a list of msps, find the min or max s coords, according to the given flag. */
int findMspListSExtent(GList *mspList, const gboolean findMin)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const*)(mspItem->data);
      
      if (first)
	{
	  result = findMin ? msp->sRange.min : msp->sRange.max;
	  first = FALSE;
	}
      else if (findMin && msp->sRange.min < result)
	{
	  result = msp->sRange.min;
	}
      else if (!findMin && msp->sRange.max > result)
	{
	  result = msp->sRange.max;
	}
    }
    
  return result;
}


/* Given a list of msps, find the min or max q coords, according to the given flag. */
int findMspListQExtent(GList *mspList, const gboolean findMin)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const*)(mspItem->data);
      
      if (first)
	{
	  result = findMin ? msp->qRange.min : msp->qRange.max;
	  first = FALSE;
	}
      else if (findMin && msp->qRange.min < result)
	{
	  result = msp->qRange.min;
	}
      else if (!findMin && msp->qRange.max > result)
	{
	  result = msp->qRange.max;
	}
    }
    
  return result;
}


/* The parsed gene name generally contains extra info that we're not interested in. This function
 * tries to extract just the part of the info that we're interested in. */
static void processGeneName(BlxSequence *blxSeq)
{
  /* Data is usually in the format: "Name=SMARCA2; Synonyms=BAF190B, BRM, SNF2A, SNF2L2;"
   * Therefore, look for the Name tag and extract everything up to the ; or end of line.
   * If there is no name tag, just include everything */
  
  if (blxSeq && blxSeq->geneName && blxSeq->geneName->str)
    {
      char *startPtr = strstr(blxSeq->geneName->str, "Name=");
      
      if (!startPtr)
        {
          startPtr = strstr(blxSeq->geneName->str, "name=");
        }
        
      if (startPtr)
        {
          startPtr += 5;
          char *endPtr = strchr(startPtr, ';');
          const int numChars = endPtr ? endPtr - startPtr : strlen(startPtr);
          
          char *result = g_malloc((numChars + 1) * sizeof(char));
          
          g_utf8_strncpy(result, startPtr, numChars);
          result[numChars] = 0;
          
          g_string_free(blxSeq->geneName, TRUE);
          blxSeq->geneName = g_string_new(result);
        }
    }
}


/* Find the first/last base in an entire sequence, if not already set. For transcripts, the
 * range should already be set from the parent transcript item. For matches, get the start/end
 * of the first/last MSP in the sequence */
static void findSequenceExtents(BlxSequence *blxSeq)
{
  blxSeq->qRange.min = findMspListQExtent(blxSeq->mspList, TRUE);
  blxSeq->qRange.max = findMspListQExtent(blxSeq->mspList, FALSE);
}


/* Adjust the MSP's q coords (as parsed from the input file) by the given offset i.e. to convert
 * them to real coords */
static void adjustMspCoordsByOffset(MSP *msp, const int offset)
{
  /* Convert the input coords (which are 1-based within the ref sequence section
   * that we're dealing with) to "real" coords (i.e. coords that the user will see). */
  msp->qRange.min += offset;
  msp->qRange.max += offset;

  /* Gap coords are also 1-based, so convert those too */
  GSList *rangeItem = msp->gaps;

  for ( ; rangeItem; rangeItem = rangeItem->next)
    {
      CoordRange *curRange = (CoordRange*)(rangeItem->data);
      curRange->qStart += offset;
      curRange->qEnd += offset;
    }
}


/* Should be called after all parsed data has been added to a BlxSequence. Calculates summary
 * data and the introns etc. */
static void finaliseBlxSequences(MSP **mspList, GList **seqList, char *opts, const int offset)
{
  /* Loop through all MSPs and adjust their coords by the offest. Remember the last in the list */
  MSP *msp = *mspList;
  MSP *lastMsp = msp;
  
  while (msp)
    {
      adjustMspCoordsByOffset(msp, offset);

      msp = msp->next;
      if (msp)
	lastMsp = msp;
    }

  /* Loop through all BlxSequences */
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      
      /* So far we only have the forward strand version of each sequence. We must complement any 
       * that need the reverse strand */
      if (blxSeq && blxSeq->strand == BLXSTRAND_REVERSE && blxSeq->sequence && blxSeq->sequence->str)
        {
          blxComplement(blxSeq->sequence->str);
        }
      
      findSequenceExtents(blxSeq);
      constructTranscriptData(blxSeq, &lastMsp, mspList, seqList, opts);
    }
}


/* Returns true if this msp is part of a feature series */
gboolean mspHasFs(const MSP *msp)
{
  gboolean result = (msp->type == BLXMSP_FS_SEG || msp->type == BLXMSP_XY_PLOT);
  return result;
}


/* Add the given sequence data to a BlxSequence. Validates that the existing sequence data is
 * either null or is the same as the new data; sets the given error if not. We claim ownership
 * of the given sequence data (either the BlxSequence owns it, or we delete it if it is not 
 * required). The given sequence should always be the forward strand; we complement it ourselves
 * here if this BlxSequence requires the reverse strand. */
void addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error)
{
  if (!sequence)
    {
      return;
    }
  
  gboolean sequenceUsed = FALSE;
  
  if (blxSeq && blxSeq->sequenceReqd)
    {
      if (!blxSeq->sequence)
        {
          /* Sequence does not yet exist, so add it */
          blxSeq->sequence = g_string_new(sequence);
          sequenceUsed = TRUE;
          }
      else if (error && *error)
        {
          /* Sequence already exists. Validate that it's the same as the existing one. */
          if (!stringsEqual(sequence, blxSeq->sequence->str, FALSE))
            {
              g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_DATA_MISMATCH, "Sequence data for '%s' does not match previously-found data.\n", blxSeq->fullName);
            }
        }
    }      
  
  if (!sequenceUsed)
    {
      g_free(sequence);
    }
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
		  char *source,
                  const gdouble score,
                  const int phase,
		  char *idTag,
                  char *qName,
                  const int qStart,
                  const int qEnd,
                  const BlxStrand qStrand,
                  const int qFrame,
                  char *sName,
                  const int sStart,
                  const int sEnd,
                  BlxStrand sStrand,
                  char *sequence,
                  char *opts, 
                  GError **error)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = mspType;
  msp->score = score; 
  msp->source = source;
  msp->phase = phase;
  
  msp->qname = qName;
  
  msp->qFrame = qFrame;
  msp->qStrand = qStrand;
  
  msp->sname = sName ? g_ascii_strup(sName, -1) : NULL;
  
  intrangeSetValues(&msp->qRange, qStart, qEnd);  
  intrangeSetValues(&msp->sRange, sStart, sEnd);
  
  /* Check that a phase wasn't passed in for anything that's not an exon. (Technically it
   * should be anything that's not a CDS, but we copy the CDS phase to exons/UTRs so we 
   * can calculate which reading frame to display them in.) */
  if (!mspIsExon(msp) && msp->phase != UNSET_INT)
    {
      g_warning("Warning: MSP '%s' (%d-%d) has phase '%d' specified but only CDS types should have phase (MSP type=%d).\n", msp->sname, msp->qRange.min, msp->qRange.max, msp->phase, msp->type);
    }

  /* Alignments don't have a phase, so set it to 0 so that it gets ignored in our reading-frame calculation. */
  if (mspIsBlastMatch(msp) && msp->phase == UNSET_INT)
    {
      msp->phase = 0;
    }
  
  /* For exons and introns, the s strand is not applicable. We always want the exon
   * to be in the same direction as the ref sequence, so set the match seq strand to be 
   * the same as the ref seq strand */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      sStrand = qStrand;
    }
  
  /* Dotter still uses the old text versions of strand and frame, so create them here */
  sprintf(msp->qframe, "(%c%d)", getStrandAsChar(qStrand), qFrame);
  sprintf(msp->sframe, "(%c%d)", getStrandAsChar(sStrand), 1);
  
  /* Add/create a BlxSequence for this MSP's sequence name */
  addBlxSequence(msp->sname, idTag, sStrand, seqList, sequence, msp, error);

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
static MSP* createEmptyMsp(MSP **lastMsp, MSP **mspList)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  msp->next = NULL;
  msp->type = BLXMSP_INVALID;
  msp->score = 0.0;
  msp->id = 0.0;
  msp->phase = UNSET_INT;
  
  msp->qname = NULL;
  msp->qFrame = UNSET_INT;
  
  msp->qframe[0] = '(';
  msp->qframe[1] = '+';
  msp->qframe[2] = '1';
  msp->qframe[3] = ')';
  msp->qframe[4] = '\0';
  
  msp->qRange.min = 0;
  msp->qRange.max = 0;
  
  msp->sSequence = NULL;
  msp->sname = NULL;
  
  msp->sRange.min = 0;
  msp->sRange.max = 0;
  
  msp->sframe[0] = '(';
  msp->sframe[1] = '+';
  msp->sframe[2] = '1';
  msp->sframe[3] = ')';
  msp->sframe[4] = '\0';
  
  msp->desc = NULL;
  msp->source = NULL;

  msp->style = NULL;

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
  if (msp->qname)
    {
      g_free(msp->qname);
      msp->qname = NULL;
    }
    
  if (msp->sname)
    {
      g_free(msp->sname);
      msp->sname = NULL;
    }
    
  if (msp->desc)
    {
      g_free(msp->desc);
      msp->desc = NULL;
    }
    
  if (msp->gaps)
    {
      g_slist_free(msp->gaps);
      msp->gaps = NULL;
    }
  
  if (msp->xy)
    {
      g_array_free(msp->xy, TRUE);
      msp->xy = NULL;
    }
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


/* Reset any global variables / singleton instances */
void blviewResetGlobals()
{
  blixemWindow = NULL;
  destroyMessageList();
}


/* Checks the list of sequences for blixem to display to see which ones need populating with
 * sequence data (if getSequenceData is true) and/or the optional-columns' data (if getOptionalData
 * is true. Returns a list of all the sequences that require the requested data. */
GList* getSeqsToPopulate(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData)
{
  GList *resultList = NULL;

  /* Loop through the input list */
  GList *inputItem = inputList;
  
  for ( ; inputItem; inputItem = inputItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(inputItem->data);
      
      gboolean getSeq = FALSE;
      
      /* We only want to get data for blast matches, which have the sequenceReqd flag set. */
      if (blxSeq->sequenceReqd)
        {
          /* Check if sequence data was requested and is not already set. */
          getSeq = (getSequenceData && blxSeq->sequence == NULL);

          /* Check if optional data was requested and is not already set. We can assume that
           * if any of the data fields is set then the parsing has been done for all of them
           * (and any remaining empty fields just don't have that data available) */
          getSeq |= (getOptionalData && 
                     !blxSeq->organism &&
                     !blxSeq->geneName &&
                     !blxSeq->tissueType &&
                     !blxSeq->strain);
          
          if (getSeq)
            {
              resultList = g_list_prepend(resultList, blxSeq);
            }
          }
    }

  return resultList;
}



/* Create a BlxStyle */
BlxStyle* createBlxStyle(const char *styleName, 
			 const char *fillColor, 
			 const char *fillColorSelected, 
			 const char *fillColorPrint, 
			 const char *fillColorPrintSelected, 
			 const char *lineColor, 
			 const char *lineColorSelected, 
			 const char *lineColorPrint, 
			 const char *lineColorPrintSelected, 
			 GError **error)
{
  BlxStyle *style = g_malloc(sizeof(BlxStyle));
  GError *tmpError = NULL;

  if (!styleName)
    {
      g_set_error(error, BLX_ERROR, 1, "Style name is NULL.\n");
    }
    
  if (!tmpError)
    {
      style->styleName = g_strdup(styleName);
    }
  
  if (!tmpError && fillColor)
    {
      getColorFromString(fillColor, &style->fillColor.normal, &tmpError);
      getSelectionColor(&style->fillColor.normal, &style->fillColor.selected); /* default if selection color not given */
      convertToGrayscale(&style->fillColor.normal, &style->fillColor.print); /* default if print color not given */
      getSelectionColor(&style->fillColor.print, &style->fillColor.printSelected); /* default if selected-print color not given */
    }
  
  if (!tmpError && lineColor)
    {
      getColorFromString(lineColor, &style->lineColor.normal, &tmpError);
      getSelectionColor(&style->lineColor.normal, &style->lineColor.selected); /* default if selection color not given */
      convertToGrayscale(&style->lineColor.normal, &style->lineColor.print); /* default if print color not given */
      getSelectionColor(&style->lineColor.print, &style->lineColor.printSelected); /* default if selected-print color not given */
    }

  if (!tmpError && fillColorSelected)
    {
      getColorFromString(fillColorSelected, &style->fillColor.selected, &tmpError);
    }
    
    if (!tmpError && lineColorSelected)
    {
      getColorFromString(lineColorSelected, &style->lineColor.selected, &tmpError);
    }

  if (!tmpError && fillColorPrint)
    {
      getColorFromString(fillColorPrint, &style->fillColor.print, &tmpError);
    }

  if (!tmpError && lineColorPrint)
    {
      getColorFromString(lineColorPrint, &style->lineColor.print, &tmpError);
    }


  if (!tmpError && fillColorPrintSelected)
    {
      getColorFromString(fillColorPrintSelected, &style->fillColor.printSelected, &tmpError);
    }

  if (!tmpError && lineColorPrintSelected)
    {
      getColorFromString(lineColorPrintSelected, &style->lineColor.printSelected, &tmpError);
    }

  if (tmpError)
    {
      g_free(style);
      style = NULL;
      prefixError(tmpError, "Error creating style '%s'. ", styleName);
      g_propagate_error(error, tmpError);
    }
  
  return style;
}


/* Destroy a BlxStyle */
void destroyBlxStyle(BlxStyle *style)
{
  if (style)
    {
      g_free(style->styleName);
    }
}


/* Get the BlxStyle with the given name. Returns null and sets the error if 
 * we expected to find the name but didn't. The special-case input string "."
 * means "value not set", so in this case we return NULL and don't set the error. */
BlxStyle* getBlxStyle(const char *styleName, GSList *styles, GError **error)
{
  BlxStyle *result = NULL;

  if (styleName && strcmp(styleName, ".") != 0)
    {
      GSList *item = styles;
      
      for ( ; item; item = item->next)
	{
	  BlxStyle *style = (BlxStyle*)(item->data);
	  if (style && stringsEqual(style->styleName, styleName, FALSE))
	    {
	      result = style;
	      break;
	    }
	}
	
      if (!result)
	{
	  g_set_error(error, BLX_ERROR, 1, "Requested style '%s' does not exist.\n", styleName);
	}
    }

  return result;
}



/* Return a pointer to the BlxColor with the given id */
BlxColor* getBlxColor(GArray *defaultColors, const BlxColorId colorId)
{
  BlxColor *result = &g_array_index(defaultColors, BlxColor, colorId);

  if (result && result->transparent)
    {
      /* return the background color instead */
      result = &g_array_index(defaultColors, BlxColor, BLXCOLOR_BACKGROUND);
    }
  else if (!result)
    {
      printf("Warning: color ID %d not found in the default colors array.\n", colorId);
    }
  
  return result;
}


/* Lookup the BlxColor with the given id in the given hash table and return a 
 * pointer to the gdk color with the given properties */
GdkColor* getGdkColor(const BlxColorId colorId, GArray *defaultColors, const gboolean selected, const gboolean usePrintColors)
{
  GdkColor *result = NULL;
  
  BlxColor *blxColor = getBlxColor(defaultColors, colorId);
  
  if (blxColor)
    {
      if (usePrintColors)
	{
	  result = selected ? &blxColor->printSelected : &blxColor->print;
	}
      else
	{
	  result = selected ? &blxColor->selected : &blxColor->normal;
	}
    }
  else
    {
      printf("Error: requested invalid color ID %d", colorId);
    }
  
  return result;
}


/***********************************************************
 *                     Message handlers
 ***********************************************************/

/* Default handler for GLib log messages (e.g. from g_message() etc.) */
void defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  /* print to console */
  printf("%s", (char*)message);
  
  /* also display the message in the status bar (unless it's debug output) */
  if (log_level != G_LOG_LEVEL_DEBUG)
    {
      printMessageToStatusbar(message, data);
    }
}


/* Default handler for critical GLib log messages (i.e. from g_error and g_critical.) */
void popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  /* Display message to screen */
  printf("***%s", (char*)message);
  
  /* also display in the status bar */
  printMessageToStatusbar(message, data);

  /* Record each message in a list */
  GSList **messageList = getMessageList();
  
  BlxMessage *blxMessage = createBlxMessage(message, log_level);
  *messageList = g_slist_append(*messageList, blxMessage);
  
  /* Display as popup or in a scrolled list, whichever method is active */
  if (getUseScrolledMessages())
    {
      displayMessageAsList(*messageList, FALSE, data);
    }
  else
    {
      displayMessageAsPopup(message, log_level, data);
    }
}


/* Print the given message to the main window's statusbar. The main window is passed
 * as the data (which may be null if the window has not been created yet). */
static void printMessageToStatusbar(const gchar *message, gpointer data)
{
  if (data)
    {
      /* Remove any newline chars */
      char *displayText = g_strdup(message);
      char *cp = strchr(displayText, '\n');
      
      while (cp)
        {
          *cp = ' ';
          cp = strchr(displayText, '\n');
        }
    
      GtkWidget *blxWindow = GTK_WIDGET(data);
      BlxViewContext *bc = blxWindowGetContext(blxWindow);
      
      if (bc && bc->statusBar)
        {
          guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(bc->statusBar), "defaultMessages");
          gtk_statusbar_pop(GTK_STATUSBAR(bc->statusBar), contextId);
          gtk_statusbar_push(GTK_STATUSBAR(bc->statusBar), contextId, displayText);
        }
        
      g_free(displayText);
    }
}


/* Display a warning/error message as a popup */
static void displayMessageAsPopup(const gchar *message, GLogLevelFlags log_level, gpointer data)
{
  GtkWindow *parent = data ? GTK_WINDOW(data) : NULL;
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Error", 
						  parent, 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 12);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), vbox, TRUE, TRUE, 0);
  
  GtkWidget *hbox = gtk_hbox_new(FALSE, 12);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
  
  GtkWidget *image = gtk_image_new_from_stock(getDialogIcon(log_level), GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start(GTK_BOX(hbox), image, TRUE, TRUE, 0);
  
  GtkWidget *label = gtk_label_new(message);
  gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
  
  GtkWidget *button = gtk_check_button_new_with_mnemonic("Switch to _scrolled message window");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), getUseScrolledMessages());
  widgetSetCallbackData(button, onSetUseScrolledMessages, NULL);
  gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  
  gtk_widget_show_all(vbox);
  
  /* Block until the user clears the dialog */
  gtk_dialog_run(GTK_DIALOG(dialog));
}


/* Display a warning/error message in the scrolled message list */
static void displayMessageAsList(GSList *messageList, const gboolean bringToFront, gpointer data)
{
  static GtkWidget *dialog = NULL;
  static GtkWidget *button = NULL;
  static GtkTextView *textView = NULL;
  static GtkTextBuffer *textBuffer = NULL;
  static GtkTextTag *normalTag = NULL;       /* for normal text */
  static GtkTextTag *highlightTag = NULL;    /* for highlighting text */
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons("Errors", 
                                           NULL, 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);
      
      /* Hide the widget instead of destroying it when it is closed */
      g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
      
      GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
      gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), vbox, TRUE, TRUE, 0);
      
      /* Create the text view */
      int width = 600;
      int height = 180;
      
      GtkWidget *child = createScrollableTextView(NULL, FALSE, NULL, TRUE, &height, &textView);
      gtk_box_pack_start(GTK_BOX(vbox), child, TRUE, TRUE, 0);
      gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
      textBuffer = gtk_text_view_get_buffer(textView);
      
      /* Always show the horizontal scrollbar, otherwise it can cover the last line in the list
       * when it magically appears */
      gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(child), GTK_POLICY_ALWAYS, GTK_POLICY_AUTOMATIC);
      
      /* Add some tags for normal/highlighted text */
      normalTag = gtk_text_buffer_create_tag(textBuffer, NULL, "foreground", "#000000", NULL);
      highlightTag = gtk_text_buffer_create_tag(textBuffer, NULL, "foreground", "#ff0000", NULL);
      
      /* Create a button to allow user to switch back to popup messages */
      button = gtk_check_button_new_with_mnemonic("Switch to _popup messages");
      widgetSetCallbackData(button, onSetUsePopupMessages, NULL);
      gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
    }
  else
    {
      /* Delete contents of text buffer */
      GtkTextIter startIter;
      GtkTextIter endIter;
      gtk_text_buffer_get_start_iter(textBuffer, &startIter);
      gtk_text_buffer_get_end_iter(textBuffer, &endIter);
      
      gtk_text_buffer_delete(textBuffer, &startIter, &endIter);
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), !getUseScrolledMessages());
  
  addMessagesToBuffer(messageList, textBuffer, normalTag, highlightTag);
  
  /* Don't bring the dialog to the front if it was already open, unless specifically requested. */
  gtk_widget_show_all(dialog);
  
  if (bringToFront || !GTK_WIDGET_VISIBLE(dialog))
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
  
  /* Scroll to the last line in the text view. We must make sure the text view height has
   * been calculated first or gtk_text_view_scroll_to_iter might not do the correct thing. */
  while(gtk_events_pending()) 
    {
      gtk_main_iteration();
    }
  
  GtkTextIter textIter;
  gtk_text_buffer_get_end_iter(textBuffer, &textIter);
  gtk_text_view_scroll_to_iter(textView, &textIter, 0.0, FALSE, 0.0, 0.0);
}


/* Create a BlxMessage. Result should be destroyed with destroyBlxMessage. Timestamps the
 * message with the current system time. */
static BlxMessage* createBlxMessage(const char *text, const GLogLevelFlags logLevel)
{
  BlxMessage *blxMessage = g_malloc(sizeof(BlxMessage));
  
  blxMessage->text = g_strdup(text);
  blxMessage->timestamp = time(NULL);
  blxMessage->logLevel = logLevel;
  
  return blxMessage;
}

/* free the memory used by a blxMessage */
static void destroyBlxMessage(BlxMessage **blxMessage)
{
  if (blxMessage && *blxMessage)
    {
      if ((*blxMessage)->text)
	{
	  g_free((*blxMessage)->text);
	}
      
      g_free(*blxMessage);
    }
}


/* Get the message list. Creates the list the first time it is requested. */
static GSList** getMessageList()
{
  static GSList *messageList = NULL;
  
  return &messageList;
}


/* Clear the message list, freeing any memory that it uses */
void destroyMessageList()
{
  GSList **messageList = getMessageList();
  
  if (messageList && *messageList)
    {
      GSList *msgItem = *messageList;
      for ( ; msgItem; msgItem = msgItem->next)
	{
	  BlxMessage *blxMessage = (BlxMessage*)(msgItem->data);
	  destroyBlxMessage(&blxMessage);
	}
      
      g_slist_free(*messageList);
      *messageList = NULL;
    }
}


/* Get the display text for a message, including the timestamp if 'incTimestamp' is true.
 * The result should be free'd with g_free */
static char* blxMessageGetDisplayText(const BlxMessage *msg, const gboolean incTimestamp)
{
  /* Use ctime to format the timestamp but cut off the newline char at the end. */
  char *timeText = g_strdup(ctime(&msg->timestamp));
  char *cutPoint = strchr(timeText, '\n');
  if (cutPoint)
    {
      *cutPoint = '\0';
    }
  
  char separatorText[] = " - "; /* separates timestamp from message body */
  
  char *result = g_malloc(strlen(msg->text) + strlen(separatorText) + strlen(timeText) + 1);
  
  sprintf(result, "%s%s%s", timeText, separatorText, msg->text);
  
  return result;
}


/* Add each BlxMessage in the given list to the given text buffer */
static void addMessagesToBuffer(GSList *messageList, 
                                GtkTextBuffer *textBuffer, 
                                GtkTextTag *normalTag, 
                                GtkTextTag *highlightTag)
{
  GSList *msgItem = messageList;  
  
  for ( ; msgItem; msgItem = msgItem->next)
    {
      const BlxMessage *msg = (const BlxMessage*)(msgItem->data);
      
      /* If it's the last in the list, highlight it. */
      GtkTextTag *tag = (msgItem->next == NULL ? highlightTag : normalTag);
      
      char *displayText = blxMessageGetDisplayText(msg, TRUE);
      
      GtkTextIter endIter;
      gtk_text_buffer_get_end_iter(textBuffer, &endIter);
      gtk_text_buffer_insert_with_tags(textBuffer, &endIter, displayText, -1, tag, NULL);
      
      g_free(displayText);
    }
}


/* Internal function to get a pointer to the useScrolledMessages flag. Creates the flag the first 
 * time it is requested.  */
static gboolean* getUseScrolledMessagesPtr()
{
  static gboolean useScrolledMessages = FALSE;
  return &useScrolledMessages;
}


/* Access function to set the 'use scrolled messages' flag */
static gboolean getUseScrolledMessages()
{
  return *getUseScrolledMessagesPtr();
}

/* Access function to set the 'use scrolled messages' flag */
static void setUseScrolledMessages(const gboolean newValue)
{
  gboolean *useScrolledMessages = getUseScrolledMessagesPtr();
  *useScrolledMessages = newValue;
}


/* Called when user requests to view messages as a list */
static void onSetUseScrolledMessages(GtkWidget *button, const gint responseId, gpointer data)
{
  const gboolean useScrolledMessages = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  setUseScrolledMessages(useScrolledMessages);
}


/* Called when user requests to view messages as popups */
static void onSetUsePopupMessages(GtkWidget *button, const gint responseId, gpointer data)
{
  const gboolean usePopupMessages = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  setUseScrolledMessages(!usePopupMessages);
}


/* Utility to return the stock icon to use for a dialog for a given message severity */
static char* getDialogIcon(GLogLevelFlags log_level)
{
  char *result = NULL;
  
  switch (log_level)
    {
      case G_LOG_LEVEL_ERROR:
        result = GTK_STOCK_DIALOG_ERROR;
        break;
        
      case G_LOG_LEVEL_CRITICAL:
      case G_LOG_LEVEL_WARNING:
        result = GTK_STOCK_DIALOG_WARNING;
        break;
        
      default:
        result = GTK_STOCK_DIALOG_INFO;
    }
    
  return result;
}

/***************** end of file ***********************/
