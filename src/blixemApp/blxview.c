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
 * Exported functions: See blxview.h
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
 * CVS info:   $Id: blxview.c,v 1.87 2010-11-08 18:41:29 gb10 Exp $
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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <gdk/gdkkeysyms.h>

#include <blixemApp/blixem_.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blxdotter.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/utilities.h>


#define MAXALIGNLEN                   10000
#define ORGANISM_PREFIX_SEPARATOR     "  "



static void            blviewCreate(char *align_types, const char *paddingSeq, GList* featureLists[], GList *seqList, GSList *supportedTypes, CommandLineOptions *options, const char *net_id, int port, const gboolean External) ;
static void            finaliseBlxSequences(GList *featureLists[], MSP **mspList, GList **seqList, const int offset);
static void            processGeneName(BlxSequence *blxSeq);
static void            processOrganism(BlxSequence *blxSeq);


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


/* Check we have everything we need */
static void validateInput(CommandLineOptions *options)
{
  if (!options->refSeq)
    {
      g_error("Error: reference sequence not specified.\n");
    }
  
  /* if we were given the blast mode instead of seq type, determine the display sequence type */
  if (options->seqType == BLXSEQ_INVALID && options->blastMode != BLXMODE_UNSET)
    {
      if (options->blastMode == BLXMODE_BLASTN)
        options->seqType = BLXSEQ_DNA;
      else
        options->seqType = BLXSEQ_PEPTIDE;
    }
  
  /* Check we have a valid seq type */
  if (options->seqType == BLXSEQ_INVALID)
    {
      printf("\nNo sequence type specified. Detected ");
      
      if (determineSeqType(options->refSeq) == BLXSEQ_PEPTIDE)
	{
	  printf("protein sequence. Will try to run Blixem in protein mode.\n");
	  options->seqType = BLXSEQ_PEPTIDE;
	}
      else
	{
	  printf("nucleotide sequence. Will try to run Blixem in nucelotide mode.\n");
	  options->seqType = BLXSEQ_DNA;
	}
    }
  
  /* Ideally we'd get rid of blast mode and just deal with real sequence types (ref seq type, match
   * seq type and display type), but it's too in-built for now so make sure it's set */
  if (options->blastMode == BLXMODE_UNSET)
    options->blastMode = (options->seqType == BLXSEQ_DNA ? BLXMODE_BLASTN : BLXMODE_BLASTX);
  
  /* Set the number of reading frames */
  if (options->seqType == BLXSEQ_PEPTIDE)
    options->numFrames = 3;
  else
    options->numFrames = 1;
  
  /* Now we've got the ref seq, we can work out the zoom, if zooming to view the whole sequence */
  if (options->zoomWhole)
    {
      options->bigPictZoom = strlen(options->refSeq);
    }
  else
    {
      if (options->seqType == BLXSEQ_DNA)
        options->bigPictZoom = 30;
      else
        options->bigPictZoom = 10;
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
      processOrganism(blxSeq);
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
 *    pfetch:  if non-NULL, then we use pfetch instead of efetch for
 *             _all_ sequence fetching (use the node/port info. in the
 *             pfetch struct to locate the pfetch server).
 *
 */
gboolean blxview(CommandLineOptions *options,
                 GList* featureLists[],
                 GList *seqList,
                 GSList *supportedTypes,
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

  validateInput(options);
  
  const char *net_id = NULL;
  int port = UNSET_INT;
  gboolean status = blxviewFetchSequences(pfetch, External, options->parseFullEmblInfo, options->seqType, seqList, &options->fetchMode, &net_id, &port);
  
  if (status)
    {
      /* Construct missing data and do any other required processing now we have all the sequence data */
      finaliseBlxSequences(featureLists, &options->mspList, &seqList, options->refSeqOffset);
    }

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(align_types, padseq, featureLists, seqList, supportedTypes, options, net_id, port, External) ;
    }
  
  return status;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *align_types, 
			 const char *paddingSeq,
                         GList* featureLists[],
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
      blixemWindow = createBlxWindow(options, paddingSeq, featureLists, seqList, supportedTypes, net_id, port, External);

      gboolean pep_nuc_align = (options->blastMode == BLXMODE_BLASTX || options->blastMode == BLXMODE_BLASTN);
      
      char *title = blxprintf("Blixem %s%s%s:   %s",
                              (pep_nuc_align ? "  (" : ""),
                              (align_types ? align_types : (options->seqType == BLXSEQ_PEPTIDE ? "peptide" :
                                                            (options->seqType == BLXSEQ_DNA ? "nucleotide" : ""))),
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
      char *dotterSSeq = getDotterSSeq(blixemWindow, &error);
      
      if (!error)
        {
          callDotter(blixemWindow, FALSE, dotterSSeq, &error);
        }
        
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  if (options->startNextMatch)
    {
      /* Set the start coord to be the start of the next MSP on from the default start coord */
      nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
    }
}

/***********************************************************
 ***********************************************************/

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
                                    BlxSequence *blxSeq, GList* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, 
                                    GError **error)
{
  BlxMspType newType = BLXMSP_INVALID;
  MSP **ptrToUpdate = NULL;
  int newStart = UNSET_INT;
  int newEnd = UNSET_INT;
  int newPhase = 0;
  BlxStyle *newStyle = NULL;

  if (!*exon && (*cds || *utr))
    {
      /* The exon spans both the cds and utr */
      if (*cds && *utr)
        {
          newStart = min((*cds)->qRange.min, (*utr)->qRange.min);
          newEnd = max((*cds)->qRange.max, (*utr)->qRange.max);
          newPhase = (*cds)->phase;
          newStyle = (*cds)->style;
          newType = BLXMSP_EXON;
          ptrToUpdate = exon;
        }
      else if (*cds)
        {
          newStart = (*cds)->qRange.min;
          newEnd = (*cds)->qRange.max;
          newPhase = (*cds)->phase;
          newStyle = (*cds)->style;
          newType = BLXMSP_EXON;
          ptrToUpdate = exon;
        }
      else
        {
          newStart = (*utr)->qRange.min;
          newEnd = (*utr)->qRange.max;
          newPhase = (*utr)->phase;
          newStyle = (*utr)->style;
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
          newStyle = (*exon)->style;
          newType = BLXMSP_CDS;
          ptrToUpdate = cds;
        }
      else if ((*exon)->qRange.min < (*utr)->qRange.min)
        {
          newStart = (*exon)->qRange.min;
          newEnd = (*utr)->qRange.min - 1;
          newPhase = (*exon)->phase;
          newStyle = (*exon)->style;
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
          newStyle = (*cds)->style;
          newType = BLXMSP_UTR;
          ptrToUpdate = utr;
        }
      else if ((*cds)->qRange.max < (*exon)->qRange.max)
        {
          newStart = (*cds)->qRange.max + 1;
          newEnd = (*exon)->qRange.max;
          newPhase = (*cds)->phase;
          newStyle = (*cds)->style;
          newType = BLXMSP_UTR;
          ptrToUpdate = utr;
        }
    }
  else if (*exon)
    {
      /* We just have an exon. Assume it is all utr. */
      newStart = (*exon)->qRange.min;
      newEnd = (*exon)->qRange.max;
      newPhase = 0;
      newStyle = (*exon)->style;
      newType = BLXMSP_UTR;
      ptrToUpdate = utr;
    }
  
  if (newType != BLXMSP_INVALID)
    {
      /* Create the new exon/cds/utr */
      DEBUG_OUT("Creating MSP for transcript '%s' of type %d.\n", blxSeq->fullName, newType);
      
      GError *tmpError = NULL;
      
      *ptrToUpdate = createNewMsp(featureLists, lastMsp, mspList, seqList, newType, NULL, UNSET_INT, UNSET_INT, newPhase, NULL, blxSeq->idTag,
                                  NULL, newStart, newEnd, blxSeq->strand, UNSET_INT, blxSeq->fullName,
                                  UNSET_INT, UNSET_INT, blxSeq->strand, NULL, &tmpError);
      
      (*ptrToUpdate)->style = newStyle;
      
      if (tmpError)
        {
          prefixError(tmpError, "Error constructing missing exon/cds/utr [type='%d']", newType);
          g_propagate_error(error, tmpError);
        }
    }
    
  /* We should now have all the bits. Set the relationship data. */
  if (*exon && (*cds || *utr))
    {
      if (*cds)
        (*exon)->childMsps = g_list_append((*exon)->childMsps, *cds);
        
      if (*utr)
        (*exon)->childMsps = g_list_append((*exon)->childMsps, *utr);
    }
}


/* Construct any missing transcript data, i.e.
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructTranscriptData(BlxSequence *blxSeq, GList* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList)
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
              ((prevMsp && mspIsExon(prevMsp) && !rangesOverlap(&prevMsp->qRange, &msp->qRange)) ||
               (!prevMsp && blxSeq->qRange.min < msp->qRange.min)))
            {
              foundGap = TRUE;
            }

          if (foundGap || msp == NULL)
            {
              /* We've found a gap between exons, or reached the end. First, see if the current exon/cds or utr
               * is missing and construct it if possible. Also do this if we're at the last MSP. */
              createMissingExonCdsUtr(&curExon, &curCds, &curUtr, blxSeq, featureLists, lastMsp, mspList, seqList, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
              copyCdsReadingFrame(curExon, curCds, curUtr);

              IntRange newRange = {UNSET_INT, UNSET_INT};
              
              if (prevExon && curExon && !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron to span the gap */
                  newRange.min = prevExon->qRange.max + 1;
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (!prevExon && curExon && blxSeq->qRange.min < curExon->qRange.min && !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the start */
                  newRange.min = blxSeq->qRange.min;
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (msp == NULL && curExon && blxSeq->qRange.max > curExon->qRange.max && !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the end */
                  newRange.min = curExon->qRange.max + 1;
                  newRange.max = blxSeq->qRange.max;
                }
              
              if (newRange.min != UNSET_INT && newRange.max != UNSET_INT)
                {
                  createNewMsp(featureLists, lastMsp, mspList, seqList, BLXMSP_INTRON, NULL, UNSET_INT, UNSET_INT, UNSET_INT, NULL, blxSeq->idTag,
                              NULL, newRange.min, newRange.max, blxSeq->strand, UNSET_INT, blxSeq->fullName,
                               UNSET_INT, UNSET_INT, blxSeq->strand, NULL, &tmpError);
                  
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


/* Add a prefix to the given BlxSequence's organism field, e.g. change "Homo sapiens (human)" 
 * to "Hs - Homo sapiens (human)". This is so that we can have a very narrow column that just 
 * displays the prefix part of the field. */
static void processOrganism(BlxSequence *blxSeq)
{
  if (blxSeq && blxSeq->organism && blxSeq->organism->str)
    {
      /* To create the alias, we'll take the first letter of each word until we hit a non-alpha
       * character, e.g. the alias for "Homo sapiens (human)" would be "Hs" */
      int srcIdx = 0;
      int srcLen = strlen(blxSeq->organism->str);
      
      int destIdx = 0;
      int destLen = 5; /* limit the length of the alias */
      char alias[destLen + 1];
      
      gboolean startWord = TRUE;
      
      for ( ; srcIdx < srcLen && destIdx < destLen; ++srcIdx)
        {
          if (blxSeq->organism->str[srcIdx] == ' ')
            {
              startWord = TRUE;
            }
          else if (g_ascii_isalpha(blxSeq->organism->str[srcIdx]))
            {
              if (startWord)
                {
                  alias[destIdx] = blxSeq->organism->str[srcIdx];
                  ++destIdx;
                  startWord = FALSE;
                }
            }
          else
            {
              /* Finish if we've hit a char that's not an alpha char or space. */
              break;
            }
        }

      alias[destIdx] = '\0';
      
      if (destIdx > 0)
        {
          g_string_prepend(blxSeq->organism, ORGANISM_PREFIX_SEPARATOR);
          g_string_prepend(blxSeq->organism, alias);
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
static void finaliseBlxSequences(GList* featureLists[], MSP **mspList, GList **seqList, const int offset)
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
      if (blxSeq && blxSeq->type == BLXSEQUENCE_MATCH && blxSeq->strand == BLXSTRAND_REVERSE && blxSeq->sequence && blxSeq->sequence->str)
        {
          blxComplement(blxSeq->sequence->str);
        }
      
      findSequenceExtents(blxSeq);
      constructTranscriptData(blxSeq, featureLists, &lastMsp, mspList, seqList);
    }
}


/* Returns true if this msp is part of a feature series */
gboolean mspHasFs(const MSP *msp)
{
  gboolean result = (msp->type == BLXMSP_FS_SEG || msp->type == BLXMSP_XY_PLOT);
  return result;
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

  if (msp->url)
    {
      g_free(msp->url);
      msp->url = NULL;
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


/* Get the full range of the given MSP that we want to display, in s coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
void mspGetFullSRange(const MSP const *msp, 
                      const gboolean seqSelected,
                      const gboolean *flags,
                      const int numUnalignedBases, 
                      const GList const *polyASiteList,
                      IntRange *result)
{
  /* Normally we just display the part of the sequence in the alignment */
  result->min = msp->sRange.min;
  result->max = msp->sRange.max;
  
  if (flags[BLXFLAG_SHOW_UNALIGNED] && mspGetMatchSeq(msp) && (!flags[BLXFLAG_SHOW_UNALIGNED_SELECTED] || seqSelected))
    {
      /* We're displaying additional unaligned sequence outside the alignment range. Get 
       * the full range of the match sequence */
      result->min = 1;
      result->max = mspGetMatchSeqLen(msp);
      
      if (flags[BLXFLAG_LIMIT_UNALIGNED_BASES])
        {
          /* Only include up to 'numUnalignedBases' each side of the MSP range (still limited
           * to the range we found above, though). */
          result->min = max(result->min, msp->sRange.min - numUnalignedBases);
          result->max = min(result->max, msp->sRange.max + numUnalignedBases);
        }
    }
  
  if (flags[BLXFLAG_SHOW_POLYA_SITE] && mspHasPolyATail(msp, polyASiteList) && (!flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED] || seqSelected))
    {
      /* We're displaying polyA tails, so override the 3' end coord with the full extent of
       * the s sequence if there is a polyA site here. The 3' end is the min q coord if the
       * match is on forward ref seq strand or the max coord if on the reverse. */
      const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
      
      if (sameDirection)
        {
          result->max = mspGetMatchSeqLen(msp);
        }
      else
        {
          result->min = 1;
        }
    }
}


/* Get the full range of the given MSP that we want to display, in q coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
void mspGetFullQRange(const MSP const *msp, 
                      const gboolean seqSelected,
                      const gboolean *flags,
                      const int numUnalignedBases, 
                      const GList const *polyASiteList,
                      const int numFrames,
                      IntRange *result)
{
  /* Default to the alignment range so we can exit quickly if there are no special cases */
  result->min = msp->qRange.min;
  result->max = msp->qRange.max;
  
  if (flags[BLXFLAG_SHOW_UNALIGNED] || flags[BLXFLAG_SHOW_POLYA_SITE])
    {
      /* Get the full display range of the MSP including any unaligned portions of the match sequence. */
      IntRange fullSRange;
      mspGetFullSRange(msp, seqSelected, flags, numUnalignedBases, polyASiteList, &fullSRange);
      
      /* Find the offset of the start and end of the full range compared to the alignment range and
       * offset the ref seq range by the same amount. We need to multiply by the number of reading
       * frames because q coords are in nucleotides and s coords are in peptides. */
      const int startOffset = (msp->sRange.min - fullSRange.min) * numFrames;
      const int endOffset = (fullSRange.max - msp->sRange.max) * numFrames;
      
      const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
      result->min -= sameDirection ? startOffset : endOffset;
      result->max += sameDirection ? endOffset : startOffset;
    }
}


/* Given a base index on the query sequence, find the corresonding base 
 * in subject sequence. The return value is always UNSET_INT if there is not
 * a corresponding base at this position. However, in this case the start/end
 * of the sequence (or nearest gap) is returned in the nearestIdx argument. If
 * the base exists then the return idx and nearestIdx will be the same. 
 * nearestIdx can be null if not required. */
int gapCoord(const MSP *msp, 
	     const int qIdx, 
	     const int numFrames, 
	     const BlxStrand strand, 
	     const gboolean displayRev,
             const gboolean seqSelected,
             const int numUnalignedBases,
             gboolean *flags,
             const GList const *polyASiteList)
{
  int result = UNSET_INT;
  
  if (mspIsBlastMatch(msp) || mspIsExon(msp))
    {
      const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sameDirection = (qForward == sForward);
      const gboolean inGapsRange = (qIdx >= msp->qRange.min && qIdx <= msp->qRange.max);
      
      if (msp->gaps && g_slist_length(msp->gaps) >= 1 && inGapsRange)
        {
          /* Gapped alignment. Look to see if x lies inside one of the "gaps" ranges. */
          GSList *rangeItem = msp->gaps;
          
          for ( ; rangeItem ; rangeItem = rangeItem->next)
            {
              CoordRange *curRange = (CoordRange*)(rangeItem->data);
              
              int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
              getCoordRangeExtents(curRange, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
              
              /* We've "found" the value if it's in or before this range. Note that the
               * the range values are in decreasing order if the q strand is reversed. */
              gboolean found = qForward ? qIdx <= qRangeMax : qIdx >= qRangeMin;
              
              if (found)
                {
                  gboolean inRange = (qForward ? qIdx >= qRangeMin : qIdx <= qRangeMax);
                  
                  if (inRange)
                    {
                      /* It's inside this range. Calculate the actual index. */
                      int offset = (qIdx - qRangeMin) / numFrames;
                      result = sameDirection ? sRangeMin + offset : sRangeMax - offset;
                    }
                  
                  break;
                }
            }
        }
      else if (!inGapsRange && mspIsBlastMatch(msp))
        {
          /* The q index is outside the alignment range but the option to show unaligned sequence 
           * is enabled, so check if the q index is within the full display range including 
           * any unaligned portions of the sequence that we're displaying. */
          
          if (qIdx < msp->qRange.min)
            {
              /* Find the offset backwards from the low coord and subtract it from the low end
               * of the s coord range (or add it to the high end, if the directions are opposite). */
              const int offset = (msp->qRange.min - qIdx) / numFrames;
              result = sameDirection ? msp->sRange.min - offset : msp->sRange.max + offset;
            }
          else
            {
              /* Find the offset forwards from the high coord and add it to the high end of the
               * s coord range (or subtract it from the low end, if the directions are opposite). */
              const int offset = (qIdx - msp->qRange.max) / numFrames;
              result = sameDirection ? msp->sRange.max + offset : msp->sRange.min - offset;
            }
          
          /* Get the full display range of the match sequence. If the result is still out of range
           * then there's nothing to show for this qIdx. */
          IntRange fullSRange;
          mspGetFullSRange(msp, seqSelected, flags, numUnalignedBases, polyASiteList, &fullSRange);
          
          if (!valueWithinRange(result, &fullSRange))
            {
              result = UNSET_INT;
            }
        }
      else
        {
          /* If strands are in the same direction, find the offset from qRange.min and add it to 
           * sRange.min. If strands are in opposite directions, find the offset from qRange.min and
           * subtract it from sRange.max. Note that the offset could be negative if we're outside
           * the alignment range. */
          int offset = (qIdx - msp->qRange.min)/numFrames ;
          result = (sameDirection) ? msp->sRange.min + offset : msp->sRange.max - offset ;
          
          if (result < msp->sRange.min || result > msp->sRange.max)
            {
              result = UNSET_INT;
            }
          
        }
    }
  
  return result;
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
      
      /* Check if sequence data was requested and is not already set. */
      gboolean getSeq = (blxSequenceRequiresSeqData(blxSeq) && getSequenceData && blxSeq->sequence == NULL);

      /* Check if optional data was requested and is not already set. We can assume that
       * if any of the data fields is set then the parsing has been done for all of them
       * (and any remaining empty fields just don't have that data available) */
      getSeq |= (blxSequenceRequiresOptionalData(blxSeq) &&
                 getOptionalData && 
                     !blxSeq->organism &&
                     !blxSeq->geneName &&
                     !blxSeq->tissueType &&
                     !blxSeq->strain);
          
      if (getSeq)
        {
          resultList = g_list_prepend(resultList, blxSeq);
        }
    }

  return resultList;
}


/* Set the given BlxColor elements from the given color string(s). Works out some good defaults
 * for blxColor.selected blxcolor.print etc if these colors are not given */
static void setBlxColorValues(const char *normal, const char *selected, BlxColor *blxColor, GError **error)
{
  GError *tmpError = NULL;
  
  getColorFromString(normal, &blxColor->normal, &tmpError);
  getSelectionColor(&blxColor->normal, &blxColor->selected); /* will use this if selection color is not given/has error */

  /* Calculate print coords as a greyscale version of normal colors */
  convertToGrayscale(&blxColor->normal, &blxColor->print);            /* calculate print colors, because they are not given */
  getSelectionColor(&blxColor->print, &blxColor->printSelected);
  
  /* Use the selected-color string, if given */
  if (!tmpError && selected)
    {
      getColorFromString(selected, &blxColor->selected, &tmpError);
    }
  
  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Create a BlxStyle. For transcripts, CDS and UTR features can have different colors, but they
 * are from the same source and hence have the same style object. We therefore use the default
 * fillColor, lineColor etc. variables for CDS features and provide additional options to supply 
 * UTR colors as well. */
BlxStyle* createBlxStyle(const char *styleName, 
			 const char *fillColor, 
			 const char *fillColorSelected, 
			 const char *lineColor, 
			 const char *lineColorSelected, 
			 const char *fillColorUtr, 
			 const char *fillColorUtrSelected, 
			 const char *lineColorUtr, 
			 const char *lineColorUtrSelected, 
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
      setBlxColorValues(fillColor, fillColorSelected, &style->fillColor, &tmpError);
      setBlxColorValues(fillColor, fillColorSelected, &style->fillColorUtr, &tmpError); /* default UTR to same as CDS */
    }
    
  if (!tmpError && lineColor)
    {
      setBlxColorValues(lineColor, lineColorSelected, &style->lineColor, &tmpError);
      setBlxColorValues(lineColor, lineColorSelected, &style->lineColorUtr, &tmpError); /* default UTR to same as CDS */
    }
    
  if (!tmpError && fillColorUtr)
    {
      setBlxColorValues(fillColorUtr, fillColorUtrSelected, &style->fillColorUtr, &tmpError);
    }
  
  if (!tmpError && lineColorUtr)
    {
      setBlxColorValues(lineColorUtr, lineColorUtrSelected, &style->lineColorUtr, &tmpError);
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

/***************** end of file ***********************/
