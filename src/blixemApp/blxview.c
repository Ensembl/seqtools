/*  File: blxview.c
 *  Author: Erik Sonnhammer, 1993-02-20
 *  Copyright (c) 2009 - 2010 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
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
 * Description: Setup and miscellaneous functions for Blixem.
 *               
 *              This file is a bit of a mish-mash. It originally contained all
 *              of the graphical stuff for Blixem, but that has been moved to
 *              other files. This file now mostly contains setup functions,
 *              but it also contains other functions that don't really belong
 *              anywhere else.
 *----------------------------------------------------------------------------
 */


/*
MSP score codes (for obsolete exblx file format):
----------------
-1 exon                  -> Big picture + Alignment
-2 intron                -> Big picture + Alignment
-3 Any colored segment   -> Big picture
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



static void            blviewCreate(char *align_types, const char *paddingSeq, GArray* featureLists[], GList *seqList, GSList *supportedTypes, CommandLineOptions *options, const char *net_id, int port, const gboolean External) ;
static void            processGeneName(BlxSequence *blxSeq);
static void            processOrganism(BlxSequence *blxSeq);
static GList*          getSeqsToPopulate(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData);
static GHashTable*     getSeqsToPopulateByMode(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData, const char *defaultFetchMode);


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
static gboolean fetchSequences(GList *seqsToFetch, 
                               GList *seqList,
                               const gchar *fetchMode, 
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
      else if (strcmp(fetchMode, BLX_FETCH_DB) == 0)
        {
          g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", fetchMode);
        }
      else if (strcmp(fetchMode, BLX_FETCH_REGION) == 0)
        {
          g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", fetchMode);
        }
      else if (fetchMode)
        {
          g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", fetchMode);
        }
      else
        {
          g_set_error(error, BLX_ERROR, 1, "Cannot fetch sequences - fetch mode not specified.\n");
        }
    }
    
  
  if (parseOptionalData)
    {
      /* Once we've fetched all the sequences we need to do some post-processing. Loop 
       * twice: once to modify any fields in our own custom manner, and once more to
       * see if any with missing data can copy it from their parent. (Need to do these in
       * separate loops or we don't know if the data we're copying is processed or not.) */
      GList *seqItem = seqList;
      for ( ; seqItem; seqItem = seqItem->next)
        {
          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
          processGeneName(blxSeq);
          processOrganism(blxSeq);
        }

      for (seqItem = seqList; seqItem; seqItem = seqItem->next)
        {
          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
          populateMissingDataFromParent(blxSeq, seqList);
        }
    }
  
  return success;
}


/* Find out if we need to fetch any sequences (they may all be contained in the input
 * files), if we do need to, then fetch them by the appropriate method. */
gboolean blxviewFetchSequences(gboolean External, 
                               const gboolean parseFullEmblInfo,
                               const gboolean parseSequenceData,
                               const BlxSeqType seqType,
                               GList *seqList, /* list of BlxSequence structs for all required sequences */
                               char *bulkFetchMode,
                               const char *net_id,
                               const int port)
{
  gboolean success = FALSE; /* will get set to true if any of the fetch methods succeed */
  
  /* Fetch any sequences that do not have their sequence data already populated (or
   * optional data too, if requested). */
  GHashTable *seqsTable = getSeqsToPopulateByMode(seqList, TRUE, parseFullEmblInfo, bulkFetchMode);

  /* Loop through each fetch mode */
  GHashTableIter iter;
  g_hash_table_iter_init(&iter, seqsTable);
  gpointer key, value;
  
  while (g_hash_table_iter_next(&iter, &key, &value))
    {
      GQuark dataTypeQuark = GPOINTER_TO_INT(key);
      const gchar *fetchMode = g_quark_to_string(dataTypeQuark);
      GList *seqsToFetch = (GList*)value;
      
      GError *error = NULL;
      
      if (fetchSequences(seqsToFetch, seqList, fetchMode, seqType, net_id, port, parseFullEmblInfo, parseSequenceData, External, &error))
        {
          success = TRUE;
        }
      else
        {
          prefixError(error, "Error fetching sequences by '%s' mode. ");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
      
      /* We're done with this list now, so free the memory. Don't delete it from
       * the table yet, though, because that will invalidate the iterators. */
      g_list_free(seqsToFetch);
    }
  
  /* We're done with the hash table now */
  g_hash_table_unref(seqsTable);
  
  if (!success)
    {
      g_critical("Error fetching sequences.\n");
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
                 GArray* featureLists[],
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
  setupFetchModes(pfetch, &options->bulkFetchMode, &options->userFetchMode, &net_id, &port);
  
  gboolean status = blxviewFetchSequences(External, options->parseFullEmblInfo, TRUE, options->seqType, seqList, options->bulkFetchMode, net_id, port);
  
  if (status)
    {
      /* Construct missing data and do any other required processing now we have all the sequence data */
      finaliseBlxSequences(featureLists, &options->mspList, &seqList, options->refSeqOffset, options->seqType, options->numFrames, &options->refSeqRange, TRUE);

      /* Sort msp arrays by start coord (only applicable to msp types that
       * appear in the detail-view because the order is only applicable when
       * filtering detail-view rows) */
      int typeId = 0;
      for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
        {
          if (typeShownInDetailView(typeId))
            g_array_sort(featureLists[typeId], compareFuncMspArray);
        }
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
                         GArray* featureLists[],
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
      (mspIsBlastMatch(options->mspList) || options->mspList->type == BLXMSP_HSP || options->mspList->type == BLXMSP_GSP))
    {
      /* We must select the sequence before calling callDotter. Get the first
       * sequence to the right of the start coord. */
      MSP *msp = nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
      
      if (msp)
        {
          blxWindowSelectSeq(blixemWindow, msp->sSequence);
          
          GError *error = NULL;
          char *dotterSSeq = getDotterSSeq(blixemWindow, &error);
          
          if (!error)
            {
              callDotter(blixemWindow, FALSE, dotterSSeq, NULL);
            }
            
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
      else 
        {
          g_error("Could not run Dotter on first sequence; no sequences found to the right of the start coord.\n");
        }
    }

  if (options->startNextMatch)
    {
      /* Set the start coord to be the start of the next MSP on from the default start coord */
      nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
    }
}

/***********************************************************
 ***********************************************************/

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


/* Returns true if this msp is part of a feature series */
gboolean mspHasFs(const MSP *msp)
{
  gboolean result = (msp->type == BLXMSP_FS_SEG || msp->type == BLXMSP_XY_PLOT);
  return result;
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


/* Return the (cached) full extent of the match that we're showing in match seq coords */
const IntRange* mspGetFullSRange(const MSP const *msp, const gboolean seqSelected, const BlxViewContext const *bc)
{
  const IntRange *result = NULL;
  
  if (seqSelected || !bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED] || !bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED])
    result = &msp->fullSRange;
  else
    result = &msp->sRange;
  
  return result;
}

/* Return the (cached) full extent of the match on the ref sequence in display coords
 * (including any portions of unaligned sequence that we're showing). Depending on the
 * options, the range may depend on whether the sequence is selected or not. */
const IntRange* mspGetFullDisplayRange(const MSP const *msp, const gboolean seqSelected, const BlxViewContext const *bc)
{
  const IntRange *result = NULL;
  
  if (seqSelected || !bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED] || !bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED])
    result = &msp->fullRange;
  else
    result = &msp->displayRange;
  
  return result;
}

/* Return the (cached) extent of the alignment that we're showing in display coords
 * (excluding unaligned sequence) */
const IntRange* mspGetDisplayRange(const MSP const *msp)
{
  return &msp->displayRange;
}


/* Get the full range of the given MSP that we want to display, in s coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
static void mspCalcFullSRange(const MSP const *msp, 
			      const gboolean *flags,
			      const int numUnalignedBases, 
			      const GArray const *polyASiteList,
			      IntRange *result)
{
  /* Normally we just display the part of the sequence in the alignment */
  result->min = msp->sRange.min;
  result->max = msp->sRange.max;
  
  if (mspIsBlastMatch(msp))
    {
      if (flags[BLXFLAG_SHOW_UNALIGNED] && mspGetMatchSeq(msp))
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
      
      if (flags[BLXFLAG_SHOW_POLYA_SITE] && mspHasPolyATail(msp, polyASiteList))
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
}


/* Get the full range of the given MSP that we want to display, in q coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
static void mspCalcFullQRange(const MSP const *msp, 
			      const gboolean *flags,
			      const int numUnalignedBases, 
			      const GArray const *polyASiteList,
			      const int numFrames,
			      const IntRange const *fullSRange,
			      IntRange *result)
{
  /* Default to the alignment range so we can exit quickly if there are no special cases */
  result->min = msp->qRange.min;
  result->max = msp->qRange.max;
  
  if (mspIsBlastMatch(msp) && (flags[BLXFLAG_SHOW_UNALIGNED] || flags[BLXFLAG_SHOW_POLYA_SITE]))
    {
      /* Find the offset of the start and end of the full range compared to the alignment range and
       * offset the ref seq range by the same amount. We need to multiply by the number of reading
       * frames because q coords are in nucleotides and s coords are in peptides. */
      const int startOffset = (msp->sRange.min - fullSRange->min) * numFrames;
      const int endOffset = (fullSRange->max - msp->sRange.max) * numFrames;
      
      const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
      result->min -= sameDirection ? startOffset : endOffset;
      result->max += sameDirection ? endOffset : startOffset;
    }
}


/* Calculate the full extent of the match sequence to display, in display coords,
 * and cache the result in the msp. Includes any portions of unaligned sequence that we're
 * displaying */
void mspCalculateFullExtents(MSP *msp, const BlxViewContext const *bc, const int numUnalignedBases)
{
  mspCalcFullSRange(msp, bc->flags, numUnalignedBases, bc->featureLists[BLXMSP_POLYA_SITE], &msp->fullSRange);
  mspCalcFullQRange(msp, bc->flags, numUnalignedBases, bc->featureLists[BLXMSP_POLYA_SITE], bc->numFrames, &msp->fullSRange, &msp->fullRange);
 
  /* convert the Q range to display coords */
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int coord1 = convertDnaIdxToDisplayIdx(msp->fullRange.min, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->fullRange.max, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  intrangeSetValues(&msp->fullRange, coord1, coord2);
  
  /* Remember the max len of all the MSPs in the detail-view */
  if (typeShownInDetailView(msp->type) && getRangeLength(&msp->fullRange) > getMaxMspLen())
    {
      setMaxMspLen(getRangeLength(&msp->fullRange));
    }
}


/* Convert the ref-seq range of the given msp in display coords and cache it in the msp */
static void mspCalculateDisplayRange(MSP *msp, const BlxViewContext const *bc)
{
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  intrangeSetValues(&msp->displayRange, coord1, coord2);  
}


/* This caches the display range (in display coords rather than dna coords,
 * and inverted if the display is inverted) for each MSP */
void cacheMspDisplayRanges(const BlxViewContext const *bc, const int numUnalignedBases)
{
  /* This also calculates the max msp len */
  setMaxMspLen(0);
  
  MSP *msp = bc->mspList;
  for ( ; msp; msp = msp->next)
    {
      mspCalculateDisplayRange(msp, bc);
      mspCalculateFullExtents(msp, bc, numUnalignedBases);
    }
}


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the MSP is a gapped MSP and the ref-seq coord is known to lie within the
 * MSP's alignment range. */
static int mspGetGappedAlignmentCoord(const MSP *msp, const int qIdx, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));

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
              int offset = (qIdx - qRangeMin) / bc->numFrames;
              result = sameDirection ? sRangeMin + offset : sRangeMax - offset;
            }
          
          break;
        }
    }
  
  return result;
}


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the MSP is an ungapped MSP and the ref-seq coord is known to lie within the
 * MSP's alignment range. */
static int mspGetUngappedAlignmentCoord(const MSP *msp, const int qIdx, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  /* If strands are in the same direction, find the offset from qRange.min and add it to 
   * sRange.min. If strands are in opposite directions, find the offset from qRange.min and
   * subtract it from sRange.max. Note that the offset could be negative if we're outside
   * the alignment range. */
  int offset = (qIdx - msp->qRange.min) / bc->numFrames ;
  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));

  result = (sameDirection) ? msp->sRange.min + offset : msp->sRange.max - offset ;
  
  if (result < msp->sRange.min || result > msp->sRange.max)
    {
      result = UNSET_INT;
    }
  
  return result;
}  


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the ref-seq coord is known to lie outside the MSP's alignment range. The 
 * result will be unset unless the option to display unaligned portions of 
 * sequence is enabled. */
static int mspGetUnalignedCoord(const MSP *msp, const int qIdx, const gboolean seqSelected, const int numUnalignedBases, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  /* First convert to display coords */
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int displayIdx = convertDnaIdxToDisplayIdx(qIdx, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const IntRange const *mspRange = mspGetDisplayRange(msp);

  /* Note that because we have converted to display coords the ref seq coords are
   * now in the direction of the display, regardless of the ref seq strand. */
  const gboolean sameDirection = (mspGetMatchStrand(msp) == mspGetRefStrand(msp));
  const gboolean sForward = (sameDirection != bc->displayRev);
  
  if (displayIdx < mspRange->min)
    {
      /* Find the offset backwards from the low coord and subtract it from the low end
       * of the s coord range (or add it to the high end, if the directions are opposite).
       * We're working in display coords here. */
      const int offset = mspRange->min - displayIdx;
      result = sForward ? msp->sRange.min - offset : msp->sRange.max + offset;
    }
  else
    {
      /* Find the offset forwards from the high coord and add it to the high end of the
       * s coord range (or subtract it from the low end, if the directions are opposite). */
      const int offset = displayIdx - mspRange->max;
      result = sForward ? msp->sRange.max + offset : msp->sRange.min - offset;
    }
  
  /* Get the full display range of the match sequence. If the result is still out of range
   * then there's nothing to show for this qIdx. */
  const IntRange *fullSRange = mspGetFullSRange(msp, seqSelected, bc);
  
  if (!valueWithinRange(result, fullSRange))
    {
      result = UNSET_INT;
    }
  
  return result;
}



/* Given a base index on the reference sequence, find the corresonding base 
 * in the match sequence. The return value is always UNSET_INT if there is not
 * a corresponding base at this position. */
int mspGetMatchCoord(const MSP *msp, 
                     const int qIdx, 
                     const gboolean seqSelected,
                     const int numUnalignedBases,
                     BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  if (mspIsBlastMatch(msp) || mspIsExon(msp))
    {
      const gboolean inMspRange = valueWithinRange(qIdx, &msp->qRange);
      
      if (msp->gaps && g_slist_length(msp->gaps) >= 1 && inMspRange)
        {
          result = mspGetGappedAlignmentCoord(msp, qIdx, bc);
        }
      else if (!inMspRange && mspIsBlastMatch(msp))
        {
          /* The q index is outside the alignment range but if the option to show
           * unaligned sequence is enabled we may still have a valie result. */
          result = mspGetUnalignedCoord(msp, qIdx, seqSelected, numUnalignedBases, bc);
        }
      else
        {
          result = mspGetUngappedAlignmentCoord(msp, qIdx, bc);
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
 * is true. 
 * Returns a list of all the sequences that require the requested data. If you require a breakdown of
 * sequences-to-fetch by the fetch-mode that they should use, use getSeqsToPopulateByMode instead. */
static GList* getSeqsToPopulate(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData)
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


/* Checks the list of sequences for blixem to display to see which ones need populating with
 * sequence data (if getSequenceData is true) and/or the optional-columns' data (if getOptionalData
 * is true. 
 * Returns lists of all the sequences that require the requested data, categorised by the fetch 
 * mode that should be used to fetch them. The return value is a map of a GQuark (representing the
 * fetch-mode string) to the GList of sequences. */
static GHashTable* getSeqsToPopulateByMode(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData, const char *defaultFetchMode)
{
  GHashTable *resultTable = g_hash_table_new(g_direct_hash, g_direct_equal);
  
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
          GQuark dataTypeQuark = 0;
          
          if (blxSeq->dataType && blxSeq->dataType->bulkFetch)
            dataTypeQuark = g_quark_from_string(blxSeq->dataType->bulkFetch);
          else if (defaultFetchMode)
            dataTypeQuark = g_quark_from_string(defaultFetchMode);
          
          if (dataTypeQuark)
            {
              /* Get the result list for this fetch method. It's ok if it is 
               * null because the list will be created by g_list_prepend. */
              GList *resultList = g_hash_table_lookup(resultTable, GINT_TO_POINTER(dataTypeQuark));
              resultList = g_list_prepend(resultList, blxSeq);
              
              /* Update the existing (or insert the new) list */
              g_hash_table_insert(resultTable, GINT_TO_POINTER(dataTypeQuark), resultList);
            }
        }
    }
  
  return resultTable;
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
