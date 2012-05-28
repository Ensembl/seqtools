/*  File: blxview.c
 *  Author: Erik Sonnhammer, 1993-02-20
 *  Copyright (c) 2009 - 2012 Genome Research Ltd
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
#include <unistd.h>
#include <ctype.h>
#include <gdk/gdkkeysyms.h>

#include <blixemApp/blixem_.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blxdotter.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/blxparser.h>
#include <seqtoolsUtils/utilities.h>


#define MAXALIGNLEN			      10000
#define ORGANISM_PREFIX_SEPARATOR	      "  "
#define MKSTEMP_CONST_CHARS		      "BLIXEM_gff"	  /* the prefix to use when creating a temp file name */
#define MKSTEMP_REPLACEMENT_CHARS	      "XXXXXX"		  /* the required string that will be replaced by unique chars when creating a temp file name */
#define GFF3_VERSION_HEADER		      "##gff-version 3"	  /* the header line of a GFF v3 file */
#define GFF3_SEQUENCE_REGION_HEADER	      "##sequence-region" /* the start comment of the sequence-region comment line in a GFF v3 file */
#define MIN_GAP_HIGHLIGHT_WIDTH		      5			  /* minimum width of assembly gaps markers */


static void            blviewCreate(char *align_types, const char *paddingSeq, GArray* featureLists[], GList *seqList, GSList *supportedTypes, CommandLineOptions *options, const gboolean External) ;
static void            processGeneName(BlxSequence *blxSeq);
static void            processOrganism(BlxSequence *blxSeq);
static GHashTable*     getSeqsToPopulate(GList *inputList, const GArray *defaultFetchMethods, const int attempt);


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
      
      GError *error = NULL;
      BlxSeqType seqType = determineSeqType(options->refSeq, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
      
      if (seqType == BLXSEQ_PEPTIDE)
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
static gboolean socketFetchList(GList *seqsToFetch, 
                                BlxFetchMethod *fetchMethod,
                                GList *seqList, 
                                const gboolean External,
                                const BlxSeqType seqType,
                                GError **error)
{
  gboolean success = FALSE;

  if (fetchMethod->output == BLXFETCH_OUTPUT_EMBL)
    {
      success = populateFullDataPfetch(seqsToFetch,
                                       fetchMethod,
                                       External,
                                       seqType,
                                       error) ;
    }
  else if (fetchMethod->output == BLXFETCH_OUTPUT_FASTA)
    {
      success = populateFastaDataPfetch(seqsToFetch,
                                        fetchMethod,
                                        External,
                                        seqType,
                                        error) ;
    }
  else 
    {
      g_set_error(error, BLX_ERROR, 1, "Invalid output format for fetch method %s (expected '%s' or '%s')\n", 
                  g_quark_to_string(fetchMethod->name),
                  FETCH_OUTPUT_FASTA,
                  FETCH_OUTPUT_EMBL);
    }

  return success;   
}


#ifdef PFETCH_HTML
static gboolean httpFetchList(GList *seqsToFetch, 
                              BlxFetchMethod *fetchMethod,
                              GList *seqList, 
                              const BlxSeqType seqType,
                              GError **error)
{
  gboolean success = FALSE;

  success = populateSequenceDataHtml(seqsToFetch, seqType, fetchMethod);

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


void appendNewSequences(MSP *newMsps, GList *newSeqs, MSP **mspList, GList **seqList)
{
  /* Append new MSPs to MSP list */
  MSP *lastMsp = *mspList;
  
  while (lastMsp->next)
    lastMsp = lastMsp->next;
  
  lastMsp->next = newMsps;
  
  /* Append new sequences to sequence list */
  GList *seqItem = newSeqs;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      *seqList = g_list_prepend(*seqList, seqItem->data);
    }
}


/* Utility function to load the contents of the given file into blixem. The 
 * new features are appended onto the existing sequence and MSP lists. */
void loadGffFile(const char *fileName,
                 GKeyFile *keyFile,
                 BlxBlastMode *blastMode,
                 GArray* featureLists[],
                 GSList *supportedTypes, 
                 GSList *styles,
                 MSP **newMsps,
                 GList **newSeqs)
{
  if (!fileName)
    return;
  
  FILE *inputFile = fopen(fileName, "r");

  if (!inputFile)
    {
      g_critical("Failed to open file.\n");
      return;
    }
  
  char *dummyseq1 = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname1[FULLNAMESIZE+1] = "";
  char *dummyseq2 = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname2[FULLNAMESIZE+1] = "";
  
  parseFS(newMsps, inputFile, blastMode, featureLists, newSeqs, supportedTypes, styles,
          &dummyseq1, dummyseqname1, NULL, &dummyseq2, dummyseqname2, keyFile) ;
  
  fclose(inputFile);
}


/* Called by regionFetchList. Fetches the sequences for a specific region.
 *   tmpDir is the directory in which to place the temporary files
 *   script is the script to call to do the fetch
 *   dataset will be passed as the -dataset argument to the script if it is not null */
static void regionFetchFeature(const MSP const *msp, 
                               const BlxSequence const *blxSeq,
                               BlxFetchMethod *fetchMethod,
                               const char *script,
                               const char *dataset,
                               const char *tmpDir,
                               const int refSeqOffset,
                               BlxBlastMode *blastMode,
                               GList **seqList,
                               MSP **mspListIn,
                               GArray* featureLists[],
                               GSList *supportedTypes, 
                               GSList *styles,
                               const gboolean saveTempFiles,
                               const IntRange const *refSeqRange,
                               GError **error)
{
  /* Create a temp file for the results */
  char *fileName = blxprintf("%s/%s_%s", tmpDir, MKSTEMP_CONST_CHARS, MKSTEMP_REPLACEMENT_CHARS);
  int fileDesc = mkstemp(fileName);
  
  if (fileDesc == -1)
    {
      g_set_error(error, BLX_ERROR, 1, "Error creating temp file for region-fetch results (filename=%s)\n", fileName);
      g_free(fileName);
      return;
    }
  
  close(fileDesc);
  
  /* Get the command string, including the args */
  GString *command = getFetchCommand(fetchMethod, NULL, msp, mspGetRefName(msp), refSeqOffset, refSeqRange, dataset);
  
  /* Send the output to the temp file */
  if (fileName)
    g_string_append_printf(command, " > %s", fileName);

  FILE *outputFile = fopen(fileName, "w");

  g_debug("region-fetch command:\n%s\n", command->str);

  g_message_info("Calling region-fetch script...\n");
  const gboolean success = (system(command->str) == 0);
  
  fclose(outputFile);
  g_string_free(command, TRUE);
  
  if (success)
    {
      /* Parse the results */
      g_message_info("Parsing region-fetch results...");
      MSP *newMsps = NULL;
      GList *newSeqs = NULL;

      GKeyFile *keyFile = blxGetConfig();
      loadGffFile(fileName, keyFile, blastMode, featureLists, supportedTypes, styles, &newMsps, &newSeqs);
      appendNewSequences(newMsps, newSeqs, mspListIn, seqList);

      g_message_info(" complete.\n");
    }
  else
    {
      g_message_info("... failed.\n");
      g_critical("Failed to fetch sequences for region [%d, %d].\n", msp->qRange.min, msp->qRange.max);
    }
  
  /* Delete the temp file (unless the 'save temp files' option is on) */
  if (!saveTempFiles)
    {
      if (unlink(fileName) != 0)
        g_warning("Error removing temp file '%s'.\n", fileName);
    }
  
  g_free(fileName);
}


/* Fetch sequences for a given region. This uses the config file to find the
 * script and arguments to call to fetch the sequences.
 * The input GList contains a list of BlxSequences that are parent objects for
 * MSPs that identify regions. For each region, the script is called to fetch
 * all sequences that lie within that region, and the results are placed in 
 * a temporary GFF file, which is then parsed to get the results. The GFF file
 * is deleted when finished, unless the saveTempFiles argument is true. */
static void regionFetchList(GList *regionsToFetch, 
                            GList **seqList, 
                            BlxFetchMethod *fetchMethod, 
                            MSP **mspListIn,
                            BlxBlastMode *blastMode,
                            GArray* featureLists[],
                            GSList *supportedTypes, 
                            GSList *styles,
                            gboolean External,
                            const gboolean saveTempFiles,
                            const BlxSeqType seqType,
                            const int refSeqOffset,
                            const char *dataset,
                            const IntRange const *refSeqRange,
                            GError **error)
{
  /* Get the command to run */
  const gchar *script = fetchMethod->location;

  if (!script)
    {
      g_set_error(error, BLX_ERROR, 1, "Error fetching sequences; no command given for fetch-method [%s].\n", g_quark_to_string(fetchMethod->name));
      return;
    }

  /* Get the temp directory. Some systems seem to have a trailing slash, 
   * some not, so if it has one then remove it... */
  char *tmpDir = g_strdup(g_get_tmp_dir());
  
  if (tmpDir[strlen(tmpDir) - 1] == '/')
    tmpDir[strlen(tmpDir) - 1] = '\0';

  /* Loop through each region, creating a GFF file with the results for each region */
  GList *regionItem = regionsToFetch;
  GError *tmpError = NULL;

  for ( ; regionItem && !tmpError; regionItem = regionItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(regionItem->data);
      GList *mspItem = blxSeq->mspList;
    
      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP const *msp = (const MSP const*)(mspItem->data);

          /* Only fetch regions that are at least partly inside our display range */
          if (!rangesOverlap(&msp->qRange, refSeqRange))
            continue;
          
          regionFetchFeature(msp, blxSeq, fetchMethod, script, dataset, tmpDir, refSeqOffset,
                             blastMode, seqList, mspListIn, featureLists, supportedTypes,
                             styles, saveTempFiles, refSeqRange, &tmpError);
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* This function actually performs the work of fetching the 
 * given list of features. The first argument is the list 
 * of features to fetch and the second is the list of all
 * features. */
static gboolean fetchList(GList *seqsToFetch, 
                          GList **seqList,
                          BlxFetchMethod *fetchMethod,
                          const BlxSeqType seqType,
                          const gboolean saveTempFiles,
                          const gboolean External,
                          MSP **mspList,
                          BlxBlastMode *blastMode,
                          GArray* featureLists[],
                          GSList *supportedTypes, 
                          GSList *styles,
                          const int refSeqOffset,
                          const IntRange const *refSeqRange,
                          const char *dataset,
                          GError **error)
{
  gboolean success = TRUE;
  
  if (g_list_length(seqsToFetch) > 0)
    {
      if (fetchMethod)
        {
          if (fetchMethod->mode == BLXFETCH_MODE_SOCKET)
            {  
              success = socketFetchList(seqsToFetch,
                                        fetchMethod,
                                        *seqList,
                                        External,
                                        seqType,
                                        error);
            }
#ifdef PFETCH_HTML 
          else if (fetchMethod->mode == BLXFETCH_MODE_HTTP || fetchMethod->mode == BLXFETCH_MODE_PIPE)
            {
              success = httpFetchList(seqsToFetch,
                                      fetchMethod,
                                      *seqList,
                                      seqType,
                                      error);
            }
#endif
          else if (fetchMethod->mode == BLXFETCH_MODE_DB)
            {
              g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", g_quark_to_string(fetchMethod->name));
            }
          else if (fetchMethod->mode == BLXFETCH_MODE_COMMAND)
            {
              regionFetchList(seqsToFetch, 
                              seqList, 
                              fetchMethod, 
                              mspList,
                              blastMode,
                              featureLists,
                              supportedTypes,
                              styles,
                              External,
                              saveTempFiles,
                              seqType,
                              refSeqOffset,
                              dataset,
                              refSeqRange,
                              error);
            }
          else if (fetchMethod->mode == BLXFETCH_MODE_NONE)
            {
              /* do nothing */
            }
          else
            {
              g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", g_quark_to_string(fetchMethod->name));
            }
        }
      else
        {
          g_set_error(error, BLX_ERROR, 1, "Fetch mode not specified.\n");
        }
    }
  
  return success;
}


/* Once we've fetched all the sequences we need to do some post-processing. Loop 
 * twice: once to modify any fields in our own custom manner, and once more to
 * see if any with missing data can copy it from their parent. (Need to do these in
 * separate loops or we don't know if the data we're copying is processed or not.) */
void finaliseFetch(GList *seqList)
{
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


/* Find out if we need to fetch any sequences (they may all be 
 * contained in the input files so there might not be anything to 
 * fetch). If we do need to, then fetch them by the preferred method.
 * If the preferred fetch method fails, recusively try any other
 * fetch methods set up for each sequence until we have either fetched 
 * everything or run out of fetch methods to try.
 * 'attempt' should be passed as 0 for the first call. */
gboolean bulkFetchSequences(const int attempt,
                            gboolean External, 
                            const gboolean saveTempFiles,
                            const BlxSeqType seqType,
                            GList **seqList, /* list of BlxSequence structs for all required sequences */
                            const GArray *defaultFetchMethods,
                            GHashTable *fetchMethods,
                            MSP **mspList,
                            BlxBlastMode *blastMode,
                            GArray* featureLists[],
                            GSList *supportedTypes, 
                            GSList *styles,
                            const int refSeqOffset,
                            const IntRange const *refSeqRange,
                            const char *dataset)
{
  gboolean success = FALSE; /* will get set to true if any of the fetch methods succeed */
  
  /* Fetch any sequences that do not have their sequence data
   * already populated. If this is a re-try attempt, then use
   * a secondary fetch method, if one is given; otherwise, exclude
   * from the list (i.e. when we run out of fetch methods or everything
   * has been successfully fetched, then this table will be empty). */
  GHashTable *seqsTable = getSeqsToPopulate(*seqList, defaultFetchMethods, attempt);

  if (g_hash_table_size(seqsTable) < 1)
    {
      /* Nothing to fetch */
      success = TRUE;
    }
  else
    {
      /* Loop through each fetch mode */
      GHashTableIter iter;
      g_hash_table_iter_init(&iter, seqsTable);
      gpointer key, value;
      GError *error = NULL;
      
      while (g_hash_table_iter_next(&iter, &key, &value))
        {
          GQuark fetchMethodQuark = GPOINTER_TO_INT(key);
          BlxFetchMethod *fetchMethod = (BlxFetchMethod*)g_hash_table_lookup(fetchMethods, GINT_TO_POINTER(fetchMethodQuark));

          if (!fetchMethod)
            continue;

          GList *seqsToFetch = (GList*)value;          
          GError *tmpError = NULL;
          
          g_message_info("Fetching %d sequences (attempt %d, method=%s)\n", g_list_length(seqsToFetch), attempt + 1, g_quark_to_string(fetchMethod->name));
          
          if (fetchList(seqsToFetch, seqList, fetchMethod, seqType, saveTempFiles, External, mspList, blastMode, featureLists, supportedTypes, styles, refSeqOffset, refSeqRange, dataset, &tmpError))
            {
              success = TRUE;
              
              /* Compile all errors into a single error */
              if (error)
                {
                  prefixError(error, tmpError->message);
                  g_error_free(tmpError);
                  tmpError = NULL;
                }
              else
                {
                  error = tmpError;
                }
            }
          
          /* We're done with this list now, so free the memory. Don't delete it from
           * the table yet, though, because that will invalidate the iterators. */
          g_list_free(seqsToFetch);
        }
      
      if (success && error)
        {
          /* Some fetches succeeded, so just issue a warning */
          prefixError(error, "Error fetching sequences:\n");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
      else if (error)
        {
          /* All failed, so issue a critical warning */
          prefixError(error, "Error fetching sequences:\n");
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }

      /* Recurse to re-try any that failed. This returns straight
       * away everything was fetched successfully, or if there are
       * no more fetch methods to try. */
      success = bulkFetchSequences(attempt + 1, External, saveTempFiles, seqType, seqList,
                                   defaultFetchMethods, fetchMethods, mspList,
                                   blastMode, featureLists, supportedTypes, 
                                   styles, refSeqOffset, refSeqRange, dataset);
    }

  /* Clean up */
  g_hash_table_unref(seqsTable);

  return success;
}


/* Find any gaps in the reference sequence */
static void findAssemblyGaps(const char *refSeq, GArray *featureLists[], MSP **mspList, const IntRange const *refSeqRange)
{
  /* Find the last msp in the list so we can append new ones to the end */
  MSP *lastMsp = *mspList;
  while (lastMsp && lastMsp->next)
    lastMsp = lastMsp->next;
  
  /* Scan for the gap character */
  const char *cp = strchr(refSeq, SEQUENCE_CHAR_GAP);

  while (cp && *cp)
    {
      /* Found the start of a gap; remember the start coord */
      const int startCoord = cp - refSeq + refSeqRange->min;
    
      /* Loop until we find a non-gap character (or the end of the string) */
      while (cp && *cp == SEQUENCE_CHAR_GAP)
	++cp;
    
      const int endCoord = cp - refSeq - 1 + refSeqRange->min;
    
      MSP *msp = createEmptyMsp(&lastMsp, mspList);
      msp->type = BLXMSP_GAP;
      intrangeSetValues(&msp->qRange, startCoord, endCoord);
      
      featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);
    
      /* Continue looking for more gaps */
      cp = strchr(cp, SEQUENCE_CHAR_GAP);
    }
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
      if (blxGetConfig())
        {
          GError *error = NULL;
          blxInitConfig(NULL, options, &error);

          if (error)
            g_error("Config File Error: %s\n", error->message);
        }
    }
  
  validateInput(options);
  
  /* Find any assembly gaps (i.e. gaps in the reference sequence) */
  findAssemblyGaps(options->refSeq, featureLists, &options->mspList, &options->refSeqRange);
  
  gboolean status = bulkFetchSequences(
    0, External, options->saveTempFiles, options->seqType, &seqList, 
    options->bulkFetchDefault, options->fetchMethods, &options->mspList, &options->blastMode, 
    featureLists, supportedTypes, NULL, 0, &options->refSeqRange, 
    options->dataset); /* offset has not been applied yet, so pass offset=0 */

  if (status)
    {
      finaliseFetch(seqList);

      /* Construct missing data and do any other required processing now we have all the sequence data */
      finaliseBlxSequences(featureLists, &options->mspList, &seqList, 
                           options->refSeqOffset, options->seqType, 
                           options->numFrames, &options->refSeqRange, TRUE);

    }

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(align_types, padseq, featureLists, seqList, supportedTypes, options, External) ;
    }

  return status;
}


/* Find all of the different alignment types (i.e. source names) and 
 * compile them into a single string, separated by the given separator string */
static char* findAlignTypes(GArray* featureLists[], const char *separatorStr)
{
  GSList *sourceList = NULL;
  GString *resultStr = g_string_new(NULL);
  
  /* Loop through all MSPs of type 'match' */
  const GArray const *mspList = featureLists[BLXMSP_MATCH];
  
  int i = 0;
  const MSP *msp = mspArrayIdx(mspList, i);
  
  for ( ; msp; msp = mspArrayIdx(mspList, ++i))
    {
      const char *source = mspGetSource(msp);
      
      if (source)
        {
          /* See if we've already seen this source. If not, add it to the list,
           * and append it to the string. */
          GQuark quark = g_quark_from_string(source);
          
          if (!g_slist_find(sourceList, GINT_TO_POINTER(quark)))
            {
              sourceList = g_slist_prepend(sourceList, GINT_TO_POINTER(quark));
              
              if (resultStr->len)
                g_string_append(resultStr, separatorStr);
              
              g_string_append(resultStr, source);
            }
        }
    }
  
  char *result = resultStr->str;
  g_string_free(resultStr, FALSE);
  
  return result;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *align_types, 
			 const char *paddingSeq,
                         GArray* featureLists[],
                         GList *seqList,
                         GSList *supportedTypes,
			 CommandLineOptions *options,
                         const gboolean External)
{
  if (!blixemWindow)
    {
      /* Create the window */
      blixemWindow = createBlxWindow(options, paddingSeq, featureLists, seqList, supportedTypes, External);

      /* Set the window title. Get a description of all the alignment types 
       * (unless already supplied) */
      if (!align_types)
        align_types = findAlignTypes(featureLists, ", ");

      /* If no alignment description was set, create a generic description */
      if (!align_types)
        align_types = blxprintf("%s", options->seqType == BLXSEQ_PEPTIDE ? "peptide alignment" : "nucleotide alignment");
      
      char *title = blxprintf("Blixem (%s):   %s %s",
                              align_types,
                              (options->dataset ? options->dataset : ""),
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

GtkWidget* getBlixemWindow()
{
  return blixemWindow;
}


/* Reset any global variables / singleton instances */
void blviewResetGlobals()
{
  blixemWindow = NULL;
  destroyMessageList();
}


/* Checks the list of sequences for blixem to display to see which ones
 * need fetching.
 * Returns lists of all the sequences to be fetched, categorised by the fetch 
 * mode that should be used to fetch them. The return value is a map of
 * a GQuark (representing the fetch-mode string) to the GList of sequences
 * to be fetched.
 * This function can be called multiple times on the same sequences to
 * re-try fetching sequences with different fetch methods if the original
 * fetch method fails. Pass 'attempt' as 0 for the first try, 1 for 
 * the second etc. */
static GHashTable* getSeqsToPopulate(GList *inputList, 
                                     const GArray *defaultFetchMethods,
                                     const int attempt)
{
  GHashTable *resultTable = g_hash_table_new(g_direct_hash, g_direct_equal);
  
  /* Loop through the input list */
  GList *inputItem = inputList;
  
  for ( ; inputItem; inputItem = inputItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(inputItem->data);
      
      /* Check if sequence data was requested and is not already set. */
      gboolean getSeq = (blxSequenceRequiresSeqData(blxSeq) && blxSeq->sequence == NULL);
      
      /* Check if optional data was requested and is not already set. We can assume that
       * if any of the data fields is set then the parsing has been done for all of them
       * (and any remaining empty fields just don't have that data available) */
      getSeq |= (blxSequenceRequiresOptionalData(blxSeq) &&
                 !blxSeq->organism &&
                 !blxSeq->geneName &&
                 !blxSeq->tissueType &&
                 !blxSeq->strain);
      
      if (getSeq)
        {
          GQuark fetchMethodQuark = blxSequenceGetFetchMethod(blxSeq, TRUE, attempt, defaultFetchMethods);

          if (fetchMethodQuark)
            {
              /* Get the result list for this fetch method. It's ok if it is 
               * null because the list will be created by g_list_prepend. */
              GList *resultList = g_hash_table_lookup(resultTable, GINT_TO_POINTER(fetchMethodQuark));
              resultList = g_list_prepend(resultList, blxSeq);
              
              /* Update the existing (or insert the new) list */
              g_hash_table_insert(resultTable, GINT_TO_POINTER(fetchMethodQuark), resultList);
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


/* This function highlights any regions within the given display range that
 * are assembly gaps (which are BLXMSP_GAP type MSPs in the given MSP array).
 * The given rect defines the drawing area on the given drawable that corresponds
 * to the given display range. Gaps that lie within the display range are drawn
 * within the rect at the appropriate positions. */
void drawAssemblyGaps(GtkWidget *widget,
                      GdkDrawable *drawable,
                      GdkColor *color,
                      const gboolean displayRev,
                      GdkRectangle *rect, 
                      const IntRange const *dnaRange,
                      const GArray *mspArray)
{
  /* See if any gaps lie within the display range. */
  int i = 0;
  MSP *gap = mspArrayIdx(mspArray, i);
  
  for ( ; gap; gap = mspArrayIdx(mspArray, ++i))
    {
      if (rangesOverlap(&gap->qRange, dnaRange))
        {
	  /* Draw to max coord plus one (or min coord minus one if reversed) to be inclusive */
	  const int x1 = convertBaseIdxToRectPos(displayRev ? gap->qRange.min - 1 : gap->qRange.min, rect, dnaRange, TRUE, displayRev, TRUE);
	  const int x2 = convertBaseIdxToRectPos(gap->qRange.max + 1, rect, dnaRange, TRUE, displayRev, TRUE);
          
	  const int width = max(MIN_GAP_HIGHLIGHT_WIDTH, abs(x2 - x1));
	
          cairo_t *cr = gdk_cairo_create(drawable);
          gdk_cairo_set_source_color(cr, color);
          cairo_rectangle(cr, min(x1, x2), 0, width, widget->allocation.height);
          cairo_clip(cr);
          cairo_paint_with_alpha(cr, 0.1);
          cairo_destroy(cr);
        }
    }
  
}


/***************** end of file ***********************/
