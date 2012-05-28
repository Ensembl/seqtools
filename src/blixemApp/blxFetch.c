/*  File: blxFetch.c
 *  Author: Ed Griffiths, 2008-06-17
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
 * Description: Blixem functions for control of sequence fetching and
 *              display.
 *              
 *              Compiling with -DPFETCH_HTML includes code to issue
 *              pfetch requests to a proxy pfetch server using html.
 *              This requires libs libpfetch and libcurlobj
 *              which in turn require libcurl.
 *----------------------------------------------------------------------------
 */


#include <sys/socket.h> /* for socket(), connect(), send(), and recv() */
#include <netinet/in.h>
#include <arpa/inet.h>  /* for sockaddr_in and inet_addr() */
#include <netdb.h>					    /* for gethostbyname() */
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>

#include <seqtoolsUtils/utilities.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blixem_.h>

#ifdef PFETCH_HTML 
#include <libpfetch/libpfetch.h>
#endif



/* Blixem config error domain */
#define BLX_CONFIG_ERROR g_quark_from_string("Blixem config")

/* Error codes */
typedef enum
  {
    BLX_CONFIG_ERROR_NO_GROUPS,             /* no groups in config file */
    BLX_CONFIG_ERROR_INVALID_KEY_TYPE,      /* invalid key type given in config file */
    BLX_CONFIG_ERROR_MISSING_KEY,           /* a required key is missing */
    BLX_CONFIG_ERROR_INVALID_FETCH_MODE,    /* invalid fetch mode specified */
    BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT, /* invalid output format specified for fetch mode */
    BLX_CONFIG_ERROR_WARNINGS               /* warnings found while reading config file */
  } BlxConfigError;



/* States for parsing EMBL files */
typedef enum
  {
    PARSING_NEWLINE,             /* parser is at the start of a new line */
    PARSING_ID,                  /* parsing the 2-letter id at the start of a line */
    PARSING_SEQUENCE_HEADER,     /* first line of the SQ section (does not contain sequence data) */
    PARSING_SEQUENCE,            /* main body of the SQ section (contains sequence data on multiple lines) */
    PARSING_ORGANISM,            /* parsing an OS section */
    PARSING_GENE_NAME,           /* parsing a GN section */
    PARSING_FEATURE_TYPE,        /* parsing an FT section */
    PARSING_FT_TAG_NAME,         /* parsing a tag name within a FT section */
    PARSING_FT_UNQUOTED_SECTION, /* parsing the tag data of an FT tag (but not entered the quoted section yet) */
    PARSING_FT_TISSUE_TYPE,      /* parsing the /tissue_type tag data within a FT section */
    PARSING_FT_STRAIN,           /* parsing the /strain tag data within a FT section */
    PARSING_IGNORE,              /* we're parsing an area we can ignore */
    PARSING_FINISHED,            /* finished parsing all records */
    PARSING_CANCELLED            /* cancelled by user */
  } BlxEmblParserState;


/* Used to draw/update a progress meter, consider as private. */
typedef struct
{
  GtkWidget *top_level ;
  GtkWidget *progress ;
  GtkWidget *label ;
  GdkColor blue_bar_fg, red_bar_fg ;

  gulong widget_destroy_handler_id;

  gboolean cancelled ;
  int seq_total ;
} ProgressBarStruct, *ProgressBar ;



#ifdef PFETCH_HTML 
#define PFETCH_READ_SIZE 80	/* about a line */
#define PFETCH_FAILED_PREFIX "PFetch failed:"


typedef struct
{
  GtkWidget *blxWindow;
  GtkWidget *dialog;
  GtkTextBuffer *text_buffer;
  char *title;
  const BlxSequence *blxSeq;

  gulong widget_destroy_handler_id;
  PFetchHandle pfetch;
  gboolean got_response;
  int attempt;
} PFetchDataStruct, *PFetchData;


typedef struct
{
  gboolean connection_closed;                               /* Gets set to true when pfetch connection is closed */
  BlxEmblParserState parser_state;                          /* Current state of the parser */
  gboolean status ;                                         /* Gets set to false if any errors */
  char *err_txt ;					    /* if !status then err message here. */

  int seq_total ;                                           /* Number of sequences that were requested */
  int num_fetched;                                          /* Number of seqeunces fetched */
  int num_succeeded;                                        /* Number of sequences fetched successfully */

  gboolean stats ;					    /* TRUE means record and output stats. */
  int min_bytes, max_bytes, total_bytes, total_reads ;	    /* Stats. */

  GList *seqList ;                                          /* List of sequences to fetch */
  BlxSequence *currentSeq ;                                 /* Keeps track of which BlxSequence we're currently fetching */
  GList *currentSeqItem ;                                   /* Keeps track of which BlxSequence list item we're currently fetching */
  BlxSeqType seq_type ;                                     /* Whether sequences are nucleotide or peptide */
  const BlxFetchMethod *fetchMethod ;                       /* Details about this fetch method */
  
  ProgressBar bar ;                                         /* Provides graphical feedback about how many sequences have been fetched */
  PFetchHandle pfetch ;
  
  /* Additional fields for parsing full EMBL entries */
  char section_id[3];                                       /* 2-letter ID at the start of the line indicating the current section */
  gboolean found_end_quote;                                 /* Set to true when we're parsing a quoted section and we come across what looks like an end quote */
  GString *tag_name;                                        /* When parsing the FT section, the current tag name is stored this field */
} PFetchSequenceStruct, *PFetchSequence ;

#endif


#ifdef PFETCH_HTML 
static PFetchStatus pfetch_reader_func(PFetchHandle *handle,
				       char         *text,
				       guint        *actual_read,
				       GError       *error,
				       gpointer      user_data) ;
static void handle_dialog_close(GtkWidget *dialog, gpointer user_data);
static PFetchStatus pfetch_closed_func(gpointer user_data) ;


static PFetchStatus sequence_pfetch_reader(PFetchHandle *handle,
				       char         *text,
				       guint        *actual_read,
				       GError       *error,
				    gpointer      user_data) ;
static PFetchStatus sequence_pfetch_closed(PFetchHandle *handle, gpointer user_data) ;
static void sequence_dialog_closed(GtkWidget *dialog, gpointer user_data) ;
static gboolean parsePfetchHtmlBuffer(char *read_text, int length, PFetchSequence fetch_data) ;

static void httpFetchSequence(const BlxSequence *blxSeq, BlxFetchMethod *fetchMethod, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer *text_buffer, const int attempt, char **result_out);
#endif

static int socketConstruct(const char *ipAddress, int port, gboolean External) ;
static gboolean socketSend(int sock, const char *text) ;

static ProgressBar makeProgressBar(int seq_total) ;
static void updateProgressBar(ProgressBar bar, const char *sequence, int numFetched, gboolean fetch_ok) ;
static gboolean isCancelledProgressBar(ProgressBar bar) ;
static void destroyProgressBar(ProgressBar bar) ;
static void destroyProgressCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */
static void cancelCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */

static void readConfigFile(GKeyFile *key_file, char *config_file, CommandLineOptions *options, GError **error) ;

static void socketFetchSequence(const BlxSequence *blxSeq, BlxFetchMethod *fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, char **result_out);
static void commandFetchSequence(const BlxSequence *blxSeq, BlxFetchMethod *fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, char **result_out);
static void wwwFetchSequence(const BlxSequence *blxSeq, BlxFetchMethod *fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, char **result_out);

/* Pfetch local functions */
static void                     appendCharToString(const char curChar, GString **result);
static void                     appendCharToQuotedString(const char curChar, gboolean *foundEndQuote, GString **result);

static gboolean                 socketFetchInit(char *pfetchOptions, GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, int *sock);
static void                     checkProgressBar(ProgressBar bar, BlxEmblParserState *parserState, gboolean *status);

static int                      pfetchReceiveBuffer(char *buffer, const int bufferSize, const int sock, 
                                                    BlxEmblParserState *parserState, gboolean *status);

static void                     pfetchGetNextSequence(BlxSequence **currentSeq, GList **currentSeqItem, ProgressBar bar, 
                                                      const int numRequested, const int numFetched, const gboolean pfetch_ok, 
                                                      gboolean *status, BlxEmblParserState *parserState);

static void                     pfetchParseSequenceFileBuffer(const char *buffer, const int lenReceived, BlxSequence **currentSeq, GList **currentSeqItem,
                                                              ProgressBar bar, const int numRequested,int *numFetched, int *numSucceeded, 
                                                              const BlxSeqType seqType, BlxEmblParserState *parserState, gboolean *status, GError **error);

static gboolean                 pfetchGetParserStateFromId(const char *sectionId, BlxSequence *currentSeq, GString *tagName, BlxEmblParserState *parserState);

static void                     pfetchParseEmblFileBuffer(const char *buffer, const int lenReceived, BlxSequence **currentSeq, GList **currentSeqItem,
                                                          ProgressBar bar, const int numRequested, int *numFetched, int *numSucceeded, char *sectionId,
                                                          GString *tagName, gboolean *foundEndQuote,
                                                          const BlxSeqType seqType, BlxEmblParserState *parserState, gboolean *status);

static gboolean                 pfetchFinishSequence(BlxSequence *currentSeq, const BlxSeqType seqType, const int numRequested, int *numFetched, int *numSucceeded, 
                                                     BlxEmblParserState *parserState);

static void                     pfetchFinishEmblFile(BlxSequence **currentSeq, GList **currentSeqItem, const int numRequested, int *numFetched, 
                                                     int *numSucceeded, ProgressBar bar, const BlxSeqType seqType, BlxEmblParserState *parserState, gboolean *status);

static void                     pfetchGetParserStateFromTagName(GString *tagName, BlxEmblParserState *parserState);



/* global configuration object for blixem. */
static GKeyFile *blx_config_G = NULL ;


/* Fetch the given sequence and optionally display the results. Optionally
 * return the sequence; if the result_out is supplied, its contents must be
 * freed by the caller using g_free */
void fetchSequence(const BlxSequence *blxSeq, 
                   const gboolean displayResults,
                   const int attempt,
                   GtkWidget *blxWindow,
                   char **result_out)
{
  g_assert(blxSeq);
  
  /* Look up the fetch method for this sequence */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GQuark fetchMethodQuark = blxSequenceGetFetchMethod(blxSeq, FALSE, attempt, bc->userFetchDefault);
  BlxFetchMethod *fetchMethod = (BlxFetchMethod*)g_hash_table_lookup(bc->fetchMethods, GINT_TO_POINTER(fetchMethodQuark));

  if (!fetchMethod)
    {
      /* If this is the first attempt then we should have a fetch method; 
       * therefore give a warning if it was not found */
      if (attempt == 0)
        {
          if (!fetchMethodQuark)
            g_warning("Failed to find fetch method for sequence '%s'\n", blxSequenceGetFullName(blxSeq));
          else
            g_warning("Error fetching sequence '%s'; could not find details for fetch method '%s'\n", blxSequenceGetFullName(blxSeq), g_quark_to_string(fetchMethodQuark));
        }

      return;
    }

  if (fetchMethod->mode == BLXFETCH_MODE_NONE)
    {
      g_message("Fetch method for '%s' is '%s'\n", blxSequenceGetFullName(blxSeq), FETCH_MODE_NONE_STR);
      return;
    }
  
  g_message("Fetching '%s' using method '%s' (attempt %d)\n", blxSequenceGetFullName(blxSeq), g_quark_to_string(fetchMethodQuark), attempt + 1);

  
  if (fetchMethod->mode == BLXFETCH_MODE_SOCKET)
    {
      socketFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, result_out);
    }
#ifdef PFETCH_HTML 
  else if (fetchMethod->mode == BLXFETCH_MODE_HTTP || fetchMethod->mode == BLXFETCH_MODE_PIPE)
    {
      httpFetchSequence(blxSeq, fetchMethod, blxWindow, NULL, NULL, attempt, result_out);
    }
#endif
  else if (fetchMethod->mode == BLXFETCH_MODE_COMMAND)
    {
      commandFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, result_out);
    }
  else if (fetchMethod->mode == BLXFETCH_MODE_WWW)
    {
      wwwFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, result_out);
    }
  else
    {
      /* Invalid fetch method. Try again with the next fetch method, if one is specified */
      g_warning("Unknown fetch method: %s\n", g_quark_to_string(fetchMethod->name));
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, result_out);
    }
}



/* Get the command to call for the given fetch method. Parses
 * the command and arguments and substitutes the following
 * characters with values from the given sequence/msp. Note that 
 * either blxSeq or msp must be given, but not both.
 * 
 *   %p:      program name
 *   %h:      host name
 *   %u:      user name
 *   %m:      match sequence name
 *   %r:      reference sequence name
 *   %s:      start coord
 *   %e:      end coord
 *   %d:      dataset
 * 
 * Returns the command and args compiled into a single string.
 * The caller must free the result with g_string_free.
 * Returns an empty string if the command/args are empty.
*/
GString* getFetchCommand(BlxFetchMethod *fetchMethod,
                         const BlxSequence *blxSeq,
                         const MSP* const msp,
                         const char *refSeqName,
                         const int refSeqOffset,
                         const IntRange* const refSeqRange,
                         const char *dataset)
{
  g_assert(blxSeq || msp);

  GString *result = g_string_new("");
  
  /* Compile the command and args into a single string */
  char *command = NULL;
  
  if (fetchMethod->location && fetchMethod->args)
    command = g_strdup_printf("%s %s", fetchMethod->location, fetchMethod->args);
  else if (fetchMethod->location)
    command = g_strdup(fetchMethod->location);
  else if (fetchMethod->args)
    command = g_strdup(fetchMethod->args);
  else
    return result;

  /* Extract info about the sequence / feature */
  const char *name = blxSeq ? blxSequenceGetFullName(blxSeq) : mspGetSName(msp);
  int startCoord = blxSeq ? blxSequenceGetStart(blxSeq, blxSeq->strand) : mspGetQStart(msp);
  int endCoord = blxSeq ? blxSequenceGetEnd(blxSeq, blxSeq->strand) : mspGetQEnd(msp);
  startCoord += refSeqOffset;
  endCoord += refSeqOffset;
  boundsLimitValue(&startCoord, refSeqRange);
  boundsLimitValue(&endCoord, refSeqRange);

  /* Loop through the command and substitute any special chars with the relevant info  */
  char *c = command;
      
  for ( ; c && *c; ++c)
    {
      /* If it's preceded by the special char, substitute it for the real value */
      if (*c == '%')
        {
          /* Move to the next char, which should tell us what type of substitution to make */
          ++c;
          
          if (c && *c)
            {
              switch (*c) {
              case 'p': 
                g_string_append(result, g_get_prgname());
                break;
              case 'h': 
                g_string_append(result, g_get_host_name());
                break;
              case 'u': 
                g_string_append(result, g_get_user_name());
                break;
              case 'm': 
                if (name)
                  g_string_append(result, name);
                break;
              case 'r': 
                if (refSeqName)
                  g_string_append(result, refSeqName);
                break;
              case 's': 
                g_string_append_printf(result, "%d", startCoord);
                break;
              case 'e': 
                g_string_append_printf(result, "%d", endCoord);
                break;
              case 'd':
                if (dataset)
                  g_string_append(result, dataset);
              default:
                break;
              };
            }
        }
      else
        {
          /* Normal char; just append to the result */
          g_string_append_c(result, *c);
        }
    }

  /* Clean up */
  g_free(command);

  return result;
}


/* Use the www-fetch method to fetch an entry and optionally display
 * the results in a dialog.
 * Opens a browser to display the results. Does nothing if 
 * not displaying results! */
static void wwwFetchSequence(const BlxSequence *blxSeq,
                             BlxFetchMethod *fetchMethod, 
                             const gboolean displayResults, 
                             const int attempt,
                             GtkWidget *blxWindow,
                             char **result_out)
{
  if (displayResults)
    {
      BlxViewContext *bc = blxWindowGetContext(blxWindow);

      GString *url = getFetchCommand(fetchMethod, 
                                     blxSeq, 
                                     NULL, 
                                     bc->refSeqName, 
                                     bc->refSeqOffset,
                                     &bc->refSeqRange,
                                     bc->dataset);

      GError *error = NULL;
      seqtoolsLaunchWebBrowser(url->str, &error);

      g_string_free(url, TRUE);

      /* If failed, re-try with the next-preferred fetch method, if there is one */
      if (error)
        fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, result_out);
    }
}


/* Use the command-fetch method to fetch an entry and optionally display
 * the results in a dialog. */
static void commandFetchSequence(const BlxSequence *blxSeq,
                                 BlxFetchMethod *fetchMethod, 
                                 const gboolean displayResults, 
                                 const int attempt,
                                 GtkWidget *blxWindow,
                                 char **result_out)
{
  const char *seqName = blxSequenceGetFullName(blxSeq);
  char *command = blxprintf("efetch '%s' &", seqName);
  
  GError *error = NULL;
  GString *resultText = getExternalCommandOutput(command, &error);

  if (resultText && resultText->str)
    {
      if (displayResults && !error)
        {
          char *title = blxprintf("%s - %s", g_get_prgname(), command);
          displayFetchResults(title, resultText->str, blxWindow, NULL);
          g_free(title);
        }

      if (result_out)
        {
          *result_out = resultText->str;
          g_string_free(resultText, FALSE);
        }
      else
        {
          g_string_free(resultText, TRUE);
        }
    }
  else
    {
      /* Try again with the next-preferred fetch method, if there is one */
      g_string_free(resultText, TRUE);
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, result_out);
    }
  
  g_free(command);
}


/* Use the given socket-fetch method to fetch an entry and optionally display the results. */
static void socketFetchSequence(const BlxSequence *blxSeq, 
                                BlxFetchMethod *fetchMethod, 
                                const gboolean displayResults, 
                                const int attempt,
                                GtkWidget *blxWindow,
                                char **result_out)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Run the fetch command */
  GError *error = NULL;
  GString *command = getFetchCommand(fetchMethod, blxSeq, NULL, bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, bc->dataset);
  GString *resultText = getExternalCommandOutput(command->str, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);

  if (resultText->len)  /* Success */
    {
      if (displayResults)
        {
          char *title = blxprintf("Blixem - %s", command->str);
          displayFetchResults(title, resultText->str, blxWindow, NULL);
          g_free(title);
        }

      if (result_out)
        {
          /* caller takes ownership of str */
          *result_out = resultText->str;
          g_string_free(resultText, FALSE);
        }
      else
        {
          g_string_free(resultText, TRUE);
        }
    }
  else   /* Failed */
    {
      g_string_free(resultText, TRUE);

      /* Try again with the next fetch method, if there is one set */
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, result_out);
    }

  g_string_free(command, TRUE);
}


#ifdef PFETCH_HTML 
/* Gets all the sequences needed by blixem but from http proxy server instead of from
 * the pfetch server direct, this enables blixem to be run and get sequences from
 * anywhere that can see the http proxy server. If fullEntry is true then the full
 * EMBL entry is fetched (if it exists) and data for the optional columns is parsed (e.g.
 * organism and tissue-type) */
gboolean populateSequenceDataHtml(GList *seqsToFetch, const BlxSeqType seqType, const BlxFetchMethod* const fetchMethod)
{
  gboolean status = FALSE ;
  gboolean debug_pfetch = FALSE ;

  GType pfetch_type = PFETCH_TYPE_HTTP_HANDLE ;
  PFetchSequenceStruct fetch_data = {FALSE} ;

  fetch_data.connection_closed = FALSE;
  fetch_data.parser_state = PARSING_NEWLINE;
  fetch_data.status = TRUE ;
  fetch_data.err_txt = NULL;
  
  fetch_data.seq_total = g_list_length(seqsToFetch);
  fetch_data.num_fetched = 0;
  fetch_data.num_succeeded = 0;
  
  fetch_data.stats = FALSE ;
  fetch_data.min_bytes = INT_MAX ;
  fetch_data.max_bytes = 0 ;
  fetch_data.total_bytes = 0 ;
  fetch_data.total_reads = 0 ;

  fetch_data.seqList = seqsToFetch;
  fetch_data.currentSeqItem = seqsToFetch;
  fetch_data.currentSeq = seqsToFetch ? (BlxSequence*)(seqsToFetch->data) : NULL;
  fetch_data.seq_type = seqType ;
  fetch_data.fetchMethod = fetchMethod;

  fetch_data.bar = makeProgressBar(fetch_data.seq_total) ;
  fetch_data.pfetch = PFetchHandleNew(pfetch_type);

  fetch_data.section_id[0] = ' ';
  fetch_data.section_id[1] = ' ';
  fetch_data.section_id[2] = '\0';
  fetch_data.found_end_quote = FALSE;
  fetch_data.tag_name = g_string_new("");
  
  g_signal_connect(G_OBJECT(fetch_data.bar->top_level), "destroy",
    	       G_CALLBACK(sequence_dialog_closed), &fetch_data) ;
  
  const gboolean fullEntry = fetchMethod->output == BLXFETCH_OUTPUT_EMBL; /* fetching full embl entry */

  PFetchHandleSettings(fetch_data.pfetch, 
                       "full",       fullEntry,
                       "port",       fetchMethod->port,
                       "debug",      debug_pfetch,
                       "pfetch",     fetchMethod->location,
                       "cookie-jar", fetchMethod->cookie_jar,
                       NULL);

 if (!fullEntry)
  {
    PFetchHandleSettings(fetch_data.pfetch,
                         "blixem-seqs",  TRUE, /* turns the next two on */
                         "one-per-line", TRUE, /* one sequence per line */
                         "case-by-type", TRUE, /* lowercase for DNA, uppercase for protein */
                         NULL);
  }
  
  g_signal_connect(G_OBJECT(fetch_data.pfetch), "reader", G_CALLBACK(sequence_pfetch_reader), &fetch_data) ;

  g_signal_connect(G_OBJECT(fetch_data.pfetch), "closed", G_CALLBACK(sequence_pfetch_closed), &fetch_data) ;


  /* Build up a string containing all the sequence names. */
  GString *seq_string = g_string_sized_new(1000) ;
  GList *seqItem = seqsToFetch;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      g_string_append_printf(seq_string, "%s ", blxSequenceGetFullName(blxSeq));
    }

  /* Set up pfetch/curl connection routines, this is non-blocking so if connection
   * is successful we block using our own flag. */
  if (PFetchHandleFetch(fetch_data.pfetch, seq_string->str) == PFETCH_STATUS_OK)
    {
      status = TRUE ;

      while (!(fetch_data.connection_closed))
        {
          gtk_main_iteration() ;
        }

      status = fetch_data.status ;
      
      if (!status)
        {
          if (fetch_data.parser_state != PARSING_CANCELLED)
            {
              g_critical("Sequence fetch from http server failed: %s\n", fetch_data.err_txt) ;
              
              if (fetch_data.err_txt)
                { 
                  g_free(fetch_data.err_txt) ;
                }
            }
        }
    }
  else
    {
      status = FALSE ;
    }
  
  g_string_free(seq_string, TRUE);

  destroyProgressBar(fetch_data.bar) ;

  return status ;
}


/* Use the http proxy to pfetch an entry */
/* Note that this uses a callback to update the display 
 * window; it cannot return the result immediately and the 
 * code is not currently structured to allow the callback to
 * do anything other than update the display window, so if
 * we're just requesting the sequence, this currently just 
 * returns null. */
static void httpFetchSequence(const BlxSequence *blxSeq,
                              BlxFetchMethod *fetchMethod,
                              GtkWidget *blxWindow, 
                              GtkWidget *dialog, 
                              GtkTextBuffer *text_buffer, 
                              const int attempt,
                              char **result_out)
{
  const char *sequence_name = blxSequenceGetFullName(blxSeq);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  gboolean debug_pfetch = FALSE ;
  
  if (fetchMethod->location != NULL)
    {
      PFetchData pfetch_data ;
      PFetchHandle pfetch = NULL ;
      GType pfetch_type = PFETCH_TYPE_HTTP_HANDLE ;

      if (fetchMethod->mode == BLXFETCH_MODE_PIPE)
	pfetch_type = PFETCH_TYPE_PIPE_HANDLE ;
	
      pfetch_data = g_new0(PFetchDataStruct, 1);

      pfetch_data->pfetch = pfetch = PFetchHandleNew(pfetch_type);
      
      pfetch_data->blxWindow = blxWindow;
      pfetch_data->blxSeq = blxSeq;
      pfetch_data->attempt = attempt;

      GString *command = getFetchCommand(fetchMethod, blxSeq, NULL, bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, bc->dataset);
      pfetch_data->title = g_strdup_printf("%s - %s", g_get_prgname(), command->str);
      
      if (pfetch_data->title)
	{
          if (dialog && text_buffer)
            {
              gtk_window_set_title(GTK_WINDOW(dialog), pfetch_data->title);
              pfetch_data->dialog = dialog;
              pfetch_data->text_buffer = text_buffer;
            }
          else
            {
              pfetch_data->dialog = displayFetchResults(pfetch_data->title, "pfetching...\n", blxWindow, &pfetch_data->text_buffer);
            }

	  pfetch_data->widget_destroy_handler_id = 
	    g_signal_connect(G_OBJECT(pfetch_data->dialog), "destroy", 
			     G_CALLBACK(handle_dialog_close), pfetch_data); 
	}
      
      const gboolean fetchFullEntry = fetchMethod->output == BLXFETCH_OUTPUT_EMBL; /* fetching full embl entry */

      if (PFETCH_IS_HTTP_HANDLE(pfetch))
        {
            PFetchHandleSettings(pfetch, 
                                 "full",       fetchFullEntry,
                                 "port",       fetchMethod->port,
                                 "debug",      debug_pfetch,
                                 "pfetch",     fetchMethod->location,
                                 "cookie-jar", fetchMethod->cookie_jar,
                                 NULL);
        }
      else
        {
            PFetchHandleSettings(pfetch, 
                                 "full",       fetchFullEntry,
                                 "pfetch",     fetchMethod->location,
                                 NULL);
        }
      
      g_signal_connect(G_OBJECT(pfetch), "reader", G_CALLBACK(pfetch_reader_func), pfetch_data);
      
      g_signal_connect(G_OBJECT(pfetch), "closed", G_CALLBACK(pfetch_closed_func), pfetch_data);
      
      PFetchHandleFetch(pfetch, (char*)sequence_name) ;
    }
  else
    {
      g_critical("%s", "Failed to obtain preferences specifying how to pfetch.\n") ;
    }
}


#endif



/*  getsseqsPfetch() adapted from Tony Cox's code pfetch.c
 *
 *  - this version incorporates a progress monitor as a window,
 *    much easier for user to control + has a cancel button.
 */
gboolean populateFastaDataPfetch(GList *seqsToFetch, 
                                 BlxFetchMethod *fetchMethod,
                                 gboolean External, 
                                 const BlxSeqType seqType, 
                                 GError **error)
{
  /* Initialise and send the requests */
  int sock;
  
  gboolean status = socketFetchInit(fetchMethod->args, seqsToFetch, fetchMethod->node, fetchMethod->port, External, &sock);
  
  /* Get the sequences back. They will be returned in the same order that we asked for them, i.e. 
   * in the order they are in our list. */
  GList *currentSeqItem = seqsToFetch;
  BlxSequence *currentSeq = (BlxSequence*)(currentSeqItem->data);
  
  if (status && currentSeq)
    {
      int numRequested = g_list_length(seqsToFetch); /* total number of sequences requested */
      ProgressBar bar = makeProgressBar(numRequested);
      
      enum {RCVBUFSIZE = 256} ;               /* size of receive buffer */
      char buffer[RCVBUFSIZE + 1] ;           /* receive buffer */
      
      int numFetched = 0;
      int numSucceeded = 0;

      GError *tmpError = NULL;
      BlxEmblParserState parserState = PARSING_NEWLINE;

      while (status && parserState != PARSING_CANCELLED && parserState != PARSING_FINISHED)
	{
          /* Receive and parse the next buffer */
          checkProgressBar(bar, &parserState, &status);
          int lenReceived = pfetchReceiveBuffer(buffer, RCVBUFSIZE, sock, &parserState, &status);
          
          pfetchParseSequenceFileBuffer(buffer, 
                                        lenReceived, 
                                        &currentSeq, 
                                        &currentSeqItem, 
                                        bar, 
                                        numRequested, 
                                        &numFetched, 
                                        &numSucceeded, 
                                        seqType,
                                        &parserState, 
                                        &status,
                                        &tmpError);
	}
      
      /* Finish up */
      shutdown(sock, SHUT_RDWR);
      destroyProgressBar(bar);
      bar = NULL ;
      
      if (status && numSucceeded != numRequested)
	{
	  double proportionOk = (float)numSucceeded / (float)numRequested;
          
	  /* We don't display an error message unless lots of sequences don't get fetched
	   * because users find it annoying as most of the time they don't mind if the
	   * odd sequence isn't fetched successfully. */
	  if (proportionOk < 0.5)
            {
              g_critical("pfetch sent back %d when %d requested\n", numSucceeded, numRequested) ;
            }
	  else
            {
              g_message("pfetch sent back %d when %d requested\n", numSucceeded, numRequested) ;
            }
	}

      if (tmpError)
        {
          g_propagate_error(error, tmpError);
        }
    }
  
  return status ;
}


/*  getsseqsPfetch() adapted from Tony Cox's code pfetch.c
 *
 *  - this version incorporates a progress monitor as a window,
 *    much easier for user to control + has a cancel button.
 *  - same as populateFastaDataPfetch but fetches the full EMBL entry
 *    in order to parse additional info such as organism. This
 *    may fail for protein variants though so blxGetSseqsPfetch
 *    should subsequently be run on any that still don't have their
 *    sequence filled in.
 *  - can be called after the fasta sequence data is already populated:
 *    in that case it will ignore the sequence data and just populate
 *    the additional data.
 *  - sequence data will also be ignored for sequences that do not 
 *    require sequence data
 */
gboolean populateFullDataPfetch(GList *seqsToFetch, 
                                BlxFetchMethod *fetchMethod,
                                gboolean External, 
                                const BlxSeqType seqType, 
                                GError **error)
{
  int sock;
  gboolean status = socketFetchInit(fetchMethod->args, seqsToFetch, fetchMethod->node, fetchMethod->port, External, &sock);

  /* Get the sequences back. They will be returned in the same order that we asked for them, i.e. 
   * in the order they are in our list. */
  GList *currentSeqItem = seqsToFetch;
  BlxSequence *currentSeq = (BlxSequence*)(currentSeqItem->data);

  if (status && currentSeq)
    {
      int numRequested = g_list_length(seqsToFetch); /* total number of sequences requested */
      ProgressBar bar = makeProgressBar(numRequested);
      
      enum {RCVBUFSIZE = 256} ;               /* size of receive buffer */
      char buffer[RCVBUFSIZE + 1] ;           /* receive buffer */

      int numFetched = 0;
      int numSucceeded = 0;

      /* All lines start with a two-letter identifier, which will be parsed into this string */
      char sectionId[3] = "  ";

      /* In the FT (feature type) section, we look for tags of the format /tagname="value". Use
       * the following string to parse the tag name into */
      GString *tagName = g_string_new("");
      
      /* When we're in a quoted section, this is used to flag that we've found a second quote that
       * we think is the end of the quoted section. However, a double quote means an escaped quote, so
       * if the next char is also a quote, we know we need to include it in the text and carry on parsing. */
      gboolean foundEndQuote = FALSE;

      
      BlxEmblParserState parserState = PARSING_NEWLINE;
      
      while (status && parserState != PARSING_CANCELLED && parserState != PARSING_FINISHED)
	{
          /* Receive and parse the next buffer */
          checkProgressBar(bar, &parserState, &status);
          const int lenReceived = pfetchReceiveBuffer(buffer, RCVBUFSIZE, sock, &parserState, &status);
          
          pfetchParseEmblFileBuffer(buffer, 
                                    lenReceived,
                                    &currentSeq, 
                                    &currentSeqItem, 
                                    bar, 
                                    numRequested, 
                                    &numFetched, 
                                    &numSucceeded, 
                                    sectionId,
                                    tagName,
                                    &foundEndQuote,
                                    seqType,
                                    &parserState, 
                                    &status);
	}
      
      /* Finish up */
      shutdown(sock, SHUT_RDWR);
      destroyProgressBar(bar);
      bar = NULL ;
      
      g_string_free(tagName, TRUE);

      if (status && numSucceeded != numRequested)
	{
	  double proportionOk = (float)numSucceeded / (float)numRequested;

	  /* We don't display a critical error message when fetching the full EMBL file because we're
           * going to re-try fetching just the fasta data anyway. Display a warning, or just an info
           * message if a small proportion failed */
	  if (proportionOk < 0.5)
            {
              g_warning("pfetch of full EMBL entries sent back %d when %d requested\n", numSucceeded, numRequested) ;
            }
	  else
            {
              g_message("pfetch of full EMBL entries sent back %d when %d requested\n", numSucceeded, numRequested) ;
            }
	}
    }

  return status ;
}


/* Set/Get global config, necessary because we don't have some blixem context pointer....
 * To do: we do have a context now, so this should be moved to there. 
 * Sets the error if there were any problems. Note that the error is not set if the
 * config file does not exist or is empty  */
void blxInitConfig(char *config_file, CommandLineOptions *options, GError **error)
{
  g_assert(!blx_config_G) ;

  blx_config_G = g_key_file_new() ;

  /* First, get any persistent settings from the settings file */
  char *settings_file = blxprintf("%s/%s", g_get_home_dir(), BLIXEM_SETTINGS_FILE);
  gchar *content1 = NULL;
  gsize len1 = 0;

  if (g_file_test(settings_file, G_FILE_TEST_EXISTS))
    {
      /* This will become the global config file if no other file is given. If
       * another is given, this will be overwritten, so save the content of the settings file. */
      g_key_file_load_from_file(blx_config_G, settings_file, G_KEY_FILE_NONE, NULL);
      content1 = g_key_file_to_data(blx_config_G, &len1, NULL);
    }
  
  /* Load the given config file, if any */
  if (config_file)
    {
      readConfigFile(blx_config_G, config_file, options, error);

      if (content1)
        {
          /* Merge both file contents. First, get the new content as a string and concatenate */
          gsize len2 = 0;
          gchar *content2 = g_key_file_to_data(blx_config_G, &len2, NULL);

          gsize len = len1 + len2 + 2; /* need 2 extra chars for the sprintf below; newline and terminating nul */
          gchar *content = g_malloc(len);
          sprintf(content, "%s\n%s", content1, content2);

          /* Load the concatenated contents into a new key file */
          GKeyFile *key_file = g_key_file_new();
          if (g_key_file_load_from_data(key_file, content, len, G_KEY_FILE_NONE, NULL))
            {
              /* Delete the original config key file and replace it with the new one */
              g_key_file_free(blx_config_G);
              blx_config_G = key_file;
            }
          else
            {
              g_key_file_free(key_file);
            }

          g_free(content);
          g_free(content2);
        }
    }
  else if (!content1)
    {
      /* Neither the config file or the settings file had any content, so
       * delete the key file. */
      g_key_file_free(blx_config_G) ;
      blx_config_G = NULL ;
    }

  if (content1)
    g_free(content1);
}



GKeyFile *blxGetConfig(void)
{
  return blx_config_G ;
}


/* 
 *                        Internal functions.
 */


static int socketConstruct (const char *ipAddress, int port, gboolean External)
{
  int sock ;                       /* socket descriptor */
  struct sockaddr_in *servAddr ;   /* echo server address */
  struct hostent *hp ;


  /* Create a reliable, stream socket using TCP */
  if ((sock = socket(PF_INET, SOCK_STREAM, 0)) < 0)
    {
      g_critical ("socket() failed\n") ;
      return -1 ;
    }


  /* Construct the server address structure */
  servAddr = (struct sockaddr_in *) g_malloc (sizeof (struct sockaddr_in)) ;
  hp = gethostbyname(ipAddress) ;
  if (!hp)
  {
    if (External)
    {
        g_critical("Failed to start external blixem: unknown host \"%s\"\n", ipAddress);
        return -1;
    }
    else
    {
        g_critical("Failed to start internal blixem: unknown host \"%s\"\n", ipAddress);
        return -1;
    }
  }

  servAddr->sin_family = AF_INET ;			    /* Internet address family */
  bcopy((char*)hp->h_addr, (char*) &(servAddr->sin_addr.s_addr), hp->h_length) ;
							    /* Server IP address */
  servAddr->sin_port = htons ((unsigned short) port) ;	    /* Server port */


  /* Establish the connection to the server */
  if (connect(sock, (struct sockaddr *) servAddr, sizeof(struct sockaddr_in)) < 0)
    {
      g_critical ("socket connect() to BLIXEM_PFETCH = %s failed\n", ipAddress) ;
      sock = -1 ;
    }

  g_free (servAddr) ;

  return sock ;
}


static gboolean socketSend (int sock, const char *text)
{
  gboolean status = TRUE ;
  int len, bytes_to_send, bytes_written ;
  char *tmp ;
  struct sigaction oursigpipe, oldsigpipe ;

  /* The adding of 0x20 to the end looks wierd but I think it's because the  */
  /* server may not hold strings in the way C does (i.e. terminating '\0'),  */
  /* so 0x20 is added to mark the end of the string and the C string term-   */
  /* inator is moved up one and is not actually sent.                        */
  len = strlen(text) ;
  tmp = g_malloc(len + 2) ;
  strcpy(tmp, text) ;
  tmp[len] = 0x20 ;					    /* Add new string terminator. */
  tmp[len + 1] = 0 ;
  bytes_to_send = len + 1 ;

  /* send() can deliver a SIGPIPE if the socket has been disconnected, by    */
  /* ignoring it we will receive -1 and can look for EPIPE as the errno.     */
  oursigpipe.sa_handler = SIG_IGN ;
  sigemptyset(&oursigpipe.sa_mask) ;
  oursigpipe.sa_flags = 0 ;
  if (sigaction(SIGPIPE, &oursigpipe, &oldsigpipe) < 0)
    g_error("Cannot set SIG_IGN for SIGPIPE for socket write operations.\n") ;

  bytes_written = send(sock, tmp, bytes_to_send, 0) ;
  if (bytes_written == -1)
    {
      if (errno == EPIPE || errno == ECONNRESET || errno == ENOTCONN)
	{
	  status = FALSE ;
          char *msg = getSystemErrorText();
	  g_critical("Socket connection to pfetch server has failed, error was: %s\n", msg) ;
          g_free(msg);
	}
      else
      {
        char *msg = getSystemErrorText();
	g_error("Fatal error on socket connection to pfetch server, error was: %s\n", msg) ;
        g_free(msg);
      }
    }
  else if (bytes_written != bytes_to_send)
    {
      g_error("send() call should have written %d bytes, but actually wrote %d.\n", bytes_to_send, bytes_written) ;
    }

  /* Reset the old signal handler.                                           */
  if (sigaction(SIGPIPE, &oldsigpipe, NULL) < 0)
    {
      g_error("Cannot reset previous signal handler for signal SIGPIPE for socket write operations.\n") ;
    }

  g_free(tmp) ;

  return status ;
}


#ifdef PFETCH_HTML 

static PFetchStatus pfetch_reader_func(PFetchHandle *handle,
				       char         *text,
				       guint        *actual_read,
				       GError       *error,
				       gpointer      user_data)
{
  PFetchData pfetch_data = (PFetchData)user_data;
  PFetchStatus status    = PFETCH_STATUS_OK;

  if (actual_read && *actual_read > 0 && pfetch_data)
    {
      GtkTextBuffer *text_buffer = GTK_TEXT_BUFFER(pfetch_data->text_buffer);

      /* clear the buffer the first time... */
      if(pfetch_data->got_response == FALSE)
        {
          gtk_text_buffer_set_text(text_buffer, "", 0);
        }

      /* If we tried fetching the full entry and failed, try again
       * with the next fetch method, if there is one */
      if (!strncasecmp(text, "no match", 8))
        {
          fetchSequence(pfetch_data->blxSeq, TRUE, pfetch_data->attempt + 1, pfetch_data->blxWindow, NULL);
        }
      else
        {
          gtk_text_buffer_insert_at_cursor(text_buffer, text, *actual_read);
        }
        
      pfetch_data->got_response = TRUE;
    }

  return status;
}



static PFetchStatus pfetch_closed_func(gpointer user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK;

#ifdef DEBUG_ONLY
  printf("pfetch closed\n");
#endif	/* DEBUG_ONLY */

  return status;
}


/* GtkObject destroy signal handler */
static void handle_dialog_close(GtkWidget *dialog, gpointer user_data)
{
  PFetchData pfetch_data   = (PFetchData)user_data;
  pfetch_data->text_buffer = NULL;
  pfetch_data->widget_destroy_handler_id = 0; /* can we get this more than once? */

  if(pfetch_data->pfetch)
    pfetch_data->pfetch = PFetchHandleDestroy(pfetch_data->pfetch);

  return ;
}



/* 
 *             Multiple sequence fetch callbacks.
 */

static PFetchStatus sequence_pfetch_reader(PFetchHandle *handle,
					   char         *text,
					   guint        *actual_read,
					   GError       *error,
					   gpointer      user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK ;
  PFetchSequence fetch_data = (PFetchSequence)user_data ;
  ProgressBar bar = fetch_data->bar ;

  if (fetch_data->parser_state != PARSING_FINISHED && fetch_data->parser_state != PARSING_CANCELLED)
    {
      if (!(*text) || *actual_read <= 0)
        {
          fetch_data->parser_state = PARSING_FINISHED ;
          fetch_data->status = FALSE ;
          fetch_data->err_txt = g_strdup("No data returned by pfetch http proxy server.") ;
        }
      else if (*actual_read > 0)
        {
          if (isCancelledProgressBar(bar))
            {
              fetch_data->status = FALSE ;
              fetch_data->parser_state = PARSING_CANCELLED ;
      
              status = PFETCH_STATUS_FAILED ;
            }
          else
	    {
	      const int len = min(strlen(text), strlen("not authorised"));
	    
	      if (strncasecmp(text, "not authorised", len) == 0 || strncasecmp(text, "not authorized", len) == 0)
		{
		  fetch_data->parser_state = PARSING_FINISHED ;
		  fetch_data->status = FALSE ;
		  fetch_data->err_txt = g_strdup("Not authorised to access pfetch proxy server.") ;
		}
	      else
		{
		  if (!parsePfetchHtmlBuffer(text, *actual_read, fetch_data))
		    {
		      status = PFETCH_STATUS_FAILED ;
		    }
		}
	    }
        }
    }

  return status ;
}

static PFetchStatus sequence_pfetch_closed(PFetchHandle *handle, gpointer user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK;
  PFetchSequence fetch_data = (PFetchSequence)user_data ;


#ifdef DEBUG_ONLY
  printf("pfetch closed\n");
#endif	/* DEBUG_ONLY */


  fetch_data->parser_state = PARSING_FINISHED ;
  fetch_data->connection_closed = TRUE;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  destroyProgressBar(fetch_data->bar) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

  if (fetch_data->stats)
    {
      printf("Stats for %d reads:\tmin = %d\tmax = %d\tmean = %d\n",
	     fetch_data->total_reads,
	     fetch_data->min_bytes, fetch_data->max_bytes,
	     (fetch_data->total_bytes / fetch_data->total_reads)) ;
    }

  return status ;
}


/* Progress meter has been destroyed so free all pfetch stuff. */
static void sequence_dialog_closed(GtkWidget *dialog, gpointer user_data)
{
  PFetchSequence fetch_data = (PFetchSequence)user_data ;

  fetch_data->pfetch = PFetchHandleDestroy(fetch_data->pfetch) ;

  return ;
}


/* Utility to record the stats for the current read buffer, if stats are enabled */
static void pfetchHtmlRecordStats(const char *read_text, const int length, PFetchSequence fetch_data)
{
  if (fetch_data->stats)
    {
      if (length < fetch_data->min_bytes)
	fetch_data->min_bytes = length ;

      if (length > fetch_data->max_bytes)
	fetch_data->max_bytes = length ;

      fetch_data->total_bytes += length ;

      fetch_data->total_reads++ ;
  
      if (length < 100)
	{
	  char *str ;

	  str = g_strndup(read_text, length) ;
	  printf("Read %d: \"%s\"\n", fetch_data->total_reads, str) ;
	  g_free(str) ;
	}
    }
}


/* Parse the buffer sent back by proxy server. The pfetch server sends back data separated
 * by newlines. The data is either a valid IUPAC dna sequence or a valid IUPAC peptide
 * sequence or the text "no match". The problem here is that the data is not returned to
 * this function in complete lines so we have to reconstruct the lines as best we can. It's
 * even possible for very long sequences that they may span several buffers. Note also
 * that the buffer is _not_ null terminated, we have to use length to know when to stop reading.
 */
static gboolean parsePfetchHtmlBuffer(char *read_text, int length, PFetchSequence fetch_data)
{
  gboolean status = TRUE ;

  /* Validate input and record stats, if requested */
  g_assert(fetch_data && fetch_data->currentSeqItem && fetch_data->currentSeqItem->data);
  pfetchHtmlRecordStats(read_text, length, fetch_data);
  
  GError *error = NULL;
  
  if (fetch_data->fetchMethod->output == BLXFETCH_OUTPUT_EMBL)
    {
      /* We're fetching the full EMBL entries */
      pfetchParseEmblFileBuffer(read_text, 
                                length,
                                &fetch_data->currentSeq, 
                                &fetch_data->currentSeqItem, 
                                fetch_data->bar,
                                fetch_data->seq_total,
                                &fetch_data->num_fetched,
                                &fetch_data->num_succeeded, 
                                fetch_data->section_id,
                                fetch_data->tag_name,
                                &fetch_data->found_end_quote,
                                fetch_data->seq_type,
                                &fetch_data->parser_state,
                                &fetch_data->status);
    }
  else if (fetch_data->fetchMethod->output == BLXFETCH_OUTPUT_FASTA)
    {
      /* The fetched entries just contain the FASTA sequence */
      pfetchParseSequenceFileBuffer(read_text, 
                                    length, 
                                    &fetch_data->currentSeq, 
                                    &fetch_data->currentSeqItem, 
                                    fetch_data->bar, 
                                    fetch_data->seq_total, 
                                    &fetch_data->num_fetched, 
                                    &fetch_data->num_succeeded, 
                                    fetch_data->seq_type,
                                    &fetch_data->parser_state, 
                                    &fetch_data->status,
                                    &error);
    }
  else 
    {
      g_set_error(&error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT, 
                  "Invalid output format specified for fetch method '%s'; expected '%s' or '%s'\n",
                  g_quark_to_string(fetch_data->fetchMethod->name), FETCH_OUTPUT_EMBL, FETCH_OUTPUT_FASTA);
    }
  
  if (error)
    {
      fetch_data->err_txt = error->message;
      g_error_free(error);
      error = NULL;
    }
  
  status = fetch_data->status ;
  
  return status;
}



#endif /* PFETCH_HTML */



/* Functions to display, update, cancel and remove a progress meter. */

static ProgressBar makeProgressBar(int seq_total)
{
  ProgressBar bar = g_new0(ProgressBarStruct, 1) ;

  bar->seq_total = seq_total ;
  bar->cancelled = FALSE;
  bar->widget_destroy_handler_id = 0;
  
  gdk_color_parse("blue", &(bar->blue_bar_fg)) ;
  gdk_color_parse("red", &(bar->red_bar_fg)) ;

  bar->top_level = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  char *title = g_strdup_printf("Blixem - pfetching %d sequences...", seq_total) ;
  gtk_window_set_title(GTK_WINDOW(bar->top_level), title) ;
  g_free(title) ;
  g_signal_connect(G_OBJECT(bar->top_level), "destroy", G_CALLBACK(destroyProgressCB), bar) ;


  gtk_window_set_default_size(GTK_WINDOW(bar->top_level), 350, -1) ;

  GtkWidget *frame = gtk_frame_new(NULL) ;
  gtk_container_add(GTK_CONTAINER(bar->top_level), frame) ;

  GtkWidget *vbox = gtk_vbox_new(FALSE, 5) ;
  gtk_container_add(GTK_CONTAINER(frame), vbox) ;

  bar->progress = gtk_progress_bar_new() ; 

  gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->blue_bar_fg)) ;

  gtk_box_pack_start(GTK_BOX(vbox), bar->progress, TRUE, TRUE, 0);

  bar->label = gtk_label_new("") ;
  gtk_box_pack_start(GTK_BOX(vbox), bar->label, TRUE, TRUE, 0);

  GtkWidget *hbox = gtk_hbox_new(FALSE, 0) ;
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
  gtk_container_border_width(GTK_CONTAINER(hbox), 5);

  GtkWidget *cancel_button = gtk_button_new_with_label("Cancel") ;
  gtk_box_pack_start(GTK_BOX(hbox), cancel_button, TRUE, TRUE, 0) ;

  gtk_signal_connect(GTK_OBJECT(cancel_button), "clicked", GTK_SIGNAL_FUNC(cancelCB), (gpointer)bar) ;

  gtk_widget_show_all(bar->top_level) ;
  
  while (gtk_events_pending())
    {
      gtk_main_iteration() ;
    }

  return bar ;
}

static void updateProgressBar(ProgressBar bar, const char *sequence, int numFetched, gboolean fetch_ok)
{
  char *label_text ;

  if (fetch_ok)
    gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->blue_bar_fg)) ;
  else
    gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->red_bar_fg)) ;

  label_text = g_strdup_printf("%s - %s", sequence, fetch_ok ? "Fetched" : "Not found") ;
  gtk_label_set_text(GTK_LABEL(bar->label), label_text) ;
  g_free(label_text) ;
  
  gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(bar->progress),
				(double)((double)numFetched / (double)(bar->seq_total))) ;

  while (gtk_events_pending())
    gtk_main_iteration();

  return ;
}


static gboolean isCancelledProgressBar(ProgressBar bar)
{
  return bar->cancelled ;
}


static void destroyProgressBar(ProgressBar bar)
{
  gtk_widget_destroy(bar->top_level) ;
  
  return ;
}


/* called from progress meter when being destroyed, do not call directly. */
static void destroyProgressCB(GtkWidget *widget, gpointer cb_data)
{
  ProgressBar bar = (ProgressBar)cb_data ;

  g_free(bar) ;

  return ;
}


/* Called from progress meter. */
static void cancelCB(GtkWidget *widget, gpointer cb_data)
{
  ProgressBar bar = (ProgressBar)cb_data ;

  bar->cancelled = TRUE ;

  return ;
}





/* 
 *            Configuration file reading.
 */

/* Get default blixem options from the "blixem" stanza */
static void readBlixemConfigGroup(GKeyFile *key_file,
                                  const char *group, 
                                  CommandLineOptions *options,
                                  GError **error)
{
  /* Get the comma-separated list of possible fetch methods */
  options->bulkFetchDefault = keyFileGetCsv(key_file, group, SEQTOOLS_BULK_FETCH);
  options->userFetchDefault = keyFileGetCsv(key_file, group, SEQTOOLS_USER_FETCH);
}


/* Create a fetch method struct */
static BlxFetchMethod* createBlxFetchMethod(const char *fetchName, 
                                            const GQuark fetchMode)
{
  BlxFetchMethod *result = g_malloc(sizeof *result);

  result->name = g_quark_from_string(fetchName);
  result->mode = BLXFETCH_MODE_NONE;

  result->location = NULL;
  result->node = NULL;
  result->port = 0;
  result->cookie_jar = NULL;
  result->args = NULL;
  result->output = 0;

  return result;
}

static BlxFetchOutput getFetchOutputFormat(GKeyFile *key_file, const char *group)
{
  BlxFetchOutput result = BLXFETCH_OUTPUT_INVALID;
  
  char *output = g_key_file_get_string(key_file, group, FETCH_OUTPUT, NULL);

  if (stringsEqual(output, FETCH_OUTPUT_FASTA, FALSE))
    result = BLXFETCH_OUTPUT_FASTA;
  else if (stringsEqual(output, FETCH_OUTPUT_EMBL, FALSE))
    result = BLXFETCH_OUTPUT_EMBL;
  if (stringsEqual(output, FETCH_OUTPUT_GFF, FALSE))
    result = BLXFETCH_OUTPUT_GFF;
  
  g_free(output);

  return result;
}


/* Get details about the given fetch method stanza and add it to 
 * the list of fetch methods in 'options' */
static void readFetchConfigGroup(GKeyFile *key_file, 
                                 const char *group, 
                                 const char *fetchMode,
                                 CommandLineOptions *options, 
                                 GError **error)
{
  BlxFetchMethod *result = createBlxFetchMethod(group, g_quark_from_string(fetchMode));

  /* Set the relevant properties for this type of fetch method */
#ifdef PFETCH_HTML
  if (stringsEqual(fetchMode, FETCH_MODE_HTTP_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_HTTP;
      result->location = g_key_file_get_string(key_file, group, HTTP_FETCH_LOCATION, NULL);
      result->port = g_key_file_get_integer(key_file, group, HTTP_FETCH_PORT, NULL);
      result->cookie_jar = g_key_file_get_string(key_file, group, HTTP_FETCH_COOKIE_JAR, NULL);
      result->args = g_key_file_get_string(key_file, group, FETCH_ARGS, NULL);
      result->output = getFetchOutputFormat(key_file, group);
    }
  else if (stringsEqual(fetchMode, FETCH_MODE_PIPE_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_PIPE;
      result->location = g_key_file_get_string(key_file, group, HTTP_FETCH_LOCATION, NULL);
      result->port = g_key_file_get_integer(key_file, group, HTTP_FETCH_PORT, NULL);
      result->cookie_jar = g_key_file_get_string(key_file, group, HTTP_FETCH_COOKIE_JAR, NULL);
      result->args = g_key_file_get_string(key_file, group, FETCH_ARGS, NULL);
      result->output = getFetchOutputFormat(key_file, group);
    }
#endif
  if (stringsEqual(fetchMode, FETCH_MODE_SOCKET_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_SOCKET;
      result->location = g_key_file_get_string(key_file, group, SOCKET_FETCH_LOCATION, NULL);
      result->node = g_key_file_get_string(key_file, group, SOCKET_FETCH_NODE, NULL);
      result->port = g_key_file_get_integer(key_file, group, SOCKET_FETCH_PORT, NULL);
      result->args = g_key_file_get_string(key_file, group, FETCH_ARGS, NULL);
      result->output = getFetchOutputFormat(key_file, group);
    }
  else if (stringsEqual(fetchMode, FETCH_MODE_WWW_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_WWW;
      result->location = g_key_file_get_string(key_file, group, WWW_FETCH_LOCATION, NULL);
    }
  else if (stringsEqual(fetchMode, FETCH_MODE_DB_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_DB;
      /* to do: not implemented */
    }
  else if (stringsEqual(fetchMode, FETCH_MODE_COMMAND_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_COMMAND;
      result->location = g_key_file_get_string(key_file, group, COMMAND_FETCH_SCRIPT, NULL);
      result->args = g_key_file_get_string(key_file, group, FETCH_ARGS, NULL);
      result->output = getFetchOutputFormat(key_file, group);
    }
  else if (stringsEqual(fetchMode, FETCH_MODE_NONE_STR, FALSE))
    {
      result->mode = BLXFETCH_MODE_NONE;
    }
  else
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_FETCH_MODE, "Unrecognised fetch mode '%s'\n", fetchMode);
      result->mode = BLXFETCH_MODE_NONE;
    }

  /* Add to list */
  GQuark fetchName = g_quark_from_string(group);
  g_hash_table_insert(options->fetchMethods, GINT_TO_POINTER(fetchName), result);
}


/* Tries to load values for this group if it is a standard, recognised group or 
 * group type. If any required values are missing, sets the error. */
static void loadConfig(GKeyFile *key_file, const char *group, CommandLineOptions *options, GError **error)
{
  GError *tmpError = NULL;
  
  /* Check for known group names first */
  if (stringsEqual(group, BLIXEM_GROUP, FALSE))
    {
      readBlixemConfigGroup(key_file, group, options, error);
    }
  else
    {
      /* Check if this is a fetch method; if it is it the group will have a fetch-mode key */
      char *fetchMode = g_key_file_get_string(key_file, group, FETCH_MODE_KEY, &tmpError);
      
      if (fetchMode && !tmpError)
        readFetchConfigGroup(key_file, group, fetchMode, options, error);
    }
}


/* Read standard stanzas from the given config file. Sets the error if any problems */
static void readConfigFile(GKeyFile *key_file, char *config_file, CommandLineOptions *options, GError **error)
{
  GKeyFileFlags flags = G_KEY_FILE_NONE ;
  
  if (g_key_file_load_from_file(key_file, config_file, flags, error))
    {
      gsize num_groups ;

      char **groups = g_key_file_get_groups(key_file, (gsize*)(&num_groups)) ;

      char **group = groups;
      int i = 0;
      
      for ( ; i < num_groups ; i++, group++)
	{
          GError *tmpError = NULL;
          loadConfig(key_file, *group, options, &tmpError) ;
          
          if (tmpError && error)
            {
              /* Compile all errors into one message */
              if (*error == NULL)
                {
                  g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_WARNINGS,
                              "Errors found while reading config file '%s':\n", config_file);
                }

              postfixError(*error, "  Error reading [%s] stanza: %s", *group, tmpError->message);
            }
	}

      g_strfreev(groups);
    }
}


/*******************************************************************
 *                      Pfetch utility functions                   *
 *******************************************************************/


/* Initialise a pfetch connection with the given options and calls pfetch on each of the 
 * sequences in the given list. */
static gboolean socketFetchInit(char *fetchOptions, 
                                GList *seqsToFetch, 
                                const char *node, 
                                int port, 
                                gboolean External, 
                                int *sock)
{
  gboolean status = TRUE;
  
  /* open socket connection */
  if (status)
    {
      *sock = socketConstruct (node, port, External) ;
      if (*sock < 0)			/* we can't connect to the server */
	{
	  status = FALSE ;
	}
    }
  
  /* send the command/names to the server */
  if (status)
    {
      /* Send '-q' to get back one line per sequence, '-C' for lowercase DNA and uppercase protein. */
      status = socketSend (*sock, fetchOptions);
    }
  
  if (status)
    {
      /* For each sequence, send a command to fetch that sequence, in the order that they are in our list */
      GList *seqItem = seqsToFetch;
      
      for ( ; seqItem && status ; seqItem = seqItem->next)
	{
          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
	  status = socketSend(*sock, blxSequenceGetFullName(blxSeq));
	}
    }
  
  if (status)
    {
      /* send a final newline to flush the socket */
      if (send(*sock, "\n", 1, 0) != 1)
	{
	  g_critical("failed to send final \\n to socket\n") ;
	  status = FALSE ;
	}
    }
  
  return status;
}


/* Check the progress bar to see if the user cancelled. If so, sets the state to cancelled
 * and the status to false (meaning don't continue). */
static void checkProgressBar(ProgressBar bar, BlxEmblParserState *parserState, gboolean *status)
{
  if (isCancelledProgressBar(bar))
    {
      *status = FALSE;
      *parserState = PARSING_CANCELLED;
    }
}


/* Receive the next buffer back from the pfetch server. Only does anything if the status is ok
 * (i.e. true). Checks the length etc and sets the state to finished if there is no more data
 * to receive. Sets 'status' to false if there was an error. Returns the number of chars received. */
static int pfetchReceiveBuffer(char *buffer, const int bufferSize, const int sock, BlxEmblParserState *parserState, gboolean *status)
{
  int lenReceived = 0;
  
  if (*status == FALSE)
    {
      return lenReceived;
    }
  
  /* Ask for the next chunk of data to be put into our buffer */
  lenReceived = recv(sock, buffer, bufferSize, 0);
  
  if (lenReceived < 0)
    {
      /* Problem with this one - skip to the next */
      status = FALSE;
      char *msg = getSystemErrorText();
      g_critical("Could not retrieve sequence data from pfetch server, error was: %s\n", msg) ;
      g_free(msg);
    }
  else if (lenReceived == 0)
    {
      /* No more data, so quit out of the loop */
      *parserState = PARSING_FINISHED;
    }
  
  return lenReceived;
}


/* Update the given sequence and sequence item to the next item in the list (unless the state is
 * 'finished' in which case there is nothing to do). Sets status to false if there was an
 * error (i.e. if we ran out of sequences). */
static void pfetchGetNextSequence(BlxSequence **currentSeq, 
                                  GList **currentSeqItem, 
                                  ProgressBar bar, 
                                  const int numRequested,
                                  const int numFetched, 
                                  const gboolean pfetch_ok, 
                                  gboolean *status, 
                                  BlxEmblParserState *parserState)
{
  if (*parserState != PARSING_FINISHED)
    {
      /* Check that another BlxSequence item exists in the list */
      if ((*currentSeqItem)->next)
        {
          updateProgressBar(bar, blxSequenceGetFullName(*currentSeq), numFetched, pfetch_ok) ;
          
          /* Move to the next BlxSequence */
          *currentSeqItem = (*currentSeqItem)->next;
          *currentSeq = (BlxSequence*)((*currentSeqItem)->data);
        }
      else
        {
          *status = FALSE ;
          g_critical("Unexpected data from pfetch server: received too many lines. (%d sequences were requested.)\n", numRequested);
        }
    }
}


/* Parse the given buffer that contains a section of data from returned sequence file(s) */
static void pfetchParseSequenceFileBuffer(const char *buffer,
                                          const int lenReceived, 
                                          BlxSequence **currentSeq, 
                                          GList **currentSeqItem,
                                          ProgressBar bar,
                                          const int numRequested,
                                          int *numFetched,
                                          int *numSucceeded,
                                          const BlxSeqType seqType,
                                          BlxEmblParserState *parserState, 
                                          gboolean *status,
                                          GError **error)
{
  if (*status == FALSE || *parserState == PARSING_FINISHED || *parserState == PARSING_CANCELLED)
    {
      return;
    }
  
  /* Loop through each character in the buffer */
  int i = 0;
  
  for ( ; i < lenReceived && *status; ++i)
    {
      /* Check for user cancellation again */
      checkProgressBar(bar, parserState, status);
      
      if (*parserState == PARSING_CANCELLED)
        {
          break;
        }
      
      if (buffer[i] == '\n')
        {
          /* finish up this sequence and move to the next one */
          gboolean pfetch_ok = pfetchFinishSequence(*currentSeq, seqType, numRequested, numFetched, numSucceeded, parserState);
          pfetchGetNextSequence(currentSeq, currentSeqItem, bar, numRequested, *numFetched, pfetch_ok, status, parserState);
        }
      else
        {
          /* First time round we need to create the sequence string */
          if ((*currentSeq)->sequence == NULL)
            {
              (*currentSeq)->sequence = g_string_new("");
            }
          
          /* Append this character to the sequence string */
          g_string_append_c((*currentSeq)->sequence, buffer[i]);
        }
    }
}


/* Get the parser state from the given two-letter at at the start of an EMBL file line.
 * Returns true if we need to move to the next sequence. */
static gboolean pfetchGetParserStateFromId(const char *sectionId, 
                                           BlxSequence *currentSeq, 
                                           GString *tagName,
                                           BlxEmblParserState *parserState)
{
  gboolean finishSequence = FALSE;
  
  if (stringsEqual(sectionId, "SQ", TRUE))
    {
      /* Skip the sequence section if the sequence is already populated or not required. */
      if (currentSeq->sequence && (currentSeq->sequence->str || !blxSequenceRequiresSeqData(currentSeq)))
        {
          *parserState = PARSING_IGNORE;
        }
      else
        {
          *parserState = PARSING_SEQUENCE_HEADER;
        }
    }
  else if (stringsEqual(sectionId, "OS", TRUE) && blxSequenceRequiresOptionalData(currentSeq))
    {
      *parserState = PARSING_ORGANISM;
    }
  else if (stringsEqual(sectionId, "GN", TRUE) && blxSequenceRequiresOptionalData(currentSeq))
    {
      *parserState = PARSING_GENE_NAME;
    }
  else if (stringsEqual(sectionId, "FT", TRUE) && blxSequenceRequiresOptionalData(currentSeq))
    {
      if (tagName && tagName->len)
        {
          /* If the tag name is already set we must be in the middle of parsing a multi-line tag, 
           * so continue with that tag */
          pfetchGetParserStateFromTagName(tagName, parserState);
        }
      else
        {
          /* We'll parse the general FT section looking for a tag name */
          *parserState = PARSING_FEATURE_TYPE;
        }
    }
  else if (stringsEqual(sectionId, "//", TRUE))
    {
      /* This indicates the end of the embl file, so we've finished the current sequence. Set
       * the state to NA so that if any more chars exist they are ignored till we find the next
       * newline. */
      finishSequence = TRUE;
      *parserState = PARSING_IGNORE;
    }
  else if (stringsEqual(sectionId, "no", TRUE))
    {
      /* This must be the start of the text "no match" because otherwise all new lines should
       * start with an uppercase 2-letter ID, whitespace or "//". We can't do anything with the
       * current sequence so finish up. Set the state to something that will be ignored till
       * the next newline. */
      finishSequence = TRUE;
      *parserState = PARSING_IGNORE;
    }
  else
    {
      /* Any other section we're not interested in */
      *parserState = PARSING_IGNORE;
    }
  
  return finishSequence;
}


/* Process the current character in the current buffer of EMBL data. The char might be read into
 * various data locations depending on the current parser state, or might cause the parser state
 * to change. */
static void pfetchProcessEmblBufferChar(const char curChar, 
                                        const int numRequested,
                                        int *numFetched,
                                        int *numSucceeded,
                                        ProgressBar bar,
                                        BlxSequence **currentSeq,
                                        GList **currentSeqItem,
                                        char *sectionId,
                                        GString *tagName,
                                        gboolean *foundEndQuote,
                                        const BlxSeqType seqType,
                                        BlxEmblParserState *parserState, 
                                        gboolean *status)
{
  switch (*parserState)
    {
      case PARSING_NEWLINE:
      {
        /* First char should be the start of the two-letter ID */
        sectionId[0] = curChar;
        *parserState = PARSING_ID;
        break;
      }
        
      case PARSING_ID:
      {
        /* Only one more char to read from the ID. */
        sectionId[1] = curChar;
        
        if (pfetchGetParserStateFromId(sectionId, *currentSeq, tagName, parserState))
          {
            /* finish the current seq and move to the next */
            pfetchFinishEmblFile(currentSeq, currentSeqItem, numRequested, numFetched, numSucceeded, bar, seqType, parserState, status);
          }
        
        break;
      }
        
      case PARSING_SEQUENCE:
      {
        /* Look out for the terminating "//" characters. We ignore newlines in the sequence
         * body so won't catch these by the normal method. We can assume that if we see a single
         * slash that it's the terminator, because there should never see one here otherwise */
        if (curChar == '/')
          {
            sectionId[0] = curChar;
            *parserState = PARSING_ID; /* continue parsing to read the other '/' */
          }
        else
          {
            /* Add the character to the end of the sequence, but only if it is a valid IUPAC character
             * (because we want to ignore whitespace, newlines and digits). */
            if (isValidIupacChar(curChar, seqType))
              {
                appendCharToString(curChar, &((*currentSeq)->sequence));
              }
          }
        
        break;
      }
        
      case PARSING_ORGANISM:
      {
        appendCharToString(curChar, &((*currentSeq)->organism));
        break;
      }
        
      case PARSING_GENE_NAME:
      {
        appendCharToString(curChar, &((*currentSeq)->geneName));
        break;
      }
        
      case PARSING_FEATURE_TYPE:
      {
        /* If we find a forward slash, it's the start of a tag, so we'll parse the tag name next */
        if (curChar == '/')
          {
            *parserState = PARSING_FT_TAG_NAME;
          }
        
        break;
      }
        
      case PARSING_FT_TAG_NAME:
      {
        if (curChar == '=') /* signals the end of the tag */
          {
            pfetchGetParserStateFromTagName(tagName, parserState);
          }
        else
          {
            g_string_append_c(tagName, curChar); /* read in tag name */
          }
        
        break;
      }
        
      case PARSING_FT_UNQUOTED_SECTION:
      {
        /* Look for the start quote before we start parsing the tag properly */
        if (curChar == '"')
          {
            pfetchGetParserStateFromTagName(tagName, parserState);
          }
        
        break;
      }
        
      case PARSING_FT_TISSUE_TYPE:
      {
        appendCharToQuotedString(curChar, foundEndQuote, &((*currentSeq)->tissueType));
        break;
      }
        
      case PARSING_FT_STRAIN:
      {
        appendCharToQuotedString(curChar, foundEndQuote, &((*currentSeq)->strain));
        break;
      }
        
      default:
        break;
        
    };
}


/* Parse the given buffer that contains a section of data from returned full EMBL file(s) */
static void pfetchParseEmblFileBuffer(const char *buffer,
                                      const int lenReceived,
                                      BlxSequence **currentSeq, 
                                      GList **currentSeqItem,
                                      ProgressBar bar,
                                      const int numRequested,
                                      int *numFetched,
                                      int *numSucceeded,
                                      char *sectionId,
                                      GString *tagName,
                                      gboolean *foundEndQuote,
                                      const BlxSeqType seqType,
                                      BlxEmblParserState *parserState, 
                                      gboolean *status)
{
  if (*status == FALSE || *parserState == PARSING_FINISHED || *parserState == PARSING_CANCELLED)
    {
      return;
    }

  /* Loop through each character in the buffer */
  int i = 0;
  
  for ( ; i < lenReceived && *status; ++i)
    {
      checkProgressBar(bar, parserState, status);

      if (*parserState == PARSING_CANCELLED)
        {
          break;
        }

      const char curChar = buffer[i];

      /* Special treatment if we've previously found an end quote: if this char is NOT also
       * a quote, then it means it genuinely was an end quote, so we finish the FT tag that
       * we were reading (by setting the tag name to null). We return to parsing the general FT
       * section in case there are any more tags. */
      if (*foundEndQuote && curChar != '"')
        {
          g_string_truncate(tagName, 0);
          *foundEndQuote = FALSE;
          *parserState = PARSING_FEATURE_TYPE;
        }
      

      /* Newline character means we need to re-read the 2-letter ID at the start of the
       * line to find out what section we're in. Ignore newlines in the sequence body, though,
       * because it spans several lines and does not have an ID on each line. */
      if (curChar == '\n' && *parserState != PARSING_SEQUENCE)
        {
          /* If we were previously in the sequence header, the next line is the sequence body */
          if (*parserState == PARSING_SEQUENCE_HEADER)
            {
              *parserState = PARSING_SEQUENCE;
            }
          else
            {
              *parserState = PARSING_NEWLINE;
            }
        }
      else
        {
          pfetchProcessEmblBufferChar(curChar, 
                                      numRequested, 
                                      numFetched, 
                                      numSucceeded, 
                                      bar, 
                                      currentSeq, 
                                      currentSeqItem, 
                                      sectionId, 
                                      tagName, 
                                      foundEndQuote,
                                      seqType,
                                      parserState, 
                                      status);
        }
    }
}  


/* Check the given tag name to see what type of section it indicates, and update
 * the parser state accordingly. Sets the parser state to PARSING_IGNORE if this
 * is not a tag we care about. Returns true if we've entered a free text section
 * where newline characters should be ignored. */
static void pfetchGetParserStateFromTagName(GString *tagName, BlxEmblParserState *parserState)
{
  if (stringsEqual(tagName->str, "tissue_type", TRUE))
    {
      if (*parserState == PARSING_FT_TAG_NAME)
        {
          *parserState = PARSING_FT_UNQUOTED_SECTION; /* look for quoted section before we start parsing the actual data */
        }
      else
        {
          *parserState = PARSING_FT_TISSUE_TYPE;
        }
    }
  else if (stringsEqual(tagName->str, "strain", TRUE))
    {
      if (*parserState == PARSING_FT_TAG_NAME)
        {
          *parserState = PARSING_FT_UNQUOTED_SECTION; /* look for quoted section before we start parsing the actual data */
        }
      else
        {
          *parserState = PARSING_FT_STRAIN;
        }
    }
  else
    {
      /* Ignore this tag. Reset the tag name to zero-length so we know to start again. */
      *parserState = PARSING_IGNORE;
      g_string_truncate(tagName, 0);
    }
}


/* Append the given char to the given GString, but only if it is not an end quote. Assumes the
 * start quote has already been seen. Treats double quotes ("") as a single escaped quote and
 * includes it in the text - foundEndQuote is set for the first quote found and if the next char
 * this is called with is also a quote the quote is included and foundEndQuote reset.
 * To do: this should really identify newlines and ignore the FT identifier and whitespace at the
 * start of the line. Since none of our tags currently have multi-line data that's not worth doing just yet... */
static void appendCharToQuotedString(const char curChar, gboolean *foundEndQuote, GString **result)
{
  if (curChar == '"')
    {
      if (*foundEndQuote)
        {
          /* The previous char was also a quote, which means this is an 
           * escaped quote, so we do include the char in the text. */
          appendCharToString(curChar, result);
          *foundEndQuote = FALSE;
        }
      else
        {
          /* This looks like an end quote but we won't know until the next loop, so just flag it. */
          *foundEndQuote = TRUE;
        }
    }
  else
    {
      appendCharToString(curChar, result);
    }
}


/* Append the given char to the given GString. Create the GString if it is null. Doesn't
 * add whitespace chars to the start of the string. */
static void appendCharToString(const char curChar, GString **result)
{
  /* First time round we will need to create the string for the sequence (but don't bother
   * until we get to non-whitespace characters) */
  if (*result == NULL && !isWhitespaceChar(curChar))
    {
      *result = g_string_new("");
    }
  
  if (*result)
    {
      /* Don't add multiple consecutive whitespace characters */
      const int len = (*result)->len;
      const char *str = (*result)->str;

      if (!isWhitespaceChar(curChar) || len < 0 || !isWhitespaceChar(str[len - 1]))
        {
          g_string_append_c(*result, curChar);
        }
    }
}


/* This is called when we've finished parsing a given sequence. It checks that the sequence
 * data is valid (i.e. not 'no match') and complements it if necessary. It updates the parser
 * state to finished if we've got all the sequences we requested. Returns true if the pfetch
 * was successful, false if not */
static gboolean pfetchFinishSequence(BlxSequence *currentSeq, const BlxSeqType seqType, const int numRequested, int *numFetched, int *numSucceeded, BlxEmblParserState *parserState)
{
  *numFetched += 1;
  
  if (*numFetched >= numRequested)
    {
      /* We've fetched all of the sequences we requested. */
      *parserState = PARSING_FINISHED;
    }
  
  /* The pfetch failed if our sequence is null or equal to "no match". */
  gboolean pfetch_ok = FALSE;
  
  if (currentSeq && currentSeq->sequence && currentSeq->sequence->str)
    { 
      int len = min(strlen(currentSeq->sequence->str), 8);
      if (strncasecmp(currentSeq->sequence->str, "no match", len))
	{
	  pfetch_ok = TRUE;
	}
    }
  
  if (pfetch_ok)
    {
      *numSucceeded += 1;
    }
  else if (currentSeq && currentSeq->sequence)
    {
      /* Failed. Clear the sequence string. */
      g_string_free(currentSeq->sequence, TRUE);
      currentSeq->sequence = NULL;
    }
  
  return pfetch_ok;
}


/* This is called when we've finished parsing the embl file data for the given BlxSequence. It
 * checks the data returned was valid and complements the sequence if necessary. It moves to
 * the next sequence. */
static void pfetchFinishEmblFile(BlxSequence **currentSeq, 
                                 GList **currentSeqItem,
                                 const int numRequested, 
                                 int *numFetched, 
                                 int *numSucceeded, 
                                 ProgressBar bar,
                                 const BlxSeqType seqType,
                                 BlxEmblParserState *parserState,
                                 gboolean *status)
{
  gboolean pfetch_ok = pfetchFinishSequence(*currentSeq, seqType, numRequested, numFetched, numSucceeded, parserState);
  
  pfetchGetNextSequence(currentSeq, currentSeqItem, bar, numRequested, *numFetched, pfetch_ok, status, parserState);
}




