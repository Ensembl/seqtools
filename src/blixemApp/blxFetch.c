/*  File: blxFetch.c
 *  Author: Ed Griffiths, 2008-06-17
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

#define DEFAULT_PFETCH_WINDOW_WIDTH_CHARS	      85


/* Blixem config error domain */
#define BLX_CONFIG_ERROR g_quark_from_string("Blixem config")

/* Error codes */
typedef enum
  {
    BLX_CONFIG_ERROR_NO_GROUPS,             /* no groups in config file */
    BLX_CONFIG_ERROR_INVALID_KEY_TYPE,      /* invalid key type given in config file */
    BLX_CONFIG_ERROR_MISSING_KEY,           /* a required key is missing */
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

/* 
 * config file structs/keywords
 */

typedef enum {KEY_TYPE_INVALID, KEY_TYPE_BOOL, KEY_TYPE_INT, KEY_TYPE_DOUBLE, KEY_TYPE_STRING} ConfigKeyType ;

typedef struct
{
  char *name ;
  ConfigKeyType type ;
} ConfigKeyValueStruct, *ConfigKeyValue ;

typedef struct
{
  char *name ;
  ConfigKeyValue key_values ;
} ConfigGroupStruct, *ConfigGroup ;



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
#define PFETCH_TITLE_FORMAT_F "blixem - pfetch -F \"%s\""
#define PFETCH_TITLE_FORMAT "blixem - pfetch \"%s\""


typedef struct
{
  GtkWidget *blxWindow;
  GtkWidget *dialog;
  GtkTextBuffer *text_buffer;
  char *title;
  char *sequence_name;

  gulong widget_destroy_handler_id;
  PFetchHandle pfetch;
  gboolean got_response;
  gboolean requested_full_entry;
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
  gboolean full_entry;                                      /* Whether fetching the full EMBL entry or just the sequence data */
  
  ProgressBar bar ;                                         /* Provides graphical feedback about how many sequences have been fetched */
  PFetchHandle pfetch ;
  
  /* Additional fields for parsing full EMBL entries */
  char section_id[3];                                       /* 2-letter ID at the start of the line indicating the current section */
  gboolean found_end_quote;                                 /* Set to true when we're parsing a quoted section and we come across what looks like an end quote */
  GString *tag_name;                                        /* When parsing the FT section, the current tag name is stored this field */
} PFetchSequenceStruct, *PFetchSequence ;


typedef struct
{
  char    *location;
  char    *cookie_jar;
  char    *mode;
  int      port;
} PFetchUserPrefsStruct ;
#endif


#ifdef PFETCH_HTML 
static gboolean getPFetchUserPrefs(PFetchUserPrefsStruct *pfetch) ;
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
static void pfetchHttpEntry(char *sequence_name, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer *text_buffer, const gboolean fetchFullEntry);
#endif

static int socketConstruct(const char *ipAddress, int port, gboolean External) ;
static gboolean socketSend(int sock, char *text) ;

static ProgressBar makeProgressBar(int seq_total) ;
static void updateProgressBar(ProgressBar bar, char *sequence, int numFetched, gboolean fetch_ok) ;
static gboolean isCancelledProgressBar(ProgressBar bar) ;
static void destroyProgressBar(ProgressBar bar) ;
static void destroyProgressCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */
static void cancelCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */

static gboolean readConfigFile(GKeyFile *key_file, char *config_file, GError **error) ;
ConfigGroup getConfig(char *config_name) ;
static gboolean loadConfig(GKeyFile *key_file, ConfigGroup group, GError **error) ;

static char *blxConfigGetDefaultFetchMode(const gboolean bulk);
static void pfetchEntry(char *seqName, GtkWidget *blxWindow, const gboolean displayResults, GString **result_out);

/* Pfetch local functions */
static void                     appendCharToString(const char curChar, GString **result);
static void                     appendCharToQuotedString(const char curChar, gboolean *foundEndQuote, GString **result);

static gboolean                 pfetchInit(char *pfetchOptions, GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, int *sock);
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


/* Some local globals.... */
static const char *URL = NULL ;


/* global configuration object for blixem. */
static GKeyFile *blx_config_G = NULL ;


/* Execute the given external command and return the output from the command as a GString. 
 * The result should be free'd with g_string_free. */
static GString* getExternalCommandOutput(const char *command)
{
  GString *resultText = g_string_new(NULL) ;

  char lineText[MAXLINE+1];
  
  FILE *pipe = popen (command, "r") ;
  
  while (!feof (pipe))
    { 
      if (!fgets (lineText, MAXLINE, pipe))
        {
          break;
        }
      
      int len = strlen(lineText);
      
      if (len > 0)
	{ 
	  if (lineText[len-1] == '\n') 
            {
              lineText[len-1] = '\0';
            }
          
          g_string_append_printf(resultText, "%s\n", lineText) ;
	}
    }
  
  pclose (pipe);
  
  return resultText;
}


/* Display a message dialog showing the given display text. This utility functions sets
 * things like the font and default width based on properties of the main blixem window.
 * Returns a pointer to the dialog, and optionally sets a pointer to the text buffer in the
 * textBuffer return argument. */
static GtkWidget* displayFetchResults(const char *title, const char *displayText, GtkWidget *blxWindow, GtkTextBuffer **textBuffer)
{
  /* Use the same fixed-width font that the detail view uses, but* don't use the
   * detail view's font size, because it may be zoomed in/out. */
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  PangoFontDescription *fontDesc = pango_font_description_copy(detailViewGetFontDesc(detailView));
  pango_font_description_set_size(fontDesc, pango_font_description_get_size(blxWindow->style->font_desc));

  /* Set the initial width based on the default number of characters wide */
  PangoContext *context = gtk_widget_get_pango_context(detailView);
  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
  gdouble charWidth = pango_font_metrics_get_approximate_digit_width(metrics) / PANGO_SCALE;

  const int initWidth = DEFAULT_PFETCH_WINDOW_WIDTH_CHARS * charWidth;
  const int maxHeight = blxWindow->allocation.height * 0.5;

  GtkTextView *textView = NULL;
  GtkWidget *result = showMessageDialog(title, displayText, NULL, initWidth, maxHeight, FALSE, FALSE, fontDesc, &textView);
  
  if (textBuffer && textView)
    {
      *textBuffer = gtk_text_view_get_buffer(textView);
    }

  /* Clean up */
  pango_font_metrics_unref(metrics);
  pango_font_description_free(fontDesc);
  
  return result;
}


/* SHOULD BE MERGED INTO libfree.a */
/* call an external shell command and print output in a text_scroll window
 *
 * This is a replacement for the old graph based text window, it has the advantage
 * that it uses gtk directly and provides cut/paste/scrolling but...it has the
 * disadvantage that it will use more memory as it collects all the output into
 * one string and then this is _copied_ into the text widget.
 * 
 * If this proves to be a problem I expect there is a way to feed the text to the
 * text widget a line a time. */
static void externalCommand (char *command, GtkWidget *blxWindow)
{
#if !defined(MACINTOSH)

  GString *resultText = getExternalCommandOutput(command);
  char *title = blxprintf("Blixem - %s", command);
  displayFetchResults(title, resultText->str, blxWindow, NULL);
  
  g_free(title);
  g_string_free(resultText, TRUE);

#endif
  return ;
}


/* Display the embl entry for a sequence via pfetch, efetch or whatever. */
void fetchAndDisplaySequence(char *seqName, GtkWidget *blxWindow)
{
  const char *fetchMode = blxWindowGetDefaultFetchMode(blxWindow, FALSE);
  
  GError *error = NULL;
  
  if (!strcmp(fetchMode, BLX_FETCH_PFETCH))
    {
      pfetchEntry(seqName, blxWindow, TRUE, NULL);
    }
#ifdef PFETCH_HTML 
  else if (!strcmp(fetchMode, BLX_FETCH_PFETCH_HTML))
    {
      pfetchHttpEntry(seqName, blxWindow, NULL, NULL, TRUE) ;
    }
#endif
  else if (!strcmp(fetchMode, BLX_FETCH_EFETCH))
    {
      char *command = blxprintf("efetch '%s' &", seqName);
      externalCommand(command, blxWindow);
      g_free(command);
    }
  else if (!strcmp(fetchMode, BLX_FETCH_WWW_EFETCH))
    {
      char *link = blxprintf("%s%s", URL, seqName);
      seqtoolsLaunchWebBrowser(link, &error);
      g_free(link);
    }
  else
    {
      g_critical("Unknown fetchMode: %s", fetchMode);
    }

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  return ;
}


/* Needs a better name really, sets the fetch mode by looking at env. vars etc. */
static void blxFindDefaultFetchMode(char *fetchMode)
{
  char *tmp_mode ;

  /* Check env. vars to see how to fetch EMBL entries for sequences.         */
  if ((g_getenv("BLIXEM_FETCH_PFETCH")))
    {
      strcpy(fetchMode, BLX_FETCH_PFETCH);
    }
  else if (g_getenv("BLIXEM_FETCH_WWW"))
    {
      strcpy(fetchMode, BLX_FETCH_WWW_EFETCH);
    }
  else if (g_getenv("BLIXEM_FETCH_EFETCH"))
    {
      strcpy(fetchMode, BLX_FETCH_EFETCH);
    }
  else if (g_getenv("BLIXEM_FETCH_DB"))
    {
      strcpy(fetchMode, BLX_FETCH_DB);
    }
  else if (g_getenv("BLIXEM_FETCH_REGION"))
    {
      strcpy(fetchMode, BLX_FETCH_REGION);
    }
  else if ((tmp_mode = blxConfigGetDefaultFetchMode(TRUE)))
    {
      strcpy(fetchMode, tmp_mode) ;
    }
}


/* Set program to be used for fetching sequences depending on setting
 * of fetchmode: if fetchmode is pfetch we use pfetch, _otherwise_ efetch. */
char *blxGetFetchProg(const char *fetchMode)
{
  char *fetch_prog = NULL ;

  if (!strcmp(fetchMode, BLX_FETCH_PFETCH))
    fetch_prog = "pfetch" ;
  else
    fetch_prog = "efetch" ;

  return fetch_prog ;
}



/* Use the 'pfetch' command to fetch an entry and optionally display the results.
 * If the given result argument is given it is populated with the result text and the
 * caller takes ownership; otherwise the result text is cleared up internally. */
static void pfetchEntry(char *seqName, GtkWidget *blxWindow, const gboolean displayResults, GString **result_out)
{
  /* --client gives logging information to pfetch server;
   * -F requests the full sequence entry record. For protein variants, pfetch does not
   * return anything with the -F option, so if it fails we need to re-try without -F so
   * that at least we can display the fasta sequence. */
  GString *command = g_string_sized_new(100);
  g_string_append_printf(command, "pfetch --client=%s_%s_%s -F '%s' &", g_get_prgname(), g_get_host_name(), g_get_user_name(), seqName);
  
  GString *resultText = getExternalCommandOutput(command->str);

  if (!strncasecmp(resultText->str, "no match", 8))
    {
      g_string_truncate(command, 0);
      g_string_truncate(resultText, 0);

      g_string_append_printf(command, "pfetch --client=%s_%s_%s -C '%s' &", g_get_prgname(), g_get_host_name(), g_get_user_name(), seqName);
      resultText = getExternalCommandOutput(command->str);
    }
  
  if (displayResults)
    {
      char *title = blxprintf("Blixem - %s", command->str);
      displayFetchResults(title, resultText->str, blxWindow, NULL);
      g_free(title);
    }

  if (result_out)
    {
      *result_out = resultText;
    }
  else
    {
      g_string_free(resultText, TRUE);
    }
  
  g_string_free(command, TRUE);
}


#ifdef PFETCH_HTML 
/* Gets all the sequences needed by blixem but from http proxy server instead of from
 * the pfetch server direct, this enables blixem to be run and get sequences from
 * anywhere that can see the http proxy server. If loadOptionalData is true then the full
 * EMBL entry is fetched (if it exists) and data for the optional columns is parsed (e.g.
 * organism and tissue-type) */
gboolean populateSequenceDataHtml(GList *seqsToFetch, const BlxSeqType seqType, const gboolean loadOptionalData)
{
  gboolean status = FALSE ;
  PFetchUserPrefsStruct prefs = {NULL} ;
  gboolean debug_pfetch = FALSE ;
  

  if (getPFetchUserPrefs(&prefs))
    {
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
      fetch_data.full_entry = loadOptionalData;

      fetch_data.bar = makeProgressBar(fetch_data.seq_total) ;
      fetch_data.pfetch = PFetchHandleNew(pfetch_type);

      fetch_data.section_id[0] = ' ';
      fetch_data.section_id[1] = ' ';
      fetch_data.section_id[2] = '\0';
      fetch_data.found_end_quote = FALSE;
      fetch_data.tag_name = g_string_new("");
      
      g_signal_connect(G_OBJECT(fetch_data.bar->top_level), "destroy",
		       G_CALLBACK(sequence_dialog_closed), &fetch_data) ;
      
      PFetchHandleSettings(fetch_data.pfetch, 
                           "full",       loadOptionalData,
			   "port",       prefs.port,
			   "debug",      debug_pfetch,
			   "pfetch",     prefs.location,
			   "cookie-jar", prefs.cookie_jar,
			   NULL);

     if (!loadOptionalData)
      {
        PFetchHandleSettings(fetch_data.pfetch,
                             "blixem-seqs",  TRUE, /* turns the next two on */
                             "one-per-line", TRUE, /* one sequence per line */
                             "case-by-type", TRUE, /* lowercase for DNA, uppercase for protein */
                             NULL);
      }
      
      g_free(prefs.location);
      g_free(prefs.cookie_jar);

      g_signal_connect(G_OBJECT(fetch_data.pfetch), "reader", G_CALLBACK(sequence_pfetch_reader), &fetch_data) ;

      g_signal_connect(G_OBJECT(fetch_data.pfetch), "closed", G_CALLBACK(sequence_pfetch_closed), &fetch_data) ;


      /* Build up a string containing all the sequence names. */
      GString *seq_string = g_string_sized_new(1000) ;
      GList *seqItem = seqsToFetch;
      
      for ( ; seqItem; seqItem = seqItem->next)
	{
          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
	  g_string_append_printf(seq_string, "%s ", blxSeq->fullName);
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
    }
  else
    {
      g_critical("%s", "Failed to obtain preferences specifying how to pfetch.\n") ;
      status = FALSE ;
    }

  return status ;
}


/* Use the http proxy to pfetch an entry */
static void pfetchHttpEntry(char *sequence_name,
                            GtkWidget *blxWindow, 
                            GtkWidget *dialog, 
                            GtkTextBuffer *text_buffer, 
                            const gboolean fetchFullEntry)
{
  PFetchUserPrefsStruct prefs = {NULL} ;
  gboolean debug_pfetch = FALSE ;
  
  if ((getPFetchUserPrefs(&prefs)) && (prefs.location != NULL))
    {
      PFetchData pfetch_data ;
      PFetchHandle pfetch = NULL ;
      GType pfetch_type = PFETCH_TYPE_HTTP_HANDLE ;

      if (prefs.mode && g_ascii_strncasecmp(prefs.mode, "pipe", 4) == 0)
	pfetch_type = PFETCH_TYPE_PIPE_HANDLE ;
	
      pfetch_data = g_new0(PFetchDataStruct, 1);

      pfetch_data->pfetch = pfetch = PFetchHandleNew(pfetch_type);
      
      pfetch_data->blxWindow = blxWindow;
      pfetch_data->sequence_name = sequence_name;
      pfetch_data->requested_full_entry = fetchFullEntry;
      
      if (fetchFullEntry)
        {
          pfetch_data->title = g_strdup_printf(PFETCH_TITLE_FORMAT_F, sequence_name);
        }
      else
        {
          pfetch_data->title = g_strdup_printf(PFETCH_TITLE_FORMAT, sequence_name);
        }
      
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
      
      if (PFETCH_IS_HTTP_HANDLE(pfetch))
        {
            PFetchHandleSettings(pfetch, 
                                 "full",       fetchFullEntry,
                                 "port",       prefs.port,
                                 "debug",      debug_pfetch,
                                 "pfetch",     prefs.location,
                                 "cookie-jar", prefs.cookie_jar,
                                 NULL);
        }
      else
        {
            PFetchHandleSettings(pfetch, 
                                 "full",       fetchFullEntry,
                                 "pfetch",     prefs.location,
                                 NULL);
        }
      
      g_free(prefs.location);
      g_free(prefs.cookie_jar);
      
      g_signal_connect(G_OBJECT(pfetch), "reader", G_CALLBACK(pfetch_reader_func), pfetch_data);
      
      g_signal_connect(G_OBJECT(pfetch), "closed", G_CALLBACK(pfetch_closed_func), pfetch_data);
      
      PFetchHandleFetch(pfetch, sequence_name) ;
    }
  else
    {
      g_critical("%s", "Failed to obtain preferences specifying how to pfetch.\n") ;
    }

  return ;
}


#endif



/*  getsseqsPfetch() adapted from Tony Cox's code pfetch.c
 *
 *  - this version incorporates a progress monitor as a window,
 *    much easier for user to control + has a cancel button.
 */
gboolean populateFastaDataPfetch(GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, const BlxSeqType seqType, GError **error)
{
  /* Initialise and send the requests */
  int sock;
  
  gboolean status = pfetchInit("-q -C", seqsToFetch, pfetchIP, port, External, &sock);
  
  /* Get the sequences back. They will be returned in the same order that we asked for them, i.e. 
   * in the order they are in our list. */
  GList *currentSeqItem = seqsToFetch;
  BlxSequence *currentSeq = (BlxSequence*)(currentSeqItem->data);
  BlxEmblParserState parserState = PARSING_NEWLINE;
  
  if (status && currentSeq != NULL)
    {
      int numRequested = g_list_length(seqsToFetch); /* total number of sequences requested */
      ProgressBar bar = makeProgressBar(numRequested);
      
      enum {RCVBUFSIZE = 256} ;               /* size of receive buffer */
      char buffer[RCVBUFSIZE + 1] ;           /* receive buffer */
      
      int numFetched = 0;
      int numSucceeded = 0;

      GError *tmpError = NULL;
      
      while (status && parserState != PARSING_FINISHED && parserState != PARSING_CANCELLED)
	{
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
gboolean populateFullDataPfetch(GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, const BlxSeqType seqType, GError **error)
{
  int numRequested = g_list_length(seqsToFetch); /* total number of sequences requested */

  int sock;
  gboolean status = pfetchInit("-q -C -F", seqsToFetch, pfetchIP, port, External, &sock);

  
  /* Get the sequences back. They will be returned in the same order that we asked for them, i.e. 
   * in the order they are in our list. */
  GList *currentSeqItem = seqsToFetch;
  BlxSequence *currentSeq = (BlxSequence*)(currentSeqItem->data);

  if (status && currentSeq)
    {
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


/* Set/Get global config, necessary because we don't have some blixem context pointer.... */
gboolean blxInitConfig(char *config_file, GError **error)
{
  gboolean result = FALSE ;

  g_assert(!blx_config_G) ;

  blx_config_G = g_key_file_new() ;

  if (!config_file)
    {
      result = TRUE ;
    }
  else
    {
      if (!(result = readConfigFile(blx_config_G, config_file, error)))
	{
	  g_key_file_free(blx_config_G) ;
	  blx_config_G = NULL ;
	}
    }

  return result ;
}



GKeyFile *blxGetConfig(void)
{
  return blx_config_G ;
}


//gboolean blxConfigSetFetchMode(char *mode)
//{
//  gboolean result = FALSE ;
//  GKeyFile *key_file ;
//
//  key_file = blxGetConfig() ;
//  g_assert(key_file) ;
//
//  g_key_file_set_string(key_file, BLIXEM_GROUP, BLIXEM_DEFAULT_FETCH_MODE, mode) ;
//
//  return result ;
//}


static char *blxConfigGetDefaultFetchMode(const gboolean bulk)
{
  char *fetch_mode = NULL ;
  GKeyFile *key_file ;
  GError *error = NULL ;

  key_file = blxGetConfig() ;
  g_assert(key_file) ;

  if (bulk)
    {
      fetch_mode = g_key_file_get_string(key_file, BLIXEM_GROUP, SEQTOOLS_BULK_FETCH, &error) ;
      
      if (error)
        {
          /* For backwards compatibility with old config files, try the old key name... */
          fetch_mode = g_key_file_get_string(key_file, BLIXEM_GROUP, BLIXEM_OLD_BULK_FETCH, NULL) ;
        }
    }
  else
    {
      fetch_mode = g_key_file_get_string(key_file, BLIXEM_GROUP, SEQTOOLS_USER_FETCH, &error) ;
    }
  
  return fetch_mode ;
}



gboolean blxConfigGetPFetchSocketPrefs(const char **node, int *port)
{
  gboolean result = FALSE ;
  GKeyFile *key_file ;
  GError *error = NULL ;

  if ((key_file = blxGetConfig()) && g_key_file_has_group(key_file, PFETCH_SOCKET_GROUP))
    {
      *node = g_key_file_get_string(key_file, PFETCH_SOCKET_GROUP, PFETCH_SOCKET_NODE, &error) ;
      *port = g_key_file_get_integer(key_file, PFETCH_SOCKET_GROUP, PFETCH_SOCKET_PORT, &error) ;

      /* By the time we get here _all_ the fields should have been filled in the config. */
      g_assert(!error) ;

      result = TRUE ;
    }

  return result ;
}

/* Set the preferences for pfetch mode from the config file */
gboolean blxConfigSetPFetchSocketPrefs(char *node, int port)
{
  gboolean result = TRUE ;				    /* Can't fail. */
  GKeyFile *key_file ;

  key_file = blxGetConfig() ;
  g_assert(key_file) ;

  /* Note that these calls will create the group and key if they don't exist but will
   * overwrite any existing values. */
  g_key_file_set_string(key_file, PFETCH_SOCKET_GROUP, PFETCH_SOCKET_NODE, node) ;
  g_key_file_set_integer(key_file, PFETCH_SOCKET_GROUP, PFETCH_SOCKET_PORT, port) ;

  return result ;
}


/* Set the preferences for pfetch-www mode from the config file */
gboolean blxConfigGetPFetchWWWPrefs()
{
  gboolean result = TRUE ;				    /* Can't fail. */

  /* Try the environment var first */
  URL = g_getenv("BLIXEM_FETCH_WWW");

  if (!URL)
    {
      GKeyFile *key_file = blxGetConfig() ;

      if (key_file && g_key_file_has_group(key_file, PFETCH_PROXY_GROUP))
        {
          /* Note that this call will create the URL if it doesn't exist but will
           * overwrite any existing value. */
          GError *error = NULL;

          char *tmpUrl = g_malloc(256);
          tmpUrl = g_key_file_get_string(key_file, PFETCH_PROXY_GROUP, PFETCH_PROXY_LOCATION, &error) ;
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

          URL = blxprintf("%s?request=", tmpUrl);
          g_free(tmpUrl);
        }

      if (!URL || *URL == '\0')
        {
          char *tmpUrl = g_malloc(256);
          strcpy(tmpUrl, "http://www.sanger.ac.uk/cgi-bin/seq-query?");
          URL = tmpUrl;
        }
    }


  return result ;
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


static gboolean socketSend (int sock, char *text)
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


static gboolean getPFetchUserPrefs(PFetchUserPrefsStruct *pfetch)
{
  gboolean result = FALSE ;
  GKeyFile *key_file ;
  GError *error = NULL ;

  if ((key_file = blxGetConfig()) && g_key_file_has_group(key_file, PFETCH_PROXY_GROUP))
    {
      pfetch->location = g_key_file_get_string(key_file, PFETCH_PROXY_GROUP, PFETCH_PROXY_LOCATION, &error) ;
      pfetch->cookie_jar = g_key_file_get_string(key_file, PFETCH_PROXY_GROUP, PFETCH_PROXY_COOKIE_JAR, &error) ;
      pfetch->mode = g_key_file_get_string(key_file, PFETCH_PROXY_GROUP, PFETCH_PROXY_MODE, &error) ;
      pfetch->port = g_key_file_get_integer(key_file, PFETCH_PROXY_GROUP, PFETCH_PROXY_PORT, &error) ;

      /* By the time we get here _all_ the fields should have been filled in the config. */
      g_assert(!error) ;

      result = TRUE ;
    }

  return result ;
}



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

      /* If we tried fetching the full entry and got 'no match', try again
       * without the "full" option. (It is a "feature" of pfetch that it does 
       * not return anything with the -F option for protein variants.) */
      if (pfetch_data->requested_full_entry && !strncasecmp(text, "no match", 8))
        {
          pfetchHttpEntry(pfetch_data->sequence_name, pfetch_data->blxWindow, pfetch_data->dialog, pfetch_data->text_buffer, FALSE);
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
  
  if (fetch_data->full_entry)
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
  else
    {
      /* The fetched entries just contain the sequence */
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

static void updateProgressBar(ProgressBar bar, char *sequence, int numFetched, gboolean fetch_ok)
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

gboolean readConfigFile(GKeyFile *key_file, char *config_file, GError **error)
{
  gboolean result = FALSE ;
  GKeyFileFlags flags = G_KEY_FILE_NONE ;
  
  if ((result = g_key_file_load_from_file(key_file, config_file, flags, error)))
    {
      ConfigGroup config ;
      char **groups, **group ;
      gsize num_groups ;
      int i ;
      gboolean config_loaded = FALSE ;

      groups = g_key_file_get_groups(key_file, (gsize*)(&num_groups)) ;

      for (i = 0, group = groups ; result && i < num_groups ; i++, group++)
	{
	  if ((config = getConfig(*group)))
	    {
	      config_loaded = TRUE ;
	      result = loadConfig(key_file, config, error) ;
	    }
	}

      if (!config_loaded)
	{
          g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_NO_GROUPS, "No groups found in config file.\n") ;
	  result = FALSE ;
	}
        
      g_strfreev(groups);
    }

  
  return result ;
}


/* Returns the config keyvalues pairs of the group config_name or NULL if that group cannot be found. */
ConfigGroup getConfig(char *config_name)
{
  ConfigGroup config = NULL ;
  static ConfigKeyValueStruct pfetch_http[] = {{PFETCH_PROXY_LOCATION, KEY_TYPE_STRING},
					       {PFETCH_PROXY_COOKIE_JAR, KEY_TYPE_STRING},
					       {PFETCH_PROXY_MODE, KEY_TYPE_STRING},
					       {PFETCH_PROXY_PORT, KEY_TYPE_INT},
					       {NULL, KEY_TYPE_INVALID}} ;
  
  static ConfigKeyValueStruct pfetch_socket[] = {{PFETCH_SOCKET_NODE, KEY_TYPE_STRING},
						 {PFETCH_SOCKET_PORT, KEY_TYPE_INT},
						 {NULL, KEY_TYPE_INVALID}} ;
  
  static ConfigKeyValueStruct db_fetch[] = {{DB_FETCH_DATABASE, KEY_TYPE_STRING},
                                            {NULL, KEY_TYPE_INVALID}} ;

  static ConfigKeyValueStruct region_fetch[] = {{REGION_FETCH_SCRIPT, KEY_TYPE_STRING},
                                                {NULL, KEY_TYPE_INVALID}} ;
  
  static ConfigGroupStruct groups[] = {{PFETCH_PROXY_GROUP, pfetch_http},
				       {PFETCH_SOCKET_GROUP, pfetch_socket},
                                       {DB_FETCH_GROUP, db_fetch},
                                       {REGION_FETCH_GROUP, region_fetch},
				       {NULL, NULL}} ;
  ConfigGroup tmp ;

  tmp = &(groups[0]) ;
  while (tmp->name)
    {
      if (g_ascii_strcasecmp(config_name, tmp->name) == 0)
	{
	  config = tmp ;
	  break ;
	}

      tmp++ ;
    }

  return config ;
}


/* Tries to load values for all keyvalues given in ConfigGroup fails if any are missing. */
static gboolean loadConfig(GKeyFile *key_file, ConfigGroup group, GError **error)
{
  gboolean result = TRUE;
  GError *tmpError = NULL;
  
  ConfigKeyValue key_value = group->key_values ;
  
  while (key_value->name && result)
    {
      switch(key_value->type)
	{
	case KEY_TYPE_BOOL:
	  {
	    g_key_file_get_boolean(key_file, group->name, key_value->name, &tmpError);
	    break ;
	  }
	case KEY_TYPE_INT:
	  {
	    g_key_file_get_integer(key_file, group->name, key_value->name, &tmpError);
	    break ;
	  }
	case KEY_TYPE_DOUBLE:
	  {
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
	    /* We can do this by steam using our own funcs in fact.... */

	    /* Needs 2.12.... */

	    g_key_file_get_double(key_file, group->name, key_value->name, &tmpError);

#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

	    break ;
	  }
	case KEY_TYPE_STRING:
	  {
	    g_key_file_get_string(key_file, group->name, key_value->name, &tmpError);
	    break ;
	  }
	default:
          {
            g_set_error(&tmpError, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_KEY_TYPE, "Invalid key type in configuration file.\n") ;
            break ;
          }
        }

      /* If we got an error, record what group we were reading and then quit the loop */
      if (tmpError)
        {
          g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_MISSING_KEY, "Group [%s] does not have the required key '%s': ", group->name, key_value->name);
          g_error_free(tmpError);
          tmpError = NULL;
          
          result = FALSE;
        }
      
      key_value++ ;
    }

  return result ;
}


/* Set up the net id and port for the pfetch mode, if these are available */
static gboolean setupPfetchMode(PfetchParams *pfetch, const char *fetchMode, const char **net_id, int *port, GError **error)
{
  gboolean success = TRUE;
  
  /* three ways to determine this:
   *    1) call blixem main with "-P node:port" commandline option
   *    2) setenv BLIXEM_PFETCH to a dotted quad IP address for the
   *       pfetch server, e.g. "193.62.206.200" = Plato's IP address
   *       or to the name of the machine, e.g. "plato"
   *       and optionally setenv BLIXEM_PORT to the port number for
   *       the pfetch server.
   *    3) call blixem with the -c option to specify a config file*/
  enum {PFETCH_PORT = 22100} ;			    /* default port to connect on */
  *port = PFETCH_PORT ;
  
  if (pfetch)
    {
      *net_id = pfetch->net_id ;
      if (!(*port = pfetch->port))
        *port = PFETCH_PORT ;
    }
  else if ((*net_id = g_getenv("BLIXEM_PFETCH")))
    {
      const char *port_str = g_getenv("BLIXEM_PORT");
      
      *port = port_str ? atoi(port_str) : 0 ;
      
      if (*port <= 0)
        {
          *port = PFETCH_PORT ;
        }
    }
  else
    {
      /* Lastly try for a config file. */
      blxConfigGetPFetchSocketPrefs(net_id, port) ;
    }

  if (strcmp(fetchMode, BLX_FETCH_PFETCH) == 0 && 
      (*net_id == NULL || *port == UNSET_INT))
    {
      g_set_error(error, BLX_ERROR, 1, "Network ID and port number were not found. These must be specified in the Blixem command line options, in the Blixem config file or in the BLIXEM_PFETCH and BLIXEM_PORT enviroment variables.");
      success = FALSE;
    }
  
  return success;
}


/* Set up the fetch modes. Sets required properties for each fetch mode, e.g.
 * net_id and port for pfetch-socket etc. */
void setupFetchModes(PfetchParams *pfetch, char **bulkFetchMode, char **userFetchMode, const char **net_id, int *port)
{
  *bulkFetchMode = g_malloc(32 * sizeof(char));
  *userFetchMode = g_malloc(32 * sizeof(char));
  
  /* Set up the pfetch and www-pfetch config from the file / env vars etc. We need
   * to set them up even if we're not using that mode initially, because the user can 
   * change the mode */
  blxConfigGetPFetchWWWPrefs();
      
  /* Set the fetch mode for bulk fetch. */
  if (pfetch && blxConfigSetPFetchSocketPrefs(pfetch->net_id, pfetch->port))
    {
      *bulkFetchMode = BLX_FETCH_PFETCH;
    }
  else
    {
      blxFindDefaultFetchMode(*bulkFetchMode);
    }

  /* For pfetch mode, find the net_id and port */
  setupPfetchMode(pfetch, *bulkFetchMode, net_id, port, NULL);

  /* Get the user fetch mode. This is only supplied via the config file and is
   * optional - if it is not set, use the same as the bulk fetch mode. */
  *userFetchMode = blxConfigGetDefaultFetchMode(FALSE);

  if (*userFetchMode[0] == 0)
    strcpy(*userFetchMode, *bulkFetchMode) ;
}


/* Callback called when the user has changed the fetch mode */
static gboolean onUserFetchModeChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  BlxViewContext *bc = (BlxViewContext*)data;
  GtkComboBox *combo = GTK_COMBO_BOX(widget);
  
  GtkTreeIter iter;
  
  if (gtk_combo_box_get_active_iter(combo, &iter))
    {
      GtkTreeModel *model = gtk_combo_box_get_model(combo);
      
      GValue val = {0};
      gtk_tree_model_get_value(model, &iter, 0, &val);
      
      const char *userFetchMode = g_value_get_string(&val);
      
      g_free(bc->userFetchMode);
      bc->userFetchMode = g_strdup(userFetchMode);
    }
    
  return TRUE;
}


/* Add an item to the drop-down list for selecting a user fetch mode */
static void addUserFetchItem(GtkTreeStore *store, GtkTreeIter *parent, const char *itemName, const char *currentFetchMode, GtkComboBox *combo)
{
  GtkTreeIter iter;
  gtk_tree_store_append(store, &iter, parent);
  
  gtk_tree_store_set(store, &iter, 0, itemName, -1);
  
  if (stringsEqual(itemName, currentFetchMode, FALSE))
    {
      gtk_combo_box_set_active_iter(combo, &iter);
    }
}


/* Create a drop-down box for selecting which pfetch mode blixem should use */
void createPfetchDropDownBox(GtkBox *box, GtkWidget *blxWindow)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(box, hbox, FALSE, FALSE, 0);
  
  /* Label */
  GtkWidget *label = gtk_label_new("   Select the fetch mode:    ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  
  /* Create a tree store for the list, and create the combo box itself */
  GtkTreeStore *store = gtk_tree_store_new(1, G_TYPE_STRING);
  GtkComboBox *combo = GTK_COMBO_BOX(gtk_combo_box_new_with_model(GTK_TREE_MODEL(store)));
  g_object_unref(store);
  gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo), FALSE, FALSE, 0);
  
  /* Create a cell renderer to display the list text. */
  GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
  gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, FALSE);
  gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", 0, NULL);
  
  /* Set the callback function */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  widgetSetCallbackData(GTK_WIDGET(combo), onUserFetchModeChanged, bc);
  
  /* Add the list items */
  addUserFetchItem(store, NULL, BLX_FETCH_PFETCH, bc->userFetchMode, combo);
  
#ifdef PFETCH_HTML
  addUserFetchItem(store, NULL, BLX_FETCH_PFETCH_HTML, bc->userFetchMode, combo);
#endif
  
  addUserFetchItem(store, NULL, BLX_FETCH_EFETCH, bc->userFetchMode, combo);
  addUserFetchItem(store, NULL, BLX_FETCH_WWW_EFETCH, bc->userFetchMode, combo);
}


/*******************************************************************
 *                      Pfetch utility functions                   *
 *******************************************************************/


/* Initialise a pfetch connection with the given options and calls pfetch on each of the 
 * sequences in the given list. */
static gboolean pfetchInit(char *pfetchOptions, GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, int *sock)
{
  gboolean status = TRUE;
  
  /* open socket connection */
  if (status)
    {
      *sock = socketConstruct (pfetchIP, port, External) ;
      if (*sock < 0)			/* we can't connect to the server */
	{
	  status = FALSE ;
	}
    }
  
  /* send the command/names to the server */
  if (status)
    {
      /* Send '-q' to get back one line per sequence, '-C' for lowercase DNA and uppercase protein. */
      status = socketSend (*sock, pfetchOptions);
    }
  
  if (status)
    {
      /* For each sequence, send a command to fetch that sequence, in the order that they are in our list */
      GList *seqItem = seqsToFetch;
      
      for ( ; seqItem && status ; seqItem = seqItem->next)
	{
          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
	  status = socketSend(*sock, blxSeq->fullName);
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
          updateProgressBar(bar, (*currentSeq)->fullName, numFetched, pfetch_ok) ;
          
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




