/*  File: blxmain.c
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: BLXMAIN - Standalone calling program for blxview.
 * Exported functions: 
 * HISTORY:
 * Last edited: May 26 17:13 2009 (edgrif)
 * * Aug 26 16:57 1999 (fw): added this header
 * Created: Thu Aug 26 16:56:45 1999 (fw)
 * CVS info:   $Id: blxmain.c,v 1.34 2010-11-11 15:51:48 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blxparser.h>
#include <SeqTools/blxGff3Parser.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>


/* would be good to get rid of this.... */
#define FULLNAMESIZE               255


/* Some globals.... */


/* global debug flag for blixem, set TRUE to get debugging output.           */
gboolean blixem_debug_G = FALSE ;


static char *usageText ="\n\
 Blixem - display multiple alignments against a reference sequence.\n\
\n\
 Reference:  Sonnhammer ELL & Durbin R (1994). A workbench for Large Scale\n\
 Sequence Homology Analysis. Comput. Applic. Biosci. 10:301-307.\n\
\n\
 Copyright (c) 2009-2010: Genome Research Ltd.\n\
\n\
\n\
 Usage: blixem [options] [<sequencefile>] <datafile> [X options] \n\
\n\
   <sequencefile> contains the reference sequence in FASTA format.\n\
   <datafile> is a GFF v3 file containing alignments and other features.\n\
   If <sequencefile> is ommitted, <datafile> should contain the reference\n\
   sequence in FASTA format, below a comment line that reads ##FASTA.\n\
\n\
   Both <sequencefile> and <datafile> can be substituted by \"-\"\n\
   for reading from stdin (pipe).  If <sequencefile> is piped, the first\n\
   line should contain the sequence name and the second the sequence itself.\n\
\n\
 Options:\n\
  -m <mode>, --display-mode  MANDATORY\n\
    Whether to display sequences in nucleotide or protein mode. Valid values for <mode> are:\n\
      N = nucleotide\n\
      P = protein\n\
\n\
  -a <names>, --alignment-names=<names>\n\
    Specify a string giving the names of the alignments, e.g. \"EST_mouse EST_human\" etc.\n\
\n\
  -b, --disable-big-picture\n\
    Hide the big picture section.\n\
\n\
  -c <file>, --config-file=<file>\n\
    Read configuration options from 'file'.\n\
\n\
  -h, --help\n\
    More detailed usage information.\n\
\n\
  -I, --invert-sort\n\
    Invert sorting order\n\
\n\
  -k <file>, --key-file=<file>\n\
    Read color options from a key-value file. Use --help option to see details.\n\
\n\
  -O <n>, --offset=<n>\n\
    Offset the reference sequence coordinate system by n.\n\
\n\
  -P <nodeid:port>\n\
    Causes Blixem to get sequences from a pfetch server at machine nodeid on the given\n\
    port (default 22100).\n\
\n\
  -r, --remove-input-files\n\
    Remove input files after parsing.\n\
\n\
  -s <mode>, --sort-mode=<mode>\n\
    Default sort mode. Use --help option to see details.\n\
\n\
  -S <n>, --start-coord=<n>\n\
    Start with the display centred on coordinate n.\n\
\n\
  -z, --zoom-whole\n\
    Start with the big picture zoomed out to view the full reference sequence range.\n\
\n\
  --start-next-match\n\
    Start with the display centred on the first match to the right of the default start coord.\n\
\n\
  --dotter-first\n\
    Call Dotter on the first match to the right of the default start coord.\n\
\n\
  --highlight-diffs\n\
    Enable 'highlight differences' mode, where mismatches (rather than matches) are highlighted.\n\
\n\
  --hide-inactive\n\
    Hide the inactive strand (i.e. the reverse strand, or the forward strand if the -R option\n\
    is used).\n\
\n\
  --optional-data\n\
    Parse additional data such as organism and tissue-type on start-up.\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\
\n\
"BLIXEM_AUTHOR_TEXT"\n\
 Version %s\n\%s\n" ;

static char *help_string = "\n\
FEATURES\n\
  The prefered file format for <datafile> is GFF v3. (However, Blixem is still compatible with\n\
  older file formats such as exblx and seqbl, as used by MSPcrunch).\n\
\n\
  Blixem is mainly aimed at displaying alignments, but can also show other features such as\n\
  transcripts, variations and polyA tails. The supported SO terms are:\n%s\
\n\
SORT MODE\n\
  The sort mode is specified with the -s <mode> or --sort-mode=<mode> argument, where <mode> is\n\
  one of the following:\n\
    s = by Score\n\
    i = by Identity\n\
    n = by Name\n\
    p = by Position\n\
\n\
  If optional data is loaded on start-up using the --optional-data argument, then the following\n\
  sort modes are also valid:\n\
    t = by Tissue type\n\
    m = by Strain\n\
    g = by Gene name\n\
    o = by Organism\n\
\n\
COLOR KEY FILE\n\
  The color key file is specified with the -k <file> or --key-file=<file> argument. This is a .ini-\n\
  like file that specifies attributes such as fill and line colors for features from particular \n\
  sources (say EST_Human or polya_signal). The file should contain one or more source stanzas followed\n\
  by one or more key=value pairs, i.e. \n\
\n\
    [<source>]\n\
      <key>=<value>\n\
      ...\n\
\n\
  <value> is a hex color-string, e.g. #ff0000\n\
\n\
  Possible keys are:\n\
    fill_color                    (default fill color)\n\
    fill_color_selected           (fill color when selected) \n\
    line_color                    (default line color)\n\
    line_color_selected           (line color when selected)\n\
    fill_color_utr                (default fill color for UTR regions)\n\
    fill_color_utr_selected       (fill color for UTR regions when selected) \n\
    line_color_utr                (default line color for UTR regions)\n\
    line_color_utr_selected       (line color for UTR regions when selected)\n\
  Only fill_color and line_color are mandatory; the selection colors will be calculated automatically\n\
  if not specified (a darker shade of the default color will be used when the feature is selected).\n\
  For transcripts, the fill_color/line_color/etc items are used for CDS regions and different colors\n\
  can be specified for UTR regions using fill_color_utr, line_color_utr etc.\n\
\n\
MSPcrunch\n\
To make the datafile from blast output, run MSPcrunch with option -q.\n\
\n\
 o To pipe MSPcrunch output directly to Blixem, use \"-\"\n\
   as the second parameter ([datafile]).  Example:\n\
\n\
   MSPcrunch -q <my.blast_output> | blixem <my.seq> -\n\
\n\
 o The BLAST program (blastp, blastn, blastx, tblastn, tblastx)\n\
   is automatically detected from the Blast output header by MSPcrunch\n\
   and is passed on to Blixem in the seqbl format (-q)." ;


/* set default values for command lines options */
static void initCommandLineOptions(CommandLineOptions *options, char *fetchMode, char *refSeqName)
{
  options->refSeq = NULL;
  options->refSeqName = refSeqName;
  options->refSeqRange.min = UNSET_INT;
  options->refSeqRange.max = UNSET_INT;
  options->refSeqOffset = 0;
  options->startCoord = 1;
  options->mspList = NULL;
  options->geneticCode = stdcode1;
  options->activeStrand = BLXSTRAND_FORWARD;
  options->zoomWhole = FALSE;
  options->bigPictZoom = 10;          
  options->bigPictON = TRUE;          
  options->hideInactive = FALSE;         
  options->initSortColumn = BLXCOL_ID;
  options->sortInverted = FALSE;	
  options->highlightDiffs = FALSE;   
  options->dotterFirst = FALSE;	
  options->startNextMatch = FALSE;
  options->parseFullEmblInfo = FALSE;
  options->blastMode = BLXMODE_UNSET;
  options->seqType = BLXSEQ_INVALID;
  options->numFrames = 1;
  options->fetchMode = fetchMode;
}
  

/* Utility to extract a color string with the key name 'key' from the given group in the
 * given key file. Does nothing if the given error is already set. */
static char* getColorFromKeyFile(GKeyFile *keyFile, const char *group, const char *key, GError **error)
{
  char *result = NULL;
  
  if (error == NULL || *error == NULL)
    {
      result = g_key_file_get_value(keyFile, group, key, error);

      if (error && *error)
	{
	  prefixError(*error, "Required key not found. ");
	  postfixError(*error, "\n");
	}
    }
    
  return result;
}


/* Read in the key file, which contains style information. Returns a list of
 * style structs for each style found. */
static GSList* blxReadStylesFile(char *keyFileName, GError **error)
{
  GSList *result = NULL;
  
  if (!keyFileName)
    {
      return result;
    }
  
  /* Load the key file */
  GKeyFile *keyFile = g_key_file_new();
  GKeyFileFlags flags = G_KEY_FILE_NONE ;

  if (g_key_file_load_from_file(keyFile, keyFileName, flags, error))
    {
      /* Get all the groups (i.e. style names) from the file */
      gsize num_groups;
      char **groups = g_key_file_get_groups(keyFile, &num_groups) ;

      /* Loop through each style */
      char **group;
      int i;
      GError *tmpError = NULL;
      
      for (i = 0, group = groups ; i < num_groups && !tmpError ; i++, group++)
	{
          /* Look for the keys corresponding to the style values we require */
	  char *fillColor = getColorFromKeyFile(keyFile, *group, "fill_color", &tmpError);
	  char *lineColor = getColorFromKeyFile(keyFile, *group, "line_color", &tmpError);
	  
	  /* Look for optional keys (passing error as null means we don't care if it's not found) */
	  char *fillColorSelected = getColorFromKeyFile(keyFile, *group, "fill_color_selected", NULL);
	  char *lineColorSelected = getColorFromKeyFile(keyFile, *group, "line_color_selected", NULL);
	  char *fillColorUtr = getColorFromKeyFile(keyFile, *group, "fill_color_utr", NULL);
	  char *lineColorUtr = getColorFromKeyFile(keyFile, *group, "line_color_utr", NULL);
	  char *fillColorUtrSelected = getColorFromKeyFile(keyFile, *group, "fill_color_utr_selected", NULL);
	  char *lineColorUtrSelected = getColorFromKeyFile(keyFile, *group, "line_color_utr_selected", NULL);

	  /* If there was an error, skip this group. Otherwise, go ahead and create the style */
	  if (!tmpError)
	    {
	      BlxStyle *style = createBlxStyle(*group, fillColor, fillColorSelected, lineColor, lineColorSelected, 
                                               fillColorUtr, fillColorUtrSelected, lineColorUtr, lineColorUtrSelected, &tmpError);
	      result = g_slist_append(result, style);
	    }
        }
      
      if (tmpError)
        {
          g_propagate_error(error, tmpError);
        }
      
      g_strfreev(groups);
    }

  g_key_file_free(keyFile) ;
  keyFile = NULL ;
  
  if (error && *error)
    {
      prefixError(*error, "Error reading key file '%s'. ", keyFileName);
      postfixError(*error, "\n");
    }
  
  return result;
}


/* Determine the sequence type (nucleotide or peptide) from the given char */
static BlxSeqType getSeqTypeFromChar(char seqChar)
{
  BlxSeqType result = BLXSEQ_INVALID;
  
  if (seqChar == 'n' || seqChar == 'N')
    result = BLXSEQ_DNA;
  else if (seqChar == 'p' || seqChar == 'P')
    result = BLXSEQ_PEPTIDE;
  else
    g_error("Bad display mode '%c'\n", seqChar);
  
  return result;
}


/* For backwards compatibility only: get the blast mode from one of the chars X,P,N,L,T */
static BlxBlastMode getBlastModeFromChar(char modeChar)
{
  BlxBlastMode result = BLXMODE_UNSET;
  
  switch (modeChar)
    {
      case 'N': result = BLXMODE_BLASTN;        break;
      case 'X': result = BLXMODE_BLASTX;        break;
      case 'P': result = BLXMODE_BLASTP;        break;
      case 'T': result = BLXMODE_TBLASTN;       break;
      case 'L': result = BLXMODE_TBLASTX;       break;
      default: break;
    };
  
  return result;
}


/* Get the sort mode from a char representing that mode */
static BlxColumnId getSortModeFromChar(char sortChar)
{
  BlxColumnId result;
  
  switch (sortChar)
  {
    case 's':
      result = BLXCOL_SCORE;
      break;
    case 'i':
      result = BLXCOL_ID;
      break;
    case 'n':
      result = BLXCOL_SEQNAME;
      break;
    case 'p':
      result = BLXCOL_START;
      break;
    case 't':
      result = BLXCOL_TISSUE_TYPE;
      break;
    case 'm':
      result = BLXCOL_STRAIN;
      break;
    case 'g':
      result = BLXCOL_GENE_NAME;
      break;
    case 'o':
      result = BLXCOL_ORGANISM;
      break;

    default:
      g_error("Bad sort mode: %c\n", sortChar); 
  }
  
  return result;
}


static char* getSupportedTypesAsString(GSList *supportedTypes)
{
  GString *resultStr = g_string_new(NULL);
  GSList *item = supportedTypes;
  
  for ( ; item; item = item->next)
    {
      BlxGffType *gffType = (BlxGffType*)(item->data);
      g_string_append_printf(resultStr, "    %s\n", gffType->name);
    }
  
  char *result = resultStr->str;
  g_string_free(resultStr, FALSE);
  return result;
}


/* Entry point for blixem standalone program, you should be aware when
 * altering this that blxview.c is also compiled into acedb and that
 * blxview() is called directly by acedb code.
 */
int main(int argc, char **argv)
{
  FILE *seqfile, *FSfile;
  char seqfilename[1000] = {'\0'};
  char FSfilename[1000] = {'\0'};
  char line[MAXLINE+1];
  
  int install = 1;
  
  gboolean rm_input_files = FALSE ; /* whether to remove input files once we're done with them */
  PfetchParams *pfetch = NULL ;
  gboolean xtra_data = FALSE ;      /* whether we have an extra data file to parse */
  FILE *xtra_file = NULL ;          /* the extra data file */
  char xtra_filename[1000] = {'\0'} ;
  char *align_types = NULL ;        /* string containing alignment types, to display in the title */
  char *config_file = NULL ;        /* optional blixem config file (usually "blixemrc") */
  char *key_file = NULL ;           /* optional keyword file for passing style information */
  GError *error = NULL ;
 
  char fetchMode[32] = BLX_FETCH_EFETCH;
  char refSeqName[FULLNAMESIZE+1] = "";

  static CommandLineOptions options;
  initCommandLineOptions(&options, fetchMode, refSeqName);
 
  /* Set up the GLib message handlers
   *
   * For now, also set up the acedb message package. I'm trying to phase out use of the acedb
   * message package so we can get rid of dependencies on acedb: for now, the GLib message
   * handlers just pass the messages to the acedb package. This is so that we can gradually switch
   * everything over to using GLib calls instead of calling the acedb package directly. Then when
   * we're ready we can just change the GLib handlers to use something else instead of acedb. 
   * 
   * There are two handlers: the default one for all non-critical messages, which will just log
   * output to the console, and one for critical messages and errors, which will display a 
   * pop-up message (the idea being that we don't bother the user unless it's something serious).
   * So, to get a pop-up message use g_critical, and to log a message or warning use g_message, 
   * g_warning, g_debug etc. Note that g_error is always fatal.
   */
  g_log_set_default_handler(defaultMessageHandler, NULL);
  g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL, popupMessageHandler, NULL);

  /* Stick version info. into usage string. */
  char usage[strlen(usageText) + strlen(blixemVersion) + 10];
  sprintf(usage, usageText, blixemVersion, "") ;

  /* Get the list of supported GFF types, in case we need to print them out in the usage text */
  GSList* supportedTypes = blxCreateSupportedGffTypeList();
  char *supported_types_string = getSupportedTypesAsString(supportedTypes);

  /* Get the input args. We allow long args, so we need to create a long_options array */
  static struct option long_options[] =
    {
      {"start-next-match",      no_argument,        &options.startNextMatch, 1},
      {"dotter-first",          no_argument,        &options.dotterFirst, 1},
      {"highlight-diffs",       no_argument,        &options.highlightDiffs, 1},
      {"hide-inactive",         no_argument,        &options.hideInactive, 1},
      {"optional-data",         no_argument,        &options.parseFullEmblInfo, 1},

      {"alignment-names",       required_argument,  0, 'a'},
      {"disable-big-picture",   no_argument,        0, 'b'},
      {"config-file",           required_argument,  0, 'c'},
      {"single-intput-file",    required_argument,  0, 'F'}, /* obsolete */
      {"help",                  no_argument,        0, 'h'},
      {"disable-install",       no_argument,        0, 'i'}, /* "secret" option (hide from user) */
      {"invert-sort",           no_argument,        0, 'I'},
      {"key-file",              required_argument,  0, 'k'},
      {"tblastx",               no_argument,        0, 'l'}, /* obsolete */
      {"display-mode",          required_argument,  0, 'm'},
      {"blastn",                no_argument,        0, 'n'}, /* obsolete */
      {"negate-coords",         no_argument,        0, 'N'}, /* "secret" option (hide from user) */
      {"offset",                required_argument,  0, 'O'},
      {"blastp",                no_argument,        0, 'p'}, /* obsolete */
      {"pfetch-server",         required_argument,  0, 'P'},
      {"remove-input-files",    no_argument,        0, 'r'},
      {"reverse-strand",        no_argument,        0, 'R'},
      {"start-coord",           required_argument,  0, 'S'},
      {"sort-mode",             required_argument,  0, 's'},
      {"tblastn",               no_argument,        0, 't'}, /* obsolete */
      {"extra-file",            required_argument,  0, 'x'},
      {"zoom-whole",            no_argument,        0, 'z'},
      {0, 0, 0, 0}
   };

  char        *optstring="a:bc:F:hiIk:lm:nNo:O:pP:rRS:s:tx:z";
  extern int   optind;
  extern char *optarg;
  int          optionIndex; /* getopt_long stores the index into the option struct here */
  int          optc;        /* the current option gets stored here */

  while ((optc = getopt_long(argc, argv, optstring, long_options, &optionIndex)) != EOF)
    {
      switch (optc)
	{
        case 0:
          break; /* we get here if getopt_long set a flag; nothing else to do */
          
        case '?':
          break; /* getopt_long already printed an error message */
          
	case 'a':
	  align_types = blxprintf("%s", optarg) ;
	  break;
	case 'b':
          options.bigPictON = FALSE;
	  break;
	case 'c': 
	  config_file = g_strdup(optarg) ;
	  break;
	case 'F': 
	  strcpy(FSfilename, optarg);
	  break;
	case 'h': 
          {
            char helpText[strlen(help_string) + strlen(supported_types_string) + 10];
            sprintf(helpText, help_string, supported_types_string);
            fprintf(stderr, usageText, blixemVersion, helpText) ;
            exit(EXIT_FAILURE) ;
            break;
          }
	case 'i':
	  install = 0;
	  break;
	case 'I':
	  options.sortInverted = TRUE;
	  break;
	case 'k': 
	  key_file = g_strdup(optarg) ;
          break;
	case 'l':
	  options.blastMode = BLXMODE_TBLASTX;
	  break;
	case 'm':
	  options.seqType = getSeqTypeFromChar(*optarg);
	  break;
	case 'n':
	  options.blastMode = BLXMODE_BLASTN;
	  break;
        case 'N':
          options.negateCoords = TRUE;
          break;
        case 'o':
          /* obsolete: but for backwards compatibility with ZMap, get the mode from the first char in the opts string */
          options.blastMode = getBlastModeFromChar(*optarg);
          break;
	case 'O':
          options.refSeqOffset = convertStringToInt(optarg);
	  break;
	case 'p':
	  options.blastMode = BLXMODE_BLASTP;
	  break;
	case 'P':
	  pfetch = g_malloc(sizeof(PfetchParams)) ;
	  pfetch->net_id = strtok(optarg, ":") ;
	  pfetch->port = atoi(strtok(NULL, ":")) ;
	  break;
	case 'r':
	  rm_input_files = TRUE ;
	  break ;
        case 'R':
          options.activeStrand = BLXSTRAND_REVERSE;
          break ;
        case 'S': 
	  options.startCoord = atoi(optarg);
	  break;
	case 's': 
	  options.initSortColumn = getSortModeFromChar(*optarg);
	  break;
	case 't':
	  options.blastMode = BLXMODE_TBLASTN;
	  break;
	case 'x': 
	  xtra_data = TRUE ;
	  strcpy(xtra_filename, optarg);
	  break;
        case 'z': 
          options.zoomWhole = TRUE;
          break;
            
	default : g_error("Illegal option\n");
	}
    }

  /* We expect one or two input files, or 0 input files if the FSfilename was already set with the -F option */
  const int numFiles = argc - optind;
  if (!(numFiles == 1 || numFiles == 2 || (numFiles == 0 && *FSfilename != '\0')))
    {
      fprintf(stderr, "%s\n", usage); 
      exit(EXIT_FAILURE);
    }


  /* Add -install for private colormaps */
  if (install)
    argvAdd(&argc, &argv, "-install");

  gtk_init(&argc, &argv);

  /* Set up program configuration. */
  if (!blxInitConfig(config_file, &error))
    {
      g_error("Config File Error: %s\n", error->message) ;
    }

  /* Read in the key file, if there is one */
  GSList *styles = blxReadStylesFile(key_file, &error);
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  if (numFiles == 1)
    {
      /* We have a single file containing both the aligments and the ref seq */
      strcpy(FSfilename, argv[optind]);
    }
  else if (numFiles == 2)
    {
      /* The ref seq is in a separate file. Read it in now. */
      strcpy(seqfilename, argv[optind]);
      strcpy(FSfilename, argv[optind+1]);

      if (!strcmp(seqfilename, "-"))
	{
	  seqfile = stdin;
	}
      else if(!(seqfile = fopen(seqfilename, "r")))
	{
	  g_error("Cannot open %s\n", seqfilename);
	}
	
      /* Read in query sequence */
      if (seqfile == stdin)
	{
	  /* Same file for query and FS data - sequence on first two lines */
	    
	  if (!fgets(line, MAXLINE, seqfile))
	    {
	      g_error("Error reading seqFile.\n");
	    }

	  sscanf(line, "%s", options.refSeqName);

	  /* Scan the input stream and place the chars in an extendable string */
	  GString *resultStr = g_string_sized_new(5000);
	  char charFromStream = fgetc(seqfile);
	  
	  while (charFromStream != '\n')
	    {
	      g_string_append_c(resultStr, charFromStream);
	      charFromStream = fgetc(seqfile);
	    }
	  
	  /* Set the reference sequence and free the GString (but not its data, which we've used) */
	  options.refSeq = resultStr->str;
	  g_string_free(resultStr, FALSE);
	}
      else
	{
	  fseek(seqfile, 0, SEEK_END);
	  options.refSeq = g_malloc(ftell(seqfile)+1);
	  char *refSeqChar = options.refSeq;
	  fseek(seqfile, 0, SEEK_SET);
	  
	  while (!feof(seqfile))
	    {
	      if (!fgets(line, MAXLINE, seqfile))
		{
		  break;
		}
		  
	      char *newLineCharPtr = (char *)strchr(line, '\n');
	      if (newLineCharPtr)
		{
		  *newLineCharPtr = 0;
		}
		    
	      if (*line == '>')
		{
		  strncpy(options.refSeqName, line+1, 255);
		  options.refSeqName[255]=0;
		  
		  char *spaceCharPtr = (char *)strchr(options.refSeqName, ' ');
		  if (spaceCharPtr)
		    {
		      *spaceCharPtr = 0;
		    }
		  
		  continue;
		}

	      char *lineChar = NULL;
	      for (lineChar = line; *lineChar; lineChar++) 
		{
		  *refSeqChar++ = *lineChar;
		}
		
	      *refSeqChar = 0;
	    }
	  fclose(seqfile);
	}
    }

  /* Set the blast mode from the sequence type, or vice versa, depending on which was supplied.
   * Ideally we'll get rid of blast mode eventually. */
  if (options.blastMode == BLXMODE_UNSET && options.seqType != BLXSEQ_INVALID)
    options.blastMode = (options.seqType == BLXSEQ_PEPTIDE ? BLXMODE_BLASTX : BLXMODE_BLASTN);
  else if (options.seqType == BLXSEQ_INVALID && options.blastMode != BLXMODE_UNSET)
    options.seqType = (options.blastMode == BLXMODE_BLASTN ? BLXSEQ_DNA : BLXSEQ_PEPTIDE);

  /* Parse the data file containing the homol descriptions.                */
  if (!strcmp(FSfilename, "-"))
    {
      FSfile = stdin;
    }
  else if(!(FSfile = fopen(FSfilename, "r")))
    {
      g_error("Cannot open %s\n", FSfilename);
    }
  
  /* Parser compiles lists of MSPs per type into the following array. Initialise each GList in the array to NULL */
  GList* featureLists[BLXMSP_NUM_TYPES];
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    featureLists[typeId] = NULL;
  
  GList *seqList = NULL; /* parser compiles a list of BlxSequences into this list */

  char *dummyseq = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname[FULLNAMESIZE+1] = "";
  
  parseFS(&options.mspList, FSfile, &options.blastMode, featureLists, &seqList, supportedTypes, styles,
          &options.refSeq, options.refSeqName, &options.refSeqRange, &dummyseq, dummyseqname) ;
  
  if (FSfile != stdin)
    {
      fclose(FSfile) ;
    }

  /* There may an additional file containing homol data in an alternative format. */
  if (xtra_data)
    {
      if(!(xtra_file = fopen(xtra_filename, "r")))
	{
	  g_error("Cannot open %s\n", xtra_filename) ;
	}
      
      parseFS(&options.mspList, xtra_file, &options.blastMode, featureLists, &seqList, supportedTypes, styles, &options.refSeq, options.refSeqName, &options.refSeqRange, &dummyseq, dummyseqname) ;
      fclose(xtra_file) ;
    }

  /* Remove the input files if requested.                                  */
  if (rm_input_files)
    {
      if(seqfilename[0] != '\0' && unlink(seqfilename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of sequence input file \"%s\" failed: %s\n", seqfilename, msg) ;
          g_free(msg);
        }
        
      if(FSfilename[0] != '\0' && unlink(FSfilename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of MSP input file \"%s\" failed: %s\n", FSfilename, msg) ;
          g_free(msg);
        }
        
      if (xtra_filename[0] != '\0' && unlink(xtra_filename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of extra MSP sequence input file \"%s\" failed: %s\n", xtra_filename, msg) ;
          g_free(msg);
        }
    }

  
  /* Now display the alignments, this call does not return. (Note that
   * TRUE signals blxview() that it is being called from this standalone
   * blixem program instead of as part of acedb. */
  if (blxview(&options, featureLists, seqList, supportedTypes, pfetch, align_types, TRUE))
    {
      gtk_main();
    }
    
  /* We should not get here.... */
  return (EXIT_FAILURE) ;
}



