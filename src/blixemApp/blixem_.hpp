/*  File: blixem_.h
 *  Author: Ed Griffiths, 2001-11-29
 *  Copyright [2018-2024] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
 * Description: Internal header for Blixem package-wide code
 *----------------------------------------------------------------------------
 */

#ifndef DEF_BLIXEM_P_H
#define DEF_BLIXEM_P_H

#include <gtk/gtk.h>
#include <config.h>

#include <seqtoolsUtils/utilities.hpp>
#include <seqtoolsUtils/blxmsp.hpp>
#include <seqtoolsUtils/version.hpp>
#include <seqtoolsUtils/seqtoolsFetch.hpp>
#include <gbtools/gbtools.hpp>


class BlxContext;



/*            blixem program version and information.                        */
#define BLIXEM_TITLE   "Blixem"
#define BLIXEM_PREFIX  BLIXEM_TITLE" - "  /* The prefix to use for window titles etc. */
#define BLIXEM_PREFIX_ABBREV  "B: "       /* An abbreviated version of the prefix to save space e.g. good for displaying in the taskbar */
#define BLIXEM_DESC    "Multiple alignment visualisation tool."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define BLIXEM_VERSION_STRING      SEQTOOLS_VERSION
#define BLIXEM_PACKAGE_VERSION     UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, SEQTOOLS_VERSION)
#define BLIXEM_TITLE_STRING        UT_MAKE_TITLE_STRING(BLIXEM_TITLE, BLIXEM_VERSION_STRING)
#define BLIXEM_VERSION_COMPILE     BLIXEM_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define BLIXEM_COPYRIGHT_STRING    UT_MAKE_COPYRIGHT_STRING("2009-2015")
#define BLIXEM_WEBSITE_STRING      "http://www.sanger.ac.uk/resources/software/seqtools/"
#define BLIXEM_LICENSE_STRING      UT_MAKE_LICENCE_STRING(BLIXEM_TITLE)



#define COLUMN_WIDTHS_GROUP             "column-widths"   /* group name in the config file */
#define COLUMN_SUMMARY_GROUP            "summary-columns" /* group name in the config file for
                                                             columns to include in the summary info */

/* Define the columns' default widths and titles. */
#define BLXCOL_DEFAULT_WIDTH            50    /* default width for generic columns */
#define BLXCOL_SEQNAME_WIDTH            120   /* default width for the name column */
#define BLXCOL_SCORE_WIDTH              40    /* default width for the score column */
#define BLXCOL_ID_WIDTH                 45    /* default width for the ID column */
#define BLXCOL_SOURCE_WIDTH             85    /* default width for source column  */
#define BLXCOL_GROUP_WIDTH              58    /* default width for group column  */
#define BLXCOL_START_WIDTH              50    /* default width for the start coord column */
#define BLXCOL_END_WIDTH                80    /* default width for end coord column (bigger because it also spans the scrollbar) */
#define BLXCOL_SEQUENCE_WIDTH           40    /* default width for sequence column */
#define BLXCOL_ORGANISM_WIDTH           40    /* default width for organism column */
#define BLXCOL_GENE_NAME_WIDTH          58    /* default width for gene-name column  */
#define BLXCOL_STRAIN_WIDTH             85    /* default width for strain column  */
#define BLXCOL_TISSUE_TYPE_WIDTH        100   /* default width for tissue-type column  */


/*
 * config file groups/keywords, these should not be changed willy nilly as they
 * are used external programs and users when constructing config files.
 * If you change things or add additional groups or values, then you should
 * update the getConfig function.
 */

/* For overall blixem settings. */
#define BLIXEM_OLD_BULK_FETCH      "default-fetch-mode"      /* for compatibility with old config files (new config files use SEQTOOLS_BULK_FETCH) */

/* Fetch settings */
#define FETCH_MODE_KEY             "fetch-mode"  /* any group with this key is a fetch method, and this specifies what type of fetch to do */


/* Required keys for http-and pipe-fetch groups */
#ifdef PFETCH_HTML
#define HTTP_FETCH_LOCATION       "url"
#define HTTP_FETCH_PORT           "port"
#define HTTP_FETCH_ARGS           "request"
#define HTTP_FETCH_COOKIE_JAR     "cookie-jar"
#define HTTP_FETCH_PROXY          "proxy"
#define HTTP_FETCH_IPRESOLVE      "ipresolve"
#define HTTP_FETCH_CAINFO         "cainfo"

#define PIPE_FETCH_LOCATION       "command"
#define PIPE_FETCH_ARGS           "args"
#endif

/* Required keys for socket-fetch groups */
#define SOCKET_FETCH_LOCATION     "command"
#define SOCKET_FETCH_NODE         "node"
#define SOCKET_FETCH_PORT         "port"
#define SOCKET_FETCH_ARGS         "args"

/* Required keys for www-fetch groups */
#define WWW_FETCH_LOCATION       "url"
#define WWW_FETCH_ARGS           "request"

/* Required keys for db-fetch groups */
#define DB_FETCH_LOCATION        "location"
#define DB_FETCH_QUERY           "query"
#define DB_FETCH_COLUMNS         "columns"

/* Required keys for command-fetch groups */
#define COMMAND_FETCH_SCRIPT        "command"
#define COMMAND_FETCH_ARGS          "args"


/* Generic fetch keys (applicable to more than one method) */
#define FETCH_OUTPUT       "output"      /* output format */
#define FETCH_SEPARATOR    "separator"   /* separator, when combining multiple sequences */
#define FETCH_ERRORS       "errors"      /* list of messages that indicate errors */
#define FETCH_DEBUG        "curl-debug"  /* enable verbose debug output */


/* For settings */
#define BLIXEM_SETTINGS_FILE             ".blixemrc"  /* default file name for saving blixem settings to */
#define SETTINGS_GROUP                   "user-settings"
#define STYLES_FILE_KEY                  "stylesfile" /* styles-file key in the [blixem] group */

#define SETTING_NAME_INVERT_SORT "invert-sort"
#define SETTING_NAME_HIGHLIGHT_DIFFS "highlight-diffs"
#define SETTING_NAME_HIGHLIGHT_VARIATIONS "highlight-variations"
#define SETTING_NAME_SHOW_VARIATION_TRACK "show-variations-track"
#define SETTING_NAME_SHOW_UNALIGNED "show-unaligned"
#define SETTING_NAME_SHOW_UNALIGNED_SELECTED "show-unaligned-selected-seq"
#define SETTING_NAME_LIMIT_UNALIGNED_BASES "limit-unaligned"
#define SETTING_NAME_SHOW_POLYA_SITE "show-polya-site"
#define SETTING_NAME_SHOW_POLYA_SITE_SELECTED "show-poly-site-selected-seq"
#define SETTING_NAME_SHOW_POLYA_SIG "show-poly-sig"
#define SETTING_NAME_SHOW_POLYA_SIG_SELECTED "show-polya-sig-selected-seq"
#define SETTING_NAME_SHOW_SPLICE_SITES "show-splice-sites"
#define SETTING_NAME_SHOW_MAYBE_CANONICAL "show-maybe-canonical"
#define SETTING_NAME_SQUASH_MATCHES "squash-matches"
#define SETTING_NAME_SHOW_COLINEARITY "show-colinearity"
#define SETTING_NAME_SHOW_COLINEARITY_SELECTED "show-colinearity-selected"


/* would be good to get rid of this.... */
#define FULLNAMESIZE               255


#define MKSTEMP_CONST_CHARS_GFF               "BLIXEM_gff"        /* the prefix to use when creating a temp file name */
#define MKSTEMP_REPLACEMENT_CHARS             "XXXXXX"            /* the required string that will be replaced by unique chars when creating a temp file name */


#define HIGHLIGHT_BOX_Y_PAD             2         /* this provides space between highlight box and the top/bottom of the grid */
#define HIGHLIGHT_BOX_MIN_WIDTH         5         /* minimum width of the highlight box */


/* Blixem config/fetch error domain */
#define BLX_CONFIG_ERROR g_quark_from_string("Blixem config")
#define BLX_FETCH_ERROR g_quark_from_string("Blixem config")

/* Error codes */
typedef enum
  {
    BLX_CONFIG_ERROR_NO_GROUPS,             /* no groups in config file */
    BLX_CONFIG_ERROR_INVALID_KEY_TYPE,      /* invalid key type given in config file */
    BLX_CONFIG_ERROR_MISSING_KEY,           /* a required key is missing */
    BLX_CONFIG_ERROR_INVALID_FETCH_MODE,    /* invalid fetch mode specified */
    BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT, /* invalid output format specified for fetch mode */
    BLX_CONFIG_ERROR_WARNINGS,              /* warnings found while reading config file */
    BLX_CONFIG_ERROR_SUBSTITUTION,          /* error with substitution string */
    BLX_CONFIG_ERROR_INVALID_FETCH_METHOD,  /* null fetch method */
    BLX_CONFIG_ERROR_NO_EXE,                /* fetch method executable does not exist */
    BLX_CONFIG_ERROR_NULL_FETCH,            /* fetch method is null */
    BLX_CONFIG_ERROR_NO_ARGS,               /* mandatory args weren't specified */
    BLX_CONFIG_ERROR_INVALID_IPRESOLVE      /* invalid value given for ipresolve config */
  } BlxConfigError;


/* Error codes */
typedef enum
  {
    BLX_FETCH_ERROR_SOCKET,                /* error creating socket */
    BLX_FETCH_ERROR_HOST,                  /* unknown host */
    BLX_FETCH_ERROR_CONNECT,               /* error connecting to host */
    BLX_FETCH_ERROR_SEND,                  /* error sending to socket */
    BLX_FETCH_ERROR_PATH                   /* executable not found in path */
  } BlxFetchError;



/* Function pointer to a function that performs a fetch of a particular sequence */
typedef void(*FetchFunc)(const char *seqName, gpointer fetchMethod, const gboolean bulk, GtkWidget *blxWindow);


/* The following are used to define default colors for certain types of features in Blixem.
 * One of several different actual colors from the BlxColor struct may be used depending
 * on state, e.g. we use a different color if "print colors" (i.e. black and
 * white mode) is on. */

typedef enum
  {
    BLXCOLOR_MIN,           /* dummy value so that we don't get a zero ID */

    BLXCOLOR_BACKGROUND,    /* background color of the widgets */
    BLXCOLOR_REF_SEQ,       /* default background color for the reference sequence */
    BLXCOLOR_MATCH,         /* background color for an exact match */
    BLXCOLOR_CONS,          /* background color for a conserved match */
    BLXCOLOR_MISMATCH,      /* background color for a mismatch */
    BLXCOLOR_INSERTION,     /* color for an insertion marker */
    BLXCOLOR_EXON_START,    /* color for the start boundary line of an exon */
    BLXCOLOR_EXON_END,      /* color for the end boundary line of an exon */
    BLXCOLOR_CODON,         /* color in which to highlight the nucleotides for the currently-selected codon */
    BLXCOLOR_MET,           /* background color for MET codons in the three frame translation */
    BLXCOLOR_STOP,          /* background color for STOP codons in the three frame translation */
    BLXCOLOR_GRID_LINE,     /* color of the gridlines in the big picture grids */
    BLXCOLOR_GRID_TEXT,     /* color of the text in the big picture grids */
    BLXCOLOR_HIGHLIGHT_BOX, /* color of the highlight box in the big picture */
    BLXCOLOR_PREVIEW_BOX,   /* color of the preview box in the big picture */
    BLXCOLOR_MSP_LINE,      /* color of the MSP lines in the big picture */
    BLXCOLOR_SNP,           /* background color for SNPs */
    BLXCOLOR_GROUP,         /* default highlight color for generic groups */
    BLXCOLOR_MATCH_SET,     /* default highlight color for the special match-set group */
    BLXCOLOR_EXON_FILL,     /* fill color for an exon in the big picture */
    BLXCOLOR_EXON_LINE,     /* line color for an exon in the big picture */
    BLXCOLOR_CDS_FILL,      /* fill color for a CDS in the big picture (coding region) */
    BLXCOLOR_CDS_LINE,      /* line color for a CDS in the big picture (coding region) */
    BLXCOLOR_UTR_FILL,      /* fill color for an UTR in the big picture (non-coding/untranslated region) */
    BLXCOLOR_UTR_LINE,      /* line color for an UTR in the big picture (non-coding/untranslated region) */
    BLXCOLOR_UNALIGNED_SEQ, /* color in which to show additional sequence in the match that is not part of the alignment */
    BLXCOLOR_CANONICAL,     /* background highlight color for canonical intron bases */
    BLXCOLOR_NON_CANONICAL, /* background highlight color for non-canonical intron bases */
    BLXCOLOR_MAYBE_CANONICAL,/* background highlight color for splice sites that are "maybe"
                              * canonical, that is, that would be canonical if on the other
                              * strand (helps to identify problems in data; common with BAM) */
    BLXCOLOR_POLYA_TAIL,    /* background color for polyA tails in the detail view */
    BLXCOLOR_POLYA_SIGNAL,  /* background color for non-annotated polyA signals in the detail view */
    BLXCOLOR_POLYA_SIGNAL_ANN,/* background color for annotated polyA signals in the detail view */
    BLXCOLOR_POLYA_SITE_ANN,/* background color for annotated polyA sites in the detail view */
    BLXCOLOR_TREE_GRID_LINES,/* color of the tree grid lines (i.e. column separator lines) */
    BLXCOLOR_CLIP_MARKER,   /* color of the marker line used to indicate a match has been clipped */
    BLXCOLOR_COVERAGE_PLOT, /* color of the coverage plot */
    BLXCOLOR_ASSEMBLY_GAP,  /* highlight color for assembly gaps */
    BLXCOLOR_SELECTION,     /* highlight color for selections */
    BLXCOLOR_PARTIAL_EXON_CROSSHATCH, /* color of cross-hatch lines for partial exons */
    BLXCOLOR_COLINEAR_PERFECT, /* color of lines joining alignment blocks that are perfectly colinear */
    BLXCOLOR_COLINEAR_IMPERFECT, /* color of lines joining alignment blocks that are imperfectly colinear */
    BLXCOLOR_COLINEAR_NOT,  /* color of lines joining alignment blocks that are not colinear */

    BLXCOL_NUM_COLORS
  } BlxColorId;


/* Utility struct to pass around a set of colors from a styles file */
typedef struct _BlxStyleColors
{
  char *fillColor;
  char *lineColor;
  char *fillColorSelected;
  char *lineColorSelected;
  char *fillColorCds;
  char *lineColorCds;
  char *fillColorCdsSelected;
  char *lineColorCdsSelected;
  gboolean normalFound;
  gboolean cdsFound;
} BlxStyleColors;


/* This enum contains a list of all the boolean options that the user can toggle on/off */
typedef enum
  {
    BLXFLAG_MIN,                    /* Start index for looping through flags */

    BLXFLAG_INVERT_SORT,            /* Inverts the default sort order */
    BLXFLAG_HIGHLIGHT_DIFFS,        /* Hides matching bases and highlights mis-matching ones */
    BLXFLAG_HIGHLIGHT_VARIATIONS,   /* Whether to highlight bases that have variations in the reference sequence */
    BLXFLAG_SHOW_VARIATION_TRACK,   /* Shows the Variations track */
    BLXFLAG_SHOW_UNALIGNED,         /* Shows additional bits of the match sequence that are not part of the aligned section */
    BLXFLAG_SHOW_UNALIGNED_SELECTED,/* Only show unaligned bits of sequence for the currently-selected sequence(s) */
    BLXFLAG_LIMIT_UNALIGNED_BASES,  /* If the show-unaligned-sequence option is on, limits how many bases from the unaligned sequence are shown */
    BLXFLAG_SHOW_POLYA_SITE,        /* Show polyA tails */
    BLXFLAG_SHOW_POLYA_SITE_SELECTED,/* Only show polyA tails for the currently-selected sequence(s) */
    BLXFLAG_SHOW_POLYA_SIG,         /* Show polyA signals in the reference sequence */
    BLXFLAG_SHOW_POLYA_SIG_SELECTED,/* Only show polyA signals for the currently-selected sequence(s) */
    BLXFLAG_SHOW_SPLICE_SITES,      /* Highlights splice sites in the reference sequence for the currently-selected MSPs */
    BLXFLAG_SHOW_MAYBE_CANONICAL,   /* Highlights "maybe canonical" splice sites */
    BLXFLAG_OPTIONAL_COLUMNS,       /* Gets set to true if the optional columns have been loaded */
    BLXFLAG_SHOW_CDS,               /* True if CDS/UTR regions should be shown; false if plain exons should be shown */
    BLXFLAG_NEGATE_COORDS,          /* True if coords should be negated when display is reversed (so coords appear to increase left-to-right when really they decrease) */
    BLXFLAG_SAVE_TEMP_FILES,        /* save any temporary files that blixem creates, e.g. the GFF file created by the region-fetch fetch mode */
    BLXFLAG_ABBREV_TITLE,           /* whether to abbreviate the window titles to save space */
    BLXFLAG_LINK_FEATURES,          /* whether features with the same name should be linked */
    BLXFLAG_SHOW_COLINEARITY,       /* whether to show colinearity lines between alignment blocks */
    BLXFLAG_SHOW_COLINEARITY_SELECTED, /* whether to show colinearity lines for selected sequence only */

    BLXFLAG_NUM_FLAGS               /* Total number of flags e.g. for creating arrays and loops etc */
  } BlxFlag;


/* Structure representing a group of sequences. This groups several BlxSequences together and
 * sets various flags so that we can hide/highlight/etc all of the sequences in a group. */
typedef struct _SequenceGroup
  {
    char *groupName;               /* user-friendly name for the group (should be unique to save confusion) */
    int groupId;                   /* unique ID number for the group */
    int order;                     /* field for sorting - lower numbers will be listed first */
    GList *seqList;                /* list of BlxSequences */
    gboolean ownsSeqNames;         /* If true, the group will free the sequence names when it is destroyed */
    gboolean hidden;               /* true if the group should be hidden from the detail view */
    gboolean isFilter;             /* true if this is a filter (i.e. hide everything not in a
                                    * filter) */
    gboolean isQuickGroup;         /* true if the group was created by a quick-group or
                                    * quick-filter option */
    gboolean highlighted;          /* true if the group should be highlighted */
    GdkColor highlightColor;       /* the color to highlight the group's sequences in (in both the big picture and the detail view) */
  } SequenceGroup;


/* This enum gives a more meaningful way of indexing the "opts" string */
typedef enum
  {
    BLXOPT_MODE = 0,              /* Blast mode: L = tblastx, N = blastn, P = blastp, T = tblastn, X = blastx */
    BLXOPT_START_NEXT_MATCH = 1,  /* 'M' means start at next match. Blank to ignore. */
    BLXOPT_SHOW_BIG_PICT = 2,     /* 'B' to show big picture, 'b' to hide it */
    BLXOPT_REV_BIG_PICT = 3,      /* 'R' to show the reverse strand in the big picture; 'r' just show the forward strand (blxparser.c) */
    BLXOPT_FULL_EMBL_INFO = 4,    /* 'F' to parse full embl info on startup; blank to just parse the sequence data (quicker) */
    BLXOPT_FULL_ZOOM = 5,         /* 'Z' to use full zoom by default; blank otherwise (blxparser.c) */
    BLXOPT_INVERT_SORT = 6,       /* 'I' to invert the default sort order; blank otherwise */
    BLXOPT_HSP_GAPS = 7,          /* 'G' for HSP gaps mode; blank otherwise (blxparser.c) */
    BLXOPT_SORT_MODE = 8,         /* Initial sort mode: i=ID, n=name, p=position, s=score, t=tissue type, m=strain, g=gene name, o=organism;
                                   * OR: set to 'd' to automatically dotter the first sequence; blank otherwise (blxselect.c) */

    BLXOPT_NUM_OPTS               /* Total number of options (must always be last in list) */
  } BlxOptsIdx ;


/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the BlxContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    BLXDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */

    BLXDIALOG_HELP,                 /* The Help dialog */
    BLXDIALOG_SETTINGS,             /* The Settings dialog */
    BLXDIALOG_SORT,                 /* The Sort dialog */
    BLXDIALOG_FIND,                 /* The Find dialog */
    BLXDIALOG_GROUPS,               /* The Groups dialog */
    BLXDIALOG_VIEW,                 /* The View dialog */
    BLXDIALOG_DOTTER,               /* The Dotter dialog */

    BLXDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } BlxDialogId;


/* This enum determines the how to find the match to call dotter on */
typedef enum
  {
    BLXDOTTER_MATCH_SELECTED,  /* call dotter on the currently-selected match sequence */
    BLXDOTTER_MATCH_ADHOC,     /* call dotter on a manually-pasted sequence */
    BLXDOTTER_MATCH_SELF       /* call dotter on the reference sequence vs itself */
  } DotterMatchType ;

/* This enum determines how to find the reference sequence range to call dotter on */
typedef enum
  {
    BLXDOTTER_REF_AUTO,        /* call dotter on an automatically-calculated ref seq range */
    BLXDOTTER_REF_MANUAL,      /* call dotter on a manually entered ref seq range */
    BLXDOTTER_REF_TRANSCRIPT   /* call dotter on a transcript range */
  } DotterRefType;


/* This enum gives access into the depthArray for the specific counter (i.e. per-base or all) */
typedef enum
  {
    DEPTHCOUNTER_NONE,

    DEPTHCOUNTER_ALL_F, // all reads at this coord on + strand
    DEPTHCOUNTER_GAP_F, // all reads at this coord on + strand
    DEPTHCOUNTER_A_F,   // all reads with 'a' at this coord on forward strand
    DEPTHCOUNTER_C_F,   // all reads with 'c' at this coord on forward strand
    DEPTHCOUNTER_G_F,   // all reads with 'g' at this coord on forward strand
    DEPTHCOUNTER_T_F,   // all reads with 't' or 'u' at this coord on forward strand
    DEPTHCOUNTER_N_F,   // all reads with 'n' at this coord on forward strand

    DEPTHCOUNTER_ALL_R, // all reads at this coord on - strand
    DEPTHCOUNTER_GAP_R, // all reads at this coord on - strand
    DEPTHCOUNTER_A_R,   // all reads with 'a' at this coord on reverse strand
    DEPTHCOUNTER_C_R,   // all reads with 'c' at this coord on reverse strand
    DEPTHCOUNTER_G_R,   // all reads with 'g' at this coord on reverse strand
    DEPTHCOUNTER_T_R,   // all reads with 't' or 'u' at this coord on reverse strand
    DEPTHCOUNTER_N_R,   // all reads with 'n' at this coord on reverse strand

    DEPTHCOUNTER_NUM_ITEMS /*must be last*/
} DepthCounter;


/* Struct to hold all the settings that come from the command line options */
class CommandLineOptions
{
public:
  char *refSeq;                   /* the section of reference sequence we're viewing */
  char *refSeqName;               /* the name of the reference sequence */
  IntRange refSeqRange;           /* the range of the reference sequence (before any offset is applied) */
  int refSeqOffset;               /* if non-zero, all parsed coords will be offset by this amount */
  int startCoord;                 /* which coord to start the initial display range at */
  gboolean startCoordSet;         /* true if the start coord has been specified on the command line */
  MSP *mspList;                   /* the list of alignments */
  GList *columnList;              /* the list of display columns */
  char **geneticCode;             /* the genetic code */

  BlxStrand activeStrand;         /* which strand will initially be the active one */
  gboolean negateCoords;          /* if this option is true, the display will show coords as negative when the reverse strand is active */
  gboolean zoomWhole;             /* whether to zoom out to view the entire big picture range on startup */
  int bigPictZoom;                /* initial zoom level for the big picture (as a multiple of the initial detail view range) */
  IntRange bigPictRange;          /* initial range for the big picture (will override bigPictZoom if both are set) */
  gboolean bigPictON;             /* whether to show the big picture by default */
  gboolean hideInactive;          /* whether to show the inactive strand in the big picture/detail view */
  BlxColumnId initSortColumn;     /* initial column to sort by */
  gboolean sortInverted;          /* whether initial sort order should be inverted */
  gboolean highlightDiffs;        /* whether the initial display should highlight mismatches rather than matches */
  gboolean dotterFirst;           /* open dotter when blixem starts */
  gboolean startNextMatch;        /* start at the coord of the next match from the default start coord */
  gboolean squashMatches;         /* start with the 'squash matches' option on */
  gboolean optionalColumns;       /* load data for optional columns on startup to populate additional info like tissue-type and strain */
  gboolean saveTempFiles;         /* save any temporary files that blixem creates */
  gboolean coverageOn;            /* show the coverage view on start-up */
  gboolean abbrevTitle;           /* if true, use a abbreviated window titles to save space */

  gboolean mspFlagDefaults[MSPFLAG_NUM_FLAGS]; /* default values for MSP flags */

  BlxBlastMode blastMode;         /* the blast match mode */
  BlxSeqType seqType;             /* whether the display shows sequences as peptides or nucleotides */
  int numFrames;                  /* the number of reading frames */
  GHashTable *fetchMethods;       /* table of fetch methods (keyed on name as a GQuark) */
  GArray *bulkFetchDefault;       /* the default method(s) for bulk fetching sequences (can be overridden by an MSPs data-type properties) */
  GArray *userFetchDefault;       /* the default method(s) for fetching individual sequences interactively */
  GArray *optionalFetchDefault;   /* the default method(s) for bulk fetching optional sequence data */
  char *dataset;                  /* the name of a dataset, e.g. 'human' */
  BlxMessageData msgData;         /* data to be passed to the message handlers */
  gboolean mapCoords;             /* whether the map-coords command-line argument was specified */
  int mapCoordsFrom;              /* the coord to map from */
  int mapCoordsTo;                /* the coord to map to */
  char *windowColor;              /* if not null, set the main window background color to this */

  bool fetch_debug;               /* whether to include verbose debug output for fetch */
  long ipresolve;                 /* whether to make curl use ipv4 or ipv6 */
  const char *cainfo;             /* location of curl cainfo file */
};


/* blixem can use either efetch (default) or a pfetch server to get
 * sequences, to use pfetch node/port information must be specified. */
typedef struct
  {
    char *net_id ;
    int port ;
  } PfetchParams ;


/* Used to pass around info about a match sequence when getting fetch command arguments */
typedef struct
{
  const char *match_name;
  const char *ref_name;
  int match_start;
  int match_end;
  const char *dataset;
  const char *source;
  const char *filename;
} MatchSequenceData ;


/* Class to perform a user-fetch operation */
class UserFetch
{

public:
  UserFetch();

  UserFetch(const BlxSequence *blxSeq,
            const gboolean displayResults,
            GtkWidget *blxWindow,
            GtkWidget *dialog,
#ifdef PFETCH_HTML
            long ipresolve,
            const char *cainfo,
#endif
            bool debug);

  void performFetch();

  GtkTextBuffer *getTextBuffer();
  void setTextBuffer(GtkTextBuffer *text_buffer);

#ifdef PFETCH_HTML
  bool httpFetchSequence(const BlxFetchMethod *fetchMethod);
#endif
  void socketFetchSequence(const BlxFetchMethod *fetchMethod);
  void commandFetchSequence(const BlxFetchMethod *fetchMethod);
  void internalFetchSequence(const BlxFetchMethod *fetchMethod);
  void wwwFetchSequence(const BlxFetchMethod *fetchMethod);
  void sqliteFetchSequence(const BlxFetchMethod *fetchMethod);

private:

  const BlxSequence *blxSeq;
  gboolean displayResults;
  int attempt;
  GtkWidget *blxWindow;
  GtkWidget *dialog;
  GtkTextBuffer *text_buffer;
  bool debug;

#ifdef PFETCH_HTML
  long ipresolve;
  const char *cainfo;
#endif
};


/* Class to perform a bulk-fetch operation */
class BulkFetch
{
public:
  BulkFetch(gboolean External,
            gboolean saveTempFiles,
            BlxSeqType seqType,
            GList **seqList,
            GList *columnList,
            GArray *defaultFetchMethods,
            GHashTable *fetchMethods,
            MSP **mspList,
            BlxBlastMode *blastMode,
            GArray* featureLists[],
            GSList *supportedTypes,
            GSList *styles,
            int refSeqOffset,
            IntRange* const refSeqRange,
            char *dataset,
            gboolean optionalColumns,
            GHashTable *lookupTable,
#ifdef PFETCH_HTML
            long ipresolve,
            const char *cainfo,
#endif
            bool debug);

  gboolean performFetch();
  gboolean fetchList(GList *seqsToFetch, const BlxFetchMethod* const fetchMethod, GError **error);
  gboolean httpFetchList(GList *seqsToFetch, const BlxFetchMethod* const fetchMethod, GError **error);
  gboolean socketFetchList(GList *seqsToFetch, const BlxFetchMethod* const fetchMethod, GError **error);
  void regionFetchList(GList *regionsToFetch, const BlxFetchMethod* const fetchMethod, GError **error);
  void commandFetchList(GList *regionsToFetch, const BlxFetchMethod* const fetchMethod, GError **error);


private:

  void regionFetchFeature(const MSP* const msp,
                          const BlxFetchMethod* const fetchMethod,
                          const char *script,
                          const char *dataset,
                          const char *tmpDir,
                          GError **error);

  int attempt;
  gboolean External;
  gboolean saveTempFiles;
  BlxSeqType seqType;
  GList **seqList; /* list of BlxSequence structs for all required sequences */
  GList *columnList;
  GArray *defaultFetchMethods;
  GHashTable *fetchMethods;
  MSP **mspList;
  BlxBlastMode *blastMode;
  GArray** featureLists;
  GSList *supportedTypes;
  GSList *styles;
  int refSeqOffset;
  IntRange* refSeqRange;
  char *dataset;
  gboolean optionalColumns;
  GHashTable *lookupTable;
  bool debug;

#ifdef PFETCH_HTML
  long ipresolve;
  const char *cainfo;
#endif
};



/* blxview.c */
/* Function to show blixem window, can be called from any application. */
gboolean                            blxview(CommandLineOptions *options,
                                            GArray* featureLists[],
                                            GList *seqList,
                                            GSList *supportedTypes,
                                            PfetchParams *pfetch,
                                            char *align_types,
                                            gboolean External,
                                            GSList *styles,
                                            GHashTable *lookupTable) ;

BlxColumnInfo*                     getColumnInfo(const GList *columnList, const BlxColumnId columnId);
int                                getColumnWidth(const GList *columnList, const BlxColumnId columnId);
const char*                        getColumnTitle(const GList *columnList, const BlxColumnId columnId);
void                               getColumnXCoords(const GList *columnList, const BlxColumnId columnId, IntRange *xRange);
void                               saveColumnWidths(GList *columnList, GKeyFile *key_file);
void                               saveSummaryColumns(GList *columnList, GKeyFile *key_file);
gboolean                           showColumn(BlxColumnInfo *columnInfo);
void                               resetColumnWidths(GList *columnList);

void                               blviewRedraw(void);
GtkWidget*                         getBlixemWindow(void);
const IntRange*                    mspGetFullSRange(const MSP* const msp, const gboolean seqSelected, const BlxContext* const bc);
const IntRange*                    mspGetDisplayRange(const MSP* const msp);
const IntRange*                    mspGetFullDisplayRange(const MSP* const msp, const gboolean seqSelected, const BlxContext* const bc);
void                               mspCalculateFullExtents(MSP *msp, const BlxContext* const bc, const int numUnalignedBases);
void                               cacheMspDisplayRanges(const BlxContext* const bc, const int numUnalignedBases);

gboolean                           mspGetMatchCoord(const MSP *msp,
                                                    const int qIdx,
                                                    const gboolean seqSelected,
                                                    const int numUnalignedBases,
                                                    BlxContext *bc,
                                                    int *result_out);


void                               drawAssemblyGaps(GtkWidget *widget,
                                                    GdkDrawable *drawable,
                                                    GdkColor *color,
                                                    const gboolean displayRev,
                                                    GdkRectangle *rect,
                                                    const IntRange* const dnaRange,
                                                    const GArray *mspArray);

GList*                             blxCreateColumns(const gboolean optionalColumns, const gboolean customSeqHeader);

GSList*                            blxReadStylesFile(const char *keyFileName_in, GError **error);

const char*                        blxGetAppName();
const char*                        blxGetTitlePrefix(const BlxContext* const bc);
const char*                        blxGetCopyrightString();
const char*                        blxGetWebSiteString();
char*                              blxGetCommentsString();
const char*                        blxGetLicenseString();
const char*                        blxGetVersionString();

/* dotter.c */
//void                               selectFeatures(void);

/* blxparser.c */
gboolean                           mspHasFs(const MSP *msp);
char*                              readFastaSeq(FILE *seqfile, char *qname, int *startCoord, int *endCoord, const BlxSeqType seqType);

/* blxFetch.c */
BlxFetchMethod*                    getFetchMethodDetails(GQuark fetchMethodQuark, GHashTable *fetchMethods);
GString*                           getFetchCommand(const BlxFetchMethod* const fetchMethod, const BlxSequence *blxSeq, const MSP* const msp, const char *refSeqName, const int refSeqOffset, const IntRange* const refSeqRange, const char *dataset, GError **error);
GString*                           doGetFetchCommand(const BlxFetchMethod* const fetchMethod,MatchSequenceData *match_data, GError **error);
GString*                           getFetchArgs(const BlxFetchMethod* const fetchMethod, const BlxSequence *blxSeq,const MSP* const msp,const char *refSeqName,const int refSeqOffset,const IntRange* const refSeqRange,const char *dataset,GError **error);
GString*                           getFetchArgsMultiple(const BlxFetchMethod* const fetchMethod, GList *seqsToFetch, GError **error);
void                               finaliseFetch(GList *seqList, GList *columnList);
void                               sendFetchOutputToFile(GString *command, GKeyFile *keyFile, BlxBlastMode *blastMode,GArray* featureLists[],GSList *supportedTypes, GSList *styles, GList **seqList, MSP **mspListIn,const char *fetchName, const gboolean saveTempFiles, MSP **newMsps, GList **newSeqs, GList *columnList, GHashTable *lookupTable, const int refSeqOffset, const IntRange* const refSeqRange, GError **error);
const char*                        outputTypeStr(const BlxFetchOutputType outputType);

void                               fetchSeqsIndividually(GList *seqsToFetch, GtkWidget *blxWindow);
gboolean                           populateSequenceDataHtml(GList *seqsToFetch, const BlxSeqType seqType, const BlxFetchMethod* const fetchMethod) ;
gboolean                           populateFastaDataPfetch(GList *seqsToFetch, BlxFetchMethod *fetchMethod, gboolean External, const BlxSeqType seqType, GError **error) ;

gboolean                           populateFullDataPfetch(GList *seqsToFetch, BlxFetchMethod *fetchMethod, gboolean External, const BlxSeqType seqType, GError **error);
void                               blxInitConfig(const char *config_file, CommandLineOptions *options, GError **error) ;
GKeyFile*                          blxGetConfig(void) ;

void                               loadNativeFile(const char *filename, const char *buffer, GKeyFile *keyFile, BlxBlastMode *blastMode, GArray* featureLists[], GSList *supportedTypes, GSList *styles, MSP **newMsps, GList **newSeqs, GList *columnList, GHashTable *lookupTable, const int refSeqOffset, const IntRange* const refSeqRange, GError **error);
void                               blxMergeFeatures(MSP *newMsps, GList *newSeqs, MSP **mspList, GList **seqList);

/* Create/destroy sequences and MSPs */
void                               blviewResetGlobals();

BlxStyle*                          createBlxStyle(const char *styleName, const char *fillColor, const char *fillColorSelected, const char *lineColor, const char *lineColorSelected, const char *fillColorUtr, const char *fillColorUtrSelected, const char *lineColorUtr, const char *lineColorUtrSelected, GError **error);
void                               destroyBlxStyle(BlxStyle *style);

void                               createPfetchDropDownBox(GtkBox *box, GtkWidget *blxWindow);


/* Dotter/Blixem Package-wide variables...........MORE GLOBALS...... */
extern char      *stdcode1[];      /* 1-letter amino acid translation code */
extern int       aa_atob[];
extern int       PAM120[23][23];
extern int       dotterGraph;
extern float     fsPlotHeight;
extern GtkWidget *blixemWindow;


/* blxFetchDb.c */

/* Function pointer for sqlite callback functions */
typedef int (*SqliteFunc)(void*,int,char**,char**);

void sqliteFetchSequences(GList *seqsToFetch, const BlxFetchMethod* const fetchMethod, GList *columnList, GError **error);
void sqliteFetchSequence(const BlxSequence* const blxSeq, const BlxFetchMethod* const fetchMethod,const gboolean displayResults, const int attempt,GtkWidget *blxWindow);
void sqliteValidateFetchMethod(const BlxFetchMethod* const fetchMethod, GError **error);
int sqliteDisplayResultsCB(void *data, int argc, char **argv, char **azColName);
void sqliteRequest(const char *database, const char *query, SqliteFunc callbackFunc, void *callbackData, GError **error);


#endif /*  !defined DEF_BLIXEM_P_H */
