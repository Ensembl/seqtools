/*  File: blixem_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2001
 *-------------------------------------------------------------------
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
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Internal header for Dotter and Blixem package-wide code
 * 
 * HISTORY:
 * Last edited: Aug 26 09:09 2009 (edgrif)
 * Created: Thu Nov 29 10:59:09 2001 (edgrif)
 * CVS info:   $Id: blixem_.h,v 1.61 2010-11-09 10:13:48 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_BLIXEM_P_H
#define DEF_BLIXEM_P_H

#include <gtk/gtk.h>
#include <config.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/version.h>


/*            blixem program version and information.                        */
#define BLIXEM_TITLE   "Blixem"
#define BLIXEM_DESC    "Multiple alignment visualisation tool."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define BLIXEM_VERSION_STRING	   PACKAGE_VERSION
#define BLIXEM_PACKAGE_VERSION	   UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, PACKAGE_VERSION)
#define BLIXEM_TITLE_STRING	   UT_MAKE_TITLE_STRING(BLIXEM_TITLE, BLIXEM_VERSION_STRING)
#define BLIXEM_VERSION_COMPILE	   BLIXEM_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define BLIXEM_COPYRIGHT_STRING	   UT_MAKE_COPYRIGHT_STRING("2009-2010")
#define BLIXEM_WEBSITE_STRING	   ""
#define BLIXEM_LICENSE_STRING	   UT_MAKE_LICENCE_STRING(BLIXEM_TITLE)

#define BLIXEM_COMMENTS_STRING()                                \
"("BLIXEM_TITLE_STRING", "					\
UT_COMPILE_PHRASE " " UT_MAKE_COMPILE_DATE() ")\n\n"            \
AUTHOR_TEXT "\n"


/* 
 * config file groups/keywords, these should not be changed willy nilly as they
 * are used external programs and users when constructing config files.
 */

/* For overall blixem settings. */
#define BLIXEM_GROUP               "blixem"
#define BLIXEM_DEFAULT_FETCH_MODE  "default-fetch-mode"


/* For http pfetch proxy fetching of sequences/entries */
#define PFETCH_PROXY_GROUP         "pfetch-http"
#define PFETCH_PROXY_LOCATION      "pfetch"
#define PFETCH_PROXY_COOKIE_JAR    "cookie-jar"
#define PFETCH_PROXY_MODE          "pfetch-mode"
#define PFETCH_PROXY_PORT          "port"

/* For direct pfetch socket fetching of sequences/entries */
#define PFETCH_SOCKET_GROUP        "pfetch-socket"
#define PFETCH_SOCKET_NODE         "node"
#define PFETCH_SOCKET_PORT         "port"


/* Fetch programs for sequence entries. */
#define BLX_FETCH_PFETCH           PFETCH_SOCKET_GROUP
#ifdef PFETCH_HTML 
#define BLX_FETCH_PFETCH_HTML      PFETCH_PROXY_GROUP
#endif
#define BLX_FETCH_EFETCH           "efetch"
#define BLX_FETCH_WWW_EFETCH       "WWW-efetch"


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
    BLXCOLOR_POLYA_TAIL,    /* background color for polyA tails in the detail view */
    BLXCOLOR_TREE_GRID_LINES,/* color of the tree grid lines (i.e. column separator lines) */
    BLXCOLOR_CLIP_MARKER,   /* color of the marker line used to indicate a match has been clipped */

    BLXCOL_NUM_COLORS
  } BlxColorId;
  
  
/* This enum contains a list of all the boolean options that the user can toggle on/off */
typedef enum
  {
    BLXFLAG_MIN,		    /* Start index for looping through flags */
  
    BLXFLAG_SQUASH_MATCHES,	    /* Puts all MSPs from the same sequence on the same row in the detail view */
    BLXFLAG_INVERT_SORT,	    /* Inverts the default sort order */
    BLXFLAG_HIGHLIGHT_DIFFS,	    /* Hides matching bases and highlights mis-matching ones */
    BLXFLAG_HIGHLIGHT_VARIATIONS,   /* Whether to highlight bases that have variations in the reference sequence */
    BLXFLAG_SHOW_VARIATION_TRACK,   /* Shows the Variations track */
    BLXFLAG_SHOW_UNALIGNED,	    /* Shows additional bits of the match sequence that are not part of the aligned section */
    BLXFLAG_SHOW_UNALIGNED_SELECTED,/* Only show unaligned bits of sequence for the currently-selected sequence(s) */
    BLXFLAG_LIMIT_UNALIGNED_BASES,  /* If the show-unaligned-sequence option is on, limits how many bases from the unaligned sequence are shown */
    BLXFLAG_SHOW_POLYA_SITE,        /* Show polyA tails */
    BLXFLAG_SHOW_POLYA_SITE_SELECTED,/* Only show polyA tails for the currently-selected sequence(s) */
    BLXFLAG_SHOW_POLYA_SIG,         /* Show polyA signals in the reference sequence */
    BLXFLAG_SHOW_POLYA_SIG_SELECTED,/* Only show polyA signals for the currently-selected sequence(s) */
    BLXFLAG_SHOW_SPLICE_SITES,	    /* Highlights splice sites in the reference sequence for the currently-selected MSPs */
    BLXFLAG_EMBL_DATA_LOADED,       /* Gets set to true if the full EMBL data is parsed and populated in the MSPs */
    BLXFLAG_SHOW_CDS,               /* True if CDS/UTR regions should be shown; false if plain exons should be shown */
    BLXFLAG_NEGATE_COORDS,          /* True if coords should be negated when display is reversed (so coords appear to increase left-to-right when really they decrease) */
    
    BLXFLAG_NUM_FLAGS		    /* Number of flags, for looping through flags or creating an array */
  } BlxFlag;


/* Structure representing a group of sequences. This groups several BlxSequences together and 
 * sets various flags so that we can hide/highlight/etc all of the sequences in a group. */
typedef struct _SequenceGroup
  {
    char *groupName;		   /* user-friendly name for the group (should be unique to save confusion) */
    int groupId;		   /* unique ID number for the group */
    int order;			   /* field for sorting - lower numbers will be listed first */
    GList *seqList;		   /* list of BlxSequences */
    gboolean ownsSeqNames;	   /* If true, the group will free the sequence names when it is destroyed */
    gboolean hidden;		   /* true if the group should be hidden from the detail view */
    gboolean highlighted;	   /* true if the group should be highlighted */
    GdkColor highlightColor;	   /* the color to highlight the group's sequences in (in both the big picture and the detail view) */
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



/* COLUMNS: To add a new column you must do the following:
 *    - add an identifier for the column to the BlxColumnId enum;
 *    - add the type to the TREE_COLUMN_TYPE_LIST definition; and
 *    - specify the data source in addSequenceStructToRow and addMspToRow;
 *    - create the column in createColumns(...)
 * and optionally:
 *    - add a custom data function in createTreeColumn;
 *    - add a custom header widget and/or header refresh function in createTreeColHeader;
 *    - specify sort behaviour in sortColumnCompareFunc. */

/* This enum declares identifiers for each column in the detail-view trees. If you add an enum
 * here you must also add its type to the TREE_COLUMN_TYPE_LIST definition below. */
typedef enum
  {
    BLXCOL_SEQNAME,             /* The match sequence's name */
    BLXCOL_SOURCE,              /* The match's source */
    BLXCOL_ORGANISM,            
    BLXCOL_GENE_NAME,            
    BLXCOL_TISSUE_TYPE,            
    BLXCOL_STRAIN,            
    BLXCOL_GROUP,               /* The group that this alignment belongs to */
    BLXCOL_SCORE,               /* The alignment's score */
    BLXCOL_ID,                  /* The alignment's %ID */
    BLXCOL_START,               /* The start coord of the alignment on the match sequence */
    BLXCOL_SEQUENCE,            /* This column will display the part of the alignment currently in the display range. */
    BLXCOL_END,                 /* The end coord of the alignment on the match sequence */
    
    BLXCOL_NUM_COLUMNS          /* The number of columns; must always be the last item in this enum */
  } BlxColumnId;


/* This defines the variable type for each detail-view-tree column. These MUST be the 
 * correct types (in the correct order) for the columns listed in the BlxColumnId enum above. */
#define TREE_COLUMN_TYPE_LIST                     \
    G_TYPE_STRING,              /* seq name */    \
    G_TYPE_STRING,              /* source */      \
    G_TYPE_STRING,              /* organism */    \
    G_TYPE_STRING,              /* gene name */   \
    G_TYPE_STRING,              /* tissue type */ \
    G_TYPE_STRING,              /* strain */      \
    G_TYPE_STRING,              /* group */       \
    G_TYPE_DOUBLE,              /* score */       \
    G_TYPE_DOUBLE,              /* id */          \
    G_TYPE_INT,                 /* start */       \
    G_TYPE_POINTER,             /* sequence */    \
    G_TYPE_INT                  /* end */


/* Struct to hold all the settings that come from the command line options */
typedef struct _CommandLineOptions
{
  char *refSeq;			  /* the section of reference sequence we're viewing */
  char *refSeqName;               /* the name of the reference sequence */
  IntRange refSeqRange;           /* the range of the reference sequence (before any offset is applied) */
  int refSeqOffset;               /* if non-zero, all parsed coords will be offset by this amount */
  int startCoord;		  /* which coord to start the initial display range at */
  MSP *mspList;			  /* the list of alignments */
  char **geneticCode;             /* the genetic code */
  
  BlxStrand activeStrand;         /* which strand will initially be the active one */
  gboolean negateCoords;          /* if this option is true, the display will show coords as negative when the reverse strand is active */
  gboolean zoomWhole;             /* whether to zoom out to view the entire big picture range on startup */
  int bigPictZoom;                /* initial zoom level for the big picture (as a multiple of the initial detail view range) */
  gboolean bigPictON;             /* whether to show the big picture by default */
  gboolean hideInactive;          /* whether to show the inactive strand in the big picture/detail view */
  BlxColumnId initSortColumn;     /* initial column to sort by */
  gboolean sortInverted;	  /* whether initial sort order should be inverted */
  gboolean highlightDiffs;        /* whether the initial display should highlight mismatches rather than matches */
  gboolean dotterFirst;		  /* open dotter when blixem starts */
  gboolean startNextMatch;	  /* start at the coord of the next match from the default start coord */
  gboolean parseFullEmblInfo;     /* parse the full EMBL files on startup to populate additional info like tissue-type */
  BlxBlastMode blastMode;         /* the blast match mode */
  BlxSeqType seqType;             /* whether the display shows sequences as peptides or nucleotides */
  int numFrames;                  /* the number of reading frames */
  char *fetchMode;		  /* the default method for fetching sequences */
} CommandLineOptions;


/* blixem can use either efetch (default) or a pfetch server to get
 * sequences, to use pfetch node/port information must be specified. */
typedef struct
  {
    char *net_id ;
    int port ;
  } PfetchParams ;


/* blxview.c */
/* Function to show blixem window, can be called from any application. */
gboolean                            blxview(CommandLineOptions *options,
                                            GList* featureLists[],
                                            GList *seqList, 
                                            GSList *supportedTypes,
	                                    PfetchParams *pfetch, 
                                            char *align_types, 
                                            gboolean External) ;

void                               blviewRedraw(void);
GList*                             getSeqsToPopulate(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData);
int				   findMspListSExtent(GList *mspList, const gboolean findMin);
int				   findMspListQExtent(GList *mspList, const gboolean findMin);
void                               mspGetFullSRange(const MSP const *msp, 
                                                    const gboolean seqSelected,
                                                    const gboolean *flags,
                                                    const int numUnalignedBases, 
                                                    const GList const *polyASiteList,
                                                    IntRange *sSeqRange);

void                                mspGetFullQRange(const MSP const *msp, 
                                                     const gboolean seqSelected,
                                                     const gboolean *flags,
                                                     const int numUnalignedBases, 
                                                     const GList const *polyASiteList,
                                                     const int numFrames, 
                                                     IntRange *sSeqRange);

int                                 gapCoord(const MSP *msp, 
                                             const int qIdx, 
                                             const int numFrames, 
                                             const BlxStrand strand, 
                                             const gboolean displayRev,
                                             const gboolean seqSelected,
                                             const int numUnalignedBases,
                                             gboolean *flags,
                                             const GList const *polyASiteList);


/* dotter.c */
//void                               selectFeatures(void);

/* blxparser.c */
gboolean                           mspHasFs(const MSP *msp);  
char*                              readFastaSeq(FILE *seqfile, char *qname);

/* blxFetch.c */
void                               fetchAndDisplaySequence(char *seqName, GtkWidget *blxWindow) ;
void                               blxFindInitialFetchMode(char *fetchMode) ;
void                               blxPfetchMenu(void) ;
char*                              blxGetFetchProg(const char *fetchMode) ;

void                               fetchSeqsIndividually(GList *seqsToFetch, GtkWidget *blxWindow);
gboolean                           populateSequenceDataHtml(GList *seqsToFetch, const BlxSeqType seqType, const gboolean loadOptionalData) ;
gboolean                           populateFastaDataPfetch(GList *seqsToFetch, const char* pfetchIP, int port, gboolean External, const BlxSeqType seqType, GError **error) ;
gboolean                           populateFullDataPfetch(GList *seqsToFetch, const char *pfetchIP, int port, gboolean External, const BlxSeqType seqType, GError **error);
gboolean                           blxInitConfig(char *config_file, GError **error) ;
GKeyFile*                          blxGetConfig(void) ;
gboolean                           blxConfigSetPFetchSocketPrefs(char *node, int port) ;
gboolean                           blxConfigGetPFetchSocketPrefs(const char **node, int *port) ;
gboolean                           blxConfigGetPFetchWWWPrefs();


/* Create/destroy sequences and MSPs */
void                               destroyMspList(MSP **mspList);
void                               destroyMspData(MSP *msp);
void                               destroyBlxSequenceList(GList **seqList);
void                               blviewResetGlobals();

BlxStyle*                          createBlxStyle(const char *styleName, const char *fillColor, const char *fillColorSelected, const char *fillColorPrint, const char *fillColorPrintSelected, const char *lineColor, const char *lineColorSelected, const char *lineColorPrint, const char *lineColorPrintSelected, GError **error);
void                               destroyBlxStyle(BlxStyle *style);

void                               createPfetchDropDownBox(GtkBox *box, GtkWidget *blxWindow);
void                               setupFetchMode(PfetchParams *pfetch, char **fetchMode, const char **net_id, int *port);

gboolean                           fetchSequences(GList *seqsToFetch, 
                                                  GList *seqList,
                                                  char *fetchMode, 
                                                  const BlxSeqType seqType,
                                                  const char *net_id, 
                                                  int port, 
                                                  const gboolean parseOptionalData,
                                                  const gboolean parseSequenceData,
                                                  const gboolean External,
                                                  GError **error);

/* Dotter/Blixem Package-wide variables...........MORE GLOBALS...... */
extern char      *blixemVersion;
extern char      *stdcode1[];      /* 1-letter amino acid translation code */
extern int       aa_atob[];
extern int       PAM120[23][23];
extern int       dotterGraph;
extern float     fsPlotHeight;
extern GtkWidget *blixemWindow;


#endif /*  !defined DEF_BLIXEM_P_H */
