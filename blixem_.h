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
 * CVS info:   $Id: blixem_.h,v 1.55 2010-10-22 11:58:58 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_BLIXEM_P_H
#define DEF_BLIXEM_P_H

#include <gtk/gtk.h>
#include <SeqTools/blxview.h>
#include <SeqTools/utilities.h>
#include <wh/version.h>


/*            blixem program version and information.                        */
#define BLIXEM_TITLE   "Blixem program"
#define BLIXEM_DESC    "Sequence alignment tool."

#define BLIXEM_VERSION 4
#define BLIXEM_RELEASE 0
#define BLIXEM_UPDATE  4
#define BLIXEM_VERSION_NUMBER	   UT_MAKE_VERSION_NUMBER(BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_VERSION_STRING	   UT_MAKE_VERSION_STRING(BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_TITLE_STRING	   UT_MAKE_TITLE_STRING(BLIXEM_TITLE, BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_VERSION_COMPILE	   BLIXEM_VERSION_STRING "  " __TIME__ " "__DATE__

#define BLIXEM_COPYRIGHT_STRING	   "Copyright (c) 2009-2010: Genome Research Ltd."
#define BLIXEM_WEBSITE_STRING	   ""
#define BLIXEM_LICENSE_STRING	   "Blixem is distributed under the GNU Public License, see http://www.gnu.org/copyleft/gpl.txt"

#define BLIXEM_AUTHOR_LIST	   " Originally written by Erik Sonnhammer, <Erik.Sonnhammer@sbc.su.se>",\
				   " Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define BLIXEM_AUTHOR_TEXT	   " Originally written by Erik Sonnhammer, <Erik.Sonnhammer@sbc.su.se>\n" \
				   " Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define BLIXEM_COMMENTS_STRING(TITLE, VERSION, RELEASE, UPDATE)	\
"("BLIXEM_TITLE_STRING", "					\
"compiled on - "__DATE__" "__TIME__")\n\n"			\
BLIXEM_AUTHOR_TEXT "\n"


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
#ifdef ACEDB
#define BLX_FETCH_ACEDB            "acedb"
#define BLX_FETCH_ACEDB_TEXT       "acedb text"
#endif


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


/* This enum is to record the type of data currently being parsed by the parser. An input file can 
 * contain multiple types of data. The start of a new section of data is indicated by a header 
 * line beginning with a hash and a known type name, e.g. "# exblx_x" or "##gff-version   3"
 *
 * FS and SFS (Feature Series) data type names have "FS" or "SFS" followed by a specific format 
 * type, e.g. "# SFS type=SEG". 
 *
 * For some data types, additional data is included in the header line(s) as well as in the 'body' 
 * section below it. For these, there are two data types in the enum, one postfixed with _HEADER 
 * and one with _BODY. When the header line is detected the initial type is set to the _HEADER enum, 
 * and when we are finished processing the header information it is set to _BODY, so that we go on
 * to process the body of the data on the next loop through. For types with no information in the 
 * header, there is only a _BODY enum.
 */
typedef enum
  {
    PARSER_START,	           /* indicates that we haven't started processing yet */
    PARSER_ERROR,	           /* indiates an error state */
    
    GFF_3_HEADER,	           /* GFF 3 header */
    GFF_3_BODY,		           /* GFF 3 data */
    
    FASTA_SEQ_HEADER,              /* FASTA sequence header */
    FASTA_SEQ_BODY,                /* Sequence data in FASTA format */

    EXBLX_BODY,                    /* Old style sequence entries. */
    SEQBL_BODY,                    /* Old style sequence entries. */
    EXBLX_X_BODY,	           /* New style sequence entries with gaps and match strand. (_X stands for eXtended.) */ 
    SEQBL_X_BODY,	           /* New style sequence entries with gaps and match strand. (_X stands for eXtended.) */

    FS_HSP_BODY,                   /* feature-series HSP data */
    
    FS_GSP_HEADER,                 /* feature-series GSP data header */
    FS_GSP_BODY,                   /* feature-series GSP data */
    
    FS_GFF_BODY,                   /* feature-series GFF data */
    
    FS_SEG_BODY,                   /* feature-series segment data */
    
    FS_XY_HEADER,                  /* feature-series XY data header */
    FS_XY_BODY,                    /* feature-series XY data */
    
    FS_SEQ_HEADER,                 /* feature-series sequence data header */
    FS_SEQ_BODY                    /* feature-series sequence data */
  } BlxParserState ;
  

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



/* blxview.c */
void                               blviewRedraw(void);
GList*                             getSeqsToPopulate(GList *inputList, const gboolean getSequenceData, const gboolean getOptionalData);
int				   findMspListSExtent(GList *mspList, const gboolean findMin);
int				   findMspListQExtent(GList *mspList, const gboolean findMin);
void				   defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void                               popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);


/* dotter.c */
void                               argvAdd(int *argc, char ***argv, char *s);
void                               loadFeatures(FILE* fil, MSP **msp);
float                              mspGetFsBottomEdge(MSP *msp, float *maxy, float height);
char                               Seqtype(char *seq);
void                               selectFeatures(void);
float                              fsTotalHeight(MSP *msplist);

/* blxparser.c */
gboolean                           mspHasFs(const MSP *msp);  
char*                              readFastaSeq(FILE *seqfile, char *qname);

/* blxFetch.c */
#ifdef ACEDB
void                               fetchAndDisplaySequence(char *seqName, const KEY key, GtkWidget *blxWindow) ;
#else
void                               fetchAndDisplaySequence(char *seqName, GtkWidget *blxWindow) ;
#endif

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
MSP*                               createNewMsp(MSP **lastMsp, MSP **mspList, GList **seqList, const BlxMspType mspType, char *source, const gdouble score, const int phase,
                                                char *url, char *idTag, char *qName, const int qStart, const int qEnd, const BlxStrand qStrand, const int qFrame, 
                                                char *sName, const int sStart, const int sEnd, const BlxStrand sStrand, char *sequence, 
                                                char *opts, GError **error);  

void                               destroyMspList(MSP **mspList);
void                               destroyMspData(MSP *msp);
void                               destroyBlxSequenceList(GList **seqList);
void                               blviewResetGlobals();

void                               addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error);
BlxSequence*                       addBlxSequence(const char *name, const char *idTag, BlxStrand strand, GList **seqList, char *sequence, MSP *msp, GError **error);

BlxStyle*                          createBlxStyle(const char *styleName, const char *fillColor, const char *fillColorSelected, const char *fillColorPrint, const char *fillColorPrintSelected, const char *lineColor, const char *lineColorSelected, const char *lineColorPrint, const char *lineColorPrintSelected, GError **error);
void                               destroyBlxStyle(BlxStyle *style);
BlxStyle*                          getBlxStyle(const char *styleName, GSList *styles, GError **error);

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
