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
 * CVS info:   $Id: blixem_.h,v 1.45 2010-07-21 14:46:45 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_BLIXEM_P_H
#define DEF_BLIXEM_P_H

#include <gtk/gtk.h>
#include <wh/version.h>
#include <wh/dict.h>
#include <wh/smap.h>
#include <SeqTools/blxview.h>


/*            blixem program version and information.                        */
#define BLIXEM_TITLE   "Blixem program"
#define BLIXEM_DESC    "Sequence alignment tool."

#define BLIXEM_VERSION 4
#define BLIXEM_RELEASE 0
#define BLIXEM_UPDATE  2
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


/* Special characters for displaying in sequences */
#define SEQUENCE_CHAR_GAP    '.'   /* represents a gap in the match sequence */
#define SEQUENCE_CHAR_PAD    '-'   /* used for padding when the sequence is unavailable */
#define SEQUENCE_CHAR_BLANK  '-'   /* used to display a blank when we're not interested in what the actual base is */
#define SEQUENCE_CHAR_STOP   '*'   /* STOP codon */
#define SEQUENCE_CHAR_MET    'M'   /* MET codon */

#define selectFeaturesStr          "Feature series selection tool"
#define XY_NOT_FILLED -1000        /* Magic value meaning "value not provided" */

/* Really the buffers that use this should be dynamic but I'm not going to do that, this
 * code is so poor that it doesn't warrant the effort.... */
#define NAMESIZE      12
#define LONG_NAMESIZE 1000
#define INITDBSEQLEN  50000        /* Initial estimate of max database sequence length */
#define MAXLINE       10000


/* Main Blixem error domain */
#define BLX_ERROR g_quark_from_string("Blixem")

/* Error codes */
typedef enum
  {
    BLX_ERROR_SEQ_SEGMENT,	      /* error finding sequence segment */
    BLX_ERROR_EMPTY_STRING,           /* error code for when user entered a zero-length string */
    BLX_ERROR_STRING_NOT_FOUND,       /* error code for when a search string is not found */
    BLX_ERROR_SEQ_NAME_NOT_FOUND,     /* the sequence name(s) being searched for were not found */
    BLX_ERROR_SEQ_DATA_MISMATCH       /* same sequence was parsed more than once and data does not match */
  } BlxError;


/* Fundamental type of sequence (DNA really means nucleotide, because it could be RNA as well). */
typedef enum
  {
    BLXSEQ_INVALID, 
    BLXSEQ_DNA, 
    BLXSEQ_PEPTIDE
  } BlxSeqType ;

/* Fundamental blast mode used */
typedef enum
  {
    BLXMODE_UNSET, 
    BLXMODE_BLASTX, 
    BLXMODE_TBLASTX, 
    BLXMODE_BLASTN, 
    BLXMODE_TBLASTN, 
    BLXMODE_BLASTP
  } BlxBlastMode ;

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
    G_TYPE_INT,                 /* score */       \
    G_TYPE_INT,                 /* id */          \
    G_TYPE_INT,                 /* start */       \
    G_TYPE_POINTER,             /* sequence */    \
    G_TYPE_INT                  /* end */



/* blxview.c */
void                               blviewRedraw(void);

/* dotter.c */
void                               argvAdd(int *argc, char ***argv, char *s);
void                               loadFeatures(FILE* fil, MSP **msp);
float                              mspGetFsBottomEdge(MSP *msp, float *maxy, float height);
char                               Seqtype(char *seq);
void                               selectFeatures(void);
float                              fsTotalHeight(MSP *msplist);

/* blxparser.c */
void                               parseFS(MSP **MSPlist, FILE *file, char *opts, GList **seqList, GSList *styles,
	                                   char **seq1, char *seq1name, char **seq2, char *seq2name, const int qOffset) ;
void                               insertFS(MSP *msp, char *series);
gboolean                           mspHasFs(const MSP *msp);  
char*                              readFastaSeq(FILE *seqfile, char *qname);

/* blxFetch.c */
void                               fetchAndDisplaySequence(char *seqName, const KEY key, GtkWidget *blxWindow) ;
void                               blxFindInitialFetchMode(char *fetchMode) ;
void                               blxPfetchMenu(void) ;
char*                              blxGetFetchProg(const char *fetchMode) ;

void                               fetchSeqsIndividually(GList *seqsToFetch, GtkWidget *blxWindow);
gboolean                           blxGetSseqsPfetchHtml(GList *seqsToFetch, BlxSeqType seqType) ;
gboolean                           blxGetSseqsPfetch(GList *seqsToFetch, char* pfetchIP, int port, BOOL External) ;
gboolean                           blxGetSseqsPfetchFull(GList *seqsToFetch, char *pfetchIP, int port, BOOL External);
BOOL                               blxInitConfig(char *config_file, GError **error) ;
GKeyFile*                          blxGetConfig(void) ;
gboolean                           blxConfigSetPFetchSocketPrefs(char *node, int port) ;
gboolean                           blxConfigGetPFetchSocketPrefs(char **node, int *port) ;

/* translate.c */
char*                              blxTranslate(char *seq, char **code);
void                               blxComplement(char *seq) ;    
char*                              revComplement(char *comp, char *seq) ;

/* Create/destroy sequences and MSPs */
MSP*                               createEmptyMsp(MSP **lastMsp, MSP **mspList);
MSP*                               createNewMsp(MSP **lastMsp, MSP **mspList, GList **seqList, const BlxMspType mspType, const int score, 
                                                const char *qName, const int qStart, const int qEnd, const BlxStrand qStrand, const int qFrame, 
                                                const char *sName, const int sStart, const int sEnd, const BlxStrand sStrand, char *sequence, 
                                                char *opts, GError **error);  

void                               destroyMspList(MSP **mspList);
void                               destroyMspData(MSP *msp);
void                               destroyBlxSequenceList(GList **seqList);
void                               blviewResetGlobals();

void                               addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error);
BlxSequence*                       addBlxSequence(MSP *msp, BlxStrand strand, GList **seqList, char *sequence, GError **error);

BlxStyle*                          createBlxStyle(const char *styleName, const char *fillColor, const char *fillColorSelected, const char *fillColorPrint, const char *fillColorPrintSelected, const char *lineColor, const char *lineColorSelected, const char *lineColorPrint, const char *lineColorPrintSelected, GError **error);
void                               destroyBlxStyle(BlxStyle *style);
BlxStyle*                          getBlxStyle(const char *styleName, GSList *styles, GError **error);

BlxColor*			   getBlxColor(GArray *defaultColors, const BlxColorId colorId);
GdkColor*			   getGdkColor(BlxColorId colorId, GArray *defaultColors, const gboolean selected, const gboolean usePrintColors);

void                               createPfetchDropDownBox(GtkBox *box, GtkWidget *blxWindow);
void                               setupFetchMode(PfetchParams *pfetch, char **fetchMode, char **net_id, int *port);

gboolean                           fetchSequences(GList *seqsToFetch, 
                                                  GList *seqList,
                                                  char *fetchMode, 
                                                  const BlxSeqType seqType,
                                                  char *net_id, 
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
extern Array     fsArr;		   /* in dotter.c */
extern int       dotterGraph;
extern float     fsPlotHeight;
extern GtkWidget *blixemWindow;


#endif /*  !defined DEF_BLIXEM_P_H */
