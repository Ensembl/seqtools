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
 * Description: Internal header for blixem code.
 * HISTORY:
 * Last edited: Aug 26 09:09 2009 (edgrif)
 * Created: Thu Nov 29 10:59:09 2001 (edgrif)
 * CVS info:   $Id: blixem_.h,v 1.27 2010-05-10 16:14:49 gb10 Exp $
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
#define BLIXEM_UPDATE  1
#define BLIXEM_VERSION_NUMBER	    UT_MAKE_VERSION_NUMBER(BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_VERSION_STRING	    UT_MAKE_VERSION_STRING(BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_TITLE_STRING	    UT_MAKE_TITLE_STRING(BLIXEM_TITLE, BLIXEM_VERSION, BLIXEM_RELEASE, BLIXEM_UPDATE)
#define BLIXEM_VERSION_COMPILE	    BLIXEM_VERSION_STRING "  " __TIME__ " "__DATE__

#define BLIXEM_COPYRIGHT_STRING	    "Copyright (c) 2009-2010: Genome Research Ltd."
#define BLIXEM_WEBSITE_STRING	    ""
#define BLIXEM_LICENSE_STRING	    "Blixem is distributed under the GNU Public License, see http://www.gnu.org/copyleft/gpl.txt"

#define BLIXEM_AUTHOR_LIST	    " Originally written by Erik Sonnhammer, <Erik.Sonnhammer@sbc.su.se>",\
				    " Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define BLIXEM_AUTHOR_TEXT	    " Originally written by Erik Sonnhammer, <Erik.Sonnhammer@sbc.su.se>\n" \
				    " Rewritten by Gemma Barson, Sanger Institute, UK <gb10@sanger.ac.uk>"

#define BLIXEM_COMMENTS_STRING(TITLE, VERSION, RELEASE, UPDATE)	\
"("BLIXEM_TITLE_STRING", "					\
"compiled on - "__DATE__" "__TIME__")\n\n"			\
BLIXEM_AUTHOR_TEXT "\n"



/* Special characters for displaying in sequences */
#define SEQUENCE_CHAR_GAP	'.'   /* represents a gap in the match sequence */
#define SEQUENCE_CHAR_PAD	'-'   /* used for padding when the sequence is unavailable */
#define SEQUENCE_CHAR_BLANK	'-'   /* used to display a blank when we're not interested in what the actual base is */
#define SEQUENCE_CHAR_RES	'*'   /* not sure what this should be called, technically - residue? */


/* MSP list is sorted by one of these criteria, currently SORTBYID is the default. */
typedef enum {SORTBYUNSORTED, SORTBYSCORE, SORTBYID, SORTBYNAME, SORTBYPOS, SORTBYGROUPORDER} SortByType ;

/* Fundamental strand direction. */
typedef enum {FORWARD_STRAND, REVERSE_STRAND} Strand ;

/* Fundamental type of sequence. */
typedef enum {BLXSEQ_INVALID, BLXSEQ_DNA, BLXSEQ_PEPTIDE} BlxSeqType ;

/* Fundamental blast match mode used */
typedef enum {BLXMODE_UNSET, BLXMODE_BLASTX, BLXMODE_TBLASTX, BLXMODE_BLASTN, BLXMODE_TBLASTN, BLXMODE_BLASTP} BlxBlastMode ;


/* Really the buffers that use this should be dynamic but I'm not going to do that, this
 * code is so poor that it doesn't warrant the effort.... */
#define NAMESIZE    12
#define LONG_NAMESIZE 1000

#define INITDBSEQLEN 50000				    /* Initial estimate of max database sequence length */

#define MAXLINE 10000


/* Types and struct to support retrieving data from the window systems clipboard. */
typedef enum
  {
    BLXPASTE_INVALID,
    BLXPASTE_MATCHSET,					    /* Data is a set of matche names. */
    BLXPASTE_MOVETO					    /* Data is a single match. */
  } BlxPasteType ;

typedef struct
{
  BlxPasteType type ;					    /* defines type of data in union. */
  union
  {
    char **match_set ;
    struct
    {
      char *match ;
      int start, end ;
    } move_to ;
  } data ;
} BlxPasteDataStruct, *BlxPasteData ;


/* Structure that groups several SequenceStructs in order to hide/highlight/sort them etc. */
typedef struct _SequenceGroup
  {
    char *groupName;		/* user-friendly name for the group (should be unique to save confusion) */
    int groupId;		/* unique ID number for the group */
    int order;			/* field for sorting - lower numbers will be listed first */
    GList *seqList;		/* list of SequenceStructs */
    gboolean ownsSeqNames;	/* If true, the group will free the sequence names when it is destroyed */
    gboolean hidden;		/* true if the group should be hidden from the detail view */
    gboolean highlighted;	/* true if the group should be highlighted */
    GdkColor highlightColour;	/* the colour to highlight the group's sequences in (in both the big picture and the detail view) */
  } SequenceGroup;


/* remove ?? */
#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))


#define selectFeaturesStr     "Feature series selection tool"
#define FS(msp) (msp->type == FSSEG || msp->type == XY)
#define XY_NOT_FILLED -1000  /* Magic value meaning "value not provided" */

/* Shapes of XY data */
enum { XY_PARTIAL, XY_INTERPOLATE, XY_BADSHAPE };

  
typedef struct featureSeries_ {
    char  *name;
    int    nr;
    int    on;
    float  x;	      /* Series offset on x axis, to bump series on the screen */
    float  y;	      /* Series offset on y axis */
    int    xy;	      /* Flag for XY plot series */
} FEATURESERIES;



/* 
 * config file groups/keywords, these should not be changed willy nilly as they
 * are used external programs and users when constructing config files.
 */

/* For overall blixem settings. */
#define BLIXEM_GROUP               "blixem"
#define BLIXEM_DEFAULT_FETCH_MODE  "default-fetch-mode"


/* For http pfetch proxy fetching of sequences/entries */
#define PFETCH_PROXY_GROUP      "pfetch-http"
#define PFETCH_PROXY_LOCATION   "pfetch"
#define PFETCH_PROXY_COOKIE_JAR "cookie-jar"
#define PFETCH_PROXY_MODE       "pfetch-mode"
#define PFETCH_PROXY_PORT       "port"

/* For direct pfetch socket fetching of sequences/entries */
#define PFETCH_SOCKET_GROUP  "pfetch-socket"
#define PFETCH_SOCKET_NODE   "node"
#define PFETCH_SOCKET_PORT   "port"






/* Fetch programs for sequence entries. */
#define BLX_FETCH_PFETCH      PFETCH_SOCKET_GROUP
#ifdef PFETCH_HTML 
#define BLX_FETCH_PFETCH_HTML PFETCH_PROXY_GROUP
#endif
#define BLX_FETCH_EFETCH      "efetch"
#define BLX_FETCH_WWW_EFETCH  "WWW-efetch"
#ifdef ACEDB
#define BLX_FETCH_ACEDB       "acedb"
#define BLX_FETCH_ACEDB_TEXT  "acedb text"
#endif



/* Dotter/Blixem Package-wide functions */
void blxreadhsp(FILE *seqfile, FILE *exblxfile, char *featurefile, char *qname, 
		 int dispstart, int qoffset, char *opts, int *argc, char **argv);
char *blxTranslate(char *seq, char **code);
char *revcomp(char *comp, char *seq);
void *compl(char *seq);
void  argvAdd(int *argc, char ***argv, char *s);
void  loadFeatures(FILE* fil, MSP **msp);
float fs2y(MSP *msp, float *maxy, float height);
char  Seqtype(char *seq);
void  blviewRedraw(void);
void  selectFeatures(void);
float fsTotalHeight(MSP *msplist);
void  parseFS(MSP **MSPlist, FILE *file, char *opts,
	      char **seq1, char *seq1name, char **seq2, char *seq2name, const int qOffset) ;
void insertFS(MSP *msp, char *series);
char *readFastaSeq(FILE *seqfile, char *qname);

void blxPfetchEntry(char *sequence_name) ;
void fetchAndDisplaySequence(char *seqName, const KEY key, GtkWidget *blxWindow) ;
void blxFindInitialFetchMode(char *fetchMode) ;
void blxPfetchMenu(void) ;
char *blxGetFetchProg(const char *fetchMode) ;

BOOL blxGetSseqsPfetchHtml(MSP *msplist, DICT *dict, BlxSeqType seqType) ;
BOOL blxGetSseqsPfetch(MSP *msplist, DICT *dict, char* pfetchIP, int port, BOOL External) ;
void blxAssignPadseq(MSP *msp, MSP *msplist) ;

BOOL blxInitConfig(char *config_file, GError **error) ;
GKeyFile *blxGetConfig(void) ;
gboolean blxConfigSetPFetchSocketPrefs(char *node, int port) ;
gboolean blxConfigGetPFetchSocketPrefs(char **node, int *port) ;
char *blxConfigGetFetchMode(void) ;

void blxSeq2MSP(MSP *msp, char *seq) ;

BOOL blxRevComplement(char *seq_in, char **revcomp_seq_out) ;
char *translate(char  *seq, char **code) ;
char *revComplement(char *comp, char *seq) ;
void blxComplement(char *seq) ;

char *getqseq(int start, int end, const char const *refSeq);

void blviewDestroy(GtkWidget *blxWindow);

/* Dotter/Blixem Package-wide variables...........MORE GLOBALS...... */
extern char *blixemVersion ;
extern char *stdcode1[];        /* 1-letter amino acid translation code */
extern int   aa_atob[];
extern int PAM120[23][23];
extern Array fsArr;		/* in dotter.c */
extern int dotterGraph;
extern float fsPlotHeight;
extern GtkWidget *blixemWindow;


#endif /*  !defined DEF_BLIXEM_P_H */
