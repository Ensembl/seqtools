/*  File: blxview.h
 *  Author: Erik Sonnhammer, 92-02-20
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
 * Description: External interface to blixem window code, this interface
 *              is used both by xace and the stand alone blixem.
 * HISTORY:
 * Last edited: Aug 21 13:57 2009 (edgrif)
 * * Aug 26 16:57 1999 (fw): added this header
 * Created: Thu Aug 26 16:57:17 1999 (fw)
 * CVS info:   $Id: blxview.h,v 1.30 2010-08-05 08:55:05 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_BLXVIEW_H
#define DEF_BLXVIEW_H

#include <gtk/gtk.h>

#ifdef ACEDB
#include <wh/regular.h>
#endif

/* Only used in pephomolcol.c, would be good to get rid of this.... */
#define FULLNAMESIZE               255



/* blixem can use either efetch (default) or a pfetch server to get
 * sequences, to use pfetch node/port information must be specified. */
typedef struct
{
  char *net_id ;
  int port ;
} PfetchParams ;


/* Blixem can read a number of file formats which are documented in:
 * 
 * The following defines give keyword strings for these file formats
 * which should be used in code that writes these files.
 * 
 *  */
#define BLX_GAPS_TAG               "Gaps"
#define BLX_DESCRIPTION_TAG        "Description"
#define BLX_SEQUENCE_TAG           "Sequence"

/* Utility struct to hold a range of integers */
typedef struct _IntRange
  {
    int min;
    int max;
  } IntRange ;
  

/* Fundamental strand direction. */
typedef enum
  {
    BLXSTRAND_NONE, 
    BLXSTRAND_FORWARD, 
    BLXSTRAND_REVERSE
  } BlxStrand ;


/* Structure that contains information about a sequence */
typedef struct _BlxSequence
{
  char *fullName;                  /* full name of the sequence and variant, including prefix characters, e.g. EM:AV274505.2 */
  char *shortName;                 /* short name of the sequence, excluding prefix and variant, e.g. AV274505 */
  char *variantName;               /* short name of the variant, excluding prefix but including variant number, e.g. AV274505.2 */
  
  GString *organism;               /* organism from the EMBL data OS line */
  GString *geneName;               /* gene name from the EMBL data GN line */
  GString *tissueType;             /* tissue type from the /tissue_type attribute from the FT lines in the EMBL file */
  GString *strain;                 /* strain from the /strain attribute from the FT lines in the EMBL file */
  
  BlxStrand strand;                /* which strand of the sequence this is */
  GString *sequence;               /* the actual sequence data */
  gboolean sequenceReqd;           /* whether the sequence data is required (e.g. it is not needed for exons/introns etc.) */
  
  GList *mspList;                  /* list of MSPs from this sequence */
} BlxSequence;


typedef struct _FeatureSeries
  {
    char        *name;             /* Name of the Feature Series */
    int         on;                /* Flag to show/hide the Feature Series */
    int         order;             /* Order number used for sorting */
    float       x;	           /* Series offset on x axis, to bump series on the screen */
    float       y;	           /* Series offset on y axis */
    int         xy;	           /* Flag for XY plot series */
  } FeatureSeries;


/* Shapes of XY curves */
typedef enum
  { 
    BLXCURVE_PARTIAL, 
    BLXCURVE_INTERPOLATE, 
    BLXCURVE_BADSHAPE
  } BlxCurveShape;


/* Supported types of MSP */
typedef enum 
  {
    BLXMSP_INVALID,                /* No valid type was set */
    
    BLXMSP_MATCH,                  /* A match (i.e. alignment) */
    BLXMSP_MATCH_SET,              /* The parent of a set of matches. Can be used to specify generic properties such as color. */
    BLXMSP_EXON_CDS,               /* Exon (coding section) */
    BLXMSP_EXON_UTR,               /* Exon (untranslated region i.e. non-coding section) */
    BLXMSP_EXON_UNK,               /* Exon (unknown type) */
    BLXMSP_INTRON,                 /* Intron */
    BLXMSP_POLYA_TAIL,		   /* polyA tail */
    
    BLXMSP_SNP,                    /* Single Nucleotide Polymorphism */
    
    BLXMSP_HSP,                    /*  */
    BLXMSP_GSP,                    /*  */

    BLXMSP_FS_SEG,                 /* Feature Series Segment */
    BLXMSP_XY_PLOT                 /* x/y coordinates - for plotting feature-series curves */
  } BlxMspType;



/* The following are used to define default colors for certain types of features in Blixem.
 * One of several different actual colors from the BlxColor struct may be used depending 
 * on state, e.g. we use a different color if "print colors" (i.e. black and 
 * white mode) is on. */

typedef enum 
  {
    BLXCOLOR_MIN,		  /* dummy value so that we don't get a zero ID */
  
    BLXCOLOR_BACKGROUND,	  /* background color of the widgets */
    BLXCOLOR_REF_SEQ,	  /* default background color for the reference sequence */  
    BLXCOLOR_MATCH,	  /* background color for an exact match */
    BLXCOLOR_CONS,	  /* background color for a conserved match */
    BLXCOLOR_MISMATCH,	  /* background color for a mismatch */
    BLXCOLOR_EXON_CDS,	  /* background color for an exon (coding region) */
    BLXCOLOR_EXON_UTR,	  /* background color for an exon (non-coding/untranslated region) */
    BLXCOLOR_INSERTION,	  /* color for an insertion marker */
    BLXCOLOR_EXON_START,	  /* color for the start boundary line of an exon */
    BLXCOLOR_EXON_END,  	  /* color for the end boundary line of an exon */
    BLXCOLOR_CODON,	  /* color in which to highlight the nucleotides for the currently-selected codon */
    BLXCOLOR_MET,		  /* background color for MET codons in the three frame translation */
    BLXCOLOR_STOP,	  /* background color for STOP codons in the three frame translation */
    BLXCOLOR_GRID_LINE,	  /* color of the gridlines in the big picture grids */
    BLXCOLOR_GRID_TEXT,	  /* color of the text in the big picture grids */
    BLXCOLOR_HIGHLIGHT_BOX, /* color of the highlight box in the big picture */
    BLXCOLOR_PREVIEW_BOX,	  /* color of the preview box in the big picture */
    BLXCOLOR_MSP_LINE,	  /* color of the MSP lines in the big picture */
    BLXCOLOR_SNP,		  /* background color for SNPs */
    BLXCOLOR_GROUP,	  /* default highlight color for generic groups */
    BLXCOLOR_MATCH_SET,	  /* default highlight color for the special match-set group */
    BLXCOLOR_EXON_FILL_CDS, /* fill color for an exon in the big picture (coding region) */
    BLXCOLOR_EXON_FILL_UTR, /* fill color for an exon in the big picture (non-coding/untranslated region) */
    BLXCOLOR_EXON_LINE_CDS, /* line color for an exon in the big picture (coding region) */
    BLXCOLOR_EXON_LINE_UTR, /* line color for an exon in the big picture (non-coding/untranslated region) */
    BLXCOLOR_UNALIGNED_SEQ, /* color in which to show additional sequence in the match that is not part of the alignment */
    BLXCOLOR_CANONICAL,     /* background highlight color for canonical intron bases */
    BLXCOLOR_NON_CANONICAL, /* background highlight color for non-canonical intron bases */
    BLXCOLOR_POLYA_TAIL,    /* background color for polyA tails in the detail view */
    BLXCOLOR_TREE_GRID_LINES,/* color of the tree grid lines (i.e. column separator lines) */
    BLXCOLOR_CLIP_MARKER,   /* color of the marker line used to indicate a match has been clipped */

    BLXCOL_NUM_COLORS
  } BlxColorId;

typedef struct _BlxColor
  {
    char *name;			  /* meaningful name for the color e.g. "Match" */
    char *desc;			  /* meaningful description for what the color is used for e.g. "background color for exact matches" */
    gboolean transparent;	  /* if this is true, the colors below are not specified and the background color should be used instead */
    
    GdkColor normal;		  /* the color in normal operation */
    GdkColor selected;		  /* the color in a selected state */
    GdkColor print;		  /* the color used for printing */
    GdkColor printSelected;	  /* the selected-state color used for printing */
  } BlxColor;


/* Define a drawing style for an MSP */
typedef struct _BlxStyle
  {
    char *styleName;
    BlxColor fillColor;
    BlxColor lineColor;
  } BlxStyle;


/* Structure holding information about an alignment */
typedef struct _MSP
{
  struct _MSP       *next;
  BlxMspType        type;          /* See enum above */
  int               score;
  int               id;

  char              *qname;        /* For Dotter, the MSP can belong to either sequence */
  char              qframe[8];     /* obsolete - use qFrame and qStrand instead */
  IntRange	    qRange;	   /* the range of coords on the ref sequence where the alignment lies */
  BlxStrand         qStrand;       /* which strand on the reference sequence the match is on */
  int               qFrame;        /* which frame on the reference sequence the match is on */

  BlxSequence       *sSequence;    /* pointer to a struct holding info about the sequence/strand this match is from */
  char              *sname;        /* sequence name (could be different to the sequence name in the blxSequence e.g. exons have a postfixed 'x') */
  char              sframe[8];     /* obsolete - use sStrand instead */
  IntRange	    sRange;	   /* the range of coords on the match sequence where the alignment lies */
  
  char              *desc;         /* Optional description text for the MSP */
  char              *source;       /* Optional source text for the MSP */
  GSList            *gaps;         /* Array of "gaps" in this homolgy (this is a bit of a misnomer because the array
                                    * gives the ranges of the bits that align rather than the ranges of the gaps in between */
                         
  BlxStyle          *style;        /* Specifies drawing style for this MSP, e.g. fill color and line color */

  FeatureSeries     *fs;           /* Feature series that this MSP belongs to */
  int               fsColor;       /* Color to draw this MSP in the feature series */
  BlxCurveShape     fsShape;       /* Shape data for drawing feature series curves, i.e. XY type PARTIAL or INTERPOLATE shapes */
  GArray            *xy;            /* For XY plot feature series */

  
#ifdef ACEDB
  KEY      key;
#endif
} MSP ;


/* Function to show blixem window, can be called from any application. */
gboolean                            blxview(char *refSeq, 
                                            char *refSeqName,
	                                    int start, 
                                            int qOffset, 
                                            MSP *msplist, 
                                            GList *seqList, 
                                            char *opts, 
	                                    PfetchParams *pfetch, 
                                            char *align_types, 
                                            gboolean External) ;


#endif /*  !defined DEF_BLXVIEW_H */
