/*  File: dotter.h
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
 * Description: 
 * Exported functions:
 *   Only 3 parameters are mandatory, the rest can be set to NULL.
 *   A minimal call would look like:
 *
 *   dotter(type, 0, 0, qseq, 0, 0, sseq, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
 *
 *   NOTE: qseq and sseq must be g_malloc'ed in the calling routine.  
 *   They are g_free'd by Dotter.
 *
 * HISTORY:
 * Last edited: Nov 14 09:19 2007 (edgrif)
 * Created: Thu Aug 26 17:16:19 1999 (fw)
 * CVS info:   $Id: dotter.h,v 1.9 2010-11-08 15:52:49 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_DOTTER_H
#define DEF_DOTTER_H

#include <seqtoolsUtils/blxmsp.h>


/* Options specifying the initial state for dotter */
typedef struct _DotterOptions
  {
    int qoffset;              /* qoffset + 1 gives the value of the first coord on the ref seq */
    int soffset;              /* soffset + 1 gives the value of the first coord on the match seq */
    gboolean selfcall;        /* whether called internally, i.e. so that features/sequences will be piped into dotter rather than read from files */
    int qlen;                 /* length of the ref seq */
    int slen;                 /* length of the match seq */
    int dotterZoom;           /* initial zoom level */
    int install : 1;          /* whether to add -install to the dotter args (for private colormaps) */
    int pixelFacset;
    int seqInSFS;             /* whether the sequences are in the features file, i.e. there are no separate sequence files */
    
    float memoryLimit;
    
    char *savefile;           /* file to save the dot-plot to */
    char *loadfile;           /* file to load a dot-plot from */
    char *FSfilename;         /* file containing features i.e. MSPs */
    char *mtxfile;            /* caller-supplied matrix file */
    
    char *winsize;            /* caller-supplied sliding-window size */
    char *xOptions;           /* x-options */
    
    char *qname;              /* reference (horizontal) sequence name */
    char *qseq;               /* reference (horizontal) sequence data */
    char *sname;              /* match (vertical) sequence name */
    char *sseq;               /* match (vertical) sequence data */
    
    gboolean mirrorImage;     /* display mirror image in self comparisons (i.e. so we only have to calculate half of the dot-plot) */
    gboolean watsonOnly;      /* only show the watson (forward) strand of the ref seq */
    gboolean crickOnly;       /* only show the crick (reverse) strand of the ref seq */
    gboolean hspsOnly;        /* only draw HSPs (i.e. don't calculate the dot-plot, just draw lines where we know HSPs should be) */
    gboolean swapGreyramp;    /* swap the default black/white points on the greyramp tool (inverts the colors) */
    gboolean fsEndLinesOn;    /* to do: not used? */
    gboolean hozScaleRev;     /* revese the horizontal scale */
    gboolean vertScaleRev;    /* revese the vertical scale */
  } DotterOptions;


void dotter(
	const BlxBlastMode blastMode, /* Mandatory, one of { BLASTP, BLASTN, BLASTX } 
			      P -> Protein-Protein
			      N -> DNA-DNA
			      X -> DNA-Protein */
	
	DotterOptions *options, /* Optional, may be NULL 
                                   Various options for display features */

	const char *queryname,   /* Optional, may be NULL 
				    Name of Horizontal sequence */
	
	char *queryseq,	   /* Mandatory, NULL terminated string
			      Horisontal sequence - g_free'd by Dotter */

	int   qoff,	   /* Optional, may be NULL
			      Coordinate offset of horisontal sequence */

	const BlxStrand refSeqStrand, /* which strand of the reference sequence was passed */

	const char *subjectname, /* Optional, may be NULL 
			      Name of vertical sequence */

	char *subjectseq,  /* Mandatory, NULL terminated string
			      vertical sequence - g_free'd by Dotter */

	int   soff,	   /* Optional, may be NULL 
			      Coordinate offset of horisontal sequence */

	const BlxStrand matchSeqStrand, /* which strand of the match sequence was passed */

	int   qcenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	int   scenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	char *savefile,	   /* Optional, may be NULL 
			      Filename to save dotplot to. Invokes batch mode */

	char *loadfile,	   /* Optional, may be NULL 
			      Filename to load dotplot from */

	char *mtxfile,	   /* Optional, may be NULL 
			      Filename to load score matrix from */

	double memoryLimit, /* Optional, may be NULL 
			      Maximum Mb allowed for dotplot */

	int   zoomFac,	   /* Optional, may be NULL
			      Compression of dotplot {1, 2, 3, ... }
			      Automatically calculated if NULL */

	MSP *MSPs,	   /* Optional, may be NULL
			      List of MSPs containing genes and blast matches */

	GList *seqList,	   /* Optional, may be NULL
			      List of all match sequences, as BlxSequences */
	    
	int   MSPoff,	   /* Optional, may be NULL
			      Coordinate offset of MSPs */

	char *winsize,	   /* String determining the window size */

	int   pixelFacset  /* Preset pixel factor */
);


#endif /*  !defined DEF_DOTTER_H */