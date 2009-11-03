/*  File: hexcode.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
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
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: hexcode.c,v 1.1 2009-11-03 18:28:23 edgrif Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 18 13:54 2000 (rd)
 * Created: Sun Aug 27 16:08:28 1995 (rd)
 *-------------------------------------------------------------------
 */

#ifdef ACEDB
/********* part of an ACEDB build **********/
#include "acedb.h"
#include "dict.h"
#include "dbpath.h"

#include "fmap.h"		/* for fMapAddGfSeg() */

/*-----------------------------------------------------------*/
#else				/* provide library functions */

/********** stand alone program ************/
#include <stdio.h>

#ifdef FALSE
  typedef int BOOL ;
#else
  typedef enum {FALSE=0,TRUE=1} BOOL ;
#endif

typedef int KEY ;

#define messalloc(z) malloc(z)
#define messfree(z) free(z)
#define messerror(x,y) fprintf(stderr,x,y)

/* simple version that just writes to stdout */
static void fMapAddGfSeg (int type, KEY key, int x1, int x2, float score);
#endif


/*-----------------------------------------------------------*/

static BOOL readTable (char *name, float **tab)
{
  int j, n = 0 ;
  FILE *fil ;
  int level ;
  char *filname ;

#ifdef ACEDB
  filname = dbPathFilName ("wgf", name, "hex", "r", 0) ;
  if (!filname)
    return FALSE ;
  if (!(fil = filopen (filname, 0, "r")))
    { messfree (filname) ; 
      return FALSE ;
    }
  level = freesetfile (fil, 0) ;
  
  if (!*tab)
    *tab = (float*) messalloc (4096*sizeof(float)) ;
  
  while (freecard (level))
    { for (j = 0 ; j < 16 && n < 4096 ; ++j)
        if (freefloat (&(*tab)[n]))
	  ++n ;
        else
	  { messerror ("short line %d in hexamer file %s", 
		       freestreamline (level), filname) ;
	    freeclose (level) ;
	    messfree (filname) ;
	    return FALSE ;
	  }
    }
#else
  filname = name ;
  if (!(fil = fopen (filname, "r")))
    return FALSE ;
  if (!*tab)
    *tab = (float*) messalloc (4096*sizeof(float)) ;
  for (j = 0 ; j < 4096 && !feof(fil) ; ++j)
    if (fscanf (fil, "%f", &(*tab)[n]))
      ++n ;
    else
      { messerror ("can't find entry %d in table file\n", j) ;
        return FALSE ;
      }
#endif

  if (n == 4096)
    {
#ifdef ACEDB
      messfree (filname) ;
#endif
      return TRUE ;
    }
  else
    { messerror ("problem reading hexamer file %s", filname) ;
#ifdef ACEDB
      messfree (filname) ;
#endif
      return FALSE ;
    }
}

/******************************/

				/* acedb base codes */
#define A_ 1
#define T_ 2
#define G_ 4
#define C_ 8

static void convertAceSequence (char *new, char *old, int len)
{
  while (len--)
    switch (*old++)
      { 
      case A_: *new++ = 0 ; break ;
      case C_: *new++ = 1 ; break ;
      case G_: *new++ = 2 ; break ;
      case T_: *new++ = 3 ; break ;
      default: *new++ = 1 ; break ; /* arbitrarily set to C */
      }
}

/*********************************************/
/********** form partial sum array ***********/

static float makePartial (char *s, int len, float *tab, int skip,
			  float *partial)
{
  register int i, j, index ;
  register float score = 0 ;

  if (len < 6)
    return 0 ;

  index = (s[0] << 10) + (s[1] << 8) + (s[2] << 6) +
		(s[3] << 4) + (s[4] << 2) + s[5] ;
  s += 6 ;
  for (i = 3 ; i <= len-3 ; i += skip)
    { score += tab[index] ;
      partial[i] = score ;
      for (j = skip ; j-- ;)
	index = (index << 2) + *s++ ;
      index &= 0xfff ;
    }

  return score ;
}

/********************************************************/
/***** find maximal segments, and add to look->segs *****/

static void processPartial (int type, KEY key, 
			    int step, float thresh, BOOL isRC,
			    int offset, float *partial, int len)
{ 
  register int i, k ;
  int loclen ;
  static int *maxes = 0, *mins = 0, slen = 0 ;

  loclen = len - offset ;
  while (loclen % step) --loclen ;
  partial += offset ;

  if (slen < len)
    { messfree (maxes) ; maxes = (int*) messalloc (len * sizeof(int)) ;
      messfree (mins) ; mins = (int*) messalloc (len * sizeof(int)) ;
      slen = len ;
    }

  k = 3 ;					/* make mins */
  for (i = 3 ; i <= loclen-3 ; i += step)
    { if (partial[i] < partial[k]) k = i ;
      mins[i] = k ;
    }
  k = loclen - 3 ;				/* make maxes */
  for (i = loclen - 3 ; i >= 3 ; i -= step)
    { if (partial[i] > partial[k]) k = i ;
      maxes[i] = k ;
    }

  for (i = 3 ; i <= loclen - 3 ; i += step)
    if (mins[maxes[i]] == i && 
	partial[maxes[i]] - partial[i] > thresh)
      if (isRC)
	fMapAddGfSeg (type, key, 
		      len-1 - maxes[i] - offset + 1, len-1 - i - offset,
		      partial[maxes[i]] - partial[i]) ;
      else
	fMapAddGfSeg (type, key, i+offset, maxes[i]+offset - 1, 
		      partial[maxes[i]] - partial[i]) ;
}

/****************************************************************/
/********* entry points: ACEDB first then stand alone ***********/

#ifdef ACEDB

static DICT *dict = 0 ;
static Array tabArray ;		/* of float* */

BOOL hexAddSegs (char *name, int type, KEY key, 
		 int step, float thresh, BOOL isRC,
		 char *dna, int len, float* partial)
{
  char *seq, c ;
  int i ;
  float *tab = 0 ;
  BOOL isLocalPartial = (partial == 0) ;

  if (!dict)
    { dict = dictCreate (8) ;
      tabArray = arrayCreate (8, float*) ;
    }

  if (dictAdd (dict, name, &i))
    { readTable (name, &tab) ;
      array(tabArray,i,float*) = tab ;
    }
  else
    tab = arr(tabArray,i,float*) ;

  if (!tab)
    return FALSE ;

  seq = (char*) messalloc (len) ;
  convertAceSequence (seq, dna, len) ;
  if (isRC)
    for (i = 0 ; i < len-1-i ; ++i)
      { c = 3 - seq[i] ;	/* NB "3 -" does complement */
        seq[i] = 3 - seq[len-1-i] ; 
	seq[len-1-i] = c ; 
      }

  if (isLocalPartial)
    partial = (float*) messalloc (sizeof(float)*len) ;

  for (i = 0 ; i < step ; ++i)
    makePartial (seq+i, len-i, tab, step, partial+i) ;
  for (i = 0 ; i < step ; ++i)
    processPartial (type, key, step, thresh, isRC, i, partial, len) ;

  messfree (seq) ;
  if (isLocalPartial) messfree (partial) ;

  return TRUE ;
}

/*-----------------------------------------------------------*/
#else  /* not ACEDB  */
/****************************************************************
** stand alone version
** writes a gff file to stdout
*****************************************************************/

#include "readseq.h"

static void usage (void)
{
  fprintf (stdout, "Usage: hexamer [opts] <tableFile> <seqFile>\n") ;
  fprintf (stdout, "options: -T <threshold>          0\n") ;
  fprintf (stdout, "         -F <feature name>       tableFile name\n") ;
  fprintf (stdout, "         -n	   flag for noncoding (no triplet frame)\n") ;
  exit (-1) ;
}

static char *seqName = "" ;
static char *tableName ;
static char *featName ;
static char strand ;
static char frame = '0' ;

/*****************************
** use type for BOOL isFrame
*****************************/

static void fMapAddGfSeg (int type, KEY key, int x1, int x2, float score)
     /* simple version for standalone program */
{
  printf ("%s\t%s\t%s\t%d\t%d\t%.4f\t%c\t%c\n", 
	  seqName, "hexamer", featName, x1+1, x2+1, 
	  score, strand, frame) ;
}

int main (int argc, char **argv)
{
  int i ;
  int thresh = 0.0 ;
  float *tab, *partial ;
  int step = 3 ;
  FILE *seqFile ;
  char *seq ;
  char c ;
  int len ;
  BOOL isRC ;
  extern float atof(char*) ;

  --argc ; ++argv ;		/* remove program name */

  while (argc > 2)
    if (!strcmp (*argv, "-T"))
      { thresh = atof (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-F"))
      { featName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-n"))
      { step = 1 ; frame = '.' ;
	argc -= 1 ; argv += 1 ;
      }
    else if (**argv == '-')
      { fprintf (stderr, "Unrecognised option %s\n", *argv) ;
	usage() ;
      }
    else
      usage() ;

  if (argc != 2)
    usage() ;

  tableName = *argv ; --argc ; ++argv ;
  if (!featName) featName = tableName ;
  if (!readTable (tableName, &tab))
    { fprintf (stderr, "Failed to open table file %s\n") ;
      usage () ;
    }

  if (!strcmp (*argv, "-"))
    seqFile = stdin ;
  else if (!(seqFile = fopen (*argv, "r")))
    { fprintf (stderr, "Failed to open sequence file %s\n") ;
      usage() ;
    }

  if (!readSequence (seqFile, dna2indexConv, 
		     &seq, &seqName, 0, &len))
    { fprintf (stderr, "Errors reading sequence %s\n") ;
      usage() ;
    }
  for (i = 0 ; i < len ; ++i)	/* remove N's - map to C for this */
    if (seq[i] > 3) seq[i] = 1 ;

  partial = (float*) messalloc (sizeof(float)*len) ;

				/* first do forward direction */
  isRC = FALSE ; strand = '+' ;
  for (i = 0 ; i < step ; ++i)
    makePartial (seq+i, len-i, tab, step, partial+i) ;

  for (i = 0 ; i < step ; ++i)
    processPartial (step, 0, step, thresh, isRC, i, partial, len) ;

				/* then reverse complement */
  isRC = TRUE ; strand = '-' ;
  for (i = 0 ; i < len-1-i ; ++i)
    { c = 3 - seq[i] ;	/* NB "3 -" does complement */
      seq[i] = 3 - seq[len-1-i] ; 
      seq[len-1-i] = c ; 
    }

  for (i = 0 ; i < step ; ++i)
    makePartial (seq+i, len-i, tab, step, partial+i) ;
  for (i = 0 ; i < step ; ++i)
    processPartial (step, 0, step, thresh, isRC, i, partial, len) ;

  messfree (seq) ;
  messfree (partial) ;
}

#endif /* not ACEDB */
/*-----------------------------------------------------------*/

/**************** end of file ****************/
