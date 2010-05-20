/*  File: blxparser.c
 *  Author: Erik Sonnhammer
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
 * Description: blxparser - parses MSPcrunch output for blixem
 * Exported functions: See blxview.h
 * HISTORY:
 * Last edited: Sep 29 17:00 2009 (edgrif)
 * * 94-01-29  Added SEQBL format support.
 * * 94-01-30  Added proper option parsing and pipe'ability.
 * * 94-03-27  Added (fudged) Tblastn support.
 * * 95-02-06  Added (fudged) Tblastx support.
 * *           Added autmatic mode detection.
 * * 95-06-21  Query sequence parsing that allows spaces and '*'.
 * * 98-02-19  Changed MSP parsing to handle all SFS formats.
 * * 99-07-29  Added support for SFS type=HSP and GFF.
 * Created: 93-05-17
 * CVS info:   $Id: blxparser.c,v 1.20 2010-05-20 15:25:59 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/graph.h>
#include <wh/gex.h>
#include <wh/smap.h>					    /* for SMapMap struct. */
#include <SeqTools/utilities.h>

static char *	    nextLine(FILE *file, GString *line_string);

static void	    parseEXBLXSEQBL(MSP *msp, const BlxMSPType msp_type, char *opts, GString *line_string);
static void	    parseEXBLXSEQBLExtended(MSP *msp, const BlxMSPType msp_type, char *opts, GString *line_string);
static BOOL	    parseGaps(char **text, MSP *msp) ;
static BOOL	    parseDescription(char **text, MSP *msp_unused) ;
static BOOL	    parseSequence(char **text, MSP *msp) ;

static int	    fsorder(void *x_in, void *y_in) ;
static int	    fsnrorder(void *x_in, void *y_in) ;
static void	    getDesc(MSP *msp, char *s1, char *s2) ;
static void	    prepSeq(MSP *msp, char *seq, char *opts) ;
static void	    parseLook(MSP *msp, char *s) ;

static gboolean	    parseNonMspData(char *line, MSP *msp, char **seq1, char *seq1name, char **seq2, char *seq2name, 
				    BlxMSPType *mspType, char **readSeq, int *readSeqLen, int *readSeqMaxLen);

static gboolean	    parseMspData(char *line, char *opts, MSP *msp, GString *line_string, char **seq1, 
				 char *seq1name, char **seq2, char *seq2name, BlxMSPType *mspType);

static void	    parseXY(char *line, char *opts, MSP *msp, char **seq1, char *seq1name, char **seq2, char *seq2name, BlxMSPType *mspType);
static void	    parseGff(char *line, char *opts, MSP *msp);
static void	    parseFsseg(char *line, char *opts, MSP *msp);
static void	    parseGsp(char *line, char *opts, MSP *msp, BlxMSPType *mspType);
static void	    parseHsp(char *line, char *opts, MSP *msp);
static gboolean	    parseHeaderInfo(char *line, char *opts, MSP *msp, BlxMSPType *mspType);
static MSP*	    createEmptyMsp();

/*
 *  Globals
 */

static char *colorNames[NUM_TRUECOLORS] =
  {
    "WHITE", 
    "BLACK", 
    "LIGHTGRAY", 
    "DARKGRAY",
    "RED", 
    "GREEN", 
    "BLUE",
    "YELLOW", 
    "CYAN", 
    "MAGENTA",
    "LIGHTRED", 
    "LIGHTGREEN", 
    "LIGHTBLUE",
    "DARKRED", 
    "DARKGREEN", 
    "DARKBLUE",
    "PALERED", 
    "PALEGREEN", 
    "PALEBLUE",
    "PALEYELLOW", 
    "PALECYAN", 
    "PALEMAGENTA",
    "BROWN", 
    "ORANGE", 
    "PALEORANGE",
    "PURPLE", 
    "VIOLET", 
    "PALEVIOLET",
    "GRAY", 
    "PALEGRAY",
    "CERISE", 
    "MIDBLUE"
};


/* Try and get rid of these.... */
static char sname[MAXLINE+1], seq[MAXLINE+1] ;


/* Function: parse a stream of SFS data
 *
 * Assumptions:
 *
 *  *seq1 and *seq2 may or may not be allocated.
 *
 *  *seq1name and *seq2name are allocated at least 255 bytes
 *
 */
void parseFS(MSP **MSPlist, FILE *file, char *opts,
	     char **seq1, char *seq1name, char **seq2, char *seq2name, const int qOffset)
{
  if (!fsArr) 
    fsArr = arrayCreate(50, FEATURESERIES);
  else
    arraySort(fsArr, fsorder);

  /* Find the last MSP in the list */
  MSP *msp = NULL;
  if (*MSPlist)
    {
      msp = *MSPlist;  
      while(msp->next)
	msp = msp->next;
    }

  /* Allocate reusable/extendable string as our buffer..*/
  GString *line_string = g_string_sized_new(MAXLINE + 1); 

  BlxMSPType mspType = FSSEG;
  
  /* Variables for reading in SEQ data */
  char **readSeq;	/* buffer to hold the sequence we're reading in */
  int readSeqMaxLen = 0;    /* current max length of the buffer */
  int readSeqLen = 0;   /* current end pos of the data in the buffer */
  
  while (!feof(file))
    {
      line_string = g_string_truncate(line_string, 0) ;	    /* Reset buffer pointer. */

      char *line = nextLine(file, line_string);
      if (!line)
	{
	  break;
	}

      const int lineLen = strlen(line);
      if (lineLen == 0)
	{
	  continue; /* empty file??? */
	}
      
      /* get rid of any trailing '\n', there may not be one if the last line of the file
       * didn't have one. */
      char *charPtr = strchr(line, '\n');
      if (charPtr)
	{
	  *charPtr = 0;
	}

      /* Header info */
      if (parseHeaderInfo(line, opts, msp, &mspType))
	{
	  continue; 
	}
	
      /* Data that don't go into an MSP */
      if (parseNonMspData(line, msp, seq1, seq1name, seq2, seq2name, &mspType, readSeq, &readSeqLen, &readSeqMaxLen))
	{
	  continue; 
	}

      /* Allocate a new MSP but only if this is a new record. */
      /* If it's a chunk of a dna sequence, leave as is.      */
      if (!*MSPlist) 
	{
	  msp = *MSPlist = createEmptyMsp();
	}
      else 
	{
	  if (strcspn(line, "acgt") != 0)	/* ie line contains chars other than acgt */
	    {
	      msp->next = createEmptyMsp();
	      msp = msp->next;
	    }
	}
	
      /* Data that do go into a new MSP */
      if (parseMspData(line, opts, msp, line_string, seq1, seq1name, seq2, seq2name, &mspType))
	{
	  continue;
	}
    }

  g_string_free(line_string, TRUE) ;			    /* free everything, buffer and all. */

  /* Sort feature segment array by number */
  arraySort(fsArr, fsnrorder);

  return ;
}




void insertFS(MSP *msp, char *series)
{
    FEATURESERIES fs;
    int i;
    static int curnr=0;

    if (!fsArr) 
	fsArr = arrayCreate(50, FEATURESERIES);

    fs.on = 1;
    fs.y = 0.0;
    fs.xy = (msp->type == XY ? 1 : 0);

    fs.name = g_malloc(strlen(series)+1);
    strcpy(fs.name, series);

    if (arrayFind(fsArr, &fs, &i, fsorder)) {
	msp->fs = arrp(fsArr, i, FEATURESERIES)->nr;
	g_free(fs.name);
    }
    else {
	msp->fs = fs.nr = curnr++;
	arrayInsert(fsArr, &fs, fsorder);
    }
}





char *readFastaSeq(FILE *seqfile, char *qname)
{
    char  line[MAXLINE+1];
    char *q, *c, *cp, *cq, ch ;

    /* Read in seqfile */
    if (seqfile == stdin) {
	Array arr = arrayCreate(5000, char);
	int i=0;

	if (!fgets(line, MAXLINE, seqfile))
	  g_error("Error reading seqFile") ;

	sscanf(line, "%s", qname);

	while ((ch = fgetc(seqfile)) != '\n') {
	    if (isalpha(ch) || ch == SEQUENCE_CHAR_GAP || ch == SEQUENCE_CHAR_STOP) 
		array(arr, i++, char) = ch;
	}
	q = g_malloc(arrayMax(arr)+1);
	cq = q;
	for (i = 0; i < arrayMax(arr);) *cq++ = arr(arr, i++, char);
	arrayDestroy(arr);
    }
    else {
	fseek(seqfile, 0, SEEK_END);
	q = g_malloc(ftell(seqfile)+1);
	cq = q;
	fseek(seqfile, 0, SEEK_SET);
	while (!feof(seqfile))
	{
	    if (!fgets(line, MAXLINE, seqfile)) break;
	    while (strchr(line, '>'))
	      {
		strncpy(qname, line+1, 255); qname[255]=0;
		if ((c = (char *)strchr(qname, ' ')))
		  *c = 0;
		if ((c = (char *)strchr(qname, '\n')))
		  *c = 0;
		if (!fgets(line, MAXLINE, seqfile))
		  break;
	    }
	    
	    for (cp = line; *cp; cp++)
	      {
		if (isalpha(*cp) || *cp == SEQUENCE_CHAR_STOP)
		  {
		    *cq++ = *cp;
		  }
	      }
	  
	    /* Allow stops in query sequence */
	    *cq = 0;
	}
    }

    return q;
}



/* 
 *               Internal functions.
 */


static int parseColor(char *s) 
{
    int i;

    for (i = 0; i < NUM_TRUECOLORS; i++) {
	if (!strcasecmp(colorNames[i], s)) break;
    }
    return i;
}

static int parseShape(char *s) 
{
    if (!strcasecmp(s, "interpolate")) 
	return XY_INTERPOLATE;
    else if (!strcasecmp(s, "partial")) 
	return XY_PARTIAL;
    else 
	return XY_BADSHAPE;
}

static void parseLook(MSP *msp, char *s)
{
    char *cp, *s2;

    /* Make copy of string to mess up */
    s2 = g_malloc(strlen(s)+1);
    strcpy(s2, s);

    cp = strtok(s2, "," );
    while (cp) {
	
	if (parseColor(cp) != NUM_TRUECOLORS) {
	    msp->color = parseColor(cp);
	}
	else if (parseShape(cp) != XY_BADSHAPE) {
	    msp->shape = parseShape(cp);
	}
	else 
	    messout("Unrecognised Look: %s", cp);
	
	cp = strtok(0, "," );
    }
    g_free(s2);
}


/* Copy 'remainder of s1 after word s2' into msp->desc  */
static void getDesc(MSP *msp, char *s1, char *s2)
{
    char *cp;
    
    if (!(cp = strstr(s1, s2))) {
	messout("Can't find back %s in %s", s2, s1);
	return;
    }

    cp += strlen(s2)+1;

    msp->desc = g_malloc(strlen(cp)+1);
    strcpy(msp->desc, cp);
}



/* Check if we have a reversed subject and, if so, if this is allowed. Throws an error if not. */
static void CheckReversedSubjectAllowed(const MSP *msp, const char *opts)
{
  if (mspGetMatchStrand(msp) == BLXSTRAND_REVERSE && *opts != 'T' && *opts != 'L' && *opts != 'N')
    {
      g_error("Reversed subjects are not allowed in modes blastp or blastx");
    }
}



static void prepSeq(MSP *msp, char *seq, char *opts)
{
  if (*opts == 'T' || *opts == 'L') 
    {
      msp->sseq = g_malloc(strlen(seq)+1) ;
      strcpy(msp->sseq, seq) ;
    }
  else 
    {
      CheckReversedSubjectAllowed(msp, opts);

      /* I think this is a mistake in the code, we get the whole sequence so shouldn't be
       * prepending the '-'s */
      msp->sseq = g_malloc(msp->sstart+strlen(seq)+1);
      memset(msp->sseq, SEQUENCE_CHAR_PAD, msp->sstart); /* Fill up with dashes */
      strcpy(msp->sseq + msp->sstart - 1, seq);
    }

  return ;
}



/* This routine parses MSP files that have either the exblx or seqbl format.
 * 
 * Format for both exblx and seqbl files is seven tab separated fields:
 * 
 * score reference_strand_and_frame reference_start reference_end match_start match_end match_name
 * 
 * For exblx this is optionally followed by:
 * 
 *             [match_description]
 * 
 * and optionally for seqbl by:
 * 
 *             [match_sequence]
 * 
 * e.g. for exblx
 * 
 *    11 (+2) 49052 49783 102 328 SW:YCF2_MARPO some description or other
 * 
 */
static void parseEXBLXSEQBL(MSP *msp, BlxMSPType msp_type, char *opts, GString *line_string)
{
  int   i, qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL ;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%d%d%s", 
	     &msp->score, msp->qframe, &msp->qstart, &msp->qend, 
	     &msp->sstart, &msp->send, sname) != 7)
    {
      g_critical("Incomplete MSP data");
      abort() ;
    }

  msp->type = msp_type ;
  
  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (*opts == 'T')
    strcpy(msp->qframe, "(+1)");
  
  msp->sname = g_malloc(strlen(sname)+1);
  strcpy(msp->sname, sname); 
      
  /* Convert to upper case (necessary?) */
  for (i=0; msp->sname[i]; i++)
    {
      msp->sname[i] = freeupper(msp->sname[i]);
    }
      
  /* Convert subject names to fetchable ones if from NCBI server 
         
  Rule 1: If there is a gi, use that.
  Rule 2: If no gi, use the first and last non-blank field as db:id.
  */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_malloc(strlen(msp->sname)+1);
      strcpy(src, msp->sname);
  	
      p = strtok(src, "|");
  	
      if (!strcasecmp(p, "GI"))
	{
	  /* Use only GI number */
	  p = strtok(0, "|");
	  strcpy(msp->sname, "gi");
	  strcat(msp->sname, ":");
	  strcat(msp->sname, p);
	}
      else
	{
	  /* Try to make a proper db:name.  Use last non-blank field */
	  char *db=p, *last;
  	    
	  p = strtok(0, "|");
	  while (p) {
	    if (*p && *p != ' ') last = p;
	    p = strtok(0, "|");
	  }
	  strcpy(msp->sname, db);
	  strcat(msp->sname, ":");
	  strcat(msp->sname, last);
	}
  	
      g_free(src);
    }
  
  if (!(cp = strstr(line, sname)))
    {
      g_error("Line does not include %s", sname);
    }
      
  seq_pos = cp + strlen(sname) ;
            
  qlen = abs(msp->qend - msp->qstart)+1;
  slen = abs(msp->send - msp->sstart)+1;


  if (*opts == ' ')
    {
      /* Guess mode from prefix or from coordinates */

      if (qlen == slen)
	{
	  /* Could be blastp, blastn or tblastx */
	  if (strchr(msp->sname, '_'))
	    *opts = 'P';
	  else if (strstr(msp->sname, "PIR:"))
	    *opts = 'P';
	  else if (strstr(msp->sname, "EM:"))
	    *opts = 'N';
	  else
	    *opts = 'N';				    /* Could be P or L just as well */
	}
      else if (qlen > slen)
	*opts = 'X';
      else if (qlen < slen)
	*opts = 'T';
    }


  if (msp_type == EXBLX)
    {
      opts[5] = ' '; /* Don't use full zoom default */
      
      /* skip over description */
      while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	seq_pos++;
      if (*seq_pos && !isdigit(*seq_pos))
	{
	  while (*seq_pos && *seq_pos != ' ' && *seq_pos != '\t')
	    seq_pos++;
	  while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	    seq_pos++;
	}
    }
  else if (msp_type == SEQBL)
    {
      CheckReversedSubjectAllowed(msp, opts);
    
      if (*opts == 'L')
	{
	  slen = (abs(msp->send - msp->sstart) + 1)/3;
	}
	

      /* Line contains chars other than sequence so get the starter data...not sure this test
       * is necessary any more now we have a better mechanism of getting a whole line 
       * from a file. All this is a horrible mixture of strtok and sscanf but what else
       * to do.... */
      if (strcspn(line, "acgt") != 0)
	{
	  char *seq ;

	  seq = g_malloc(line_string->len + 1) ;
	  
	  if (sscanf(seq_pos, "%s", seq) != 1)
	    {
	      g_error("Error parsing %s", line);
	    }
	  blxSeq2MSP(msp, seq) ;


	   while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	     seq_pos++;
	   while (*seq_pos && *seq_pos != ' ' && *seq_pos != '\t')
	     seq_pos++;
	  while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	    seq_pos++;
	}
    }
  
  if (*opts == 'N')
    opts[3] = 'R';
  
  return;
}


/* This routine parses MSP files that have either the _extended_ exblx or
 * seqbl formats which borrow from the GFF format in having a set number
 * of fixed fields followed by a variable number of attribute fields.
 * These formats are distinguished by having exblx_x or seqbl_x in the
 * file comment header:
 *
 * "# seqbl_x"  or  "# exblx_x"
 *
 * Format for both exblx_x and seqbl_x files is eight tab separated fields:
 * 
 * score reference_strand_frame reference_start reference_end match_start match_end match_strand match_name
 * 
 * For exblx_x this is optionally followed by:
 * 
 *             [gaps data] [match_description]
 * 
 * and for seqbl_x by:
 * 
 *             [gaps data] [match_sequence]
 * 
 * The format of these extra fields is "tag value(s) ;", i.e.
 * 
 * "Gaps [ref_start ref_end match_start match_end]+ ;"
 * 
 * "Description the sequence description ;"
 * 
 * "Sequence aaagggtttttcccccc ;"
 * 
 * Currently acedb & ZMap export this format but other programs could also use it.
 * 
 */
static void parseEXBLXSEQBLExtended(MSP *msp, BlxMSPType msp_type, char *opts, GString *line_string)
{
  BOOL result = FALSE ;
  int   i, qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL, *first_pos ;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%s%d%d%s", 
	     &msp->score,
	     msp->qframe, &msp->qstart, &msp->qend, 
	     msp->sframe, &msp->sstart, &msp->send, sname) != 8)
    {
      g_critical("Incomplete MSP data");
      abort() ;
    }

  msp->type = msp_type ;
  
  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (*opts == 'T')
    strcpy(msp->qframe, "(+1)");
  
  msp->sname = g_malloc(strlen(sname)+1);
  strcpy(msp->sname, sname); 
      
  /* Convert to upper case (necessary?) */
  for (i=0; msp->sname[i]; i++)
    {
      msp->sname[i] = freeupper(msp->sname[i]);
    }
      
  /* Convert subject names to fetchable ones if from NCBI server 
   *  Rule 1: If there is a gi, use that.
   *  Rule 2: If no gi, use the first and last non-blank field as db:id.
   */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_malloc(strlen(msp->sname)+1);
      strcpy(src, msp->sname);
  	
      p = strtok(src, "|");
  	
      if (!strcasecmp(p, "GI"))
	{
	  /* Use only GI number */
	  p = strtok(0, "|");
	  strcpy(msp->sname, "gi");
	  strcat(msp->sname, ":");
	  strcat(msp->sname, p);
	}
      else
	{
	  /* Try to make a proper db:name.  Use last non-blank field */
	  char *db=p, *last;
  	    
	  p = strtok(0, "|");
	  while (p) {
	    if (*p && *p != ' ') last = p;
	    p = strtok(0, "|");
	  }
	  strcpy(msp->sname, db);
	  strcat(msp->sname, ":");
	  strcat(msp->sname, last);
	}
  	
      g_free(src);
    }
  
  if (!(cp = strstr(line, sname)))
    {
      g_error("Line does not include %s", sname);
    }
      
  seq_pos = cp + strlen(sname) ;
            
  qlen = abs(msp->qend - msp->qstart)+1;
  slen = abs(msp->send - msp->sstart)+1;

  if (*opts == ' ')
    {
      /* Guess mode from prefix or from coordinates */

      if (qlen == slen)
	{
	  /* Could be blastp, blastn or tblastx */
	  if (strchr(msp->sname, '_'))
	    *opts = 'P';
	  else if (strstr(msp->sname, "PIR:"))
	    *opts = 'P';
	  else if (strstr(msp->sname, "EM:"))
	    *opts = 'N';
	  else
	    *opts = 'N';				    /* Could be P or L just as well */
	}
      else if (qlen > slen)
	*opts = 'X';
      else if (qlen < slen)
	*opts = 'T';
    }



  /* Now read attributes. */
  if (seq_pos && *seq_pos)
    {
      first_pos = seq_pos ;
      result = TRUE ;
      while (result)
	{
	  if (!(seq_pos = strtok(first_pos, " ")))
	    {
	      result = FALSE ;
	      break ;
	    }
	  else
	    {
	      first_pos = NULL ;

	      /* Get first word and then parse.... */
	      if ((strstr(seq_pos, BLX_GAPS_TAG)))
		{
		  if (!(result = parseGaps(&seq_pos, msp)))
		    g_error("Incomplete MSP gap data") ;
		}
	      else if (msp_type == EXBLX && (strstr(seq_pos, BLX_DESCRIPTION_TAG)))
		{
		  opts[5] = ' '; /* Don't use full zoom default */

		  if (!(result = parseDescription(&seq_pos, msp)))
		    g_error("Bad description data") ;
		}
	      else if ((msp_type == SEQBL_X || mspIsSnp(msp)) && (strstr(seq_pos, BLX_SEQUENCE_TAG)))
		{
                  CheckReversedSubjectAllowed(msp, opts);
                
		  if (*opts == 'L')
		    {
		      slen = (abs(msp->send - msp->sstart) + 1)/3;
		    }
	
		  if (!(result = parseSequence(&seq_pos, msp)))
		    g_error("Bad sequence data") ;
		}

	    }
	}
    }

  if (msp_type == SEQBL_X)
    {
      CheckReversedSubjectAllowed(msp, opts);
    
      if (*opts == 'L')
	{
	  slen = (abs(msp->send - msp->sstart) + 1)/3;
	}
	

    }

  
  if (*opts == 'N')
    opts[3] = 'R';
  
  return ;
}


static int fsnrorder(void *x_in, void *y_in)
{
  FEATURESERIES *x = (FEATURESERIES *)x_in, *y = (FEATURESERIES *)y_in ;

  if (x->nr < y->nr)
    return -1;
  else if (x->nr > y->nr)
    return 1;
  else return 0;
}

static int fsorder(void *x_in, void *y_in)
{
  FEATURESERIES *x = (FEATURESERIES *)x_in, *y = (FEATURESERIES *)y_in ;

  /*printf("%s - %s : %d\n", x->name, y->name,  strcmp(x->name, y->name));*/

  return strcmp(x->name, y->name);
}




/* Read a line from a file, gets the whole line no matter how big...until you run out
 * of memory.
 * Returns NULL if there was nothing to read from the file, otherwise returns the
 * line read which is actually the string held within the GString passed in.
 * Crashes if there is a problem with the file. */
static char *nextLine(FILE *file, GString *line_string)
{
  enum {BLX_BUF_SIZE = 4096} ;				    /* Vague guess at initial length. */
  char *result = NULL ;
  BOOL line_finished ;
  char buffer[BLX_BUF_SIZE] ;

  messAssert(file) ;

  line_finished = FALSE ;
  while (!line_finished)
    {

      if (!fgets(buffer, BLX_BUF_SIZE, file))
	{
	  if (feof(file))
	    line_finished = TRUE ;
	  else
	    g_error("NULL value returned on reading input file") ;
	}
      else
	{
	  if (buffer[0] == '\0')
	    {
	      line_finished = TRUE ;
	      result = line_string->str ;
	    }
	  else
	    {
	      int line_end ;

	      g_string_sprintfa(line_string, "%s", &buffer[0]) ;

	      line_end = strlen(&buffer[0]) - 1 ;

	      if (buffer[line_end] == '\n')
		{
		  line_finished = TRUE ;
		  result = line_string->str ;
		}
	    }
	}
    }

  return result ;
}


/* Expects a string in the format "Gaps [ref_start ref_end match_start match_end]+ ; more text....."
 * and parses out the coords. Only spaces allowed in the string, not tabs. Moves text
 * to first char after ';'. */
static BOOL parseGaps(char **text, MSP *msp)
{
  BOOL result = TRUE ;
  SMapMap *gap ;
  char *next_gap ;
  BOOL end ;

  msp->gaps = arrayCreate(10, SMapMap) ;	    

  end = FALSE ;
  next_gap = strtok(NULL, " ") ;
  while (result && !end && next_gap)
    {
      int i ;

      for (i = 0 ;  i < 4 ; i++)
	{
	  if (!next_gap)
	    {
	      result = FALSE ;
	      break ;
	    }
	  else if (*next_gap == ';')
	    {
	      end = TRUE ;

	      break ;
	    }
	     
	  switch (i)
	  {
	    case 0:
	    {
	      /* First value is start of subject sequence range */
	      gap = arrayp(msp->gaps, arrayMax(msp->gaps), SMapMap) ;
	      gap->s1 = atoi(next_gap);
	      break;
	    }
	      
	    case 1:
	    {
	      /* Second value is end of subject sequence range. Order values so that
	       * s1 is less than s2 if we have the forward strand or v.v. if the reverse. */
	      gap->s2 = atoi(next_gap);
	      sortValues(&gap->s1, &gap->s2, mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
	      break;
	    }
	      
	    case 2:
	    {
	      /* Third value is start of reference sequence range */
	      gap->r1 = atoi(next_gap);
	      break;
	    }
	      
	    case 3:
	    {
	      /* Fourth value is end of reference sequence range. Order values so that
	       * r1 is less than r2 if ref sequence is forward strand or v.v. if the reverse. */
	      gap->r2 = atoi(next_gap);
	      sortValues(&gap->r1, &gap->r2, mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
	    }
	  }

	  next_gap = strtok(NULL, "\t ") ; 
	}  
    }

  return result ;
}


/* Description, should just be plain text, format is:
 * 
 * "Description text ; more text....."
 * 
 * text following Description must not contain ';'. Moves text to first
 * char after ';'. */
static BOOL parseDescription(char **text, MSP *msp)
{
  BOOL result = TRUE ;
  char *cp ;

  cp = strtok(NULL, ";") ;

  if (cp && *cp)
    msp->desc = strnew(cp, 0) ;

  return result ;
}


static BOOL parseSequence(char **text, MSP *msp)
{
  BOOL result = FALSE ;
  char *cp = *text ;

  cp = strtok(NULL, "\t ") ;				    /* skip "Sequence" */

  size_t origLen = strlen(cp);
  size_t validLen = strspn(cp, "acgtnACGTN"); 
  
  if (validLen < 1 || validLen < strlen(cp))
    {
      g_error("Error parsing %s, coords are %d -> %d (%d), but valid length sequence is only %d (out of total length supplied = %d).",
		msp->sname, msp->sstart, msp->send, msp->send - msp->sstart, (int)validLen, (int)origLen) ;
    }
  else
    {
      char *seq ;

      seq = g_malloc(validLen + 1) ;
	  
      if (sscanf(cp, "%s", seq) != 1)
	{
	  g_error("Error parsing %s", cp) ;
	}

      blxSeq2MSP(msp, seq) ;

      result = TRUE ;
    }

  return result ;
}



/* Allocate memory for an MSP and initialise all its fields to relevant 'empty' values */
static MSP* createEmptyMsp()
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  msp->next = NULL;
  msp->type = BLX_MSP_INVALID;
  msp->score = 0;
  msp->id = 0;
  
  msp->qname = NULL;
  msp->qframe[0] = 0;
  msp->qstart = 0;
  msp->qend = 0;
  
  msp->sSequence = NULL;
  msp->sname = NULL;
  msp->sseq = NULL;
  
  msp->sframe[0] = 0;
  msp->slength = 0;
  msp->sstart = 0;
  msp->send = 0;
  
  msp->desc = NULL;
  msp->box = 0;
  
  msp->color = 0;
  msp->shape = 0;
  msp->fs = 0;
  
  msp->xy = NULL;
  msp->gaps = NULL;
  
#ifdef ACEDB
  msp->key = 0;
#endif
  
  return msp;
}


/* Utility called by parseFS to parse the header info of a line from a file. Returns true
 * if the line was processed, false if further processing is required */
static gboolean parseHeaderInfo(char *line, char *opts, MSP *msp, BlxMSPType *mspType)
{
  gboolean processed = FALSE;
  
  if (!strncasecmp(line, "# seqbl_x", 9))
    {
      *mspType = SEQBL_X ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx_x", 9))
    {
      *mspType = EXBLX_X ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# seqbl", 7))
    {
      /* Only for backwards compatibility */
      *mspType = SEQBL ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx", 7))
    {
      /* Only for backwards compatibility */
      *mspType = EXBLX ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=HSP", 13) || 
	   !strncasecmp(line, "# SFS type=HSP", 14))
    {
      *mspType = HSP;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GSP", 13) || 
	   !strncasecmp(line, "# SFS type=GSP", 14))
    {
      *mspType = GSP;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# dotter feature format 2", 25) ||
	   !strncasecmp(line, "# FS type=SEG", 13) ||
	   !strncasecmp(line, "# SFS type=SEG", 14))
    {
      *mspType = FSSEG;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GFF", 13) || 
	   !strncasecmp(line, "# SFS type=GFF", 14))
    {
      *mspType = GFF;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=XY", 12) ||
	   !strncasecmp(line, "# SFS type=XY", 13))
    {
      *mspType = XY;
    }
  else if (!strncasecmp(line, "# FS type=SEQ", 13) ||
	   !strncasecmp(line, "# SFS type=SEQ", 14))
    {
      *mspType = SEQ;
    }
  else if (!strncasecmp(line, "# FS type=", 10) ||
	   !strncasecmp(line, "# SFS type=", 11))
    {
      g_error("Unrecognised SFS type: %s\n", line);
    }
  else if (*line == '#')
    {
      /* Very ugly; only for backwards compatibility */
      /* Changed to soft parsing (unknown labels ignored) so that
       any comment can be used */
      if (!strncasecmp(line, "# blastp" , 8))
	*opts = 'P';
      else if (!strncasecmp(line, "# tblastn", 9))
	*opts = 'T';
      else if (!strncasecmp(line, "# tblastx", 9))
	*opts = 'L';
      else if (!strncasecmp(line, "# blastn" , 8))
	*opts = 'N';
      else if (!strncasecmp(line, "# blastx" , 8))
	*opts = 'X';
      else if (!strncasecmp(line, "# hspgaps", 9))
	{
	  opts[7] = 'G';
	}
      else if (!strncasecmp(line, "# DESC ", 7) &&
	       (*mspType == HSP || *mspType == GSP || *mspType == SEQBL))
	{
	  if (msp)
	    getDesc(msp, line, sname);
	}
      
      /* move on to the next line. */
      processed = TRUE ;
    }
  
  return processed ;
}


/* Parse data from the given line into an MSP of type HSP */
static void parseHsp(char *line, char *opts, MSP *msp)
{
  msp->type = HSP;
  char qname[MAXLINE+1];
  
  /* <score> <qname> <qframe> <qstart> <qend> <sname> <sframe> <sstart> <ssend> <sequence> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s%s%d%d%s", 
	     &msp->score, 
	     qname, msp->qframe+1, &msp->qstart, &msp->qend, 
	     sname, msp->sframe+1, &msp->sstart, &msp->send,
	     seq) != 10)
    {
      g_critical("Error parsing data, type HSP: \"%s\"\n", line);
      abort();
    }
  
  msp->qname = g_malloc(strlen(qname)+1);
  strcpy(msp->qname, qname);
  msp->sname = g_malloc(strlen(sname)+1);
  strcpy(msp->sname, sname);
  
  *msp->qframe = *msp->sframe = '(';
  msp->qframe[3] = msp->sframe[3] = ')'; /* Too lazy to change code... */
  
  prepSeq(msp, seq, opts);
  
}


/* Parse data from the given line into an MSP of type GSP */
static void parseGsp(char *line, char *opts, MSP *msp, BlxMSPType *mspType)
{
  msp->type = GSP;
  
  /* Will write this as soon as MSPcrunch generates it */
  
  /* Start parsing the GSP data next */
  *mspType = GSPdata ;
}


/* Parse data from the given line into an MSP of type FFSEG */
static void parseFsseg(char *line, char *opts, MSP *msp)
{
  msp->type = FSSEG;
  char series[MAXLINE+1];
  char qname[MAXLINE+1];
  char look[MAXLINE+1];
  
  /* <score> <sequencename> <seriesname> <start> <end> <look> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s",
	     &msp->score, qname, series, &msp->qstart, &msp->qend, look) != 6)
    {
      g_critical("Error parsing data, type FSSEG: \"%s\"\n", line);
      abort();
    }
  
  msp->sstart = msp->qstart;
  msp->send = msp->qend;
  
  msp->qname = g_malloc(strlen(qname)+1);
  strcpy(msp->qname, qname);
  
  msp->sname = g_malloc(strlen(series)+1); 
  strcpy(msp->sname, series);
  
  strcpy(msp->qframe, "(+1)");
  
  parseLook(msp, look);
  
  getDesc(msp, line, look);
  
  insertFS(msp, series);
}


/* Parse data from the given line into an MSP of type GFF */
static void parseGff(char *line, char *opts, MSP *msp)
{
  char scorestring[256];
  char series[MAXLINE+1];
  char qname[MAXLINE+1];
  char look[MAXLINE+1];
  
  msp->type = FSSEG;
  
  /* <sequencename> <seriesname> <look> <start> <end> <score> <strand> <transframe> [annotation] */
  if (sscanf(line, "%s%s%s%d%d%s%s%s",
	     qname, series, look, &msp->qstart, &msp->qend, scorestring, 
	     msp->qframe+1, msp->qframe+2) != 8)
    {
      g_critical("Error parsing data, type GFF: \"%s\"\n", line);
      abort();
    }
  
  if (!strcmp(scorestring, ".")) msp->score = 100;
  else msp->score = 50.0*atof(scorestring);
  
  msp->qframe[0] = '(';
  msp->qframe[3] = ')';
  
  msp->sstart = msp->qstart;
  msp->send = msp->qend;
  
  msp->qname = g_malloc(strlen(qname)+1);
  strcpy(msp->qname, qname);
  
  msp->sname = g_malloc(strlen(series)+1); 
  strcpy(msp->sname, series);
  
  msp->desc = g_malloc(strlen(series)+1); 
  strcpy(msp->desc, series);
  
  msp->color = 6;	/* Blue */
  
  insertFS(msp, series);
}


/* Parse data from the given line into an MSP of type XY */
static void parseXY(char *line, char *opts, MSP *msp, char **seq1, char *seq1name, char **seq2, char *seq2name, BlxMSPType *mspType)
{
  int i, seqlen;
  
  msp->type = XY;
  
  char series[MAXLINE+1];
  char qname[MAXLINE+1];
  char look[MAXLINE+1];
  
  /* # FS type=XY <sequencename> <seriesname> <look> [annotation] */
  if (sscanf(line+13, "%s%s%s", qname, series, look) != 3)
    {
      g_critical("Error parsing data, type XY: \"%s\"\n", line);
      abort();
    }
  
  if (!seq1name || !seq2name)
    g_error("Sequencenames not provided");
  
  if (!strcasecmp(qname, seq1name) || !strcmp(qname, "@1"))
    {
      if (!seq1 || !*seq1)
	g_error("Sequence for %s not provided", qname);
      seqlen = strlen(*seq1);
    }
  else if (!strcasecmp(qname, seq2name) || !strcmp(qname, "@2"))
    {
      if (!seq2 || !*seq2)
	g_error("Sequence for %s not provided", qname);
      seqlen = strlen(*seq2);
    }
  else
    g_error("Invalid sequence name: %s", qname);
  
  if (!seqlen)
    g_error("Sequence for %s is empty", qname);
  
  msp->xy = arrayCreate(seqlen, int);
  for (i = 0; i < seqlen; i++)
    array(msp->xy, i, int) = XY_NOT_FILLED;
  
  msp->shape = XY_INTERPOLATE; /* default */
  
  *mspType = XYdata; /* Start parsing XY data next */
  
  msp->qname = g_malloc(strlen(qname)+1);
  strcpy(msp->qname, qname);
  
  msp->sname = g_malloc(strlen(series)+1); 
  strcpy(msp->sname, series);
  
  strcpy(msp->qframe, "(+1)");
  
  parseLook(msp, look);
  
  getDesc(msp, line, look);
  
  insertFS(msp, series); 
}


/* Utility to parse a file line for data that goes in an MSP */
static gboolean parseMspData(char *line, char *opts, MSP *msp, GString *line_string,
			     char **seq1, char *seq1name, char **seq2, char *seq2name, 
			     BlxMSPType *mspType)
{
  gboolean processed = FALSE;
  
  switch (*mspType)
  {
    case SEQBL: /* fall through */
    case EXBLX:
      parseEXBLXSEQBL(msp, *mspType, opts, line_string) ;
      processed = TRUE ;
      break;
      
    case SEQBL_X: /* fall through */
    case EXBLX_X:
      parseEXBLXSEQBLExtended(msp, *mspType, opts, line_string) ;
      processed = TRUE ;
      break;
      
    case HSP:
      parseHsp(line, opts, msp);
      break;
      
    case GSP:
      parseGsp(line, opts, msp, mspType);
      break;
      
    case FSSEG:
      parseFsseg(line, opts, msp);
      break;
      
    case GFF:
      parseGff(line, opts, msp);
      break;
      
    case XY:
      parseXY(line, opts, msp, seq1, seq1name, seq2, seq2name, mspType);
      break;
      
    default:
      break;
  };
  
  return processed;
}


/* Utility to parse a file line for data that doesn't go in an MSP */
static gboolean parseNonMspData(char *line, MSP *msp, 
				char **seq1, char *seq1name, char **seq2, char *seq2name, 
				BlxMSPType *mspType,
				char **readSeq, int *readSeqLen, int *readSeqMaxLen)
{
  gboolean processed = FALSE ;
  
  /* Data that don't go into a new MSP */
  if (*mspType == XYdata) 
    {
      int x, y;
      if (sscanf(line, "%d%d", &x, &y) != 2) 
	{
	  g_critical("Error parsing data file, type XYdata: \"%s\"\n", line);
	  abort();
	}
      array(msp->xy, x-1, int) = y;
      
      processed = TRUE ;
    }
  else if (*mspType == SEQdata) 
    {
      
      /* Realloc if necessary */
      if (*readSeqLen + strlen(line) > *readSeqMaxLen) 
	{
	  char *tmp;
	  *readSeqMaxLen += MAXLINE + strlen(line);
	  tmp = g_malloc(*readSeqMaxLen + 1);
	  strcpy(tmp, *readSeq);
	  g_free(*readSeq);
	  *readSeq = tmp;
	}
      strcpy(*readSeq + *readSeqLen, line);
      
      *readSeqLen += strlen(line);
      
      processed = TRUE;
    }
  else if (*mspType == SEQ) 
    {
      char series[MAXLINE+1];
      char qname[MAXLINE+1];
      
      if (sscanf(line+14, "%s%s", qname, series) != 2) 
	{
	  g_critical("Error parsing data file, type SEQ: \"%s\"\n", line);
	  abort();
	}
      
      if (!strcmp(qname, "@1")) 
	{
	  readSeq = seq1;
	  strcpy(seq1name, series);
	}
      else if (!strcmp(qname, "@2")) 
	{
	  readSeq = seq2;
	  strcpy(seq2name, series);
	}
      
      *readSeqMaxLen = MAXLINE;
      *readSeq = g_malloc(*readSeqMaxLen + 1);
      *readSeqLen = 0;
      
      *mspType = SEQdata;
      
      processed = TRUE;
    }
  
  return processed;
}

