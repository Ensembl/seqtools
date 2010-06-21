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
 * CVS info:   $Id: blxparser.c,v 1.28 2010-06-21 10:52:00 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/graph.h>
#include <wh/gex.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blxGff3Parser_.h>


/* Error codes and domain */
#define BLX_PARSER_ERROR g_quark_from_string("Parser")

typedef enum {
  BLX_PARSER_ERROR_NULL_MSP,	      /* MSP is null */
  BLX_PARSER_ERROR_INVALID_TYPE,      /* MSP does not have a valid type */
  BLX_PARSER_ERROR_NO_NAME,           /* MSP does not have a valid name */
  BLX_PARSER_ERROR_INVALID_COORDS,    /* MSP does not have valid ref seq coords */
  BLX_PARSER_ERROR_NO_STRAND,         /* MSP does not have a valid ref seq strand */
  BLX_PARSER_ERROR_INVALID_S_COORDS,  /* MSP does not have valid match coords */
  BLX_PARSER_ERROR_SEQS_DIFFER,       /* parsed data for the same sequence is different on different data lines */
  BLX_PARSER_ERROR_COMPLEMENT_FAILED  /* failed to complement the sequence */
} BlxDotterError;




static char *	    nextLine(FILE *file, GString *line_string);


static gboolean	    parseHeaderLine(char *line, char *opts, MSP *msp, BlxParserState *parserState);

static void	    parseBody(char *line, const int lineNum, char *opts, MSP **msp, GString *line_string, char **seq1, 
                              char *seq1name, char **seq2, char *seq2name, BlxParserState *parserState, MSP **mspList, GList **seqList,
                              char ***readSeq, int *readSeqLen, int *readSeqMaxLen);

static void	    parseEXBLXSEQBL(MSP **lastMsp, MSP **mspList, const BlxParserState parserState, char *opts, GString *line_string, GList **seqList);
static void	    parseEXBLXSEQBLExtended(MSP **lastMsp, MSP **mspList, const BlxParserState parserState, char *opts, GString *line_string, GList **seqList);
static void	    parseFsHsp(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList);
static void	    parseFsGspHeader(char *line, char *opts, MSP **lastMsp, MSP **mspList, BlxParserState *parserState, GList **seqList);
static void	    parseFsGspData(char *line, MSP *msp);
static void	    parseFsGff(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList);
static void	    parseFsSeg(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList);
static void	    parseFsXyHeader(char *line, char *opts, MSP **lastMsp, MSP **mspList, char **seq1, char *seq1name, char **seq2, char *seq2name, BlxParserState *parserState, GList **seqList);
static void	    parseFsXyData(char *line, MSP *msp);
static void         parseFsSeqHeader(char *line, char **seq1, char *seq1name, char **seq2, char *seq2name, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, BlxParserState *parserState);
static void         parseSeqData(char *line, char ***readSeq, int *readSeqLen, int *readSeqMaxLen);

static gboolean	    parseGaps(char **text, MSP *msp, const gboolean hasGapsTag) ;
static BOOL	    parseDescription(char **text, MSP *msp_unused) ;
static char*	    parseSequence(char **text, MSP *msp) ;
static void	    parseLook(MSP *msp, char *s) ;

static BlxMspType   getMspTypeFromScore(const int score);
static void	    getDesc(MSP *msp, char *s1, char *s2) ;
static char*	    prepSeq(const int sStart, char *inputSeq, char *opts) ;

static int	    fsSortByNameCompareFunc(void *fs1_in, void *fs2_in) ;
static int	    fsSortByOrderCompareFunc(void *fs1_in, void *fs2_in) ;

static void         getStrandAndFrameFromString(char *text, BlxStrand *strand, int *frame);
static void         checkReversedSubjectAllowed(const MSP *msp, const char *opts);

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


/* Function: parse a stream of SFS data
 *
 * Assumptions:
 *
 *  *seq1 and *seq2 may or may not be allocated.
 *
 *  *seq1name and *seq2name are allocated at least 255 bytes
 *
 */
void parseFS(MSP **MSPlist, FILE *file, char *opts, GList **seqList,
	     char **seq1, char *seq1name, char **seq2, char *seq2name, const int qOffset)
{
  DEBUG_ENTER("parseFS");

  if (!fsArr) 
    fsArr = arrayCreate(50, FeatureSeries);
  else
    arraySort(fsArr, fsSortByNameCompareFunc);

  /* Find the last MSP in the list - new MSPs will be tagged on to the end of this. We also
   * want to keep a record of the last MSP that was added in case additional information needs
   * be be added to the same MSP (e.g. for XY plot information, where we have multiple lines that
   * contain data for a single MSP). */
  MSP *msp = NULL;
  if (*MSPlist)
    {
      msp = *MSPlist;  
      while(msp->next)
	msp = msp->next;
    }

  /* Allocate reusable/extendable string as our buffer..*/
  GString *line_string = g_string_sized_new(MAXLINE + 1); 
  int lineNum = 0;

  BlxParserState parserState = PARSER_START;
  
  /* Variables for reading in SEQ data */
  char **readSeq = NULL;            /* buffer to hold the sequence we're reading in */
  int readSeqMaxLen = UNSET_INT;    /* current max length of the buffer */
  int readSeqLen = UNSET_INT;       /* current end pos of the data in the buffer */
  
  while (!feof(file) && parserState != PARSER_ERROR)
    { 
      ++lineNum;
      
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

      /* Check for header info first */
      if (parseHeaderLine(line, opts, msp, &parserState))
	{
	  continue; 
	}

      parseBody(line, lineNum, opts, &msp, line_string, 
                seq1, seq1name, seq2, seq2name, &parserState, MSPlist, seqList,
                &readSeq, &readSeqLen, &readSeqMaxLen);
    }

  g_string_free(line_string, TRUE) ;			    /* free everything, buffer and all. */

  /* Sort feature segment array by number */
  arraySort(fsArr, fsSortByOrderCompareFunc);

  DEBUG_EXIT("parseFS");
  return ;
}


/* Returns true if this msp is part of a feature series */
gboolean mspHasFs(const MSP *msp)
{
  gboolean result = (msp->type == BLXMSP_FS_SEG || msp->type == BLXMSP_XY_PLOT);
  return result;
}


/* Insert the given MSP into the Feature Series of the given name. If the Feature Series
 * does not exist yet, create it. Also, if the Feature Series array does not exist yet, create
 * that too. */
void insertFS(MSP *msp, char *series)
{
  if (!fsArr) 
    {
      fsArr = arrayCreate(50, FeatureSeries);
    }

  static int orderNum = 0; /* will increment this each time we add a feature series to the array */
  
  FeatureSeries fs;
  fs.on = 1;
  fs.y = 0.0;
  fs.xy = (msp->type == BLXMSP_XY_PLOT ? 1 : 0);

  fs.name = g_strdup(series);

  int i;
  if (arrayFind(fsArr, &fs, &i, fsSortByNameCompareFunc))
    {
      msp->fs = arrp(fsArr, i, FeatureSeries);
      g_free(fs.name);
    }
  else
    {
      /* Remember the order we added them so we can sort by it again later. */
      orderNum++;
      fs.order = orderNum;

      arrayInsert(fsArr, &fs, fsSortByNameCompareFunc);
      
      /* To do: there's a bug in here where we get the wrong pointer for msp->fs if the names in 
       * the file are not in alphabetical order... Intend to change this to use GArray anyway... */
      arrayFind(fsArr, &fs, &i, fsSortByNameCompareFunc);
      msp->fs = arrp(fsArr, i, FeatureSeries);
    }
}


/* Read in a FASTA sequence from a FASTA file */
char *readFastaSeq(FILE *seqfile, char *seqName)
{
  char *resultSeq;

  char line[MAXLINE+1];

  /* Read in seqfile */
  if (seqfile == stdin) 
    {
      if (!fgets(line, MAXLINE, seqfile))
        {
          g_error("Error reading seqFile.\n") ;
        }

      sscanf(line, "%s", seqName);

      /* Loop through the input text and copy into an auto-expandable array */
      Array arr = arrayCreate(5000, char);
      int i=0;
      char currentChar;

      while ((currentChar = fgetc(seqfile)) != '\n') 
        {
          if (isalpha(currentChar) || currentChar == SEQUENCE_CHAR_GAP || currentChar == SEQUENCE_CHAR_STOP) 
            {
              array(arr, i++, char) = currentChar;
            }
        }
      
      /* Allocate enough space in the result string, and copy the array contents into it */
      resultSeq = g_malloc(arrayMax(arr)+1);
      char *resultPos = resultSeq;
      
      for (i = 0; i < arrayMax(arr);) 
        {
          *resultPos++ = arr(arr, i++, char);
        }
      
      arrayDestroy(arr);
    }
  else
    {
      fseek(seqfile, 0, SEEK_END);
      resultSeq = g_malloc(ftell(seqfile)+1);
      char *resultPos = resultSeq;
      fseek(seqfile, 0, SEEK_SET);
      
      char *namePos;   /* current position in the seqName */
      char *linePos;   /* current position in the line */

      while (!feof(seqfile))
        {
          if (!fgets(line, MAXLINE, seqfile)) 
            {
              break;
            }
          
          while (strchr(line, '>'))
            {
              /* This line contains the sequence name. Copy it into the result arg */
              strncpy(seqName, line+1, 255); 
              seqName[255]=0;

              /* If the name contains a space or newline character, terminate it there */
              if ((namePos = (char *)strchr(seqName, ' ')))
                {
                  *namePos = 0;
                }
              
              if ((namePos = (char *)strchr(seqName, '\n')))
                {
                  *namePos = 0;
                }
              
              /* Get the next line */
              if (!fgets(line, MAXLINE, seqfile))
                {
                  break;
                }
            }

          /* Ok, this must be a line containing sequence data. Loop through each char in the line. */
          for (linePos = line; *linePos; linePos++)
            {
              /* If this is a valid sequence character, copy it into the result string. Allow
               * STOPs in the sequence. */
              if (isalpha(*linePos) || *linePos == SEQUENCE_CHAR_STOP)
                {
                  *resultPos++ = *linePos;
                }
            }
        
          *resultPos = 0;
        }
    }

  return resultSeq;
}



/*************************************************
 *               Internal functions.
 *************************************************/

/* Parse a string that contains a color name from our internally-defined list of accepted
 * colors. */
static int parseColor(char *s) 
{
  int i;

  for (i = 0; i < NUM_TRUECOLORS; i++) 
    {
      if (!strcasecmp(colorNames[i], s)) 
        break;
    }
    
  return i;
}

/* Parse a line that contains shape information about a curve. */
static BlxCurveShape parseShape(char *s) 
{
    if (!strcasecmp(s, "interpolate")) 
      return BLXCURVE_INTERPOLATE;
    else if (!strcasecmp(s, "partial")) 
      return BLXCURVE_PARTIAL;
    else 
      return BLXCURVE_BADSHAPE;
}

/* Parse a string that contains information about the way the given MSP
 * should look, i.e. color and shape-interpolation. */
static void parseLook(MSP *msp, char *s)
{
    char *cp, *s2;

    /* Make copy of string to mess up */
    s2 = g_strdup(s);

    cp = strtok(s2, "," );
    while (cp) {
	
	if (parseColor(cp) != NUM_TRUECOLORS) {
	    msp->fsColor = parseColor(cp);
	}
	else if (parseShape(cp) != BLXCURVE_BADSHAPE) {
	    msp->fsShape = parseShape(cp);
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

    msp->desc = g_strdup(cp);
}



/* Called when parsing HSP data. For certain modes, prepends dashes onto the start of the match
 * sequence if it does not start at 1. However, this old comment suggests that this isn't necessary:
 * "I think this is a mistake in the code, we get the whole sequence so shouldn't be prepending the
 * '-'s". */
static char* prepSeq(const int sStart, char *inputSeq, char *opts)
{
  char *result = NULL;
  
  if (*opts == 'T' || *opts == 'L') 
    {
      result = g_strdup(inputSeq) ;
    }
  else 
    {
      result = g_malloc(sStart + strlen(inputSeq)+1);
      memset(result, SEQUENCE_CHAR_PAD,sStart); /* Fill up with dashes */
      strcpy(result + sStart - 1, inputSeq);
    }

  return result;
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
static void parseEXBLXSEQBL(MSP **lastMsp, MSP **mspList, BlxParserState parserState, char *opts, GString *line_string, GList **seqList)
{
  int   qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL ;
  
  char sName[MAXLINE+1];
  char qframe[8] = "(+1)";
  int score = UNSET_INT;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sStart = UNSET_INT;
  int sEnd = UNSET_INT;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%d%d%s", 
	     &score, qframe, &qStart, &qEnd, 
	     &sStart, &sEnd, sName) != 7)
    {
      g_critical("Incomplete MSP data in input file.\n");
      abort() ;
    }

  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (*opts == 'T')
    {
      strcpy(qframe, "(+1)");
    }
  
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);
  
  const BlxMspType mspType = getMspTypeFromScore(score);

  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(lastMsp, mspList, seqList, mspType, score, 
                          NULL, qStart, qEnd, qStrand, qFrame,
                          sName, sStart, sEnd, BLXSTRAND_FORWARD, NULL,
                          opts, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  /* Convert subject names to fetchable ones if from NCBI server 
         
  Rule 1: If there is a gi, use that.
  Rule 2: If no gi, use the first and last non-blank field as db:id.
  */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_strdup(msp->sname);
  	
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
  
  if (!(cp = strstr(line, sName)))
    {
      g_error("Line does not include %s\n", sName);
    }
      
  seq_pos = cp + strlen(sName) ;
            
  qlen = mspGetQRangeLen(msp);
  slen = mspGetSRangeLen(msp);


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
  
  /* Parse the sequence, if applicable */
  char *sequence = NULL;
  
  if (parserState == EXBLX_BODY)
    {
      opts[BLXOPT_FULL_ZOOM] = ' '; /* Don't use full zoom default */
      
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
  else if (parserState == SEQBL_BODY)
    {
      if (*opts == 'L')
	{
	  slen = (abs(sEnd - sStart) + 1)/3;
	}
      
      
      /* Line contains chars other than sequence so get the starter data...not sure this test
       * is necessary any more now we have a better mechanism of getting a whole line 
       * from a file. All this is a horrible mixture of strtok and sscanf but what else
       * to do.... */
      if (strcspn(line, "acgt") != 0)
	{
	  sequence = g_malloc(line_string->len + 1) ;
	  
	  if (sscanf(seq_pos, "%s", sequence) != 1)
	    {
	      g_error("Error parsing %s\n", line);
	    }
          
          while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
            seq_pos++;
          while (*seq_pos && *seq_pos != ' ' && *seq_pos != '\t')
            seq_pos++;
	  while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	    seq_pos++;
	}
    }
  
  if (sequence)
    {
      addBlxSequenceData(msp->sSequence, sequence, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  /* Parse gaps data */
  if (seq_pos)
    {
      if (!parseGaps(&seq_pos, msp, FALSE))
        {
          g_error("Incomplete MSP gap data for MSP '%s' [%d - %d]\n", msp->sname, msp->sRange.min, msp->sRange.max) ;
        }
    }
  
  if (*opts == 'N')
    opts[3] = 'R';

  if (parserState == SEQBL_BODY)
    { 
      checkReversedSubjectAllowed(msp, opts);
    }
  
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
static void parseEXBLXSEQBLExtended(MSP **lastMsp, MSP **mspList, BlxParserState parserState, char *opts, GString *line_string, GList **seqList)
{
  BOOL result = FALSE ;
  int  qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL;
  char *first_pos ;
  
  char sName[MAXLINE+1];
  char qframe[8] ;
  char sframe[8] ;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sStart = UNSET_INT;
  int sEnd = UNSET_INT;
  int score = UNSET_INT;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%s%d%d%s", 
	     &score,
	     qframe, &qStart, &qEnd, 
	     sframe, &sStart, &sEnd, sName) != 8)
    {
      g_critical("Incomplete MSP data in input file.\n");
      abort() ;
    }

  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (*opts == 'T')
    strcpy(qframe, "(+1)");
  
  const BlxMspType mspType = getMspTypeFromScore(score) ;
  
  /* Extract frame and strand from qframe/sframe text */
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);

  BlxStrand sStrand = BLXSTRAND_NONE;
  getStrandAndFrameFromString(sframe, &sStrand, NULL);

  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(lastMsp, mspList, seqList, mspType, score,
                          NULL, qStart, qEnd, qStrand, qFrame, 
                          sName, sStart, sEnd, sStrand, NULL,
                          opts, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  
  /* Convert subject names to fetchable ones if from NCBI server 
   *  Rule 1: If there is a gi, use that.
   *  Rule 2: If no gi, use the first and last non-blank field as db:id.
   */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_strdup(msp->sname);
  	
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
  
  if (!(cp = strstr(line, sName)))
    {
      g_error("Line does not include %s\n", sName);
    }
      
  seq_pos = cp + strlen(sName) ;
            
  qlen = mspGetQRangeLen(msp);
  slen = mspGetSRangeLen(msp);

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
  char *sequence = NULL;
  
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
		  if (!(result = parseGaps(&seq_pos, msp, TRUE)))
		    g_error("Incomplete MSP gap data\n") ;
		}
	      else if (parserState == EXBLX_BODY && (strstr(seq_pos, BLX_DESCRIPTION_TAG)))
		{
		  opts[BLXOPT_FULL_ZOOM] = ' '; /* Don't use full zoom default */
                  result = parseDescription(&seq_pos, msp);
                  
		  if (!result)
                    {
                      g_error("Bad description data\n");
                    }
		}
	      else if ((parserState == SEQBL_X_BODY || mspIsSnp(msp)) && (strstr(seq_pos, BLX_SEQUENCE_TAG)))
		{
		  if (*opts == 'L')
		    {
		      slen = mspGetSRangeLen(msp) / 3;
		    }
                  
                  sequence = parseSequence(&seq_pos, msp);
                  result = (sequence != NULL);
                  
		  if (!result)
                    {
                      g_error("Bad sequence data\n");
                    }
		}
              
	    }
	}
    }
  
  if (sequence)
    {
      addBlxSequenceData(msp->sSequence, sequence, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  if (parserState == SEQBL_X_BODY)
    {
      if (*opts == 'L')
	{
	  slen = mspGetSRangeLen(msp) / 3;
	}
	

    }

  if (*opts == 'N')
    opts[3] = 'R';
  
  if (parserState == SEQBL_X_BODY)
    { 
      checkReversedSubjectAllowed(msp, opts);
    }
  
  return ;
}


/* Comparison function to sort two Feature Series by the order number stored in the FeatureSeries
 * struct. Returns -1 if the first item is before the second, 1 if the second is first, or 0 if 
 * they are equal.  */
static int fsSortByOrderCompareFunc(void *fs1_in, void *fs2_in)
{
  int result = 0;

  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;

  if (fs1->order < fs2->order)
    result = -1;
  else if (fs1->order > fs2->order)
    result = 1;
  
  return result;
}


/* Comparison function to sort two Features Series by name. */
static int fsSortByNameCompareFunc(void *fs1_in, void *fs2_in)
{
  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;

  /*printf("%s - %s : %d\n", fs1->name, fs2->name,  strcmp(fs1->name, fs2->name));*/

  return strcmp(fs1->name, fs2->name);
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
	    g_error("NULL value returned on reading input file\n") ;
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

	      if (buffer[line_end] == '\n' || buffer[line_end] == '\r')
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
 * to first char after ';'. 
 * The hasGapsTag is a bit of a hack to also allow this function to work with the old style
 * exblx file format that does not have a "Gaps" tag on the front of the gaps data. (In this case
 * we need to do something slightly different with strtok to start reading in the correct place
 * in the string.) */
static gboolean parseGaps(char **text, MSP *msp, const gboolean hasGapsTag)
{
  gboolean result = TRUE;

  CoordRange *currentGap = NULL;
  gboolean finished = FALSE;

  char *currentGapStr = hasGapsTag ? strtok(NULL, " ") : strtok(*text, " ") ;

  while (result && !finished && currentGapStr)
    {
      int i  = 0;

      /* Loop through the four coords in this range */
      for (i = 0 ;  i < 4 ; i++)
	{
	  if (!currentGapStr)
	    {
	      result = FALSE;
	      break;
	    }
	  else if (*currentGapStr == ';')
	    {
	      finished = TRUE;
	      break;
	    }
	     
	  switch (i)
	  {
	    case 0:
	    {
	      /* First value is start of subject sequence range. Create the range struct */
              currentGap = g_malloc(sizeof(CoordRange));
              msp->gaps = g_slist_append(msp->gaps, currentGap);
	      currentGap->sStart = convertStringToInt(currentGapStr);
	      break;
	    }
	      
	    case 1:
	    {
	      /* Second value is end of subject sequence range. Order values so that
	       * sStart is less than sEnd if we have the forward strand or v.v. if the reverse. */
	      currentGap->sEnd = convertStringToInt(currentGapStr);
	      sortValues(&currentGap->sStart, &currentGap->sEnd, mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
	      break;
	    }
	      
	    case 2:
	    {
	      /* Third value is start of reference sequence range */
	      currentGap->qStart = convertStringToInt(currentGapStr);
	      break;
	    }
	      
	    case 3:
	    {
	      /* Fourth value is end of reference sequence range. Order values so that
	       * qStart is less than qEnd if ref sequence is forward strand or v.v. if the reverse. */
	      currentGap->qEnd = convertStringToInt(currentGapStr);
	      sortValues(&currentGap->qStart, &currentGap->qEnd, mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
	    }
	  }

	  currentGapStr = strtok(NULL, "\t ") ; 
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


/* Description, should just be plain text, format is:
 * 
 * "Sequence text ; "
 * 
 * text following Sequence must not contain ';'. Moves text to first
 * char after ';'. */
static char* parseSequence(char **text, MSP *msp)
{
  char *result = NULL;
  char *cp = *text ;

  cp = strtok(NULL, "\t ") ;				    /* skip "Sequence" */

  size_t origLen = strlen(cp);
  size_t validLen = strspn(cp, "acgtnACGTN"); 
  
  if (validLen < 1 || validLen < strlen(cp))
    {
      g_error("Error parsing %s, coords are %d -> %d (%d), but valid length sequence is only %d (out of total length supplied = %d).\n",
		msp->sname, msp->sRange.min, msp->sRange.max, msp->sRange.max - msp->sRange.min, (int)validLen, (int)origLen) ;
    }
  else
    {
      result = g_malloc(validLen + 1) ;
	  
      if (sscanf(cp, "%s", result) != 1)
	{
	  g_error("Error parsing %s\n", cp) ;
	}
    }

  return result ;
}


/* Utility called by parseFS to parse the header info of a line from a file. Returns true if the
 * line was completely processed, false if further processing on the same line is still required */
static gboolean parseHeaderLine(char *line, char *opts, MSP *msp, BlxParserState *parserState)
{
  DEBUG_ENTER("parseHeaderLine(parserState=%d)", *parserState);
  
  gboolean processed = FALSE;
  
  if (!strncasecmp(line, "##gff-version   3", 17))
    {
      *parserState = GFF_3_HEADER ;
      processed = TRUE;
    }
  else if (!strncasecmp(line, "##FASTA", 7))
    {
      *parserState = FASTA_SEQ_HEADER ;
      processed = TRUE;
    }
  else if (!strncasecmp(line, "# seqbl_x", 9))
    {
      *parserState = SEQBL_X_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx_x", 9))
    {
      *parserState = EXBLX_X_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# seqbl", 7))
    {
      /* Only for backwards compatibility */
      *parserState = SEQBL_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx", 7))
    {
      /* Only for backwards compatibility */
      *parserState = EXBLX_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=HSP", 13) || 
	   !strncasecmp(line, "# SFS type=HSP", 14))
    {
      *parserState = FS_HSP_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GSP", 13) || 
	   !strncasecmp(line, "# SFS type=GSP", 14))
    {
      *parserState = FS_GSP_HEADER;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# dotter feature format 2", 25) ||
	   !strncasecmp(line, "# FS type=SEG", 13) ||
	   !strncasecmp(line, "# SFS type=SEG", 14))
    {
      *parserState = FS_SEG_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GFF", 13) || 
	   !strncasecmp(line, "# SFS type=GFF", 14))
    {
      *parserState = FS_GFF_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=XY", 12) ||
	   !strncasecmp(line, "# SFS type=XY", 13))
    {
      /* More info exists in the header line for this data type, so indicate we need to parse
       * header info and return false to indicate we have not finished with this line */
      *parserState = FS_XY_HEADER;
      processed = FALSE;
    }
  else if (!strncasecmp(line, "# FS type=SEQ", 13) ||
	   !strncasecmp(line, "# SFS type=SEQ", 14))
    {
      *parserState = FS_SEQ_HEADER;
      processed = FALSE; /* more info exists in the header line for this data type */
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
	  opts[BLXOPT_HSP_GAPS] = 'G';
	}
      else if (!strncasecmp(line, "# DESC ", 7) &&
	       (*parserState == FS_HSP_BODY || *parserState == FS_GSP_HEADER || *parserState == SEQBL_BODY))
	{
	  if (msp)
	    getDesc(msp, line, msp->sname);
	}
      
      /* If we're in a GFF header we want to take a look at additional comment lines;
       * otherwise we ignore additional comments. */
      processed = *parserState == GFF_3_HEADER ? FALSE : TRUE ;
    }
  else if (*parserState == GFF_3_HEADER)
    {
      /* We were processing GFF3 header lines (which all start with '#'), but this line does not
       * start with '#', so it must be the start of the GFF3 body. */
      *parserState = GFF_3_BODY;
      processed = FALSE;
    }


  DEBUG_EXIT("parseHeaderLine returning processed = %d, (parserState = %d)", processed, *parserState);
  return processed ;
}


/* Parse data from the given line into an MSP of type FS_HSP */
static void parseFsHsp(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList)
{
  char qName[MAXLINE+1];
  char sName[MAXLINE+1];
  char seq[MAXLINE+1];
  char qframe[8] = "(+1)";
  char sframe[8] = "(+1)";
  int qStart = UNSET_INT;
  int sStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sEnd = UNSET_INT;
  int score = UNSET_INT;
  
  /* <score> <qname> <qframe> <qstart> <qend> <sname> <sframe> <sstart> <ssend> <sequence> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s%s%d%d%s", 
	     &score, 
	     qName, qframe+1, &qStart, &qEnd, 
	     sName, sframe+1, &sStart, &sEnd,
	     seq) != 10)
    {
      g_critical("Error parsing data, type FS_HSP: \"%s\"\n", line);
      abort();
    }
  
  /* Get the strand and frame from the frame strings */
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);

  BlxStrand sStrand = BLXSTRAND_NONE;
  int sFrame = UNSET_INT;
  getStrandAndFrameFromString(sframe, &sStrand, &sFrame);

  /* Get the sequence data */
  char *sSeq = prepSeq(sStart, seq, opts);

  /* Create the new MSP */
  GError *error = NULL;

  MSP *msp = createNewMsp(lastMsp, mspList, seqList, BLXMSP_HSP, score,
                           qName, qStart, qEnd, qStrand, qFrame, 
                           sName, sStart, sEnd, sStrand, sSeq,
                           opts, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  checkReversedSubjectAllowed(msp, opts);
}


/* Parse data from the given line, which contains feature-series GSP header info */
static void parseFsGspHeader(char *line, char *opts, MSP **lastMsp, MSP **mspList, BlxParserState *parserState, GList **seqList)
{
  /* Data goes into a new MSP */
//  msp->type = BLXMSP_GSP;

  /* Will write this as soon as MSPcrunch generates it */
  g_warning("Parser for GSP data not implemented\n");
  
  /* Start parsing the GSP data next */
  *parserState = FS_GSP_BODY ;
}


/* Parse a line that contains feature-series GSP data */
static void parseFsGspData(char *line, MSP *msp)
{
  g_warning("Parser for GSP data not implemented\n");
}


/* Parse data from the given line into an MSP of type FS_SEG (Feature Series Segment) */
static void parseFsSeg(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList)
{
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  int score = UNSET_INT;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  
  /* <score> <sequencename> <seriesname> <start> <end> <look> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s",
	     &score, qName, series, &qStart, &qEnd, look) != 6)
    {
      g_critical("Error parsing data, type FS_SEG: \"%s\"\n", line);
      abort();
    }
  
  
  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(lastMsp, mspList, seqList, BLXMSP_FS_SEG, UNSET_INT,
                          qName, qStart, qEnd, BLXSTRAND_NONE, 1, 
                          series, qStart, qEnd, BLXSTRAND_NONE, NULL,
                          opts, &error);

  /* Parse in additional feature-series info */
  parseLook(msp, look);
  getDesc(msp, line, look);
  insertFS(msp, series);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
}


/* Parse data from a line of text that contains feature-series data in GFF format */
static void parseFsGff(char *line, char *opts, MSP **lastMsp, MSP **mspList, GList **seqList)
{
  char scorestring[256];
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  char qframe[8] = "(+1)";
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  
  /* <sequencename> <seriesname> <look> <start> <end> <score> <strand> <transframe> [annotation] */
  if (sscanf(line, "%s%s%s%d%d%s%s%s",
	     qName, series, look, &qStart, &qEnd, scorestring, 
	     qframe+1, qframe+2) != 8)
    {
      g_critical("Error parsing data, type FS_GFF: \"%s\"\n", line);
      abort();
    }
  
  int score = UNSET_INT;
  if (!strcmp(scorestring, ".")) 
    {
      score = 100;
    }
  else 
    {
      score = 50.0*atof(scorestring);
    }
  
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);
  
  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(lastMsp, mspList, seqList, BLXMSP_FS_SEG, score,
                          qName, qStart, qEnd, qStrand, qFrame, 
                          series, qStart, qEnd, BLXSTRAND_FORWARD, NULL,
                          opts, &error);
  

  /* Parse additional feature-series information */
  msp->desc = g_strdup(series);
  parseLook(msp, look);
  insertFS(msp, series);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
}


/* Parse data from the given line, which contains header info for Feature-Series XY-plot data */
static void parseFsXyHeader(char *line, char *opts, MSP **lastMsp, MSP **mspList, char **seq1, char *seq1name, char **seq2, char *seq2name, BlxParserState *parserState, GList **seqList)
{
  int i, seqlen;
  
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  
  /* # FS type=XY <sequencename> <seriesname> <look> [annotation] */
  if (sscanf(line+13, "%s%s%s", qName, series, look) != 3)
    {
      g_critical("Error parsing data, type XY: \"%s\"\n", line);
      abort();
    }
  
  if (!seq1name || !seq2name)
    g_error("Sequencenames not provided\n");
  
  if (!strcasecmp(qName, seq1name) || !strcmp(qName, "@1"))
    {
      if (!seq1 || !*seq1)
	g_error("Sequence for %s not provided\n", qName);
      seqlen = strlen(*seq1);
    }
  else if (!strcasecmp(qName, seq2name) || !strcmp(qName, "@2"))
    {
      if (!seq2 || !*seq2)
	g_error("Sequence for %s not provided\n", qName);
      seqlen = strlen(*seq2);
    }
  else
    g_error("Invalid sequence name: %s\n", qName);
  
  if (!seqlen)
    g_error("Sequence for %s is empty\n", qName);

  /* Create an MSP to put the data in */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(lastMsp, mspList, seqList, BLXMSP_XY_PLOT, 0,
                          qName, UNSET_INT, UNSET_INT, BLXSTRAND_FORWARD, 1,
                          series, UNSET_INT, UNSET_INT, BLXSTRAND_FORWARD, NULL, 
                          opts, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  /* Add additional MSP information */
  msp->xy = arrayCreate(seqlen, int);
  for (i = 0; i < seqlen; i++)
    { 
      array(msp->xy, i, int) = XY_NOT_FILLED;
    }
  
  msp->fsShape = BLXCURVE_INTERPOLATE; /* default */
  parseLook(msp, look);
  getDesc(msp, line, look);
  insertFS(msp, series); 
  
  /* Start parsing the actual XY data next */
  *parserState = FS_XY_BODY;
}


/* Parse a line that contains XY-plot data for a feature-series. Data is populated into
 * the given MSP. */
static void parseFsXyData(char *line, MSP *msp)
{
  /* XY coordinates go in to the XY list in the current MSP. */
  int x, y;
  if (sscanf(line, "%d%d", &x, &y) != 2) 
    {
      g_critical("Error parsing data file, type FS_XY_BODY: \"%s\"\n", line);
      abort();
    }
  array(msp->xy, x-1, int) = y;
}


/* Parse a line that contains feature-series sequence header info */
static void parseFsSeqHeader(char *line, 
                             char **seq1, char *seq1name, char **seq2, char *seq2name,
                             char ***readSeq, int *readSeqLen, int *readSeqMaxLen, BlxParserState *parserState)
{
  /* We have a header for sequence data. Read in the sequence name from the header line. */
  char series[MAXLINE+1];
  char qname[MAXLINE+1];
  
  if (sscanf(line+14, "%s%s", qname, series) != 2) 
    {
      g_critical("Error parsing data file, type FS_SEQ_HEADER: \"%s\"\n", line);
      abort();
    }
  
  /* We need to allocate a string to read the seqeunce data in to. We want to allocate the releavant
   * return pointer (seq1 or seq2) depending on which string was indicated in the header line. */
  if (!strcmp(qname, "@1")) 
    {
      *readSeq = seq1;
      strcpy(seq1name, series);
    }
  else if (!strcmp(qname, "@2")) 
    {
      *readSeq = seq2;
      strcpy(seq2name, series);
    }
  
  *readSeqMaxLen = MAXLINE;
  **readSeq = g_malloc(*readSeqMaxLen + 1);
  *readSeqLen = 0;
  
  /* Update the parser state so that we proceed to parse the sequence data next */
  *parserState = FS_SEQ_BODY;
}


/* Parse a line that contains a chunk of sequence data. Concatenates the contents of the
 * line onto readSeq, extending readSeq if necessary. readSeqLen is the length of the current
 * contents of readSeq, and readSeqMaxLen is the currently allocated space for readSeq. */
static void parseSeqData(char *line, char ***readSeq, int *readSeqLen, int *readSeqMaxLen)
{
  if (*readSeqLen == UNSET_INT)
    {
      /* Not initialised for reading in sequence; do nothing */
      return;
    }
  
  /* Read in the actual sequence data. It may span several lines, so concatenate each
   * one on to the end of our result. */
  
  /* First, reealloc if necessary */
  if (*readSeqLen + strlen(line) > *readSeqMaxLen) 
    {
      char *tmp;
      *readSeqMaxLen += MAXLINE + strlen(line);
      tmp = g_malloc(*readSeqMaxLen + 1);
      strcpy(tmp, **readSeq);
      g_free(**readSeq);
      **readSeq = tmp;
    }
  
  /* Copy in the contents of the line */
  strcpy(**readSeq + *readSeqLen, line);
  
  /* Update the length counter */
  *readSeqLen += strlen(line);
}


/* Utility to call the relevant parser function to parse the current data type */
static void parseBody(char *line, const int lineNum, char *opts, MSP **lastMsp, GString *line_string,
                      char **seq1, char *seq1name, char **seq2, char *seq2name, 
                      BlxParserState *parserState, MSP **mspList, GList **seqList,
                      char ***readSeq, int *readSeqLen, int *readSeqMaxLen)
{
  DEBUG_ENTER("parseBody(parserState=%d, line=%d)", *parserState, lineNum);
  
  /* Call the relevant function for the current type of data being parsed */
  switch (*parserState)
  {
    case GFF_3_HEADER:
      parseGff3Header(lineNum, lastMsp, mspList, parserState, opts, line_string, seqList, seq1name);
      break;
      
    case GFF_3_BODY:
      parseGff3Body(lineNum, lastMsp, mspList, parserState, opts, line_string, seqList);
      break;

    case SEQBL_BODY: /* fall through */
    case EXBLX_BODY:
      parseEXBLXSEQBL(lastMsp, mspList, *parserState, opts, line_string, seqList) ;
      break;
      
    case SEQBL_X_BODY: /* fall through */
    case EXBLX_X_BODY:
      parseEXBLXSEQBLExtended(lastMsp, mspList, *parserState, opts, line_string, seqList) ;
      break;

    case FS_HSP_BODY:
      parseFsHsp(line, opts, lastMsp, mspList, seqList);
      break;

    case FS_GSP_HEADER:
      parseFsGspHeader(line, opts, lastMsp, mspList, parserState, seqList);
      break;

    case FS_GSP_BODY:
      parseFsGspData(line, *lastMsp);
      break;
      
    case FS_GFF_BODY:
      parseFsGff(line, opts, lastMsp, mspList, seqList);
      break;

    case FS_SEG_BODY:
      parseFsSeg(line, opts, lastMsp, mspList, seqList);
      break;
      
    case FS_XY_HEADER:
      parseFsXyHeader(line, opts, lastMsp, mspList, seq1, seq1name, seq2, seq2name, parserState, seqList);
      break;
      
    case FS_XY_BODY:
      parseFsXyData(line, *lastMsp);
      break;

    case FS_SEQ_HEADER:
      parseFsSeqHeader(line, seq1, seq1name, seq2, seq2name, readSeq, readSeqLen, readSeqMaxLen, parserState);
      break;

    case FASTA_SEQ_HEADER:
      parseFastaSeqHeader(line, lineNum, seq1, seq1name, readSeq, readSeqLen, readSeqMaxLen, parserState);
      break;

    case FS_SEQ_BODY: /* fall through */
    case FASTA_SEQ_BODY:
      parseSeqData(line, readSeq, readSeqLen, readSeqMaxLen);
      break;
      
    case PARSER_ERROR:
      break;

    default:
      g_warning("Unexpected state '%d' while parsing input file.\n", (int)(*parserState));
      *parserState = PARSER_ERROR;
      break;
  };
  
  DEBUG_EXIT("parseBody");
}


/* For old file types, the hacky way of specifying the type of the MSP was by using a negative score
 * for anything that is not a match. This function determines the type from the given score. */
static BlxMspType getMspTypeFromScore(const int score)
{
  BlxMspType result = BLXMSP_INVALID;
  
  if (score >= 0)
    {
      result = BLXMSP_MATCH;
    }
  else if (score == -1)
    {
      result = BLXMSP_EXON_CDS;
    }
  else if (score == -2)
    {
      result = BLXMSP_INTRON;
    }
  else if (score == -3)
    {
      result = BLXMSP_SNP;
    }

  return result;  
}



/* Add the given sequence data to a BlxSequence. Validates that the existing sequence data is
 * either null or is the same as the new data; sets the given error if not. We claim ownership
 * of the given sequence data (either the BlxSequence owns it, or we delete it if it is not 
 * required). The given sequence should always be the forward strand; we complement it ourselves
 * here if this BlxSequence requires the reverse strand. */
void addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error)
{
  if (!sequence)
    {
      return;
    }
  
  gboolean sequenceUsed = FALSE;
  
  if (blxSeq && blxSeq->sequenceReqd)
    {
      if (!blxSeq->sequence)
        {
          /* Sequence does not yet exist, so add it */
          if (blxSeq->strand == BLXSTRAND_REVERSE)
            {
              blxComplement(sequence) ;
	      
	      if (!sequence)
		{
		  g_set_error(error, BLX_PARSER_ERROR, BLX_PARSER_ERROR_COMPLEMENT_FAILED, 
			      "Failed to complement sequence '%s'", blxSeq->fullName);
		}
            }
          
          blxSeq->sequence = sequence;
          sequenceUsed = TRUE;
        }
      else if (error && *error)
        {
          /* Sequence already exists. Validate that it's the same as the existing one.
           * (Remember to complement the passed in one if the existing one is complemented.) */
          if (blxSeq->strand == BLXSTRAND_REVERSE)
            {
              blxComplement(sequence);
            }
          
          if (strcmp(sequence, blxSeq->sequence) != 0)
            {
              g_set_error(error, BLX_PARSER_ERROR, BLX_PARSER_ERROR_SEQS_DIFFER,
                          "Invalid sequence data: multiple matches were found for sequence '%s' (%s strand), but the parsed sequence data is not the same for each instance of this sequence. Alignments for this sequence may therefore not display correctly. Please check the input data.\n",
                          blxSeq->fullName, (blxSeq->strand == BLXSTRAND_FORWARD ? "forward" : "reverse"));
            }
        }

      /* If we were successful, calculate the sequence range */
      if (blxSeq->sequence)
	{
	  blxSeq->seqRange.min = 1;
	  blxSeq->seqRange.max = strlen(blxSeq->sequence);
	}
    }      
  
  if (!sequenceUsed)
    {
      g_free(sequence);
    }
}


/* Extract the strand and frame from the given string of the format "(+1)" and populate
 * the return values, if non-null. */
static void getStrandAndFrameFromString(char *text, BlxStrand *strand, int *frame)
{
  if (frame)
    {
      *frame = convertStringToInt(&text[2]);
    }

  if (strand)
    {
      if (text[1] == '+')
        {
          *strand = BLXSTRAND_FORWARD;
        }
      else if (text[1] == '-')
        {
          *strand = BLXSTRAND_REVERSE;
        }
      else
        {
          *strand = BLXSTRAND_NONE;
        }
    }
}


/* Check if we have a reversed subject and, if so, if this is allowed. Throws an error if not. */
static void checkReversedSubjectAllowed(const MSP *msp, const char *opts)
{
  if (mspGetMatchStrand(msp) == BLXSTRAND_REVERSE && *opts != 'T' && *opts != 'L' && *opts != 'N')
    {
      g_error("Reversed subjects are not allowed in modes blastp or blastx.\n");
    }
}



