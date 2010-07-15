/*  File: blxGff3parser.c
 *  Author: Gemma Barson (gb10@sanger.ac.uk)
 *  Copyright (c) 2010: Genome Research Ltd.
 *-------------------------------------------------------------------
 * Blixem is free software; you can redistribute it and/or
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
 * Author: Gemma Barson (Sanger Institute, UK) gb10@sanger.ac.uk
 *
 * Description: parser for GFF version 3 format.
 *              
 * Exported functions: 
 *-------------------------------------------------------------------
 */

#include <SeqTools/blxGff3Parser_.h>
#include <SeqTools/utilities.h>
#include <SeqTools/blixem_.h>


/* Error codes and domain */
#define BLX_GFF3_ERROR g_quark_from_string("GFF 3 parser")

typedef enum {
  BLX_GFF3_ERROR_INVALID_STRAND,	      /* invalid strand in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_TYPE,                /* invalid type in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_NUM_TOKENS,          /* invalid number of columns from a line of the input file */
  BLX_GFF3_ERROR_INVALID_TAG,                 /* invalid format for a tag/data pair */
  BLX_GFF3_ERROR_INVALID_SEQ,                 /* invalid sequence data */
  BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT,        /* invalid CIGAR format */
  BLX_GFF3_ERROR_INVALID_MSP,                 /* MSP has invalid/missing data */
  BLX_GFF3_ERROR_UNKNOWN_MODE,                /* unknown blast mode */
  BLX_GFF3_ERROR_BAD_COLOR                    /* Bad color string found when parsing color */
} BlxDotterError;


static void           parseGffColumns(GString *line_string, const int lineNum, const char *opts, MSP *msp, GList **seqList, GSList *styles, GError **error);
static void           parseAttributes(char *attributes, MSP *msp, GList **seqList, const char *opts, const int lineNum, GError **error);
static void           parseTagDataPair(char *text, MSP *msp, char **sequence, BlxStrand *strand, char **gapString, GList **seqList, const int lineNum, GError **error);
static void           parseNameTag(char *data, MSP *msp, const int lineNum, GError **error);
static void           parseTargetTag(char *data, MSP *msp, BlxStrand *strand, GList **seqList, const int lineNum, GError **error);
static void           parseGapString(char *text, const char *opts, MSP *msp, GError **error);

static BlxStrand      readStrand(char *token, GError **error);
static void           parseMspType(char *token, MSP *msp, GError **error);
static void           parseCigarStringSection(const char *text, MSP *msp, const int qDirection, const int sDirection, const int numFrames, int *q, int *s, GError **error);
static int            validateNumTokens(char **tokens, const int minReqd, const int maxReqd, GError **error);
static void           validateMsp(const MSP *msp, GError **error);


/* Parse GFF3 header information */
void parseGff3Header(const int lineNum,
                     MSP **lastMsp, 
		     MSP **mspList, 
		     BlxParserState *parserState, 
		     char *opts, 
		     GString *line_string, 
		     GList **seqList,
                     char *refSeqName)
{
  DEBUG_ENTER("parseGff3Header [line=%d]", lineNum);
  
  /* Look for the "sequence-region" comment line, which tells us info about the reference
   * sequence. The format is as follows: ##sequence-region    qname qstart qend */
  char qName[MAXLINE + 1];
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;

  if (!strncasecmp(line_string->str, "##sequence-region", 17))
    {
      if (sscanf(line_string->str, "##sequence-region%s%d%d", qName, &qStart, &qEnd) < 1)
        {
          g_critical("Error parsing data file, type GFF_3_HEADER: \"%s\"\n", line_string->str);
          abort();
        }
      
      DEBUG_OUT("Found reference sequence name=%s [start=%d, end=%d]\n", qName, qStart, qEnd);

      /* If the ref seq name is already populated, check it's the same as the one we've just read */
      if (*refSeqName != '\0')
        {
          if (!stringsEqual(refSeqName, qName, FALSE))
            {
              g_critical("Reference sequence name differs between files: GFF file name is '%s' but FASTA file name was '%s'. Ignoring GFF file name.\n", qName, refSeqName);
            }
        }
      else
        {
          strcpy(refSeqName, qName);
        }
      
      /* to do: do something with qstart and qend? We will find the ref seq coords from the
       * offset and the ref seq length, so we should check that these are consistent with those... */
    }
  
  DEBUG_EXIT("parseGff3Header");
}


/* Parse GFF3 data */
void parseGff3Body(const int lineNum,
                   MSP **lastMsp, 
		   MSP **mspList, 
		   BlxParserState *parserState, 
		   char *opts, 
		   GString *line_string, 
		   GList **seqList,
                   GSList *styles)
{
  DEBUG_ENTER("parseGff3Body [line=%d]", lineNum);
  
  /* Data goes into a new MSP */
  MSP *prevMsp = *lastMsp;
  MSP *msp = createEmptyMsp(lastMsp, mspList);

  GError *error = NULL;
  parseGffColumns(line_string, lineNum, opts, msp, seqList, styles, &error);
  
  if (!error)
    {
      validateMsp(msp, &error);
    }
  
  if (error)
    {
      /* Destroy the incomplete msp and remove it from the list */
      if (prevMsp)
        {
          prevMsp->next = NULL;
          *lastMsp = prevMsp;
        }
      else
        {
          /* List was empty */
          *mspList = NULL;
        }
      
      destroyMspData(msp);
      g_free(msp);
      msp = NULL;

      prefixError(error, "[line %d] Error parsing GFF data. ", lineNum);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  DEBUG_EXIT("parseGff3Body");
}


/* Parse header info for a FASTA sequence in a GFF 3 file. The only info there should be is the
 * sequence name, on its own line below the "##FASTA" header. This should be in the format ">name".
 * Sets the readSeq pointer to point to the sequence that needs to be populated, according to the
 * sequence name that is parsed. */
/* to do: currently this only allows the reference sequence to be in fasta format. could be extended
 * to read match sequences in fasta too, although we should probably make the way dotter and blixem
 * use sequences more consistent before doing that to save refactoring a lot of code later (ie.
 * if all sequences (including ref seq) were stored in a BlxSequence then it would make things
 * much easier). */
void parseFastaSeqHeader(char *line, const int lineNum,
                         char **refSeq, char *refSeqName,
                         char ***readSeq, int *readSeqLen, int *readSeqMaxLen,
                         BlxParserState *parserState)
{
  char seqName[MAXLINE + 1];
  
  if (sscanf(line, ">%s", seqName) != 1)
    {
      g_critical("Error parsing data file, type FASTA_SEQ_HEADER: \"%s\"\n", line);
      abort();
    }

  /* The reference sequence name should already be set */
  if (refSeqName && stringsEqual(refSeqName, seqName, FALSE))
    {
      /* Only bother if the reference sequence is not already populated */
      if (*refSeq == NULL)
        {
          *readSeq = refSeq;
          *readSeqMaxLen = MAXLINE;
          **readSeq = g_malloc(*readSeqMaxLen + 1);
          *readSeqLen = 0;
        }
      
      /* Update the parser state so that we proceed to parse the sequence data next. (Even if
       * we're not populating the ref seq, we still need to loop over these lines. Leaving the
       * readSeqLen as unset will mean that the fasta sequence parser will ignore the input.) */
      *parserState = FASTA_SEQ_BODY;
    }
  else
    {
      g_critical("[line %d] Failed to parse FASTA sequence: unrecognised sequence name '%s'\n", lineNum, seqName);
      *parserState = PARSER_ERROR;
    }
}


/*********************************************************
 *                  Internal functions
 *********************************************************/

static gboolean isNucleotideMatch(char *text)
{
  return (stringsEqual(text, "nucleotide_match", FALSE) || stringsEqual(text, "SO:0000347", FALSE));
}

static gboolean isProteinMatch(char *text)
{
  return (stringsEqual(text, "protein_match", FALSE) || stringsEqual(text, "SO:0000349", FALSE));
}

static gboolean isMatchPart(char *text)
{
  return (stringsEqual(text, "match_part", FALSE) || stringsEqual(text, "SO:0000039", FALSE));
}

static gboolean isMatchSet(char *text)
{
  return (stringsEqual(text, "match_set", FALSE) || stringsEqual(text, "SO:0000038", FALSE));
}

static gboolean isMatch(char *text)
{
  return (isNucleotideMatch(text) || 
          isProteinMatch(text) || 
          isMatchPart(text) ||
          isMatchSet(text) ||
          stringsEqual(text, "match", FALSE) || stringsEqual(text, "SO:0000343", FALSE));
}

static gboolean isCds(char *text)
{
  return (stringsEqual(text, "CDS", FALSE) || stringsEqual(text, "SO:0000316", FALSE));
}

static gboolean isUtr(char *text)
{
  return (stringsEqual(text, "UTR", FALSE) || stringsEqual(text, "SO:0000203", FALSE));
}

static gboolean isExon(char *text)
{
  return (stringsEqual(text, "exon", FALSE) || stringsEqual(text, "SO:0000147", FALSE));
}

static gboolean isIntron(char *text)
{
  return (stringsEqual(text, "intron", FALSE) || stringsEqual(text, "SO:0000188", FALSE));
}

static gboolean isSnp(char *text)
{
  return (stringsEqual(text, "SNP", FALSE) || stringsEqual(text, "SO:0000694", FALSE));
}

static gboolean isPolyATail(char *text)
{
  return (stringsEqual(text, "polyA_sequence", FALSE) || stringsEqual(text, "SO:0000610", FALSE));
}



/* Get the MSP type from a string continaing a valid GFF3 type name */
static void parseMspType(char *token, MSP *msp, GError **error)
{
  msp->type = BLXMSP_INVALID;

  if (isMatchSet(token))
    {
      msp->type = BLXMSP_MATCH_SET;
    }
  else if (isMatch(token))
    {
      msp->type = BLXMSP_MATCH;
    }
  else if (isExon(token))
    {
      msp->type = BLXMSP_EXON_UNK;
    }
  else if (isCds(token))
    {
      msp->type = BLXMSP_EXON_CDS;
    }
  else if (isUtr(token))
    {
      msp->type = BLXMSP_EXON_UTR;
    }
  else if (isIntron(token))
    {
      msp->type = BLXMSP_INTRON;
    }
  else if (isSnp(token))
    {
      msp->type = BLXMSP_SNP;
    }
  else if (isPolyATail(token))
    {
      msp->type = BLXMSP_POLYA_TAIL;
    }
  else
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_TYPE, "Unexpected MSP type '%s' in input file.", token);
    }
}


/* Get the strand enum from a string containing '+' for the forward strand, '-' for the
 * reverse strand or '.' if no strand is specified */
static BlxStrand readStrand(char *token, GError **error)
{
  BlxStrand result = BLXSTRAND_NONE;
  
  if (!strcmp(token, "+"))
    {
      result = BLXSTRAND_FORWARD;
    }
  else if (!strcmp(token, "-"))
    {
      result = BLXSTRAND_REVERSE;
    }
  else if (!strcmp(token, "."))
    {
      result = BLXSTRAND_NONE;
    }
  else
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_STRAND, "Invalid strand '%s' in input file.\n", token);
    }
  
  return result;
}


/* Parse the columns in a GFF line and populate the parsed info into the given MSP. */
static void parseGffColumns(GString *line_string, const int lineNum, const char *opts, MSP *msp, GList **seqList, GSList *styles, GError **error)
{
    /* Split the line into its tab-separated columns. We should get 9 of them */
  char **tokens = g_strsplit_set(line_string->str, "\t", -1);   /* -1 means do all tokens. */
  
  /* This error should get set if there is a fatal error reading this line. */
  GError *tmpError = NULL;

  validateNumTokens(tokens, 9, 9, &tmpError);

  if (!tmpError)
    {
      msp->qname = g_ascii_strup(tokens[0], -1);
      msp->source = g_strdup(tokens[1]);
      msp->style = getBlxStyle(msp->source, styles, &tmpError);
      
      if (tmpError)
	{
	  /* An error getting the style is not critical so report it and clear the error */
	  prefixError(tmpError, "[line %d] Error getting style for MSP (default styles will be used instead). ", lineNum);
	  reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING); 
	}
	
      parseMspType(tokens[2], msp, &tmpError);
    }
    
  if (!tmpError)
    {
      int qStart = convertStringToInt(tokens[3]);
      int qEnd = convertStringToInt(tokens[4]);
      intrangeSetValues(&msp->qRange, qStart, qEnd);
      
      msp->score = convertStringToInt(tokens[5]);
      msp->qStrand = readStrand(tokens[6], &tmpError);
    }

  if (!tmpError)
    {
      msp->qFrame = convertStringToInt(tokens[7]);
  
      /* Parse the optional attributes */
      char *attributes = tokens[8];
      parseAttributes(attributes, msp, seqList, opts, lineNum, &tmpError);
    }
    
  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
}


/* Parse the given text, which contains attributes of the format "tag=data". The data
 * can contain multiple values, separated by spaces. Space characters within the data must 
 * be escaped. Populates the match sequence into 'sequence' if found in one of the attributes. */
static void parseAttributes(char *attributes, MSP *msp, GList **seqList, const char *opts, const int lineNum, GError **error)
{
  /* Attributes are separated by semi colons */
  char **tokens = g_strsplit_set(attributes, ";", -1);   /* -1 means do all tokens. */

  /* Loop through all the tags and read their data. The Target data will be read into 
   * the MSP and a BlxSequence (either a new or existing one). The actual sequence data will be
   * read from a different tag, but it requires that the Target tag is processed first so that the
   * BlxSequence exists. Remember both the BlxSequence and the sequence data so we can consolidate 
   * them. Also, the "gaps" string cannnot be processed until the Target tag has been processed because
   * we need certain information about the match sequence before we know how to process it, so just
   * remember the gaps string on the first pass through. */
  char *sequence = NULL;
  BlxStrand strand = BLXSTRAND_NONE;
  char *gapString = NULL;
  char **token = tokens;
  GError *tmpError = NULL;
  
  while (token && *token && **token && !tmpError)
    {
      parseTagDataPair(*token, msp, &sequence, &strand, &gapString, seqList, lineNum, &tmpError);
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
      ++token;
    }

  /* Create the BlxSequence for this MSP (or add it to an existing one of the same name, if exists) */
  if (!tmpError)
    {
      addBlxSequence(msp, strand, seqList, sequence, &tmpError);
    }
  
  /* Parse the gaps string */
  if (!tmpError && gapString)
    {
      parseGapString(gapString, opts, msp, &tmpError);
    }
  else if (gapString)
    {
      g_free(gapString);
      gapString = NULL;
    }

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
}


/* Parse a tag/data pair of the format "tag=data" */
static void parseTagDataPair(char *text, 
                             MSP *msp, 
                             char **sequence, 
			     BlxStrand *strand,
                             char **gapString,
                             GList **seqList, 
                             const int lineNum,
                             GError **error)
{
  DEBUG_ENTER("parseTagDataPair(text='%s')", text);
              
  /* Split on the "=" and check that we get 3 tokens */
  char **tokens = g_strsplit_set(text, "=", -1);
  
  GError *tmpError = NULL;
  validateNumTokens(tokens, 2, 2, &tmpError);
  
  if (!tmpError)
    {
      /* Call the relevant function to parse data for this tag */
      if (!strcmp(tokens[0], "Name"))
        {
          parseNameTag(tokens[1], msp, lineNum, &tmpError);
        }
      else if (!strcmp(tokens[0], "Target"))
        {
          parseTargetTag(tokens[1], msp, strand, seqList, lineNum, &tmpError);
        }
      else if (!strcmp(tokens[0], "sequence"))
        {
          *sequence = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "Gap"))
        {
          *gapString = g_strdup(tokens[1]);
        }
      else
        {
          DEBUG_OUT("Unknown tag: ignorning.\n");
        }
    }
  
  if (tmpError)
    {
      prefixError(tmpError, "Error processing data for tag '%s'. ", text);
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
  
  DEBUG_EXIT("parseTagDataPair");
}


/* Parse the data from the 'Name' tag */
static void parseNameTag(char *data, MSP *msp, const int lineNum, GError **error)
{
  if (!msp->sname)
    {
      msp->sname = g_ascii_strup(data, -1);
    }
  else if (!stringsEqual(data, msp->sname, FALSE))
    {
      g_warning("[line %d] Warning: Name attribute '%s' differs from previously-set MSP name '%s'. Ignoring new value.\n", lineNum, data, msp->sname);
    }
}


/* Parse the data from a 'Target' tag */
static void parseTargetTag(char *data, MSP *msp, BlxStrand *strand, GList **seqList, const int lineNum, GError **error)
{
  /* Split on spaces */
  char **tokens = g_strsplit_set(data, " ", -1); /* -1 means all tokens */
  
  GError *tmpError = NULL;
  int numTokens = validateNumTokens(tokens, 3, 4, &tmpError);
  
  if (!tmpError)
    {
      if (!msp->sname)
        {
          msp->sname = g_ascii_strup(tokens[0], -1);
        }
      else if (!stringsEqual(msp->sname, tokens[0], FALSE))
        {
          g_warning("[line %d] Warning: Target name '%s' differs from previously-set MSP name '%s'. Overriding old value.\n", lineNum, tokens[0], msp->sname);

          /* It's easiest if the Target tag overrides other values, because this is where we set
           * the name in the BlxSequence. */
          g_free(msp->sname);
          msp->sname = g_ascii_strup(tokens[0], -1);
        }
      
      const int sStart = convertStringToInt(tokens[1]);
      const int sEnd = convertStringToInt(tokens[2]);
      intrangeSetValues(&msp->sRange, sStart, sEnd);
      
      if (numTokens == 4)
        {
          *strand = readStrand(tokens[3], &tmpError);
        }
      }
  
   if (tmpError)
     {
       prefixError(tmpError, "Error parsing 'Target' tag '%s'", data);
       g_propagate_error(error, tmpError);
     }
  
  g_strfreev(tokens);
}


/* Find out how many reading frames there are based on the options (1 for nucleotide matches, 
 * 3 for protein matches) */
static int getNumFrames(const char *opts, GError **error)
{
  int result = 1;
  
  switch (opts[BLXOPT_MODE])
    {
      case 'X': /* fall through */
      case 'L':
        result = 3; /* protein matches: 3 frames */
        break;
        
      case 'N': /* fall through */
      case 'T': /* fall through */
      case 'P':
        result = 1; /* nucleotide matches: 1 frame */
        break;
        
      default:
        g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_UNKNOWN_MODE, "Unknown blast mode '%c' in options. Expected X, L, N, T or P.\n", opts[BLXOPT_MODE]);
        break;
    };
  
  return result;
}


/* Parse the data from a 'Gap' tag, which uses the CIGAR format, e.g. "M8 D3 M6 I1 M6".
 * Frees the given gap string when finished with. */
static void parseGapString(char *text, const char *opts, MSP *msp, GError **error)
{
  if (!text)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_TAG, "Gaps string is empty.\n");
      return;
    }
  
  /* Split on spaces */
  char **tokens = g_strsplit_set(text, " ", -1); /* -1 means do all tokens */

  /* Find out how many reading frames there are, based on the blast mode */
  GError *tmpError = NULL;
  const int numFrames = getNumFrames(opts, &tmpError);

  /* getNumFrames() always returns something so if there was an error give a warning and try to continue */
  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
  
  /* Start at the MSP's min coord and increase values as we progress through
   * the cigar string IF the strand is forward. If it is the reverse strand, 
   * start at the max coord and decrease values. */
  const gboolean qForward = (mspGetRefStrand(msp) != BLXSTRAND_REVERSE);
  const gboolean sForward = (mspGetMatchStrand(msp) != BLXSTRAND_REVERSE);
  const int qDirection = qForward ? 1 : -1;
  const int sDirection = sForward ? 1 : -1;

  int q = qForward ? msp->qRange.min : msp->qRange.max;
  int s = sForward ? msp->sRange.min : msp->sRange.max;
  
  char **token = tokens;

  while (token && *token && **token && !tmpError)
    {
      parseCigarStringSection(*token, msp, qDirection, sDirection, numFrames, &q, &s, &tmpError);
      
      if (tmpError)
        {
          prefixError(tmpError, "Invalid CIGAR string '%s'. ", text);
        }

      ++token;
    }
  
  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
  g_free(text);
}


/* Parse a section from a CIGAR string, e.g. "M8" or "I2" or "D3" etc. and create an array of
 * ranges where bases match.
 *
 * For example, the match between the following sequences:
 *   Chr3  (reference)  1 CAAGACCTAAACTGGAT-TCCAAT  23
 *   EST23 (target)     1 CAAGACCT---CTGGATATCCAAT  21
 *
 * is represented as: "M8 D3 M6 I1 M6", and corresponds to the following ranges of matching coordinates:
 *   q = 1 to 8    matches   s = 1 to 8
 *   q = 12 to 17  matches   s = 9 to 14
 *   q = 18 to 23  matches   s = 16 to 21
 *
 * So to translate the CIGAR strings to ranges, we do the following:
 *   M8 indicates a match of 8 places, so we increase both coords by 8;
 *   I2 indicates an insertion in the Subject sequence of 2 bases, so we increase the s coord by 2;
 *   D3 indicates a deletion from the Subject sequence of 3 bases, so we increase the q coord by 3. 
 */
static void parseCigarStringSection(const char *text, 
                                    MSP *msp, 
                                    const int qDirection, 
                                    const int sDirection, 
                                    const int numFrames,
                                    int *q, 
                                    int *s,
                                    GError **error)
{
  /* Get the digit part of the string */
  const int numBases = convertStringToInt(text+1);
  
  if (text[0] == 'M' || text[0] == 'm')
    {
      /* Create a match range between the old coords and the new */
      int newQ = *q + qDirection * ( ((numBases - 1) * numFrames) + numFrames - 1);  /* convert to nucleotide coords */
      int newS = *s + sDirection *(numBases - 1);
      
      CoordRange *newRange = g_malloc(sizeof(CoordRange));
      msp->gaps = g_slist_append(msp->gaps, newRange);
      
      newRange->qStart = *q;
      newRange->qEnd = newQ;
      newRange->sStart = *s;
      newRange->sEnd = newS;
      
      *q = newQ;
      *s = newS;
    }
  else if (text[0] == 'D' || text[0] == 'd')
    {
      *q += qDirection * ( ((numBases + 1) * numFrames) - numFrames + 1);    /* convert to nucleotide coords */
      *s += sDirection;
    }
  else if (text[0] == 'I' || text[0] == 'i')
    {
      *q += qDirection * (numFrames - numFrames + 1);    /* convert to nucleotide coords */
      *s += sDirection * (numBases + 1);
    }
  else if (text[0] == 'f' || text[0] == 'F' || text[0] == 'r' || text[0] == 'R')
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Blixem does not handle Frameshift operators.\n");
    }
  else
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Unknown operator '%c'.\n", text[0]);
    }
}

                           
/* Validate the number of tokens in a given list. Checks that there are between 'min' and 'max' 
 * tokens. If not, set the error. Returns the number of tokens. */
static int validateNumTokens(char **tokens, const int minReqd, const int maxReqd, GError **error)
{
  char **token = tokens;
  int count = 0;
  
  while (token && *token)
    {
      ++count;
      ++token; 
    }
  
  if (count < minReqd || count > maxReqd)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_NUM_TOKENS, "Expected between %d and %d columns but found %d\n", minReqd, maxReqd, count);
    }
  
  return count;
}


/* Check the given MSP has been populated with valid data. Sets the error if not. */
static void validateMsp(const MSP *msp, GError **error)
{
  g_assert(error && !*error);
  
  if (!msp)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "MSP is null.\n");
    }
  else if (msp->type == BLXMSP_INVALID)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "MSP type is invalid.\n");
    }
  else if (!msp->qname)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Reference sequence for MSP is not set.\n");
    }
  else if (msp->qRange.min == UNSET_INT || msp->qRange.max == UNSET_INT)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Reference sequence coords are not set.\n");
    }
  else if (msp->qStrand == BLXSTRAND_NONE)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Reference sequence strand not specified.\n");
    }
  else if (mspHasSName(msp) && !msp->sname)
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Match sequence name is null.\n");
    }
  else if (mspHasSCoords(msp) && (msp->sRange.min == UNSET_INT || msp->sRange.max == UNSET_INT))
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Match sequence coords are not set.\n");
    }
  else if (mspHasSStrand(msp) && mspGetMatchStrand(msp) == BLXSTRAND_NONE) /* if coords are set, we also need to know the direction */
    {
      g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_MSP, "Match sequence strand is not set.\n");
    }
}
    
    
