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
#include <string.h>


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


/* Utility struct to compile GFF fields and attributes into */
typedef struct _BlxGffData 
  {
    /* standard fields */
    char *qName;	/* ref seq name */
    char *source;	/* source */
    BlxMspType mspType;	/* type */
    int qStart;		/* start coord on the ref seq */
    int qEnd;		/* end coord on the ref seq */
    gdouble score;	/* score */
    BlxStrand qStrand;	/* ref seq strand */
    int phase;		/* phase */
    
    /* Attributes */
    char *sName;	/* target name */
    BlxStrand sStrand;	/* target sequence strand */
    int sStart;		/* target start coord */
    int sEnd;		/* target end coord */
    char *idTag;	/* ID of the item */
    char *parentIdTag;	/* Parent ID of the item */
    char *sequence;	/* sequence data */
    char *gapString;	/* the gaps cigar string */
  } BlxGffData;



static void           parseGffColumns(GString *line_string, const int lineNum, const char *opts, GList **seqList, GSList *supportedTypes, GSList *styles, BlxGffData *gffData, GError **error);
static void           parseAttributes(char *attributes, GList **seqList, const char *opts, const int lineNum, BlxGffData *gffData, GError **error);
static void           parseTagDataPair(char *text, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseNameTag(char *data, char **sName, const int lineNum, GError **error);
static void           parseTargetTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseGapString(char *text, const char *opts, MSP *msp, GError **error);

static BlxStrand      readStrand(char *token, GError **error);
//static void           parseMspType(char *token, MSP *msp, GSList *supportedTypes, GError **error);
static void           parseCigarStringSection(const char *text, MSP *msp, const int qDirection, const int sDirection, const int numFrames, int *q, int *s, GError **error);
static int            validateNumTokens(char **tokens, const int minReqd, const int maxReqd, GError **error);
//static void           validateMsp(const MSP *msp, GError **error);
static void           addGffType(GSList **supportedTypes, char *name, char *soId, BlxMspType blxType);
static void           destroyGffType(BlxGffType **gffType);



/* Free all memory used by the given list of supported GFF types */
void blxDestroyGffTypeList(GSList **supportedTypes)
{
  GSList *item = *supportedTypes;
  
  for ( ; item; item = item->next)
    {
      BlxGffType *gffType = (BlxGffType*)(item->data);
      destroyGffType(&gffType);
    }
  
  g_slist_free(*supportedTypes);
  *supportedTypes = NULL;
}


/* Create the list of supported types */
GSList* blxCreateSupportedGffTypeList()
{
  GSList *supportedTypes = NULL;
  
  addGffType(&supportedTypes, "match", "SO:0000343", BLXMSP_MATCH);
  addGffType(&supportedTypes, "nucleotide_match", "SO:0000347", BLXMSP_MATCH);
  addGffType(&supportedTypes, "protein_match", "SO:0000349", BLXMSP_MATCH);
  addGffType(&supportedTypes, "match_part", "SO:0000039", BLXMSP_MATCH);
  
  addGffType(&supportedTypes, "match_set", "SO:0000038", BLXMSP_MATCH_SET);
  
  addGffType(&supportedTypes, "transcript", "SO:0000673", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "primary_transcript", "SO:0000185", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "processed_transcript", "SO:0000233", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "mRNA", "SO:0000234", BLXMSP_TRANSCRIPT);

  addGffType(&supportedTypes, "CDS", "SO:0000316", BLXMSP_CDS);
  addGffType(&supportedTypes, "UTR", "SO:0000203", BLXMSP_UTR);
  addGffType(&supportedTypes, "exon", "SO:0000147", BLXMSP_EXON);

  addGffType(&supportedTypes, "intron", "SO:0000188", BLXMSP_INTRON);
  
  addGffType(&supportedTypes, "SNP", "SO:0000694", BLXMSP_SNP);
  addGffType(&supportedTypes, "polyA_sequence", "SO:0000610", BLXMSP_POLYA_TAIL);
  
  return supportedTypes;
}


/* Get the internal blixem type from the given GFF type (as either the type name or
 * the SO id) */
static BlxMspType getBlxType(GSList *supportedTypes, const char *typeStr, GError **error)
{
  BlxMspType result = BLXMSP_INVALID;

  /* Loop through the supported types and see if the requested type matches the name or SO id */
  GSList *item = supportedTypes;
  
  for ( ; item; item = item->next)
    {
      BlxGffType *gffType = (BlxGffType*)(item->data);
      
      if (stringsEqual(typeStr, gffType->name, FALSE) || stringsEqual(typeStr, gffType->soId, FALSE))
        {
          result = gffType->blxType;
          break;
        }
    }
  
  /* Check if it was found... */
  if (result == BLXMSP_INVALID)
  {
    g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_TYPE, "Unsupported type '%s' in input file.\n", typeStr);
  }
  
  return result;
}


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
          g_error("Error parsing data file, type GFF_3_HEADER: \"%s\"\n", line_string->str);
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


/* Create a blixem object from the given parsed GFF data. Creates an MSP if the type is
 * exon or match, or a BlxSequence if the type is transcript. Does nothing for other types. */
static void createBlixemObject(BlxGffData *gffData, 
			       MSP **lastMsp, 
			       MSP **mspList, 
			       GList **seqList, 
			       GSList *styles,
			       char *opts, 
			       GError **error)
{
  if (!gffData)
    {
      return;
    }
    
  GError *tmpError = NULL;

  if (typeIsExon(gffData->mspType) || gffData->mspType == BLXMSP_MATCH)
    {
      if (!gffData->sName && !gffData->parentIdTag)
	{
	  g_set_error(error, BLX_ERROR, 1, "Target/name/parent-ID must be specified for exons and alignments.\n");
	  return;
	}
	
      /* Get the id corresponding to the BlxSequence: we want the parent if it's an exon, or the
       * ID for this item if it's a transcript or match */
      char *idTag = typeIsExon(gffData->mspType) ? gffData->parentIdTag : gffData->idTag;
      
      /* For exons and transcripts, the target strand is irrelevant - use the ref seq strand */
      if (gffData->mspType != BLXMSP_MATCH)
	{
	  gffData->sStrand = gffData->qStrand;
	}

      MSP *msp = createNewMsp(lastMsp, 
			      mspList, 
			      seqList, 
			      gffData->mspType, 
			      gffData->source,
			      gffData->score, 
                              gffData->phase,
			      idTag,
			      gffData->qName, 
			      gffData->qStart, 
			      gffData->qEnd, 
			      gffData->qStrand, 
			      UNSET_INT,
			      gffData->sName, 
			      gffData->sStart, 
			      gffData->sEnd, 
			      gffData->sStrand, 
			      gffData->sequence, 
			      opts, 
			      &tmpError);

    if (!tmpError)
	{ 
	  /* Get the style based on the source */
	  msp->style = getBlxStyle(gffData->source, styles, &tmpError);
	  
	  if (tmpError)
	    {
	      /* style errors are not critical */
	      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
	    }

	  /* populate the gaps array */
	  parseGapString(gffData->gapString, opts, msp, &tmpError);
	}    
    }
  else if (gffData->mspType == BLXMSP_TRANSCRIPT)
    {
      addBlxSequence(gffData->sName, gffData->idTag, gffData->qStrand, seqList, gffData->sequence, NULL, &tmpError);
    }
    
  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
}


/* Parse GFF3 data */
void parseGff3Body(const int lineNum,
                   MSP **lastMsp, 
		   MSP **mspList, 
		   BlxParserState *parserState, 
		   char *opts, 
		   GString *line_string, 
		   GList **seqList,
                   GSList *supportedTypes,
                   GSList *styles)
{
  DEBUG_ENTER("parseGff3Body [line=%d]", lineNum);
  
  /* Parse the data into a temporary struct */
  BlxGffData gffData = {NULL, NULL, BLXMSP_INVALID, UNSET_INT, UNSET_INT, UNSET_INT, BLXSTRAND_NONE, UNSET_INT,
			NULL, BLXSTRAND_NONE, UNSET_INT, UNSET_INT, NULL, NULL, NULL, NULL};
		      
  GError *error = NULL;
  parseGffColumns(line_string, lineNum, opts, seqList, supportedTypes, styles, &gffData, &error);
  
  /* Create a blixem object based on the parsed data */
  if (!error)
    {
      createBlixemObject(&gffData, lastMsp, mspList, seqList, styles, opts, &error);
    }
  
  if (error)
    {
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
      g_error("Error parsing data file, type FASTA_SEQ_HEADER: \"%s\"\n", line);
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
static void parseGffColumns(GString *line_string, 
                            const int lineNum, 
                            const char *opts, 
                            GList **seqList,
                            GSList *supportedTypes,
                            GSList *styles, 
			    BlxGffData *gffData,
                            GError **error)
{
    /* Split the line into its tab-separated columns. We should get 9 of them */
  char **tokens = g_strsplit_set(line_string->str, "\t", -1);   /* -1 means do all tokens. */
  
  /* This error should get set if there is a fatal error reading this line. */
  GError *tmpError = NULL;

  validateNumTokens(tokens, 9, 9, &tmpError);

  if (!tmpError)
    {
      gffData->qName = tokens[0] ? g_ascii_strup(tokens[0], -1) : NULL;
      gffData->source = tokens[1] && strcmp(tokens[1], ".") ? g_strdup(tokens[1]) : NULL;
      
      if (tmpError)
	{
	  /* An error getting the style is not critical so report it and clear the error */
	  prefixError(tmpError, "[line %d] Error getting style (default styles will be used instead). ", lineNum);
	  reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING); 
	}
	
      gffData->mspType = getBlxType(supportedTypes, tokens[2], &tmpError);
    }
    
  if (!tmpError)
    {
      gffData->qStart = convertStringToInt(tokens[3]);
      gffData->qEnd = convertStringToInt(tokens[4]);
      
      gffData->score = g_strtod(tokens[5], NULL);
      gffData->qStrand = readStrand(tokens[6], &tmpError);
    }

  if (!tmpError)
    {
      if (stringsEqual(tokens[7], ".", TRUE))
        gffData->phase = UNSET_INT;
      else
        gffData->phase = convertStringToInt(tokens[7]);
  
      /* Parse the optional attributes */
      char *attributes = tokens[8];
      parseAttributes(attributes, seqList, opts, lineNum, gffData, &tmpError);
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
static void parseAttributes(char *attributes, 
			    GList **seqList, 
			    const char *opts, 
			    const int lineNum, 
			    BlxGffData *gffData,
			    GError **error)
{
  /* Attributes are separated by semi colons */
  char **tokens = g_strsplit_set(attributes, ";", -1);   /* -1 means do all tokens. */

  /* Loop through all the tags and read their data. */
  char **token = tokens;
  GError *tmpError = NULL;
  
  while (token && *token && **token && !tmpError)
    {
      parseTagDataPair(*token, lineNum, seqList, gffData, &tmpError);
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
      ++token;
    }

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
}


/* Parse a tag/data pair of the format "tag=data" */
static void parseTagDataPair(char *text,
                             const int lineNum,
                             GList **seqList, 
			     BlxGffData *gffData, 
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
          parseNameTag(tokens[1], &gffData->sName, lineNum, &tmpError);
        }
      else if (!strcmp(tokens[0], "Target"))
        {
          parseTargetTag(tokens[1], lineNum, seqList, gffData, &tmpError);
        }
      else if (!strcmp(tokens[0], "sequence"))
        {
          gffData->sequence = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "Gap") || !strcmp(tokens[0], "Gaps")) /* to do: get rid of "Gaps" once zmap starts supporting the correct keyword "Gap" */
        {
          gffData->gapString = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "ID"))
        {
          gffData->idTag = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "Parent"))
        {
	  gffData->parentIdTag = g_strdup(tokens[1]);
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
static void parseNameTag(char *data, char **sName, const int lineNum, GError **error)
{
  if (data)
    { 
      if (*sName == NULL)
	{
	  *sName = g_ascii_strup(data, -1);
	}
      else if (!stringsEqual(data, *sName, FALSE))
	{
	  g_warning("[line %d] Warning: Name attribute '%s' differs from previously-set name '%s'. Ignoring new value.\n", lineNum, data, *sName);
	}
    }
}


/* Parse the data from a 'Target' tag */
static void parseTargetTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error)
{
  /* Split on spaces */
  char **tokens = g_strsplit_set(data, " ", -1); /* -1 means all tokens */
  
  GError *tmpError = NULL;
  int numTokens = validateNumTokens(tokens, 3, 4, &tmpError);
  
  if (!tmpError)
    {
      if (gffData->sName == NULL)
        {
          gffData->sName = tokens[0] ? g_ascii_strup(tokens[0], -1) : NULL;
        }
      else if (!stringsEqual(gffData->sName, tokens[0], FALSE))
        {
          g_warning("[line %d] Warning: Target name '%s' differs from previously-set name '%s'. Overriding old value.\n", lineNum, tokens[0], gffData->sName);

          /* It's easiest if the Target tag overrides other values, because this is where we set
           * the name in the BlxSequence. */
          g_free(gffData->sName);
          gffData->sName = tokens[0] ? g_ascii_strup(tokens[0], -1) : NULL;
        }
      
      gffData->sStart = convertStringToInt(tokens[1]);
      gffData->sEnd = convertStringToInt(tokens[2]);
      
      if (numTokens == 4)
        {
          gffData->sStrand = readStrand(tokens[3], &tmpError);
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


/* Parse the data from the "gaps" string, which uses the CIGAR format, e.g. "M8 D3 M6 I1 M6".
 * Populates the Gaps array in the given MSP. Frees the given gap string when finished with. */
static void parseGapString(char *text, const char *opts, MSP *msp, GError **error)
{
  if (!text)
    {
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

  /* Start at one beyond the edge of the range, because it will be incremented (or decremented if
   * direction is reverse) when we construct the first range. */
  int q = qForward ? msp->qRange.min - 1 : msp->qRange.max + 1;
  int s = sForward ? msp->sRange.min - 1 : msp->sRange.max + 1;
  
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
 *
 *  qDirection and sDirection are 1 if coords are in an increasing direction or -1 if decreasing.
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
  /* Get the digit part of the string, which indicates the number of peptides. */
  const int numPeptides = convertStringToInt(text+1);
  const int numNucleotides = numPeptides * numFrames;

//  /* to do: hack to make blixem work the with currently-wrong data from zmap. take this out and use the above. */
//  const int numNucleotides = convertStringToInt(text+1);
//  const int numPeptides= numNucleotides / numFrames;

  if (text[0] == 'M' || text[0] == 'm')
    {
      /* We were at the end of the previous range or gap, so move to the next coord, where our range will start */
      *q += qDirection;
      *s += sDirection;

      /* Find the coords at the end of the range. */
      int newQ = *q + (qDirection * (numNucleotides - 1));
      int newS = *s + (sDirection * (numPeptides - 1));
      
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
      /* Deletion from the subject sequence: increase the q coord by the number of nucleotides. */
      *q += qDirection * numNucleotides;
    }
  else if (text[0] == 'I' || text[0] == 'i')
    {
      /* Insertion on the subject sequence: increase the s coord by the number of peptides. */
      *s += sDirection * numPeptides;
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


/* Create a gff type with the given info and add it to the given list */
static void addGffType(GSList **supportedTypes, char *name, char *soId, BlxMspType blxType)
{
  BlxGffType *gffType = g_malloc(sizeof *gffType);
  
  gffType->name = g_strdup(name);
  gffType->soId = g_strdup(soId);
  gffType->blxType = blxType;
  
  *supportedTypes = g_slist_prepend(*supportedTypes, gffType);
}


/* Free all memory used by the given gff type */
static void destroyGffType(BlxGffType **gffType)
{
  if (gffType && *gffType)
    {
      g_free((*gffType)->name);
      g_free((*gffType)->soId);
      
      g_free(*gffType);
      *gffType = NULL;
    }
}



