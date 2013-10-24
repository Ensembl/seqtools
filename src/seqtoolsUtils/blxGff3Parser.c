/*  File: blxGff3parser.c
 *  Author: Gemma Barson
 *  Copyright (c) 2010 - 2012 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
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
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: See blxGff3parser.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/blxGff3Parser.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <ctype.h>


#define SOURCE_DATA_TYPES_GROUP "source-data-types" /* group name for stanza where default data types are specified for sources */


/* Error codes and domain */
#define BLX_GFF3_ERROR g_quark_from_string("GFF 3 parser")

typedef enum {
  BLX_GFF3_ERROR_INVALID_STRAND,	      /* invali strand in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_TYPE,                /* invalid type in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_NUM_TOKENS,          /* invalid number of columns from a line of the input file */
  BLX_GFF3_ERROR_INVALID_TAG,                 /* invalid format for a tag/data pair */
  BLX_GFF3_ERROR_INVALID_SEQ,                 /* invalid sequence data */
  BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT,        /* invalid CIGAR format */
  BLX_GFF3_ERROR_INVALID_MSP,                 /* MSP has invalid/missing data */
  BLX_GFF3_ERROR_UNKNOWN_MODE,                /* unknown blast mode */
  BLX_GFF3_ERROR_BAD_COLOR,                   /* Bad color string found when parsing color */
  BLX_GFF3_ERROR_OUT_OF_RANGE,                /* Feature is not in the reference sequence range */
  BLX_GFF3_ERROR_DATA_TYPE                   /* Error finding data type */
} BlxGff3Error;


typedef enum {
  BLX_GAP_STRING_INVALID,
  BLX_GAP_STRING_GFF3,                        /* The Gap string used in GFF3 e.g. M23 D3 M10 I1 M20 */
  BLX_GAP_STRING_BAM_CIGAR                    /* The cigar format used by SAM/BAM, e.g. 23M3D10M1I20M */
} BlxGapFormat;


/* Utility struct to compile GFF fields and attributes into */
typedef struct _BlxGffData 
  {
    /* standard fields */
    char *qName;	/* ref seq name */
    char *source;	/* source */
    BlxMspType mspType;	/* type (converted to display type) */
    int qStart;		/* start coord on the ref seq */
    int qEnd;		/* end coord on the ref seq */
    gdouble score;	/* score */
    gdouble percentId;  /* percent ID */
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
    char *gapString;	/* the gap string (cigar) */
    BlxGapFormat gapFormat;    /* the format of the gap string */
    GQuark dataType;    /* represents a string that should correspond to a data type in the config file */
    GQuark filename;    /* optional filename e.g. for fetching data from a bam file */
  } BlxGffData;


/* Data used while parsing a gap string */
typedef struct _GapStringData
{
  BlxGapFormat gapFormat; /* the type of gap string e.g. cigar_bam */
  MSP **msp;              /* the msp that the gap string is for */
  int qDirection;         /* the direction we're parsing the reference (query) sequence: 1 for forward strand, -1 for rev */
  int sDirection;         /* the direction we're parsing the match (subject) sequence: 1 for forward strand, -1 for rev */
  int resFactor;          /* residue factor (3 for a peptide sequence, 1 for dna) */
  int *q;                 /* the current reference (query) sequence coord */
  int *s;                 /* the current match (subject) sequence coord */
  GArray **featureLists;  /* the array of lists of feature (one list per feature type) */
  MSP **lastMsp;          /* the last msp in the main msp list */
  MSP **mspList;          /* the main msp list */
  GList **seqList;        /* the list of sequence structs */
  GError *error;          /* gets set if there is an error */
} GapStringData;


static void           parseGffColumns(GString *line_string, const int lineNum, GList **seqList, GSList *supportedTypes, const IntRange const *refSeqRange, BlxGffData *gffData, GError **error);
static void           parseAttributes(char *attributes, GList **seqList, const int lineNum, BlxGffData *gffData, GError **error);
static void           parseTagDataPair(char *text, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseNameTag(char *data, char **sName, const int lineNum, GError **error);
static void           parseTargetTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseSequenceTag(const char *text, const int lineNum, BlxGffData *gffData, GError **error);
static void           parseGapString(char *text, BlxGapFormat gapFormat, MSP *msp, const int resFactor, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GError **error);

static BlxStrand      readStrand(char *token, GError **error);
//static void           parseMspType(char *token, MSP *msp, GSList *supportedTypes, GError **error);
static void           parseCigarStringSection(const char *text, GapStringData *data);
static int            validateNumTokens(char **tokens, const int minReqd, const int maxReqd, GError **error);
//static void           validateMsp(const MSP *msp, GError **error);
static void           addGffType(GSList **supportedTypes, char *name, char *soId, BlxMspType blxType);
static void           destroyGffType(BlxGffType **gffType);


/* utility to free a given string and set it to null */
static void freeAndNullString(char **ptr)
{
  g_free(*ptr);
  *ptr = NULL;
}

/* free all the memory used by the given gffdata struct (but not the struct itself) */
static void freeGffData(BlxGffData *gffData)
{
  freeAndNullString(&gffData->qName);
  freeAndNullString(&gffData->sName);
  freeAndNullString(&gffData->source);
  freeAndNullString(&gffData->idTag);
  freeAndNullString(&gffData->parentIdTag);
  freeAndNullString(&gffData->gapString);
}


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
  
  addGffType(&supportedTypes, "SNP", "SO:0000694", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "copy_number_variation", "SO:0001019", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "substitution", "SO:1000002", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "insertion", "SO:0000694", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "deletion", "SO:0000694", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "sequence_alteration", "SO:0001059", BLXMSP_VARIATION);

  addGffType(&supportedTypes, "polyA_signal_sequence", "SO:0000551", BLXMSP_POLYA_SIGNAL);
  addGffType(&supportedTypes, "polyA_site", "SO:0000553", BLXMSP_POLYA_SITE);

  addGffType(&supportedTypes, "read", "SO:0000150", BLXMSP_MATCH);
  addGffType(&supportedTypes, "read_PAIR", "SO:0000007", BLXMSP_MATCH);
  addGffType(&supportedTypes, "similarity", "SO:0000150", BLXMSP_MATCH); /* not a true gff type but temp fix because it gets put in gff by bam-get script */

  addGffType(&supportedTypes, "region", "SO:0000001", BLXMSP_REGION);

  supportedTypes = g_slist_reverse(supportedTypes);
  
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
		     GString *line_string, 
		     GList **seqList,
                     char *refSeqName,
                     IntRange *refSeqRange)
{
  DEBUG_ENTER("parseGff3Header [line=%d]", lineNum);
  
  /* Look for the "sequence-region" comment line, which tells us info about the reference
   * sequence. The format is as follows: ##sequence-region    qname qstart qend */
  char qName[MAXLINE + 1];
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;

  if (!strncasecmp(line_string->str, "##sequence-region", 17))
    {
      if (sscanf(line_string->str, "##sequence-region%s%d%d", qName, &qStart, &qEnd) < 3)
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
      
      if (refSeqRange)
	intrangeSetValues(refSeqRange, qStart, qEnd);
    }
  
  DEBUG_EXIT("parseGff3Header");
}


/* Look in the given config file to see if there is a default
 * datatype specified for the given source. If so, return its
 * name as a quark. Otherwise, return 0. */
static GQuark getBlxDataTypeDefault(const char *source, GKeyFile *keyFile)
{
  GQuark dataType = 0;

  if (keyFile && g_key_file_has_group(keyFile, SOURCE_DATA_TYPES_GROUP))
    {
      char *dataTypeName = g_key_file_get_string(keyFile, SOURCE_DATA_TYPES_GROUP, source, NULL);
      
      if (dataTypeName)
        {
          dataType = g_quark_from_string(dataTypeName);
          g_free(dataTypeName);
        }
    }

  return dataType;
}


/* Get the value for the given flag for the given group, and set 
 * it in the datatype if found */
static void getMspFlag(GKeyFile *keyFile, const char *group, const MspFlag flag, BlxDataType *dataType)
{
  /* Get the config-file key to use for this flag */
  const char *key = mspFlagGetConfigKey(flag);

  if (key)
    {
      GError *tmpError = NULL;
      gboolean result = g_key_file_get_boolean(keyFile, group, key, &tmpError);
      
      /* If found, update the value in the dataType */
      if (!tmpError)
        dataType->flags[flag] = result;
    }
}


/* Get the BlxDataType with the given name. Returns null and sets the error if 
 * we expected to find the name but didn't. */
BlxDataType* getBlxDataType(GQuark dataType, const char *source, GKeyFile *keyFile, GError **error)
{
  BlxDataType *result = NULL;

  /* If no data type was specified in the gff, see if there is a default for this gff type */
  if (!dataType)
    dataType = getBlxDataTypeDefault(source, keyFile);

  /* A keyfile might not be supplied if the calling program is not interested
   * in the data-type data (i.e. the data-type data is currently only used to
   * supply values that are used in blixem, so are irrelevant to dotter). */
  if (!keyFile || !dataType)
    return result;
  
  static GHashTable *dataTypes = NULL;

  if (!dataTypes)
    dataTypes = g_hash_table_new(g_direct_hash, g_direct_equal);

  if (dataType)
    {
      /* See if we've already got a struct for this data-type */
      result = (BlxDataType*)g_hash_table_lookup(dataTypes, GINT_TO_POINTER(dataType));

      if (!result)
        {
          /* We haven't requested this datatype before; look it up in the config file
           * and if we find it then create a new BlxDataType struct for it. */
          const gchar *typeName = g_quark_to_string(dataType);

          if (g_key_file_has_group(keyFile, typeName))
            {
              result = createBlxDataType();
              result->name = dataType;

              /* Get the values. They're all optional so just ignore any errors. */
              result->bulkFetch = keyFileGetCsv(keyFile, typeName, SEQTOOLS_BULK_FETCH, NULL); 
              result->userFetch = keyFileGetCsv(keyFile, typeName, SEQTOOLS_USER_FETCH, NULL); 
              
              /* Get the flags. Again, they're all optional. These calls update the
               * flag in place if it is found, or leave it at the pre-set default otherwise. */
              MspFlag flag = MSPFLAG_MIN + 1;
              for ( ; flag < MSPFLAG_NUM_FLAGS; ++flag)
                {
                  getMspFlag(keyFile, typeName, flag, result);
                }              
              
              /* Insert it into the table of data types */
              g_hash_table_insert(dataTypes, GINT_TO_POINTER(dataType), result);
            }
          else
            {
              g_set_error(error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_DATA_TYPE, 
                          "Config file error: data type '%s' not found", typeName);
            }
        }
    }
  
  return result;
}


/* Return the filename from the gff if given, otherwise check 
 * if the filename is given in the config and return that.
 * Returns 0 if not found. */
static GQuark getFeatureFilename(BlxGffData *gffData, GKeyFile *keyFile, GError **error)
{
  GQuark result = 0;
  
  if (gffData->filename)
    {
      /* Filename was given in the gff so use that */
      result = gffData->filename;
    }
  else if (keyFile && gffData->source)
    {
      /* Check if filename is given in the keyfile for this source */
      char *filename = g_key_file_get_string(keyFile, gffData->source, SEQTOOLS_GFF_FILENAME_KEY, error);
      result = g_quark_from_string(filename);
      g_free(filename);
    }
  
  return result;
}


/* Create a blixem object from the given parsed GFF data. Creates an MSP if the type is
 * exon or match, or a BlxSequence if the type is transcript. Does nothing for other types. */
static void createBlixemObject(BlxGffData *gffData, 
                               GArray* featureLists[],
			       MSP **lastMsp, 
			       MSP **mspList, 
			       GList **seqList, 
			       GSList *styles,
                               const int resFactor,
                               GKeyFile *keyFile,
			       GError **error)
{
  if (!gffData)
    {
      return;
    }
    
  GError *tmpError = NULL;

  /* Get the data type struct */
  BlxDataType *dataType = getBlxDataType(gffData->dataType, gffData->source, keyFile, &tmpError);
  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);

  GQuark filename = getFeatureFilename(gffData, keyFile, NULL);

  if (gffData->mspType > BLXMSP_NUM_TYPES)
    {
      /* "Invalid" MSP types, i.e. don't create a real MSP from these types. */
      
      if (gffData->mspType == BLXMSP_TRANSCRIPT)
        {
          /* For transcripts, although we don't create an MSP we do create a sequence */
          addBlxSequence(gffData->sName, gffData->idTag, gffData->qStrand,
                         dataType, gffData->source, seqList, gffData->sequence, NULL, 
                         &tmpError);
        }
    }
  else
    {
      /* For all other types, create an MSP */
      
      /* Regions don't necessarily have a name or ID, but they should have a
       * source, so use that as the name */
      if (gffData->mspType == BLXMSP_REGION && !gffData->sName)
        {
          gffData->sName = g_strdup(gffData->source);
        }
      
      if (!gffData->sName && !gffData->parentIdTag && 
	  (gffData->mspType == BLXMSP_TRANSCRIPT || typeIsExon(gffData->mspType) || 
	   typeIsMatch(gffData->mspType)))
	{
	  g_set_error(error, BLX_ERROR, 1, "Target/name/parent-ID must be specified for exons and alignments.\n");
	  return;
	}
	
      /* Get the id corresponding to the BlxSequence: we want the parent if it's an intron/exon, or the
       * ID for this item if it's a transcript or match */
      char *idTag = typeIsExon(gffData->mspType) || typeIsIntron(gffData->mspType) ? gffData->parentIdTag : gffData->idTag;
      
      /* For exons and transcripts, the target strand is irrelevant - use the ref seq strand */
      if (typeIsExon(gffData->mspType) || typeIsIntron(gffData->mspType) || gffData->mspType == BLXMSP_TRANSCRIPT)
	{
	  gffData->sStrand = gffData->qStrand;
	}

      MSP *msp = createNewMsp(featureLists,
                              lastMsp, 
			      mspList, 
			      seqList, 
			      gffData->mspType,
                              dataType,
			      gffData->source,
			      gffData->score, 
			      gffData->percentId, 
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
                              filename,
			      &tmpError);

    if (!tmpError)
	{ 
	  /* Get the style based on the source */
	  msp->style = getBlxStyle(gffData->source, styles, &tmpError);
          
	  if (tmpError)
	    {
	      /* style errors are not critical */
	      //reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
              g_error_free(tmpError);
              tmpError = NULL;
	    }

	  /* populate the gaps array */
	  parseGapString(gffData->gapString, gffData->gapFormat, msp, resFactor, featureLists, lastMsp, mspList, seqList, &tmpError);
	}
    }

  freeGffData(gffData);

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
}


/* Parse GFF3 data */
void parseGff3Body(const int lineNum,
                   GArray* featureLists[],
                   MSP **lastMsp, 
		   MSP **mspList, 
		   BlxParserState *parserState, 
		   GString *line_string, 
		   GList **seqList,
                   GSList *supportedTypes,
                   GSList *styles,
                   const int resFactor, 
                   GKeyFile *keyFile,
                   const IntRange const *refSeqRange)
{
  DEBUG_ENTER("parseGff3Body [line=%d]", lineNum);
  
  /* Parse the data into a temporary struct */
  BlxGffData gffData = {NULL, NULL, BLXMSP_INVALID, UNSET_INT, UNSET_INT, UNSET_INT, UNSET_INT, BLXSTRAND_NONE, UNSET_INT,
			NULL, BLXSTRAND_NONE, UNSET_INT, UNSET_INT, NULL, NULL, NULL, NULL, BLX_GAP_STRING_INVALID, 0, 0};
		      
  GError *error = NULL;
  parseGffColumns(line_string, lineNum, seqList, supportedTypes, refSeqRange, &gffData, &error);
  
  /* Create a blixem object based on the parsed data */
  if (!error)
    {
      createBlixemObject(&gffData, featureLists, lastMsp, mspList, seqList, styles, resFactor, keyFile, &error);
    }
  
  if (error)
    {
      prefixError(error, "[line %d] Error parsing GFF data. ", lineNum);
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
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
                         char **refSeq, char *refSeqName, IntRange *refSeqRange,
                         char ***readSeq, int *readSeqLen, int *readSeqMaxLen,
                         BlxParserState *parserState)
{
  char seqName[MAXLINE + 1];
  
  /* Read the ref seq name (and optionally the coords) from the header line */
  int startCoord = UNSET_INT, endCoord = UNSET_INT;
  const int numFound = sscanf(line, ">%s %d %d", seqName, &startCoord, &endCoord);
  
  if (numFound < 1)
    {
      /* Didn't find name - this is required */
      g_error("Error parsing data file: FASTA_SEQ_HEADER line \"%s\" is the wrong format; expected '>seq_name [start_coord end_coord].\n", line);
    }
  else if (numFound == 3 && refSeqRange)
    {
      /* We found the coords too. These override anything that was set from the
       * ##sequence_header line, which is fine; however, we don't currently 
       * handle the case where there might be more than one fasta sequence in 
       * the GFF file (we will end up overwriting the coords each time...). */
      intrangeSetValues(refSeqRange, startCoord, endCoord);
    }

  /* Set the name. Again, this will overwrite anything already set. */
  if (*refSeqName == '\0')
    {
      strcpy(refSeqName, seqName);
    }
  else if (!stringsEqual(refSeqName, seqName, FALSE))
    {
      g_warning("Reference sequence name was previously set to '%s' but the name in the FASTA data is '%s'; the name will overwritten with the new value.\n", refSeqName, seqName);
      strcpy(refSeqName, seqName);
    }

  /* Now allocate memory for the sequence data (if the sequence is not already populated) */
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
  else if (!strcmp(token, ".") || !strcmp(token, "?"))
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
                            GList **seqList,
                            GSList *supportedTypes,
                            const IntRange const *refSeqRange,
			    BlxGffData *gffData,
                            GError **error)
{
    /* Split the line into its tab-separated columns. We should get 9 of them */
  char **tokens = g_strsplit_set(line_string->str, "\t", -1);   /* -1 means do all tokens. */
  
  /* This error should get set if there is a fatal error reading this line. */
  GError *tmpError = NULL;

  validateNumTokens(tokens, 8, 9, &tmpError);

  if (!tmpError)
    {
      /* Reference sequence name */
      gffData->qName = tokens[0] ? g_ascii_strup(tokens[0], -1) : NULL;
      
      /* Source (optional) */
      if (tokens[1] && strcmp(tokens[1], "."))
          {
            gffData->source = g_uri_unescape_string(tokens[1], NULL);
          }
      
      /* Type (converted to a seqtools type) */
      gffData->mspType = getBlxType(supportedTypes, tokens[2], &tmpError);
    }
    
  if (!tmpError)
    {
      /* Reference sequence coords - ignore anything not in refSeqRange.
       * Note though that we accept features that are partially in range, and
       * also we currently accept all exons/introns. This is because we may, say, 
       * be given the exons in a transcript and be expected to calculate the
       * introns ourselves, but if such an intron is not entirely within range 
       * then we need knowledge of the adjacent exon that is out of range in 
       * order to be able to calculate that intron. Rather than get into 
       * complicated filtering to include only those exons we need, we currently 
       * just include all exons and introns in the input file. Ideally we would
       * at least filter out transcripts that are entirely out of range, but we
       * don't fully know the parent/child relationship at this point, so that
       * is again getting quite tricky (and it will not cause problems if they
       * are left in). */
      gffData->qStart = convertStringToInt(tokens[3]);
      gffData->qEnd = convertStringToInt(tokens[4]);

      /* We can only check the range if the refseqrange is set... */
      if (!typeIsExon(gffData->mspType) && !typeIsIntron(gffData->mspType) && !typeIsTranscript(gffData->mspType) &&
          refSeqRange && (refSeqRange->min != UNSET_INT || refSeqRange->max != UNSET_INT))
        {
          IntRange featureRange;
          intrangeSetValues(&featureRange, gffData->qStart, gffData->qEnd); /* makes sure min < max */
          
          if (!rangesOverlap(&featureRange, refSeqRange))
            g_set_error(&tmpError, BLX_GFF3_ERROR, BLX_GFF3_ERROR_OUT_OF_RANGE, "Feature is outside the reference sequence range.\n");
        }
    }
  
  if (!tmpError)
    {
      if (!stringsEqual(tokens[5], ".", TRUE))
        {
          gffData->score = g_ascii_strtod(tokens[5], NULL);
        }
      
      gffData->qStrand = readStrand(tokens[6], &tmpError);
    }

  if (!tmpError)
    {
      if (stringsEqual(tokens[7], ".", TRUE))
        {
          gffData->phase = 0;
          
          if (gffData->mspType == BLXMSP_CDS)
            {
              g_warning("[line %d] CDS type does not have phase specified.\n", lineNum);
            }
        }
      else
        {
          gffData->phase = convertStringToInt(tokens[7]);
        }
  
      /* Parse the optional attributes */
      char *attributes = tokens[8];
      
      if (attributes)
        {
          parseAttributes(attributes, seqList, lineNum, gffData, &tmpError);
        }
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
      else if (!strcmp(tokens[0], "Gap"))
        {
          /* This might have already been set if we have more than one type of gap string */
          if (!gffData->gapString)
            {
              /*! \todo For now, override the cigar_bam string because we are experiencing
               * bugs with it. Longer term it shouldn't really matter which we use, although 
               * we may want to give preference to more informative gap strings e.g. vulgar */
              g_free(gffData->gapString);
              gffData->gapString = NULL;
              gffData->gapFormat = BLX_GAP_STRING_INVALID;
            }
          
          gffData->gapString = g_strdup(tokens[1]);
          gffData->gapFormat = BLX_GAP_STRING_GFF3; 
        }
      else if (!strcmp(tokens[0], "cigar_bam"))
        {
          /* This might have already been set if we have more than one type of gap string */
          if (!gffData->gapString)
            {
              gffData->gapString = g_strdup(tokens[1]);
              gffData->gapFormat = BLX_GAP_STRING_BAM_CIGAR; 
            }
        }
      else if (!strcmp(tokens[0], "ID"))
        {
          gffData->idTag = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "Parent"))
        {
	  gffData->parentIdTag = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "percentID"))
        {
          gffData->percentId = g_ascii_strtod(tokens[1], NULL);
        }
      else if (!strcmp(tokens[0], "sequence"))
        {
          parseSequenceTag(tokens[1], lineNum, gffData, &tmpError);
        }
      else if (!strcmp(tokens[0], "variant_sequence"))
        {
          gffData->sequence = g_strdup(tokens[1]);
        }
      else if (!strcmp(tokens[0], "dataType"))
        {
          gffData->dataType = g_quark_from_string(tokens[1]);
        }
      else if (!strcmp(tokens[0], "file"))
        {
          gffData->filename = g_quark_from_string(tokens[1]);
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


/* Parse the data from the 'sequence' tag */
static void parseSequenceTag(const char *text, const int lineNum, BlxGffData *gffData, GError **error)
{
  gffData->sequence = g_strdup(text);
}


/* Parse the data from the "gaps" string, which uses the CIGAR format, e.g. "M8 D3 M6 I1 M6".
 * Populates the Gaps array in the given MSP.*/
static void parseGapString(char *text,
                           BlxGapFormat gapFormat,
                           MSP *msp,
                           const int resFactor,
                           GArray* featureLists[],
                           MSP **lastMsp, 
                           MSP **mspList, 
                           GList **seqList, 
                           GError **error)
{
  if (!text || gapFormat == BLX_GAP_STRING_INVALID)
    {
      return;
    }
  
  /* Split on spaces */
  char **tokens = g_strsplit_set(text, " ", -1); /* -1 means do all tokens */

  /* If we have the forward strand of either sequence, start at the min coord
   * and increase values as we progress through the cigar string; if we have the
   * reverse strand, start at the max coord and decrease. */
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
  const int qDirection = (qForward ? 1 : -1);
  const int sDirection = (sForward ? 1 : -1);

  /* Start at one beyond the edge of the range, because it will be incremented (or decremented if
   * direction is reverse) when we construct the first range. */
  int q = qForward ? msp->qRange.min - 1 : msp->qRange.max + 1;
  int s = sForward ? msp->sRange.min - 1 : msp->sRange.max + 1;
  
  char **token = tokens;
  GError *tmpError = NULL;

  GapStringData gapStringData = {gapFormat, &msp, qDirection, sDirection, resFactor, &q, &s, 
                                 featureLists, lastMsp, mspList, seqList, NULL};  

  while (token && *token && **token)
    {
      parseCigarStringSection(*token, &gapStringData);
      
      if (tmpError)
        {
          prefixError(tmpError, "Error parsing gap string '%s'. ", text);
        }

      ++token;
    }
  
  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
  
  g_strfreev(tokens);
}


/* Get the length part of a gap string section, e.g. if the text is "M76"
 * then this returns 76 */
static int getCigarStringSectionLen(const char *text, BlxGapFormat gapFormat)
{
  int result = 0;

  switch (gapFormat)
    {
    case BLX_GAP_STRING_GFF3: /* e.g. M76 */
      result = convertStringToInt(text+1);
      break;
      
    case BLX_GAP_STRING_BAM_CIGAR: /* e.g. 76M */
      convertStringToInt(text); /* uses atoi, which will ignore characters after */
      break;

    default:
      g_warning("Invalid gap string format\n");
      break;
    };

  return result;
}


/* Get the operator part of a gap string section, e.g. if the text is "M76"
 * then this returns 'M' */
static char getCigarStringSectionOperator(const char *text, BlxGapFormat gapFormat)
{
  char result = 0;

  switch (gapFormat)
    {
    case BLX_GAP_STRING_GFF3:
      result = *text;
      break;
      
    case BLX_GAP_STRING_BAM_CIGAR:
      {
        const char *cp = text;
        for ( ; cp && *cp && !isalpha(*cp); ++cp); /* find first alphabetic character */
        if (cp) 
          result = *cp;
        break;
      }
      
    default:
      g_warning("Invalid gap string format\n");
      break;
    };
  
  return result;
}


static void parseCigarStringMatch(GapStringData *data, const int numNucleotides, const int numPeptides)
{
  MSP *msp = *data->msp;

  /* We were at the end of the previous range or gap, so move to the next coord, where our range will start */
  *data->q += data->qDirection;
  *data->s += data->sDirection;
  
  /* Find the coords at the end of the range. */
  int newQ = *data->q + (data->qDirection * (numNucleotides - 1));
  int newS = *data->s + (data->sDirection * (numPeptides - 1));
  
  CoordRange *newRange = g_malloc(sizeof(CoordRange));
  msp->gaps = g_slist_append(msp->gaps, newRange);
  
  newRange->qStart = *data->q;
  newRange->qEnd = newQ;
  newRange->sStart = *data->s;
  newRange->sEnd = newS;
  
  *data->q = newQ;
  *data->s = newS;
}

static void parseCigarStringIntron(GapStringData *data, const int numNucleotides, const int numPeptides)
{
  /* Intron. Create a separate msp under the same sequence. */
  MSP *msp = *data->msp;
  MSP *newMsp = copyMsp(msp, data->featureLists, data->lastMsp, data->mspList, data->seqList, &data->error);
  
  /* end current msp at the current coords */
  if (data->qDirection > 0)
    msp->qRange.max = *data->q;
  else
    msp->qRange.min = *data->q;
  
  if (data->sDirection > 0)
    msp->sRange.max = *data->s;
  else
    msp->sRange.min = *data->s;
  
  /* start new msp at new coords */
  *data->q += data->qDirection * numNucleotides;
  
  if (data->qDirection > 0)
    newMsp->qRange.min = *data->q + 1;
  else
    newMsp->qRange.max = *data->q - 1;
  
  if (data->sDirection > 0)
    newMsp->sRange.min = *data->s + 1;
  else
    newMsp->sRange.max = *data->s - 1;
  
  *data->msp = newMsp;
}

static void parseCigarStringDeletion(GapStringData *data, const int numNucleotides, const int numPeptides)
{
  /* Deletion from the subject sequence: increase the q coord by the number of nucleotides. */
  *data->q += data->qDirection * numNucleotides;
}

static void parseCigarStringInsertion(GapStringData *data, const int numNucleotides, const int numPeptides)
{
  /* Insertion on the subject sequence: increase the s coord by the number of peptides. */
  *data->s += data->sDirection * numPeptides;
}


/* Return TRUE if the given char is a valid operator for the
 * given gap string format. */
static gboolean validateCigarOperator(char operator, BlxGapFormat gapFormat)
{
  gboolean result = FALSE;
  
  switch (operator)
    {
    case 'M':    case 'm':
    case 'N':    case 'n':
    case 'D':    case 'd':
    case 'I':    case 'i':
      result = (gapFormat == BLX_GAP_STRING_GFF3 || gapFormat == BLX_GAP_STRING_BAM_CIGAR);
      break;
      
    case 'X':    case 'x':
    case 'P':    case 'p':
    case 'S':    case 's':
    case 'H':    case 'h':
      result = (gapFormat == BLX_GAP_STRING_BAM_CIGAR);
      break;

    case 'F':    case 'f':
    case 'R':    case 'r':
      result = (gapFormat == BLX_GAP_STRING_GFF3);
      break;

    default:
      result = FALSE;
      break;
    };

  return result;
}


/* Parse a section from a CIGAR string, e.g.

 "M8" or "I2" or "D3" etc. and create an array of
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
                                    GapStringData *data)
{
  /* Get the digit part of the string, which indicates the number of display coords (peptides in peptide matches,
   * nucleotides in nucelotide matches). */
  const int numPeptides = getCigarStringSectionLen(text, data->gapFormat);
  int numNucleotides = numPeptides * data->resFactor;

  char operator = getCigarStringSectionOperator(text, data->gapFormat);
  
  /*! \todo If the operator is not valid for this type of cigar string
   * then we should set the error and return. However, for historic 
   * reasons we have allowed invalid operators in the GFF Gap string,
   * so continue allowing this for now and just give a warning. */
  if (!validateCigarOperator(operator, data->gapFormat))
    g_warning("Invalid operator '%c' for gap string format '%d'", operator, data->gapFormat);

  switch (operator)
    {
    case 'M':
    case 'm':
      parseCigarStringMatch(data, numNucleotides, numPeptides);
      break;
      
    case 'N':
    case 'n':
      parseCigarStringIntron(data, numNucleotides, numPeptides);
      break;
      
    case 'D':
    case 'd':
      parseCigarStringDeletion(data, numNucleotides, numPeptides);
      break;

    case 'I':
    case 'i':
      parseCigarStringInsertion(data, numNucleotides, numPeptides);
      break;
      
    case 'X':
    case 'x':
      /* Mismatch. For now just treat it like a match. Blixem detects mistmatches
       * anyway as long as we have the reference sequence. */
      /*! \todo It might be useful to parse mistmatches in the gap string. This 
       * extra info may help e.g. to calculate percentID for matches that extend
       * beyond the reference sequence range. We could do that calculation here
       * if percentID was not given in the GFF? */
      parseCigarStringMatch(data, numNucleotides, numPeptides);
      break;

    case 'P':
    case 'p':
      /* Padding (silent deletion from the reference sequence). Not supported, but 
       * ignoring these essentially gives us the unpadded cigar. */
      break;

    case 'S':
    case 's':
      /* Soft clipping. This indicates whether the sequence aligns from the first
       * residue to the last one, but blixem detects this anyway, so we can ignore
       * this. */
      break;
      
    case 'H':
    case 'h':
      /* Hard clipping - not supported */
      /*! \todo We could probably implement this relatively easily. Hard clipping means
       * that bases are excluded from the start/end of the sequence (to stop repeating
       * the same bit of sequence when we have multiple matches from the same sequence).
       * We'd probably need to pad the sequence (making sure we merge the bits of sequence
       * from all matches for this sequence). We don't have any real examples of this yet so
       * leaving it for now. */
      g_set_error(&data->error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Blixem does not handle the hard-clipping operator.\n");
      break;

    case 'F':
    case 'f':
    case 'R':
    case 'r':
      /* Frameshift operators - not supported */
      g_set_error(&data->error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Blixem does not handle Frameshift operators.\n");
      break;
      
    default:
      g_set_error(&data->error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Invalid operator '%c' in cigar string.\n", operator);
      break;
    };
  
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



