/*  File: blxGff3parser.c
 *  Author: Gemma Barson
 *  Copyright [2018-2024] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

#include <seqtoolsUtils/blxGff3Parser.hpp>
#include <seqtoolsUtils/blxparser.hpp>
#include <seqtoolsUtils/utilities.hpp>
#include <seqtoolsUtils/blxmsp.hpp>
#include <seqtoolsUtils/seqtoolsFetch.hpp>
#include <string.h>
#include <ctype.h>
#include <string>
#include <map>

using namespace std;


/* globals */
static std::map<GQuark, BlxDataType*> g_dataTypes;


#define SOURCE_DATA_TYPES_GROUP "source-data-types" /* group name for stanza where default data types are specified for sources */
#define DATA_TYPE_TAG "dataType" /* tag name for dataType */


/* Error codes and domain */
#define BLX_GFF3_ERROR g_quark_from_string("GFF 3 parser")

typedef enum {
  BLX_GFF3_ERROR_INVALID_STRAND,	      /* invali strand in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_TYPE,                /* invalid type in GFF3 input file */
  BLX_GFF3_ERROR_INVALID_NUM_TOKENS,          /* invalid number of columns from a line of the input file */
  BLX_GFF3_ERROR_INVALID_TAG,                 /* invalid format for a tag/data pair */
  BLX_GFF3_ERROR_INVALID_SEQ,                 /* invalid sequence data */
  BLX_GFF3_ERROR_INVALID_SEQ_NAME,            /* invalid sequence name */
  BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT,        /* invalid CIGAR format */
  BLX_GFF3_ERROR_INVALID_MSP,                 /* MSP has invalid/missing data */
  BLX_GFF3_ERROR_INVALID_HEADER,              /* invalid header line */
  BLX_GFF3_ERROR_UNKNOWN_MODE,                /* unknown blast mode */
  BLX_GFF3_ERROR_BAD_COLOR,                   /* Bad color string found when parsing color */
  BLX_GFF3_ERROR_OUT_OF_RANGE,                /* Feature/file is not in the reference sequence range */
  BLX_GFF3_ERROR_DATA_TYPE                   /* Error finding data type */
} BlxGff3Error;


typedef enum {
  BLX_GAP_STRING_INVALID,
  BLX_GAP_STRING_GFF3,                        /* The Gap string used in GFF3 e.g. M23 D3 M10 I1 M20 */
  BLX_GAP_STRING_BAM_CIGAR,                   /* The cigar format used by SAM/BAM, e.g. 23M3D10M1I20M */
  BLX_GAP_STRING_ACEDB,                       /* Legacy acedb-style gaps string */
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
    char *sName_orig;   /* target name with original case. */
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
    char *fetchCommand; /* optional command which will be used for fetching sequence data */
    char *fetchArgs;    /* optional command which will be used for fetching sequence data */
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


static void           parseGffColumns(GString *line_string, const int lineNum, GList **seqList, GSList *supportedTypes, const IntRange* const refSeqRange, BlxGffData *gffData, GError **error);
static void           parseAttributes(char *attributes, GList **seqList, const int lineNum, BlxGffData *gffData, GError **error);
static void           parseTagDataPair(char *text, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseNameTag(char *data, char **sName, const int lineNum, GError **error);
static void           parseTargetTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseCommandTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error);
static void           parseSequenceTag(const char *text, const int lineNum, BlxGffData *gffData, GError **error);
static void           parseGapString(char *text, BlxGapFormat gapFormat, MSP *msp, const int resFactor, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GError **error);

static BlxStrand      readStrand(char *token, GError **error);
//static void           parseMspType(char *token, MSP *msp, GSList *supportedTypes, GError **error);
static const char*           parseCigarStringSection(const char *text, GapStringData *data);
static int            validateNumTokens(char **tokens, const int minReqd, const int maxReqd, GError **error);
//static void           validateMsp(const MSP *msp, GError **error);
static void           addGffType(GSList **supportedTypes, const char *name, const char *soId, BlxMspType blxType);
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
  freeAndNullString(&gffData->sName_orig);
  freeAndNullString(&gffData->source);
  freeAndNullString(&gffData->idTag);
  freeAndNullString(&gffData->parentIdTag);
  freeAndNullString(&gffData->gapString);
  freeAndNullString(&gffData->fetchCommand);
  freeAndNullString(&gffData->fetchArgs);
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


/* Create the list of supported types. Filter match types by the given seqType
 * (or pass BLXSEQ_NONE to include all supported types) */
GSList* blxCreateSupportedGffTypeList(const BlxSeqType seqType)
{
  GSList *supportedTypes = NULL;

  BlxMspType dnaType = BLXMSP_MATCH;
  BlxMspType proteinType = BLXMSP_MATCH;

  /* In DNA mode we need to exclude protein matches, so set them to invalid */
  if (seqType == BLXSEQ_DNA)
    proteinType = BLXMSP_INVALID;

  /* In protein mode we need to exclude nucleotide matches, so set them to invalid */
  if (seqType == BLXSEQ_PEPTIDE)
    dnaType = BLXMSP_INVALID;

  addGffType(&supportedTypes, "nucleotide_match", "SO:0000347", dnaType);
  addGffType(&supportedTypes, "primer_match", "SO:0001472", dnaType);
  addGffType(&supportedTypes, "cross_genome_match", "SO:0000177", dnaType);
  addGffType(&supportedTypes, "translated_nucleotide_match", "SO:0000181", dnaType);
  addGffType(&supportedTypes, "expressed_sequence_match", "SO:0000102", dnaType);
  addGffType(&supportedTypes, "cDNA_match", "SO:0000689", dnaType);
  addGffType(&supportedTypes, "EST_match", "SO:0000668", dnaType);
  addGffType(&supportedTypes, "UST_match", "SO:0001470", dnaType);
  addGffType(&supportedTypes, "RST_match", "SO:0001471", dnaType);

  addGffType(&supportedTypes, "protein_match", "SO:0000349", proteinType);
  addGffType(&supportedTypes, "protein_hmm_match", "SO:0001831", proteinType);

  addGffType(&supportedTypes, "match", "SO:0000343", BLXMSP_MATCH);
  addGffType(&supportedTypes, "match_part", "SO:0000039", BLXMSP_MATCH);
  addGffType(&supportedTypes, "match_set", "SO:0000038", BLXMSP_MATCH_SET);

  addGffType(&supportedTypes, "transcript", "SO:0000673", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "primary_transcript", "SO:0000185", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "processed_transcript", "SO:0000233", BLXMSP_TRANSCRIPT);
  addGffType(&supportedTypes, "mRNA", "SO:0000234", BLXMSP_TRANSCRIPT);

  /* gb10: hack to support old Sequence type from zmap which unfortunately still gets
   * passed sometimes (it's not even a real SO term) */
  addGffType(&supportedTypes, "Sequence", NULL, BLXMSP_TRANSCRIPT);

  addGffType(&supportedTypes, "CDS", "SO:0000316", BLXMSP_CDS);
  addGffType(&supportedTypes, "UTR", "SO:0000203", BLXMSP_UTR);
  addGffType(&supportedTypes, "exon", "SO:0000147", BLXMSP_EXON);
  addGffType(&supportedTypes, "intron", "SO:0000188", BLXMSP_INTRON);

  addGffType(&supportedTypes, "SNP", "SO:0000694", BLXMSP_VARIATION);
  addGffType(&supportedTypes, "SNV", "SO:0001483", BLXMSP_VARIATION);
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
  BlxMspType result = BLXMSP_NONE;

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

  /* If it wasn't found, treat as a basic feature... */
  if (result == BLXMSP_NONE)
    {
      result = BLXMSP_BASIC;
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
                     IntRange *refSeqRange,
                     GError **error)
{
  //DEBUG_ENTER("parseGff3Header [line=%d]", lineNum);

  GError *tmpError = NULL;

  /* Look for the "sequence-region" comment line, which tells us info about the reference
   * sequence. The format is as follows: ##sequence-region    qname qstart qend */
  char qName[MAXLINE + 1];
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;

  if (!strncasecmp(line_string->str, "##sequence-region", 17))
    {
      if (sscanf(line_string->str, "##sequence-region%s%d%d", qName, &qStart, &qEnd) < 3)
        {
          g_set_error(&tmpError, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_HEADER,
                      "Invalid format in sequence-region line '%s'\n", line_string->str);
        }

      //DEBUG_OUT("Found reference sequence name=%s [start=%d, end=%d]\n", qName, qStart, qEnd);

      /* If the ref seq name is already populated, check it's the same as the one we've just read */
      if (!tmpError && *refSeqName != '\0')
        {
          if (!stringsEqual(refSeqName, qName, FALSE))
            {
              g_set_error(&tmpError, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_SEQ,
                          "The sequence name '%s' in the GFF file does not match the reference sequence '%s'.\n",
                          qName, refSeqName);
            }
        }
      else if (!tmpError)
        {
          strcpy(refSeqName, qName);
        }

      if (!tmpError && refSeqRange)
        {
          if (!refSeqRange->isSet())
            {
              /* Range is currently unset, so set it */
              refSeqRange->set(qStart, qEnd);
            }
          else if (qStart > refSeqRange->max() || qEnd < refSeqRange->min())
            {
              /* GFF file range does not overlap the existing range, so we can't load this file */
              g_set_error(&tmpError, BLX_GFF3_ERROR, BLX_GFF3_ERROR_OUT_OF_RANGE,
                          "GFF file range [%d,%d] does not overlap the reference sequence range [%d,%d]",
                          qStart, qEnd, refSeqRange->min(), refSeqRange->max());
            }
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);

  //DEBUG_EXIT("parseGff3Header");
}


/* Get the dataType for the given source from the source-to-data-type mapping
 * stanza in the config file. Returns 0 if not found. */
static GQuark getBlxDataTypeFromSourceMapping(const char *source, GKeyFile *keyFile)
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


/* Get the dataType for the given source from the source stanza. Returns 0 if
 * not found. */
static GQuark getBlxDataTypeFromSource(const char *source, GKeyFile *keyFile)
{
  GQuark dataType = 0;

  if (keyFile && source && g_key_file_has_group(keyFile, source))
    {
      char *dataTypeName = g_key_file_get_string(keyFile, source, DATA_TYPE_TAG, NULL);

      if (dataTypeName)
        {
          dataType = g_quark_from_string(dataTypeName);
          g_free(dataTypeName);
        }
    }

  return dataType;
}


/* Find the default data type for this source. */
static GQuark getBlxDataTypeDefault(const char *source, GKeyFile *keyFile)
{
  GQuark dataType = 0;

  /* Check in the source stanza first */
  dataType = getBlxDataTypeFromSource(source, keyFile);

  /* If not there, check in the source-to-data-types mapping stanza */
  if (!dataType)
    dataType = getBlxDataTypeFromSourceMapping(source, keyFile);

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

  /* If no data type was specified in the gff, see if there is a default
   * data-type for this source */
  if (!dataType)
    dataType = getBlxDataTypeDefault(source, keyFile);

  /* A keyfile might not be supplied if the calling program is not interested
   * in the data-type data (i.e. the data-type data is currently only used to
   * supply values that are used in blixem, so are irrelevant to dotter). */
  if (!keyFile || !dataType)
    return result;

  if (dataType)
    {
      /* Check if it's already cached */
      std::map<GQuark, BlxDataType*>::iterator iter = g_dataTypes.find(dataType) ;

      if (iter != g_dataTypes.end())
        {
          result = iter->second ;
        }
      else
        {
          /* look it up in the config file and if we find it then create a new BlxDataType struct for it. */
          const gchar *typeName = g_quark_to_string(dataType);

          if (g_key_file_has_group(keyFile, typeName))
            {
              result = createBlxDataType();
              result->name = dataType;

              /* Get the values. They're all optional so just ignore any errors. */
              result->bulkFetch = keyFileGetCsv(keyFile, typeName, SEQTOOLS_BULK_FETCH, NULL);
              result->userFetch = keyFileGetCsv(keyFile, typeName, SEQTOOLS_USER_FETCH, NULL);
              result->optionalFetch = keyFileGetCsv(keyFile, typeName, SEQTOOLS_OPTIONAL_FETCH, NULL);

              /* Get the flags. Again, they're all optional. These calls update the
               * flag in place if it is found, or leave it at the pre-set default otherwise. */
              int flag = MSPFLAG_MIN + 1;
              for ( ; flag < MSPFLAG_NUM_FLAGS; ++flag)
                {
                  getMspFlag(keyFile, typeName, (MspFlag)flag, result);
                }

              /* Insert it into the table of data types */
              g_dataTypes[dataType] = result ;
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


BlxDataType* constructBlxDataType(BlxGffData *gffData, GHashTable *fetchMethods, GError **error)
{
  BlxDataType *result = NULL;

  /* If no data type is specified then try to construct one from the given fetch method */
  if (gffData && !gffData->dataType && gffData->fetchCommand)
    {
      /* Create a name for the new data type. Try using the source and fetchCommand */
      std::string name("");

      if (gffData->source)
        {
          name += gffData->source;
          name += "_";
        }

      name += gffData->fetchCommand;

      gffData->dataType = g_quark_from_string(name.c_str());

      /* Ok, now create the data type struct (only if it doesn't exist). */
      if (!g_dataTypes[gffData->dataType])
        {
          /* We'll set the bulk fetch method to a new fetch method based on the given
           * fetch command. */
          BlxFetchMethod *fetchMethod = createBlxFetchMethod(name.c_str());

          fetchMethod->name = gffData->dataType;  /// for now, use the same name as the datatype
          fetchMethod->mode = BLXFETCH_MODE_COMMAND;
          fetchMethod->separator = g_strdup(" ");
          fetchMethod->outputType = BLXFETCH_OUTPUT_GFF;
          fetchMethod->location = g_strdup(gffData->fetchCommand);

          if (gffData->fetchArgs)
            fetchMethod->args = g_strdup(gffData->fetchArgs);

          /* Add the fetch method to the output list of fetch methods */
          g_hash_table_insert(fetchMethods, GINT_TO_POINTER(fetchMethod->name), fetchMethod);

          /* Now create the new datatype */
          result = createBlxDataType();

          result->name = gffData->dataType;
          result->bulkFetch = g_array_sized_new(FALSE, TRUE, sizeof(GQuark), 1);
          g_array_append_val(result->bulkFetch, fetchMethod->name);

          /* Add the new datatype to the table of datatypes */
          g_dataTypes[gffData->dataType] = result;
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
                               GList *columnList,
			       GSList *styles,
                               const int resFactor,
                               GKeyFile *keyFile,
                               GHashTable *lookupTable,
                               GHashTable *fetchMethods,
			       GError **error)
{
  if (!gffData)
    {
      return;
    }

  GError *tmpError = NULL;

  /* Get the data type struct */
  BlxDataType *dataType = getBlxDataType(gffData->dataType, gffData->source, keyFile, &tmpError);
  reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);

  /* If no dataType is given then see if we can construct one from a given fetch command */
  if (!dataType && gffData->fetchCommand)
    {
      dataType = constructBlxDataType(gffData, fetchMethods, &tmpError);
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
    }

  GQuark filename = getFeatureFilename(gffData, keyFile, NULL);

  if (gffData->mspType > BLXMSP_NUM_TYPES)
    {
      /* "Invalid" MSP types, i.e. don't create a real MSP from these types. */

      if (gffData->mspType == BLXMSP_TRANSCRIPT)
        {
          /* For transcripts, although we don't create an MSP we do create a sequence */
          addBlxSequence(gffData->sName, gffData->sName_orig, gffData->idTag, gffData->qStrand,
                         dataType, gffData->source, seqList, columnList, gffData->sequence, NULL,
                         lookupTable, NULL, &tmpError);
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
                              columnList,
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
			      gffData->sName_orig,
			      gffData->sStart,
			      gffData->sEnd,
			      gffData->sStrand,
			      gffData->sequence,
                              filename,
                              lookupTable,
                              NULL,
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

          /* Parsing the gaps array may have separated the msp into multiple different msps if it
           * contained introns. We need to make sure the parent BlxSequence's msp list is sorted
           * correctly. We do this now, after it is complete, because it's more efficient so sort
           * the list once rather than on each iteration. */
          if (msp && msp->sSequence && msp->sSequence->mspList)
            msp->sSequence->mspList = g_list_sort(msp->sSequence->mspList, compareFuncMspPos);

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
                   GList *columnList,
                   GSList *supportedTypes,
                   GSList *styles,
                   const int resFactor,
                   GKeyFile *keyFile,
                   const IntRange* const refSeqRange,
                   GHashTable *lookupTable,
                   GHashTable *fetchMethods)
{
  //DEBUG_ENTER("parseGff3Body [line=%d]", lineNum);

  static int num_errors = 0 ;
  const int max_errors = 20 ; /* Limit the number of errors we report in case there are, say, thousands
                               * of lines we can't read */

  /* Parse the data into a temporary struct */
  BlxGffData gffData = {NULL, NULL, BLXMSP_NONE,
                        UNSET_INT, UNSET_INT, UNSET_INT, UNSET_INT, BLXSTRAND_NONE, UNSET_INT,
			NULL, NULL, BLXSTRAND_NONE,
                        UNSET_INT, UNSET_INT, NULL, NULL, NULL, NULL, BLX_GAP_STRING_INVALID, 0, 0, NULL, NULL};

  GError *error = NULL;
  parseGffColumns(line_string, lineNum, seqList, supportedTypes, refSeqRange, &gffData, &error);

  /* Create a blixem object based on the parsed data */
  if (!error)
    {
      createBlixemObject(&gffData, featureLists, lastMsp, mspList, seqList, columnList, styles, resFactor, keyFile, lookupTable, fetchMethods, &error);
    }

  if (error)
    {
      ++num_errors ;

      if (num_errors <= max_errors)
        {
          prefixError(error, "[line %d] Error parsing GFF data. ", lineNum);
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }
      else if (num_errors == max_errors + 1)
        {
          g_warning("Truncating error report (more than %d errors in reading GFF file)\n", max_errors);
        }
      else
        {
          g_error_free(error);
        }
    }

  //DEBUG_EXIT("parseGff3Body");
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
  gboolean status = TRUE;
  char seqName[MAXLINE + 1];

  /* Read the ref seq name (and optionally the coords) from the header line */
  int startCoord = UNSET_INT, endCoord = UNSET_INT;
  const int numFound = sscanf(line, ">%s %d %d", seqName, &startCoord, &endCoord);

  if (numFound < 1 || !seqName[0])
    {
      /* Didn't find name - this is required */
      status = FALSE;
      g_error("Error parsing data file: FASTA_SEQ_HEADER line \"%s\" is the wrong format; expected '>seq_name [start_coord end_coord].\n", line);
    }

  /* Trim out the name. The text can have additional info separated by '|' characters that we're not
   * interested in (at the moment) so just trim everything off after the first '|' char. */
  char *cp = strchr(seqName, '|');
  if (cp)
    *cp = 0;

  /* Set the name, if not already set. (gb10: we shouldn't really get here so should probably add
   * some error checking to make sure refSeqName is set so we can check we have the right
   * sequence. However for now we're flexible and if there's only one fasta sequence in the GFF
   * then we take that to be the reference sequence. If there are multiple in the GFF and no
   * reference sequence name is specified then we'll get in trouble here because we have no way of
   * telling which is the correct one.) */
  if (status && *refSeqName == 0)
    {
      strcpy(refSeqName, seqName);
    }
  else if (status && !stringsEqual(refSeqName, seqName, FALSE))
    {
      /* Not the sequence we're looking for so quit */
      status = FALSE;
    }

  /* Check if we also found coordinates in the header line. (There should be exactly three text
   * items if so) */
  if (status && numFound == 3 && refSeqRange)
    {
      refSeqRange->set(startCoord, endCoord);
    }

  /* Now allocate memory for the sequence data (if the sequence is not already populated) */
  if (status && *refSeq == NULL)
    {
      *readSeq = refSeq;
      *readSeqMaxLen = MAXLINE;
      **readSeq = (char*)g_malloc(*readSeqMaxLen + 1);
      *readSeqLen = 0;
    }

  if (status)
    {
      /* Update the parser state so that we proceed to parse the sequence data next. (Even if
       * we're not populating the ref seq, we still need to loop over these lines. Leaving the
       * readSeqLen as unset will mean that the fasta sequence parser will ignore the input.) */
      *parserState = FASTA_SEQ_BODY;
    }
  else
    {
      *parserState = FASTA_SEQ_IGNORE;
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
                            const IntRange* const refSeqRange,
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
          refSeqRange && refSeqRange->isSet())
        {
          IntRange featureRange(gffData->qStart, gffData->qEnd); /* make sure min < max */

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
  //DEBUG_ENTER("parseTagDataPair(text='%s')", text);

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
      else if (!strcmp(tokens[0], "command"))
        {
          parseCommandTag(tokens[1], lineNum, seqList, gffData, &tmpError);
        }
      else if (!strcmp(tokens[0], "Gap"))
        {
          /* This might have already been set if we have more than one type of gap string */
          if (gffData->gapString)
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
      else if (!strcmp(tokens[0], "gaps"))
        {
          /* This might have already been set if we have more than one type of gap string */
          if (!gffData->gapString)
            {
              gffData->gapString = g_strdup(tokens[1]);
              gffData->gapFormat = BLX_GAP_STRING_ACEDB;
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
      else if (!strcmp(tokens[0], DATA_TYPE_TAG))
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

  //DEBUG_EXIT("parseTagDataPair");
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
          gffData->sName = tokens[0] ? g_strdup(tokens[0]) : NULL;
          gffData->sName_orig = tokens[0] ? g_ascii_strup(gffData->sName, -1) : NULL;
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


/* A function to replace all occurances of substring 'from' with substring 'to' in the given
 * string. Taken from http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string */
void stringReplace(string& str, const string& from, const string& to)
{
  if(from.empty())
    return;

  size_t start_pos = 0;

  while((start_pos = str.find(from, start_pos)) != string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length();
    }
}


/* Unescape gff special characters from the given string. Returns a newly-allocated string which
 * should be free'd using g_free */
static char* unescapeGffString(const char *src)
{
  char *result = NULL;

  if (src)
    {
      string dest(src);

      stringReplace(dest, "%3D", "=");
      stringReplace(dest, "%3B", ";");
      stringReplace(dest, "%26", ",");
      stringReplace(dest, "%2C", "&");

      result = g_strdup(dest.c_str());
    }

  return result;
}


/* Parse the data from a 'command' tag */
static void parseCommandTag(char *data, const int lineNum, GList **seqList, BlxGffData *gffData, GError **error)
{
  if (data && !gffData->fetchCommand)
    {
      /* Text up to first space is the executable. Anything after that is the args. */
      int len = strlen(data);
      char *cp = strchr(data, ' ');

      if (cp)
        len = cp - data;

      gffData->fetchCommand = (char*)g_malloc(sizeof(char) * (len + 1));
      strncpy(gffData->fetchCommand, data, len);
      gffData->fetchCommand[len] = '\0';

      if (cp)
        {
          /* We need to unescape gff special characters */
          gffData->fetchArgs = unescapeGffString(cp);
        }
    }
}


/* Parse the data from the 'sequence' tag */
static void parseSequenceTag(const char *text, const int lineNum, BlxGffData *gffData, GError **error)
{
  gffData->sequence = g_strdup(text);
}


/* Required after reading in legacy acedb-style gaps array; new code assumes
 * the gaps are in the forward-strand order */
static void sortGapsArray(MSP *msp)
{
  if (msp && msp->gaps)
    {
      /* They should be ordered but may be in reverse order, so if the last one is before the
       * first one then just swap the order */
      CoordRange *first_range = (CoordRange*)(msp->gaps->data);
      CoordRange *last_range = (CoordRange*)g_slist_nth_data(msp->gaps, g_slist_length(msp->gaps) - 1);

      const gboolean qRev = last_range->qStart < first_range->qStart;
      const gboolean sRev = last_range->sStart < first_range->sStart;

      if (qRev != sRev)
        {
          msp->gaps = g_slist_reverse(msp->gaps);
        }
    }
}


/* Parse the data from the "Gap" string, which uses the CIGAR format, e.g. "M8 D3 M6 I1 M6".
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
  else if (gapFormat == BLX_GAP_STRING_ACEDB)
    {
      blxParseGaps(&text, msp, FALSE); /* legacy code for lecacy acedb-style gap string */
      sortGapsArray(msp);
      return;
    }

  /* If we have the forward strand of either sequence, start at the min coord
   * and increase values as we progress through the cigar string; if we have the
   * reverse strand, start at the max coord and decrease. */
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
  const int qDirection = (qForward ? 1 : -1);
  const int sDirection = (sForward ? 1 : -1);

  /* Start at one beyond the edge of the range, because it will be incremented (or decremented if
   * direction is reverse) when we construct the first range. */
  int q = qForward ? msp->qRange.min() - 1 : msp->qRange.max() + 1;
  int s = sForward ? msp->sRange.min() - 1 : msp->sRange.max() + 1;

  GError *tmpError = NULL;

  GapStringData gapStringData = {gapFormat, &msp, qDirection, sDirection, resFactor, &q, &s,
                                 featureLists, lastMsp, mspList, seqList, NULL};

  const char *cp = text;

  while (cp && *cp)
    {
      const char *cp_new = parseCigarStringSection(cp, &gapStringData);

      if (tmpError)
        {
          prefixError(tmpError, "Error parsing gap string '%s'. ", cp);
        }

      cp = cp_new;
    }

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }
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
      result = convertStringToInt(text); /* uses atoi, which will ignore characters after */
      break;

    default:
      g_warning("Invalid gap string format\n");
      break;
    };

  return result;
}


/* Get the operator part of a gap string section, e.g. if the text is "M76"
 * then this returns 'M' */
static char getCigarStringSectionOperator(const char *text, BlxGapFormat gapFormat, const char **cp_out)
{
  char result = 0;
  const char *cp = text;

  switch (gapFormat)
    {
    case BLX_GAP_STRING_GFF3:
      result = *cp;

      /* Move cp on to the start of the next section in the cigar, i.e. next alpha char */
      for (++cp ; cp && *cp && !isalpha(*cp); ++cp);

      break;

    case BLX_GAP_STRING_BAM_CIGAR:
      {
        for ( ; cp && *cp && !isalpha(*cp); ++cp); /* find first alphabetic character */

        if (cp)
          result = *cp;

        /* Move cp on to the start of the next section in the cigar, i.e. next digit */
        for (cp++ ; cp && *cp && !isdigit(*cp); ++cp);

        break;
      }

    default:
      g_warning("Invalid gap string format\n");
      break;
    };

  if (cp_out)
    *cp_out = cp;

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

  CoordRange *newRange = new CoordRange;
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
  /* Intron. Create a separate msp under the same sequence. Add the new msp to the same parent
   * BlxSequence but pass sort=FALSE because we'll sort the resulting msp list once we've finished. */
  MSP *msp = *data->msp;
  MSP *newMsp = copyMsp(msp, data->featureLists, data->lastMsp, data->mspList, TRUE, FALSE);

  /* end current msp at the current coords */
  if (data->qDirection > 0)
    msp->qRange.setMax(*data->q);
  else
    msp->qRange.setMin(*data->q);

  if (data->sDirection > 0)
    msp->sRange.setMax(*data->s);
  else
    msp->sRange.setMin(*data->s);

  /* start new msp at new coords */
  *data->q += data->qDirection * numNucleotides;

  if (data->qDirection > 0)
    newMsp->qRange.setMin(*data->q + 1);
  else
    newMsp->qRange.setMax(*data->q - 1);

  if (data->sDirection > 0)
    newMsp->sRange.setMin(*data->s + 1);
  else
    newMsp->sRange.setMax(*data->s - 1);

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
static gboolean validateCigarOperator(char op, BlxGapFormat gapFormat)
{
  gboolean result = FALSE;

  switch (op)
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
static const char* parseCigarStringSection(const char *text,
                                     GapStringData *data)
{
  /* Get the digit part of the string, which indicates the number of display coords (peptides in peptide matches,
   * nucleotides in nucelotide matches). */
  const int numPeptides = getCigarStringSectionLen(text, data->gapFormat);
  int numNucleotides = numPeptides * data->resFactor;

  const char *cp = text;
  char op = getCigarStringSectionOperator(text, data->gapFormat, &cp);

  /*! \todo If the operator is not valid for this type of cigar string
   * then we should set the error and return. However, for historic
   * reasons we have allowed invalid operators in the GFF Gap string,
   * so continue allowing this for now and just give a warning. */
  if (!validateCigarOperator(op, data->gapFormat))
    g_warning("Invalid operator '%c' for gap string format '%d'", op, data->gapFormat);

  switch (op)
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
      g_set_error(&data->error, BLX_GFF3_ERROR, BLX_GFF3_ERROR_INVALID_CIGAR_FORMAT, "Invalid operator '%c' in cigar string.\n", op);
      break;
    };

  return cp;
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
static void addGffType(GSList **supportedTypes, const char *name, const char *soId, BlxMspType blxType)
{
  BlxGffType *gffType = new BlxGffType;

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
      if ((*gffType)->name)
        g_free((*gffType)->name);

      if ((*gffType)->soId)
        g_free((*gffType)->soId);

      delete *gffType;
      *gffType = NULL;
    }
}
