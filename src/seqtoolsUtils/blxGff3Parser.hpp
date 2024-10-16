/*  File: blxGff3parser.h
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
 * Description: GFF v3 file format parser.
 *----------------------------------------------------------------------------
 */

#ifndef BLX_GFF_P_H
#define BLX_GFF_P_H

#include <gtk/gtk.h>
#include <seqtoolsUtils/blxmsp.hpp>

/* We follow glib convention in error domain naming:
 *          "The error domain is called <NAMESPACE>_<MODULE>_ERROR" */
#define BLX_GFF_ERROR "BLX_GFF_ERROR"


/* This struct contains info about a GFF type: its name and SO id, and the blixem
 * type that it corresponds to */
typedef struct _BlxGffType
  {
    char *name;         /* the type name as given in the GFF file, e.g. nucleotide_match */
    char *soId;         /* the type's SO Id, as given in the GFF file, e.g. SO:0000347 */
    BlxMspType blxType; /* the Blixem object type that the GFF type corresponds to */
  } BlxGffType;



/* This enum is to record the type of data currently being parsed by the parser. An input file can
 * contain multiple types of data. The start of a new section of data is indicated by a header
 * line beginning with a hash and a known type name, e.g. "# exblx_x" or "##gff-version   3"
 *
 * FS and SFS (Feature Series) data type names have "FS" or "SFS" followed by a specific format
 * type, e.g. "# SFS type=SEG".
 *
 * For some data types, additional data is included in the header line(s) as well as in the 'body'
 * section below it. For these, there are two data types in the enum, one postfixed with _HEADER
 * and one with _BODY. When the header line is detected the initial type is set to the _HEADER enum,
 * and when we are finished processing the header information it is set to _BODY, so that we go on
 * to process the body of the data on the next loop through. For types with no information in the
 * header, there is only a _BODY enum.
 */
typedef enum
  {
    PARSER_START,                  /* indicates that we haven't started processing yet */
    PARSER_ERROR,                  /* indiates an error state */

    GFF_3_HEADER,                  /* GFF 3 header */
    GFF_3_BODY,                    /* GFF 3 data */

    FASTA_SEQ_HEADER,              /* FASTA sequence header */
    FASTA_SEQ_BODY,                /* Sequence data in FASTA format */
    FASTA_SEQ_IGNORE,              /* A FASTA sequence we're not interested in */

    EXBLX_BODY,                    /* Old style sequence entries. */
    SEQBL_BODY,                    /* Old style sequence entries. */
    EXBLX_X_BODY,                  /* New style sequence entries with gaps and match strand. (_X stands for eXtended.) */
    SEQBL_X_BODY,                  /* New style sequence entries with gaps and match strand. (_X stands for eXtended.) */

    FS_HSP_BODY,                   /* feature-series HSP data */

    FS_GSP_HEADER,                 /* feature-series GSP data header */
    FS_GSP_BODY,                   /* feature-series GSP data */

    FS_GFF_BODY,                   /* feature-series GFF data */

    FS_SEG_BODY,                   /* feature-series segment data */

    FS_XY_HEADER,                  /* feature-series XY data header */
    FS_XY_BODY,                    /* feature-series XY data */

    FS_SEQ_HEADER,                 /* feature-series sequence data header */
    FS_SEQ_BODY                    /* feature-series sequence data */
  } BlxParserState ;



/* External functions */
void parseGff3Header(const int lineNum,
                     MSP **lastMsp,
                     MSP **mspList,
                     BlxParserState *parserState,
                     GString *line_string,
                     GList **seqList,
                     char *refSeqName,
                     IntRange *refSeqRange,
                     GError **error);

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
                   GHashTable *fetchMethods);

void parseFastaSeqHeader(char *line, const int lineNum,
                         char **refSeq, char *refSeqName, IntRange *refSeqRange,
                         char ***readSeq, int *readSeqLen, int *readSeqMaxLen,
                         BlxParserState *parserState);


GSList*                            blxCreateSupportedGffTypeList(const BlxSeqType seqType);
void                               blxDestroyGffTypeList(GSList **supportedTypes);
BlxDataType*                       getBlxDataType(GQuark dataType, const char *source, GKeyFile *keyFile, GError **error);


#endif /* BLX_GFF_P_H */
