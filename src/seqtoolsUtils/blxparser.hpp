/*  File: blxparser.h
 *  Author: Gemma Barson, 2010-09-02
 *  Copyright [2018-2022] EMBL-European Bioinformatics Institute
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
 *  Description: Parser for Blixem and Dotter feature input files.
 *
 *               GFF v3 is the preferred input file format. This parser also
 *               includes code to parse the old exblx or feature-series file
 *               types (e.g. as output from the MSPcrunch program).
 *----------------------------------------------------------------------------
 */

void           parseFS(MSP **MSPlist, FILE *file, BlxBlastMode *blastMode, GArray* featureLists[], GList **seqList, GList *columnList, GSList *supportedTypes, GSList *styles,
		       char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, GKeyFile *keyFile, GHashTable *lookupTable, GHashTable *fetchMethods, GError **error) ;

void           parseBuffer(MSP **MSPlist, const char *buffer_in, BlxBlastMode *blastMode, GArray* featureLists[], GList **seqList, GList *columnList, GSList *supportedTypes, GSList *styles,
                           char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, GKeyFile *keyFile, GHashTable *lookupTable, GError **error) ;
gboolean       blxParseGaps(char **text, MSP *msp, const gboolean hasGapsTag);
