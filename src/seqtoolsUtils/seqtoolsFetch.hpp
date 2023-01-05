/*  File: seqtoolsFetch.hpp
 *  Author: Gemma Guest, 2015-12-21
 *  Copyright [2018-2023] EMBL-European Bioinformatics Institute
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
 * Description: Utility structs for sequence-fetching code
 *----------------------------------------------------------------------------
 */

#ifndef DEF_SEQTOOLS_FETCH_HPP
#define DEF_SEQTOOLS_FETCH_HPP

#include <glib.h>

/* These are the supported fetch modes. ***If you add anything here, also add it in fetchModeStr*** */
typedef enum
  {
#ifdef PFETCH_HTML
    BLXFETCH_MODE_HTTP,
    BLXFETCH_MODE_PIPE,
#endif
    BLXFETCH_MODE_SOCKET,
    BLXFETCH_MODE_WWW,
    BLXFETCH_MODE_SQLITE,
    BLXFETCH_MODE_COMMAND,
    BLXFETCH_MODE_INTERNAL,
    BLXFETCH_MODE_NONE,

    BLXFETCH_NUM_MODES /* must be last in list */
  } BlxFetchMode;


/* output types for fetch modes. *** if you add anything here, also add it in outputTypeStr *** */
typedef enum
{
  BLXFETCH_OUTPUT_INVALID,
  BLXFETCH_OUTPUT_RAW,      /* raw sequence data, separated by newlines */
  BLXFETCH_OUTPUT_FASTA,    /* sequence data in FASTA format */
  BLXFETCH_OUTPUT_EMBL,     /* the sequence's EMBL entry */
  BLXFETCH_OUTPUT_LIST,     /* a list of named columns is returned */
  BLXFETCH_OUTPUT_GFF,      /* a new gff for re-parsing is returned */

  BLXFETCH_NUM_OUTPUT_TYPES
} BlxFetchOutputType;




/* struct to hold info about a fetch method */
typedef struct _BlxFetchMethod
{
  GQuark name;                      /* fetch method name */
  BlxFetchMode mode;                /* the type of fetch method */

  char *location;                   /* e.g. url, script, command, db location etc. */
  char *node;                       /* for socket fetch mode */
  int port;                         /* for socket and http/pipe fetch modes */
  char *cookie_jar;                 /* for http/pipe fetch mode */
  char *proxy;                      /* for http/pipe fetch mode */
  char *args;                       /* arguments/query/request */
  GArray *columns;                  /* for db-fetch, the list of columns the query will populate */

  char *separator;                  /* separator when combining multiple sequence names into a list */
  GArray *errors;                   /* array of messages (as GQuarks) that indicate that an error occurred, e.g. "no match" */
  BlxFetchOutputType outputType;    /* the output format to expect from the fetch command */
} BlxFetchMethod;

#endif /* DEF_SEQTOOLS_FETCH_HPP */


BlxFetchMethod* createBlxFetchMethod(const char *fetchName);
