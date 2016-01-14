/*  File: seqtoolsFetch.hpp
 *  Author: Gemma Guest, 2015-12-21
 *  Copyright (c) 2015 Genome Research Ltd
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


BlxFetchMethod* createBlxFetchMethod(const char *fetchName) ;
