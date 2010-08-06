/*  File: blxGff3parser_.h
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

#ifndef BLX_GFF_P_H
#define BLX_GFF_P_H

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>


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



/* External functions */
void parseGff3Header(const int lineNum,
                     MSP **lastMsp, 
		     MSP **mspList, 
		     BlxParserState *parserState, 
		     char *opts, 
		     GString *line_string, 
		     GList **seqList,
                     char *refSeqName);

void parseGff3Body(const int lineNum,
                   MSP **lastMsp, 
		   MSP **mspList, 
		   BlxParserState *parserState, 
		   char *opts, 
		   GString *line_string, 
		   GList **seqList,
                   GSList *supportedTypes,
                   GSList *styles);

void parseFastaSeqHeader(char *line, const int lineNum,
                         char **refSeq, char *refSeqName,
                         char ***readSeq, int *readSeqLen, int *readSeqMaxLen,
                         BlxParserState *parserState);

#endif /* BLX_GFF_P_H */
