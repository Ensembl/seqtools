/*  File: blxview.h
 *  Author: Erik Sonnhammer, 92-02-20
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
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
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: External interface to blixem window code, this interface
 *              is used both by xace and the stand alone blixem.
 * HISTORY:
 * Last edited: Aug 21 13:57 2009 (edgrif)
 * * Aug 26 16:57 1999 (fw): added this header
 * Created: Thu Aug 26 16:57:17 1999 (fw)
 * CVS info:   $Id: blxview.h,v 1.41 2010-11-01 15:31:01 gb10 Exp $
 *-------------------------------------------------------------------
 */
#ifndef DEF_BLXVIEW_H
#define DEF_BLXVIEW_H

#include <gtk/gtk.h>
#include <SeqTools/blxmsp.h>

#ifdef ACEDB
#include <wh/regular.h>
#endif

/* Only used in pephomolcol.c, would be good to get rid of this.... */
#define FULLNAMESIZE               255


/* blixem can use either efetch (default) or a pfetch server to get
 * sequences, to use pfetch node/port information must be specified. */
typedef struct
{
  char *net_id ;
  int port ;
} PfetchParams ;


/* Blixem can read a number of file formats which are documented in:
 * 
 * The following defines give keyword strings for these file formats
 * which should be used in code that writes these files.
 * 
 *  */
#define BLX_GAPS_TAG               "Gaps"
#define BLX_DESCRIPTION_TAG        "Description"
#define BLX_SEQUENCE_TAG           "Sequence"


/* Function to show blixem window, can be called from any application. */
gboolean                            blxview(char *refSeq, 
                                            char *refSeqName,
	                                    int start, 
                                            int qOffset, 
                                            GList* featureLists[],
                                            MSP *msplist, 
                                            GList *seqList, 
                                            GSList *supportedTypes,
                                            char *opts, 
	                                    PfetchParams *pfetch, 
                                            char *align_types, 
                                            gboolean External) ;


#endif /*  !defined DEF_BLXVIEW_H */
