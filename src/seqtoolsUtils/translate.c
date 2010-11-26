/*  File: translate.c
 *  Author: sre
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
 * Description: functions for translating nucleic acid sequence
 *              
 * Exported functions: see blixem_.h
 *              
 * HISTORY:
 * Last edited: Sep 10 16:23 2009 (edgrif)
 * Created: Tue Jan 12 11:27:29 1993 (SRE)
 * CVS info:   $Id: translate.c,v 1.9 2010-11-03 15:23:56 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <seqtoolsUtils/utilities.h>
#include <iupac.h>
#include <glib.h>
#include <ctype.h>
#include <string.h>


/* THIS FILE NEEDS RENAMING TO SOMETHING LIKE utils.c */

/* Function: Translate(char *seq, char **code)
 * 
 * Given a ptr to the start of a nucleic acid sequence, and a genetic code, translate the sequence
 * into amino acid sequence.
 * 
 * code is an array of 65 strings, representing the translations of the 64 codons, arranged in
 * order AAA, AAC, AAG, AAU, ..., UUA, UUC, UUG, UUU.  '*' or '***' is used to represent
 * termination codons, usually. The final string, code[64], is the code for an ambiguous amino
 * acid.
 *
 * Because of the way space is allocated for the amino acid sequence, the amino acid strings
 * cannot be longer than 3 letters each. (I don't foresee using anything but the single- and
 * triple- letter codes.)
 * 
 * Returns a ptr to the translation string on success, or NULL on failure.
 */
char *blxTranslate(const char *seq, char **code)
{
  char *aaseq = NULL ;					    /* RETURN: the translation */
  int   codon;						    /* index for codon         */
  char *aaptr;						    /* ptr into aaseq */
  int   i;
  
  if (seq && *seq)
    {
      aaseq = (char *)g_malloc(strlen(seq) + 1) ;

      for (aaptr = aaseq ; *seq != '\0' && *(seq+1) != '\0' && *(seq+2) != '\0'; seq += 3)
	{
	  /* calculate the lookup value for this codon */
	  codon = 0;
	  for (i = 0; i < 3; i++)
	    {
	      codon *= 4;

	      switch (*(seq + i))
		{
		case 'A': case 'a':             break;
		case 'C': case 'c': codon += 1; break;
		case 'G': case 'g': codon += 2; break;
		case 'T': case 't': codon += 3; break;
		case 'U': case 'u': codon += 3; break;
		default: codon = 64; break;
		}

	      if (codon == 64)
		break;
	    }

	  strcpy(aaptr, code[codon]);
	  aaptr += strlen(code[codon]);
	}
    }  

  return aaseq ;
}


/* All these calls need rationalising into a single function with options. */

/* revComplement.c
 * 
 * Reverse complement of a IUPAC character string
 * 
 */
/* Ratinalise this with my func. below..... */
char *revComplement(char *comp, char *seq)
{
  long  bases;
  char *bckp, *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (comp == NULL)
    return NULL;
  if (seq == NULL)
    return NULL;

  bases = strlen(seq);

  fwdp = comp;
  bckp = seq + bases -1;
  for (pos = 0; pos < bases; pos++)
    {
      c = *bckp;
      c = toupper(c);

      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);

      if (idx > IUPACSYMNUM)
	{
	  *fwdp = '\0';
	  return NULL;
	}
      else
	{
	  *fwdp = iupac[idx].symcomp;
	}

      if (islower(*bckp))
	*fwdp = tolower(*fwdp);

      fwdp++;
      bckp--;
    }
  *fwdp = '\0';

  return comp;
}
  

/* blxComplement.c
 * 
 * Just complement of a IUPAC character string
 * 
 * Note that it overwrites calling string!!!! (revcomp doesn't)
 */
void blxComplement(char *seq)
{
  char *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (seq == NULL)
    return ;

  fwdp = seq;
  for (pos = 0; pos < strlen(seq); pos++)
    {
      c = toupper(*fwdp);

      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);

      if (idx > IUPACSYMNUM)
	{
	  *fwdp = '\0';
	  return;
	}

      else c = iupac[idx].symcomp;

      if (islower(*fwdp))
	*fwdp = tolower(c);
      else
	*fwdp = c;

      fwdp++;
    }

  *fwdp = '\0';

  return ;
}
  


