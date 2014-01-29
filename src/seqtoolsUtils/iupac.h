/*  File: iupac.h
 *  Author: Fred Wobus, 1999-08-26
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
 * Description: Globally defines the IUPAC symbols for nucleic acid sequence
 *              iupac stuff from Sean Eddy's libsquid
 *----------------------------------------------------------------------------
 */

#ifndef _iupac_h_included_
#define _iupac_h_included_


int PAM120[23][23] =
  {
    /*A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X */ 
    { 3, -3, -1,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  1,  0, -1},
    {-3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -1,  0, -1},
    {-1, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  4,  1, -1},
    { 0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  5,  3, -1},
    {-3, -4, -5, -7,  9, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -4, -6, -1},
    {-1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  1,  5, -1},
    { 0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2, -1, -2, -8, -5, -3,  3,  5, -1},
    { 1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  1, -1, -1},
    {-3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3, -1, -3,  2,  2, -1},
    {-1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -2, -2, -1},
    {-3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -3, -4, -3, -3, -2,  1, -3, -2, -1},
    {-2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  1,  0, -1},
    {-2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -3, -1, -1},
    {-4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  0, -7, -1,  8, -5, -3, -4, -1,  4, -3, -4, -5, -1},
    { 1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -3, -2, -3, -5,  6,  1, -1, -7, -6, -2, -1,  0, -1},
    { 1, -1,  1,  0,  0, -2, -1,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  1,  0, -1},
    { 1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4, -1,  2,  4, -6, -3,  0,  1, -1, -1},
    {-7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6, -1, -7, -2, -6, 12, -2, -8, -5, -6, -1},
    {-4, -5, -2, -5, -1, -5, -5, -6, -1, -2, -2, -5, -4,  4, -6, -3, -3, -2,  8, -3, -2, -4, -1},
    { 0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -2, -2, -1},
    { 1, -1,  4,  5, -4,  1,  3,  1,  2, -2, -3,  1, -3, -4, -1,  1,  1, -5, -2, -2,  6,  4, -1},
    { 0,  0,  1,  3, -6,  5,  5, -1,  2, -2, -2,  0, -1, -5,  0,  0, -1, -6, -4, -2,  4,  6, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
  };

#define NA 23   /* Not A residue - equal to X */
int aa_atob[]	/* ASCII-to-binary translation table */
	= {
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,
	NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
	15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};


struct iupactype
{
  char       sym;               /* character representation */
  char       symcomp;           /* complement (regular char */
  char       code;              /* my binary rep */
  char       comp;              /* binary encoded complement */
};


/* Binary encoding of the IUPAC code for nucleotides
 *
 *    four-bit "word", permitting rapid degenerate matching
 *         A  C  G  T/U
 *         0  0  1  0
 */
#define NTA 8
#define NTC 4
#define NTG 2
#define NTT 1
#define NTU 1
#define NTN 15                  /* A|C|G|T */
#define NTR 10                  /* A|G */
#define NTY 5                   /* C|T */
#define NTM 12                  /* A|C */
#define NTK 3                   /* G|T */
#define NTS 6                   /* C|G */
#define NTW 9                   /* A|T */
#define NTH 13                  /* A|C|T */
#define NTB 7                   /* C|G|T */
#define NTV 14                  /* A|C|G */
#define NTD 11                  /* A|G|T */
#define NTGAP 16                /* GAP */
#define NTEND 0                 /* null string terminator */

                                /* ntmatch(): bitwise comparison of two nuc's */
                                /* note that it's sensitive to the order;
                                   probe may be degenerate but target should
                                   not be */

				/* IUPAC code translations */
				/* note: sequence chars are UPPER CASE */
				/* in order of occurrence frequency;
				   to speed up loops a bit */

/* Define all of the valid characters in a DNA/RNA sequence, and their
 * complements */
/* residue:                 inverse:
 * A = adenine              T
 * C = cytosine             G
 * G = guanine              C
 * T = thymine (DNA only)   A
 * U = uracil (RNA only)    A
 * R = G/A (purine)         Y = C/T
 * Y = T/C (pyrimidine)     R = G/A
 * K = G/T (keto)           M = C/A
 * M = A/C (amino)          K = T/G
 * S = G/C (strong bonds)   S = C/G
 * W = A/T (weak bonds)     W = T/A
 * B = G/T/C (all but A)    V = C/A/G
 * D = G/A/T (all but C)    H = C/T/A
 * H = A/C/T (all but G)    D = T/G/A
 * V = G/C/A (all but T)    B = C/G/T
 * N = A/G/C/T (any)        N = A/G/C/T
 */ 
#define IUPACSYMNUM 17
struct iupactype iupac[IUPACSYMNUM] =
  {
    {'A', 'T', NTA, NTT},
    {'C', 'G', NTC, NTG},
    {'G', 'C', NTG, NTC},
    {'T', 'A', NTT, NTA},
    {'U', 'A', NTU, NTA},
    {'N', 'N', NTN, NTN},
    {'-', '-', NTGAP, NTGAP},
    {'R', 'Y', NTR, NTY},
    {'Y', 'R', NTY, NTR},
    {'M', 'K', NTM, NTK},
    {'K', 'M', NTK, NTM},
    {'S', 'S', NTS, NTS},
    {'W', 'W', NTW, NTW},
    {'H', 'D', NTH, NTD},
    {'B', 'V', NTB, NTV},
    {'V', 'B', NTV, NTB},
    {'D', 'H', NTD, NTH}
  };


const char *stdcode1[65] = {
  "K",				/* AAA */
  "N",				/* AAC */
  "K",				/* AAG */
  "N",				/* AAU */
  "T",				/* ACA */
  "T",				/* ACC */
  "T",				/* ACG */
  "T",				/* ACU */
  "R",				/* AGA */
  "S",				/* AGC */
  "R",				/* AGG */
  "S",				/* AGU */
  "I",				/* AUA */
  "I",				/* AUC */
  "M",				/* AUG */
  "I",				/* AUU */
  "Q",				/* CAA */
  "H",				/* CAC */
  "Q",				/* CAG */
  "H",				/* CAU */
  "P",				/* CCA */
  "P",				/* CCC */
  "P",				/* CCG */
  "P",				/* CCU */
  "R",				/* CGA */
  "R",				/* CGC */
  "R",				/* CGG */
  "R",				/* CGU */
  "L",				/* CUA */
  "L",				/* CUC */
  "L",				/* CUG */
  "L",				/* CUU */
  "E",				/* GAA */
  "D",				/* GAC */
  "E",				/* GAG */
  "D",				/* GAU */
  "A",				/* GCA */
  "A",				/* GCC */
  "A",				/* GCG */
  "A",				/* GCU */
  "G",				/* GGA */
  "G",				/* GGC */
  "G",				/* GGG */
  "G",				/* GGU */
  "V",				/* GUA */
  "V",				/* GUC */
  "V",				/* GUG */
  "V",				/* GUU */
  "*",				/* UAA */
  "Y",				/* UAC */
  "*",				/* UAG */
  "Y",				/* UAU */
  "S",				/* UCA */
  "S",				/* UCC */
  "S",				/* UCG */
  "S",				/* UCU */
  "*",				/* UGA */
  "C",				/* UGC */
  "W",				/* UGG */
  "C",				/* UGU */
  "L",				/* UUA */
  "F",				/* UUC */
  "L",				/* UUG */
  "F",				/* UUU */
  "X",				/* unknown */
};




const char *stdcode3[65] = {
  "Lys",			/* AAA */
  "Asn",			/* AAC */
  "Lys",			/* AAG */
  "Asn",			/* AAU */
  "Thr",			/* ACA */
  "Thr",			/* ACC */
  "Thr",			/* ACG */
  "Thr",			/* ACU */
  "Arg",			/* AGA */
  "Ser",			/* AGC */
  "Arg",			/* AGG */
  "Ser",			/* AGU */
  "Ile",			/* AUA */
  "Ile",			/* AUC */
  "Met",			/* AUG */
  "Ile",			/* AUU */
  "Gln",			/* CAA */
  "His",			/* CAC */
  "Gln",			/* CAG */
  "His",			/* CAU */
  "Pro",			/* CCA */
  "Pro",			/* CCC */
  "Pro",			/* CCG */
  "Pro",			/* CCU */
  "Arg",			/* CGA */
  "Arg",			/* CGC */
  "Arg",			/* CGG */
  "Arg",			/* CGU */
  "Leu",			/* CUA */
  "Leu",			/* CUC */
  "Leu",			/* CUG */
  "Leu",			/* CUU */
  "Glu",			/* GAA */
  "Asp",			/* GAC */
  "Glu",			/* GAG */
  "Asp",			/* GAU */
  "Ala",			/* GCA */
  "Ala",			/* GCC */
  "Ala",			/* GCG */
  "Ala",			/* GCU */
  "Gly",			/* GGA */
  "Gly",			/* GGC */
  "Gly",			/* GGG */
  "Gly",			/* GGU */
  "Val",			/* GUA */
  "Val",			/* GUC */
  "Val",			/* GUG */
  "Val",			/* GUU */
  "***",			/* UAA */
  "Tyr",			/* UAC */
  "***",			/* UAG */
  "Tyr",			/* UAU */
  "Ser",			/* UCA */
  "Ser",			/* UCC */
  "Ser",			/* UCG */
  "Ser",			/* UCU */
  "***",			/* UGA */
  "Cys",			/* UGC */
  "Trp",			/* UGG */
  "Cys",			/* UGU */
  "Leu",			/* UUA */
  "Phe",			/* UUC */
  "Leu",			/* UUG */
  "Trp",			/* UUU */
  "XXX",			/* unknown */
};


#endif /* _iupac_h_included_ */
