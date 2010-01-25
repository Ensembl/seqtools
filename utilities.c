/*
 *  utilities.c
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */

#include "SeqTools/utilities.h"



/* Utility to return the centre value of the given range (rounded down if an odd number) */
int getRangeCentre(IntRange *range)
{
  return range->min + round(((double)(range->max - range->min)) / 2.0);
}


/* Utility to calculate how many digits are in an integer */
int numDigitsInInt(int val)
{
  int count = 0;
  while (val > 0)
    {
      ++count;
      val /= 10;
    }
  
  return count;
}


/* Create a GdkColor from a pixel value (e.g. GDK_YELLOW) */
GdkColor getGdkColor(gulong colour)
{
  GdkColor result = {colour};
  return result;
}


gboolean mspIsExon(const MSP const *msp)
{
  return (msp && msp->score == -1);
}

gboolean mspIsIntron(const MSP const *msp)
{
  return (msp && msp->score == -2);
}

gboolean mspIsBlastMatch(const MSP const *msp)
{
  return (msp && msp->score > 0);
}


/* Sorts the given values so that val1 is less than val2 if forwards is true,
 * or the reverse if forwards is false. */
void sortValues(int *val1, int *val2, gboolean forwards)
{
  if ((forwards && *val1 > *val2) || (!forwards && *val1 < *val2))
    {
      int temp = *val1;
      *val1 = *val2;
      *val2 = temp;
    }
}


/* Returns the upper and lower extents of the query and subject sequence ranges in 
 * the given MSP. Any of the return values can be NULL if they are not required. */
void getMspRangeExtents(const MSP *msp, int *qSeqMin, int *qSeqMax, int *sSeqMin, int *sSeqMax)
{
  if (qSeqMin)
    *qSeqMin = msp->qstart < msp->qend ? msp->qstart : msp->qend;
  
  if (qSeqMax)
    *qSeqMax = msp->qstart < msp->qend ? msp->qend : msp->qstart;
  
  if (sSeqMin)
    *sSeqMin = msp->sstart < msp->send ? msp->sstart : msp->send;
  
  if (sSeqMax)
    *sSeqMax = msp->sstart < msp->send ? msp->send : msp->sstart;
}


/* Returns the upper and lower extents of the query and subject sequence ranges in
 * the given gap range. Any of the return values can be passed as NULL if they are not required. */
void getSMapMapRangeExtents(SMapMap *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax)
{
  if (qRangeMin)
    *qRangeMin = range->r1 < range->r2 ? range->r1 : range->r2;
  
  if (qRangeMax)
    *qRangeMax = range->r1 < range->r2 ? range->r2 : range->r1;
  
  if (sRangeMin)
    *sRangeMin = range->s1 < range->s2 ? range->s1 : range->s2;
  
  if (sRangeMax)
    *sRangeMax = range->s1 < range->s2 ? range->s2 : range->s1;
}


char convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType)
{
  char result = (seqType == BLXSEQ_PEPTIDE) ? toupper(charToConvert) : tolower(charToConvert);
  return result;
}


/* Returns the base at the given index in the reference sequence. If
 * 'complement' is true, the result is complemented. */
char getRefSeqBase(char *refSeq,
		   const int qIdx, 
		   const gboolean complement, 
		   IntRange *refSeqRange,
		   const BlxSeqType seqType)
{
  char result = ' ';
  
  if (!refSeqRange || (qIdx >= refSeqRange->min && qIdx <= refSeqRange->max))
    {
      char base = refSeq[qIdx - 1];
      base = convertBaseToCorrectCase(base, seqType);
      
      if (!complement)
	{
	  result = base;
	}
      else
	{
	  switch (base)
	  {
	    case 'c':
	      result = 'g';
	      break;
	    case 'g':
	      result = 'c';
	      break;
	    case 'a':
	      result = 't';
	      break;
	    case 't':
	      result = 'a';
	      break;
	    default:
	      result = ' ';
	      break;
	  }
	}
    }
  
  return result;
}


/* Given an index into a peptide sequence and a reading frame number, return
 * the index into the DNA sequence that will give the equivalent DNA base */
int convertPeptideToDna(const int peptideIdx, const int frame, const int numFrames)
{
  return (peptideIdx * numFrames) - numFrames + frame;
}


/* Given an index into a dna sequence and a reading frame number, return
 * the index into the protein sequence that will give the equivalent peptide */
int convertDnaToPeptide(const int dnaIdx, const int numFrames)
{
  return ceil((double)dnaIdx / (double)numFrames);
}


/* Simple utility to determine whether the given index is within the given range.
 * Returns false if the given index is an unset int or the given range is null */
gboolean indexWithinRange(const int idx, const IntRange const *range)
{
  return (idx != UNSET_INT && range != NULL && idx >= range->min && idx <= range->max);
}


/* Get the start coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getStartDnaCoord(const IntRange const *displayRange, 
		     const BlxSeqType displaySeqType, 
		     const gboolean reverse, 
		     const int numReadingFrames)
{
  int result = reverse ? displayRange->max : displayRange->min;
  
  if (displaySeqType == BLXSEQ_PEPTIDE)
    {
      /* Take the first reading frame to get the first DNA coord for this peptide (or the last frame if reversed) */
      int frame = reverse ? numReadingFrames : 1;
      result = convertPeptideToDna(result, frame, numReadingFrames);
    }
  
  return result;
}


/* Get the end coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getEndDnaCoord(const IntRange const *displayRange, 
		   const BlxSeqType displaySeqType, 
		   const gboolean reverse, 
		   const int numReadingFrames)
{
  int result = reverse ? displayRange->min : displayRange->max;
  
  if (displaySeqType == BLXSEQ_PEPTIDE)
    {
      /* Take the last reading frame to get the last DNA coord for this peptide (or the first frame if reversed) */
      int frame = reverse ? 1 : numReadingFrames;
      result = convertPeptideToDna(result, frame, numReadingFrames);
    }
  
  return result;
}




