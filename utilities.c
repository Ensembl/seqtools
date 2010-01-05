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


/* Returns true if the given MSP is a "fake" MSP (i.e. one that that is used just
 * to display the reference sequence and does not represent a real match.) */
gboolean mspIsFake(const MSP const *msp)
{
  /* This is not a great way to identify fake msp's but it'll do for now. */
  return (msp && msp->score == 0);
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