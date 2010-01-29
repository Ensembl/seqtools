/*
 *  utilities.h
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>


#define UNSET_INT  -1

/* These colours are decimal version of, e.g. "ffff00" for yellow, i.e. R=255, G=255, B=0.
 * They can be used to define a GdkColor */
#define GDK_YELLOW	    16776960 
#define GDK_DARK_YELLOW	    8421376
#define GDK_CYAN	    65535
#define GDK_DARK_CYAN	    32896
#define GDK_BLUE	    255
#define GDK_DARK_BLUE	    128
#define GDK_GREY	    12500670
#define GDK_DARK_GREY	    4210752
#define GDK_BLACK	    0
#define GDK_WHITE	    16777215


GdkDrawable*	      widgetGetDrawable(GtkWidget *widget);
void		      widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable);

gboolean	      onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer data);
GtkWidget*	      createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean enableCopyPaste, const gboolean showWhenPrinting);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
GdkColor	      getGdkColor(gulong colour);

gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);

void		      getMspRangeExtents(const MSP *msp, int *qSeqMin, int *qSeqMax, int *sSeqMin, int *sSeqMax);
void		      getSMapMapRangeExtents(SMapMap *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

int		      getRangeCentre(IntRange *range);
gboolean	      indexWithinRange(const int idx, const IntRange const *range);

int		      convertPeptideToDna(const int peptideIdx, const int frame, const int numFrames);
int		      convertDnaToPeptide(const int dnaIdx, const int numFrames);
char		      convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

char		      getRefSeqBase(char *refSeq, 
				    const int qIdx, 
				    const gboolean complement, 
				    IntRange *refSeqRange,
				    const BlxSeqType seqType);

int		      getStartDnaCoord(const IntRange const *displayRange, 
				       const BlxSeqType displaySeqType, 
				       const gboolean reverse, 
				       const int numReadingFrames);

int		      getEndDnaCoord(const IntRange const *displayRange, 
				     const BlxSeqType displaySeqType, 
				     const gboolean reverse, 
				     const int numReadingFrames);


