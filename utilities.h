/*
 *  utilities.h
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */
 
#ifndef _utilities_h_included_
#define _utilities_h_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>


#define UNSET_INT  -1

/* Colour strings that can be passed to gdk_color_parse to create a GdkColor */
#define GDK_YELLOW	      "#FFFF00" 
#define GDK_DARK_YELLOW	      "#808000"
#define GDK_CYAN	      "#00FFFF"
#define GDK_DARK_CYAN	      "#008080"
#define GDK_BLUE	      "#0000FF"
#define GDK_DARK_BLUE	      "#000080"
#define GDK_LIGHT_STEEL_BLUE  "#6EB4E6"
#define GDK_STEEL_BLUE	      "#4682B4"
#define GDK_ROYAL_BLUE	      "#4169E1"
#define GDK_GREY	      "#BEBEBE"
#define GDK_DARK_GREY	      "#808080"
#define GDK_BLACK	      "#000000"
#define GDK_WHITE	      "#FFFFFF"
#define GDK_RED		      "#FF0000"
#define GDK_DARK_RED	      "#800000"
#define GDK_ORANGE	      "#FFA500"
#define GDK_GREEN	      "#00FF00"
#define GDK_DARK_GREEN	      "#00BB00"
#define GDK_TURQUOISE	      "#40E0D0"


GdkDrawable*	      widgetGetDrawable(GtkWidget *widget);
void		      widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable);
gboolean	      widgetGetHidden(GtkWidget *widget);
void		      widgetSetHidden(GtkWidget *widget, const gboolean hidden);
void		      hideUserHiddenWidget(GtkWidget *widget, gpointer data);
void		      widgetClearCachedDrawable(GtkWidget *widget);

gboolean	      onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer data);
GtkWidget*	      createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean enableCopyPaste, const gboolean showWhenPrinting);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
GdkColor	      getGdkColor(const char *colour);

gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);

void		      getMspRangeExtents(const MSP *msp, int *qSeqMin, int *qSeqMax, int *sSeqMin, int *sSeqMax);
void		      getSMapMapRangeExtents(SMapMap *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

int		      getRangeLength(const IntRange const *range);
int		      getRangeCentre(const IntRange const *range);
gboolean	      valueWithinRange(const int value, const IntRange const *range);
void		      boundsLimitValue(int *value, const IntRange const *range);

int		      convertPeptideToDna(const int peptideIdx, const int frame, const int baseNum, const int numFrames);
int		      convertDnaToPeptide(const int dnaIdx, const int frame, const int numFrames, int *baseNum);
char		      convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

int		      mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
Strand		      mspGetRefStrand(const MSP const *msp);

void		      addMspToHashTable(GHashTable *hashTable, MSP *msp, char *hashKey);

int                   roundNearest(const double val);

char		      getRefSeqBase(char *refSeq, 
				    const int qIdx, 
				    const gboolean complement, 
				    IntRange *refSeqRange,
				    const BlxSeqType seqType);

int		      getStartDnaCoord(const IntRange const *displayRange, 
				       const int frame,
				       const BlxSeqType displaySeqType, 
				       const gboolean reverse, 
				       const int numReadingFrames);

int		      getEndDnaCoord(const IntRange const *displayRange, 
				     const int frame,
				     const BlxSeqType displaySeqType, 
				     const gboolean reverse, 
				     const int numReadingFrames);

int		      getMatchIdxFromDisplayIdx(MSP *msp,
						const int displayIdx,
						const int qFrame,
						const Strand qStrand,
						const gboolean rightToLeft,
						const BlxSeqType seqType,
						const int numFrames);

int		      gapCoord(const MSP *msp, 
			       const int qIdx, 
			       const int numFrames, 
			       const Strand strand, 
			       const gboolean rightToLeft, 
			       int *nearestIdx);

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
int		      convertStringToInt(const char *inputStr);

#endif /* _utilities_h_included_ */
