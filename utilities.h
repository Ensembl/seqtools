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
#define GDK_YELLOW	      "#ffff00" 
#define GDK_DARK_YELLOW	      "#c0c000"
#define GDK_CYAN	      "#00ffff"
#define GDK_DARK_CYAN	      "#008080"
#define GDK_BLUE	      "#0000ff"
#define GDK_DARK_BLUE	      "#000080"
#define GDK_LIGHT_STEEL_BLUE  "#78b4f0"
#define GDK_STEEL_BLUE	      "#5e76bb"
#define GDK_ROYAL_BLUE	      "#4169e1"
#define GDK_GREY	      "#bebebe"
#define GDK_DARK_GREY	      "#929292"
#define GDK_BLACK	      "#000000"
#define GDK_WHITE	      "#ffffff"
#define GDK_RED		      "#ff0000"
#define GDK_DARK_RED	      "#800000"
#define GDK_ORANGE	      "#ffa500"
#define GDK_GREEN	      "#00ff00"
#define GDK_DARK_GREEN	      "#00bb00"
#define GDK_TURQUOISE	      "#40e0d0"
#define GDK_DARK_TURQUOISE    "#32aca0"


/* This struct holds generic callback information. It can be stored as a widget
 * property and called on the widget on request (e.g. by a dialog when it applies changes). */
typedef struct _CallbackData
  {
    GtkCallback func;	  /* Callback function to be called */
    gpointer data;	  /* User data to pass to the callback function */
  } CallbackData;


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
char		      convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

int		      convertDisplayIdxToDnaIdx(const int inputIdx, 
						const BlxSeqType inputIdxType,
						const int frame, 
						const int baseNum, 
						const int numFrames,
						const gboolean rightToLeft,
						const IntRange const *dnaIdxRange);

int		      convertDnaIdxToDisplayIdx(const int dnaIdx, 
						const BlxSeqType displaySeqType,
						const int frame,
						const int numFrames, 
						const gboolean rightToLeft,
						const IntRange const *dnaIdxRange,
						int *baseNum);

int		      mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
Strand		      mspGetRefStrand(const MSP const *msp);
Strand		      mspGetSubjectStrand(const MSP const *msp);

void		      addMspToHashTable(GHashTable *hashTable, MSP *msp, char *hashKey);

int                   roundNearest(const double val);
int		      roundToValue(const int inputVal, const int roundTo);

char		      getRefSeqBase(char *refSeq, 
				    const int qIdx, 
				    const gboolean complement, 
				    const IntRange const *refSeqRange,
				    const BlxSeqType seqType);

int		      getStartDnaCoord(const IntRange const *displayRange, 
				       const int frame,
				       const BlxSeqType displaySeqType, 
				       const gboolean rightToLeft, 
				       const int numFrames,
				       const IntRange const *refSeqRange);

int		      getEndDnaCoord(const IntRange const *displayRange, 
				     const int frame,
				     const BlxSeqType displaySeqType, 
				     const gboolean rightToLeft, 
				     const int numFrames,
				     const IntRange const *refSeqRange);

int		      gapCoord(const MSP *msp, 
			       const int qIdx, 
			       const int numFrames, 
			       const Strand strand, 
			       const gboolean rightToLeft, 
			       int *nearestIdx);

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
int		      convertStringToInt(const char *inputStr);
char*		      abbreviateText(const char *inputStr, const int maxLen);

void		      showMessageDialog(const char *title,  
					const char *messageText,
					GtkWidget *parent,
					const int initWidth,
					const int maxHeight,
					const gboolean wrapText,
					PangoFontDescription *fontDesc);

void		      widgetSetCallbackData(GtkWidget *widget, GtkCallback callbackFunc, gpointer callbackData);
void		      widgetCallCallback(GtkWidget *widget);
void		      widgetCallAllCallbacks(GtkWidget *widget, gpointer data);
void		      onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data);

void		      setPrimaryClipboardText(const char *text);
void		      setDefaultClipboardText(const char *text);
void		      requestPrimaryClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);
void		      requestDefaultClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);

int		      parseMatchLine(const char *inputText,
				     char **matchNameOut,
				     int *matchStartOut, 
				     int *matchEndOut, 
				     int *matchLenOut);
			
GList*		      parseMatchList(const char *inputText);

GList*		      findStringInList(GList *list, const char *seqName);


#endif /* _utilities_h_included_ */
