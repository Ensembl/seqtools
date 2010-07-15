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


#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))


/* Debug logging macros. #define DEBUG to enable debug output. */
#ifdef DEBUG
#define DEBUG_OUT(format, args...) debugLogLevel(0); printf(format, ##args);
#else
#define DEBUG_OUT(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_ENTER(format, args...) debugLogLevel(1); printf("-->"); printf(format, ##args); printf("\n");
#else
#define DEBUG_ENTER(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_EXIT(format, args...) debugLogLevel(-1); printf("<--"); printf(format, ##args); printf("\n");
#else
#define DEBUG_EXIT(format, args...)
#endif



/* Color strings that can be passed to create a GdkColor */
#define BLX_YELLOW	      "#ffff00" 
#define BLX_DARK_YELLOW	      "#c0c000"
#define BLX_LIGHT_YELLOW      "#fafad2"
#define BLX_CYAN	      "#00ffff"
#define BLX_DARK_CYAN	      "#008080"
#define BLX_LIGHT_CYAN	      "#80FFF2"
#define BLX_VERY_LIGHT_BLUE   "#c1ccdb"
#define BLX_BLUE	      "#0000ff"
#define BLX_DARK_BLUE	      "#000080"
#define BLX_PALE_STEEL_BLUE   "#78b4f0"
#define BLX_LIGHT_STEEL_BLUE  "#b0c4ff"
#define BLX_STEEL_BLUE	      "#5e76bb"
#define BLX_ROYAL_BLUE	      "#4169e1"
#define BLX_LIGHT_SKY_BLUE    "#87cefa"
#define BLX_SKY_BLUE	      "#a0b8c8"
#define BLX_GREY	      "#bebebe"
#define BLX_DARK_GREY	      "#929292"
#define BLX_LIGHT_GREY	      "#d3d3d3"
#define BLX_BLACK	      "#000000"
#define BLX_WHITE	      "#ffffff"
#define BLX_RED		      "#ff0000"
#define BLX_LIGHTER_RED	      "#ff5353"
#define BLX_LIGHT_RED	      "#ff7373"
#define BLX_DARK_RED	      "#800000"
#define BLX_VERY_DARK_RED     "#A00000"
#define BLX_ORANGE	      "#ffa500"
#define BLX_LIGHT_ORANGE      "#FFCF55"
#define BLX_DUSK_ORANGE	      "#D7BA94"
#define BLX_GREEN	      "#00ff00"
#define BLX_DARK_GREEN	      "#00bb00"
#define BLX_VERY_DARK_GREEN   "#015800"
#define BLX_TURQUOISE	      "#40e0d0"
#define BLX_DARK_TURQUOISE    "#32aca0"
#define BLX_PURPLE	      "#9370db"
#define BLX_PLUM	      "#dda0dd"
#define BLX_THISTLE	      "#d8bfd8"
#define BLX_ORANGE_RED	      "#ff4500"
#define BLX_PALE_GREEN	      "#C1FFC1" 
#define BLX_LAWN_GREEN	      "#7cfc00"
#define BLX_LIGHT_CORAL	      "#f08080"
#define BLX_LIGHT_SALMON      "#ffa07a"
#define BLX_BURLYWOOD	      "#deb887"
#define BLX_TAN		      "#d2b48c"

/* Function pointer for callback functions used by widgets on dialog boxes */
typedef void (*BlxResponseCallback)(GtkWidget *widget, const gint responseId, gpointer data);

/* This struct holds generic callback information. It can be stored as a widget
 * property and called on the widget on request (e.g. by a dialog when it applies changes). */
typedef struct _CallbackData
  {
    BlxResponseCallback func;	  /* Callback function to be called */
    gpointer data;                /* User data to pass to the callback function */
  } CallbackData;

/* Custom dialog response types, which can be used in addition to the default types specified by GtkResponseType */
typedef enum
  {
    BLX_RESPONSE_FORWARD, 
    BLX_RESPONSE_BACK
  } BlxResponseType;


/* Struct to hold a pair of coordinate ranges, one holding ref seq coords, the other match coords */
typedef struct _CoordRange
  {
    int qStart;     /* start coord on the reference (Query) seq */
    int qEnd;       /* end coord on the reference (Query) seq */
    int sStart;     /* start coord on the match (Subject) seq */
    int sEnd;       /* end coord on the match (Subject) seq */
  } CoordRange;
  
  

GdkDrawable*	      widgetGetDrawable(GtkWidget *widget);
void		      widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable);
gboolean	      widgetGetHidden(GtkWidget *widget);
void		      widgetSetHidden(GtkWidget *widget, const gboolean hidden);
void		      hideUserHiddenWidget(GtkWidget *widget, gpointer data);
void		      widgetClearCachedDrawable(GtkWidget *widget, gpointer data);

gboolean	      onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer data);
GtkWidget*	      createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean enableCopyPaste, const gboolean showWhenPrinting);
GdkDrawable*	      createBlankPixmap(GtkWidget *widget);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
gboolean              getColorFromString(const char *colorStr, GdkColor *color, GError **error);
void		      getSelectionColor(GdkColor *origColor, GdkColor *result);
void		      getDropShadowColor(GdkColor *origColor, GdkColor *result);
void		      convertToGrayscale(GdkColor *origColor, GdkColor *result);
void		      adjustColorBrightness(GdkColor *origColor, const double factor, GdkColor *result);

gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsSnp(const MSP const *msp);
gboolean	      mspIsPolyATail(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);

gboolean              mspHasSName(const MSP const *msp);
gboolean              mspHasSSeq(const MSP  const *msp);
gboolean              mspHasSCoords(const MSP const *msp);
gboolean              mspHasSStrand(const MSP const *msp);

void		      getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

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
						const gboolean displayRev,
						const IntRange const *dnaIdxRange);

int		      convertDnaIdxToDisplayIdx(const int dnaIdx, 
						const BlxSeqType displaySeqType,
						const int frame,
						const int numFrames, 
						const gboolean displayRev,
						const IntRange const *dnaIdxRange,
						int *baseNum);

int		      mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
BlxStrand	      mspGetRefStrand(const MSP const *msp);
BlxStrand	      mspGetMatchStrand(const MSP const *msp);
const char*           mspGetMatchSeq(const MSP const *msp);
char*		      mspGetSeqName(const MSP *msp);
const IntRange const* mspGetRefCoords(const MSP const *msp);
const IntRange const* mspGetMatchCoords(const MSP const *msp);
int		      mspGetQStart(const MSP const *msp);
int		      mspGetQEnd(const MSP const *msp);
int		      mspGetSStart(const MSP const *msp);
int		      mspGetSEnd(const MSP const *msp);
int		      mspGetQRangeLen(const MSP const *msp);
int		      mspGetSRangeLen(const MSP const *msp);
int		      mspGetMatchSeqLen(const MSP const *msp);
const GdkColor*       mspGetColor(const MSP const *msp, GArray *defaultColors, const BlxSequence *blxSeq, const gboolean selected, const gboolean usePrintColors, const gboolean fill);

char                  getStrandAsChar(const BlxStrand strand);

BlxSequence*          createEmptyBlxSequence(char *fullName);
BlxSequence*          findBlxSequence(GList *seqList, const char *reqdName, const BlxStrand reqdStrand);

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
				       const gboolean displayRev, 
				       const int numFrames,
				       const IntRange const *refSeqRange);

int		      getEndDnaCoord(const IntRange const *displayRange, 
				     const int frame,
				     const BlxSeqType displaySeqType, 
				     const gboolean displayRev, 
				     const int numFrames,
				     const IntRange const *refSeqRange);

int		      gapCoord(const MSP *msp, 
			       const int qIdx, 
			       const int numFrames, 
			       const BlxStrand strand, 
			       const gboolean displayRev,
			       const gboolean showUnalignedSeq,
			       const gboolean limitUnalignedBases,
			       const gboolean numUnalignedBases);

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
int		      convertStringToInt(const char *inputStr);
char*		      abbreviateText(const char *inputStr, const int maxLen);
gboolean              stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive);

GtkWidget*	      showMessageDialog(const char *title,  
					const char *messageText,
					GtkWidget *parent,
					const int initWidth,
					const int maxHeight,
					const gboolean wrapText,
                                        const gboolean useMarkup,
					PangoFontDescription *fontDesc,
                                        GtkTextBuffer **textBuffer);

GtkWidget*		createScrollableTextView(const char *messageText,
						 const gboolean wrapText,
						 PangoFontDescription *fontDesc,
                                                 const gboolean useMarkup,
						 int *height,
                                                 GtkTextBuffer **textBufferOut);
				    
void		      widgetSetCallbackData(GtkWidget *widget, BlxResponseCallback callbackFunc, gpointer callbackData);
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

gint		      runConfirmationBox(GtkWidget *blxWindow, char *title, char *messageText);

const char*	      getSeqVariantName(const char *longName);
char*		      getSeqShortName(const char *longName);
const char*	      sequenceGetFullName(const BlxSequence *seq);
const char*	      sequenceGetVariantName(const BlxSequence *seq);
const char*	      sequenceGetDisplayName(const BlxSequence *seq);
const char*	      sequenceGetShortName(const BlxSequence *seq);
const IntRange const* sequenceGetRange(const BlxSequence *seq);
int		      sequenceGetLength(const BlxSequence *seq);

void		      destroyBlxSequence(BlxSequence *seq);

void		      defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void		      popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);

void		      prefixError(GError *error, char *prefixStr, ...);
void                  postfixError(GError *error, char *formatStr, ...);
void                  reportAndClearIfError(GError **error, GLogLevelFlags log_level);

void		      intrangeSetValues(IntRange *range, const int val1, const int val2);

#ifdef DEBUG
void		      debugLogLevel(const int increaseAmt);
#endif

void                  drawHighlightBox(GdkDrawable *drawable, const GdkRectangle const *rect, const gint minWidth, GdkColor *color);



void    gtk_text_buffer_insert_markup             (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup);

void    gtk_text_buffer_insert_markup_with_tag    (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup_with_tag       (GtkTextBuffer *buffer,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup                (GtkTextBuffer *buffer,
                                                   const gchar   *markup);

#endif /* _utilities_h_included_ */
